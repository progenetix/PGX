package lib::ReadFiles;

use Data::Dumper;
use Math::Random qw(random_normal);
use lib::PlotMakeParameters;
use lib::Helpers;

require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(
	read_probefile
	read_frequencyfile
	read_segmentfile
	read_file_to_split_array
	read_webfile_to_split_array
);

################################################################################

sub read_probefile {

=pod

Expects:
  - a standard Progenetix style probe file

  ID  chro  pos log2
  cnvi0111187 17  35295593  0.0859121900
  cnvi0111188 8 65499402  -0.1438023000
  cnvi0111189 2 177061178 -0.0113166000
  cnvi0111190 5 70255894  0.0463862400
  ...

Returns:
  - a list reference of genome position / value objects:
    [
      {
        no => __integer__,          # 1 -> n
        probe_id => __string__,
        reference_name => __string__,
        position => __integer__,
        value => __long__,
      },
      {
      ...
      },
    ]

=cut

########    ####    ####    ####    ####    ####    ####    ####    ####    ####

	my $pgx = shift;
	my $probeF = shift;
	my $probeT = shift;
	$probeT ||= 'probedata';

	$pgx->{$probeT} = [];
	my @randomV;

	if (! -f $probeF) { return $pgx->{$probeT} }

	my $numfactor = 1;
	if (
		$pgx->{parameters}->{'reverse'} =~ /y/i
		&&
		$probeT !~ /frac/i
	) { $numfactor = -1 }

	if ($pgx->{parameters}->{plot_adjust_baseline} =~ /[123456789]/) {
	if ($probeT !~ /fracb/i) {
	  $pgx->{parameters}->{probebaseline} =   $pgx->{parameters}->{plot_adjust_baseline} } }
	  
	if ($pgx->{debug}) {
		print "probes from $probeF\n" }

	my $probeData = _f2l( $probeF );
	shift @$probeData;

	my $probeNo = @$probeData;
	if ($pgx->{debug}) {
		print "$probeNo probes from $probeF\n" }

	my $i = 0;

	foreach (@$probeData) {

		$i++;
		my (
			$probe_id,
			$reference_name,
			$position,
			$value,
		) = split (/\s/, $_, 5);

		$probe_id =~  s/[^\w\-\,]/_/g;
		$reference_name =~ s/[^\dxXyY]//;
		$reference_name =~ s/^23$/X/;
		$reference_name =~ s/^24$/Y/;
		$position = sprintf "%.0f", $position;  # due to some erroneous .5 in-between pos.
		$value = sprintf "%.4f", ($pgx->{parameters}->{probebaseline} + $value);

		if ($reference_name !~ /^\w\d?$/)             { next }
		if ($position       !~ /^\d{1,9}$/)           { next }
		if ($value          !~ /^\-?\d+?(\.\d+?)?$/)  { next }

		push(
			@{ $pgx->{$probeT} },
			{
				no => $i,
				probe_id => $probe_id,
				reference_name => $reference_name,
				position => $position,
				value => $numfactor * $value,
			}
		);
	}
  
	# random values
	if ($pgx->{parameters}->{simulated_probes} =~ /y/i ) {
		my @randomV = random_normal(scalar @{ $pgx->{$probeT} }, 0, 0.25);
		foreach my $n (0..$#{ $pgx->{$probeT} }) {
			$pgx->{$probeT}->[$n]->{value}  =   $randomV[$n];
		}
	}

  return $pgx;

}

################################################################################

sub read_frequencyfile {

	no warnings 'uninitialized';

=podmd

#### Expects:

* a standard tab-delimited .pgxseg frequency file
  - an additional header may exist
* header (optional)

```
# plotpars;color_var_dup_hex=#EE4500;color_var_del_hex=#09F911
...
```

#### Returns:

* a pgx "frequencysets" list of objects

```
...
{
	"index" : 2840,
	"id" : "22:4000000-5000000",
	"reference_name" : "22",
	"start" : 4000000,
	"end" : 5000000,
	"size" : 1000000,
	"gain_frequency" : 13.208,
	"loss_frequency" : 33.962
},
...
```

=cut

########    ####    ####    ####    ####    ####    ####    ####    ####    ####

	my $pgx = shift;
	my $fF = shift;
	
	if (! -f $fF) { return $pgx }
	
	$pgx->{frequencymaps} = [];
	
	my ($header, $table) = read_file_to_split_array($fF);
	my $headerValues = _objectify_header($header);	
	$pgx->{pgxfileheader} = $headerValues;
	
	my %colOrder = (
		group_id => 0,
		reference_name => 1,
		start => 2,
		end => 3,
		gain_frequency => 4,
		loss_frequency => 5,
		"index" => 6
	);
	
	if ($table->[0]->[ $colOrder{reference_name} ] !~ /^([12]\d?)|X|Y/i) {
		shift @$table }

	#first get the collation_id values
	my $f_maps_keyed = {};
	
	foreach my $segment (@$table) {
		if ($segment->[ $colOrder{start} ] =~ /\w/) {
			$f_maps_keyed->{ $segment->[ $colOrder{group_id} ] } =  {
				name => $segment->[ $colOrder{group_id} ],
				id => $segment->[ $colOrder{group_id} ],
				interval_count => 0,
				labels => [ ],
				intervals => [ ]			
			}
		}	
	}
	
	foreach my $segment (@$table) {		
		my %segVals =  ();
		foreach (keys %colOrder) {
			$segVals{$_} = $segment->[$colOrder{$_}];
			$segVals{$_} =~ s/\s//g;
		};
		
		$segVals{ size } = $segVals{ end } - $segVals{ start };
		$segVals{ id } = $segVals{ reference_name }.":".$segVals{ start }."-".$segVals{ end };
		
		push(@{ $f_maps_keyed->{ $segment->[ $colOrder{collation_id} ] }->{ intervals } }, \%segVals);
	}
	
	foreach my $f_map_k (keys %$f_maps_keyed) {
	
		if ($pgx->{parameters}->{min_group_no} > 0) {
			if ($headerValues->{groups}->{$f_map_k}->{sample_count} < $pgx->{parameters}->{min_group_no}) {
				next;			
			}
		}

		$f_maps_keyed->{ $f_map_k }->{ interval_count } = scalar @{ $f_maps_keyed->{ $f_map_k }->{ intervals } };
				
		my $label = $f_map_k;
		
		if ($headerValues->{groups}->{$f_map_k}->{label} =~ /.../) {
			$label = $headerValues->{groups}->{$f_map_k}->{label};
			$f_maps_keyed->{ $f_map_k }->{name} = $headerValues->{groups}->{$f_map_k}->{label};
		}
		if ($headerValues->{groups}->{$f_map_k}->{sample_count} > 0) {
			$label .= ' ('.$headerValues->{groups}->{$f_map_k}->{sample_count}.')';		
		}
		
		push(
			@{ $f_maps_keyed->{ $f_map_k }->{labels} },
			{
			  label_text => $label,
			  label_link => q{},
			  label_color => random_hexcolor(),
			}		
		);
		
		push(
			@{ $pgx->{frequencymaps} },
			$f_maps_keyed->{ $f_map_k }		
		);
	}
		
	return $pgx;

}




################################################################################

sub read_segmentfile {

	no warnings 'uninitialized';

=podmd

#### Expects:

* a standard tab-delimited Progenetix segments file
  - an additional header may exist
  - only first 5 columns are necessary
  - column 5 (mean) can be empty or dot, if column 7 exists and contains status value
  - undefined fields in existing columns are replaced with the "." character
* header (optional)
  - the `biosample_id` parameter is required to assign values (e.g. group labels) to samples
  - parameter=value pairs are semicolon-separated
  - supported tags
    * `biosample_id` is required and has to correspond to column 1 values
    * `group_id`
    * `group_label`

```
# plotpars;color_var_dup_hex=#EE4500;color_var_del_hex=#09F911
# biosample_id=GSM481286;group_id=NCIT:C4017;group_label=Ductal Breast Carcinoma
# biosample_id=GSM481418;group_id=NCIT:C3059;group_label=Glioma
biosample_id	chro	start	stop	mean	variant_type	probes
GSM481286	1	742429	7883881	-0.1594	699	DEL
GSM481286	2	115673158	115705254	-0.3829	8	DEL
GSM481286	3	115722621	119771659	0.167	424	DUP
GSM481286	4	119776776	162617092	0.4168	1587	DUP
GSM481418	5	162621657	165278686	.	.	DUP
GSM481418	6	165280711	167221337	.	.	DUP
GSM481418	7	167248788	168289603	0.6784	.	DUP	
...
```

#### Returns:

* a list reference of genome CNV objects inside `$pgx->{segmentdata}`
* header derived information, e.g. group labels per sample in `$pgx->{pgxfileheader}`

```
[
	{
		no : __integer__,    # 1 -> n
		callset_id : __string__,
		reference_name : __string__,
		start : __integer__,
		end : __integer__,
		variant_type : __string__,
		info :  {
			value : __long__,
			svlen : __integer__,
			probes : __integer__,
			assembly_id : __string__,
			experiment_type : __string__,
		},
	},
	{
	...
	},
]
```

=cut

########    ####    ####    ####    ####    ####    ####    ####    ####    ####

	my $pgx = shift;
	my $segmentsF = shift;
	my $segmentsT = shift;
	$segmentsT ||= 'segmentdata';
	$pgx->{$segmentsT} = [];

	if (! -f $segmentsF) { return $pgx }

	my ($header, $table) = read_file_to_split_array($segmentsF);
	my $headerValues = _objectify_header($header);	
	$pgx->{pgxfileheader} = $headerValues;

	my $numfactor = 1;
	if (
		$pgx->{parameters}->{'reverse'} =~ /y/i
		&&
		$segmentsT !~ /frac/i
	) { $numfactor = -1 }

	if ($pgx->{parameters}->{plot_adjust_baseline} =~ /[123456789]/) {
		$pgx->{parameters}->{segbaseline} = $pgx->{parameters}->{plot_adjust_baseline} }

	my %colOrder = (
		callset_id => 0,
		reference_name => 1,
		start => 2,
		end => 3,
		value => 4,
		probes => 6,
		variant_type => 5
	);

	if ($pgx->{parameters}->{format_inputfiles} =~ /tcga/i) {
		$colOrder{value} = 5;
		$colOrder{probes} = 4;
	};

	if ($pgx->{parameters}->{format_inputfiles} =~ /dgv/i) {
		$colOrder{value} = 4;
		$colOrder{probes} = 5;
	};

	if ($table->[0]->[1] !~ /^([12]\d?)|X|Y/i) {
		shift @$table }

	my $segNo = @$table;
	if ($pgx->{debug}) {
		print "$segNo segments from $segmentsF\n" }

	my $i = 0;

	foreach my $segment (@$table) {
		
		my %segVals =  ();
		foreach (keys %colOrder) {
			$segVals{$_} = $segment->[$colOrder{$_}];
			$segVals{$_} =~ s/\s//g;
		};

		$segVals{callset_id} = _norm_sample_id($segVals{callset_id});

		$segVals{reference_name} =~ s/[^\dxXyY]//g;
		$segVals{reference_name} =~ s/^23$/X/;
		$segVals{reference_name} =~ s/^24$/Y/;
		if ($segVals{reference_name} !~ /^\w\d?$/) { next }

		$segVals{start} = sprintf "%.0f", $segVals{start};
		$segVals{end} = sprintf "%.0f", $segVals{end};
		$segVals{probes} =~ s/[^\d]//g;

		if ($segVals{start} !~ /^\d{1,9}$/)           { next }
		if ($segVals{end}   !~ /^\d{1,9}$/)           { next }
		if (
			$segVals{value} !~ /^\-?\d+?(\.\d+?)?$/
			&&
			$segVals{variant_type} !~ /^(DUP)|(DEL)|(LOH)$/
		)  { next }

		if ($segVals{value} =~ /^\-?\d+?(\.\d+?)?$/) {
			$segVals{value} = sprintf "%.4f", $segVals{value} }

		my $varStatus = '_NS_';

		if ($segmentsT !~ /fracb/i) {

			# baseline adjustment
			if ($segVals{value} =~ /^\-?\d+?(\.\d+?)?$/) {
				$segVals{value}	+=   $pgx->{parameters}->{segbaseline} }

			if ($segVals{variant_type} =~ /^(DUP)|(DEL)|(LOH)$/) {
				$varStatus = $segVals{variant_type} }
			elsif ($segVals{value} >= $pgx->{parameters}->{cna_gain_threshold}) {
				$varStatus = 'DUP' }
			elsif ($segVals{value} <= $pgx->{parameters}->{cna_loss_threshold}) {
				$varStatus = 'DEL' }
			else {
				next }
		}
		

		if (
			$segVals{probes} =~ /^\d+?$/
			&&
			$segVals{probes} < $pgx->{parameters}->{segment_probecount_min}
		) { next }

		$i++;
		
		my $info = { svlen => 1 * ($segVals{end} - $segVals{start}) };
		
		if ($segVals{probes} =~ /^\d+?$/) {
			$info->{probes} = 1* $segVals{probes} }
		if ($segVals{value} =~ /^\-?\d+?(\.\d+?)?$/) {
			$info->{value} = $numfactor * $segVals{value} }
		
		push(
			@{ $pgx->{$segmentsT} },
			{
				no => $i,
				callset_id => $segVals{callset_id},
				reference_name => $segVals{reference_name},
				variant_type =>	$varStatus,
				start => 1 * $segVals{start},
				end => 1 * $segVals{end},
				info => $info,
				digest => join(':',
					$segVals{reference_name},
					join(',', $segVals{start}.'-'.$segVals{end} ),
					$varStatus
				),
			}
		);

	}

	return $pgx;

}

################################################################################

sub _f2l {

	my $file = shift;
	open FILE, "$file" or die "No file $file $!";
	local $/;
	my $fContent = <FILE>;
	close FILE;
	return [ split(/\r\n?|\n/, $fContent) ];

}

################################################################################

sub _l2t {

	my $list = shift;
	my $header = [ ];
	my $table = [ ];
	foreach (@$list) {
		if (/^\#/) {
			push( @$header, $_) }
		else {
			push( @$table, [ split("\t", $_) ] ) }
	}
	return ($header, $table);

}

################################################################################

sub _norm_sample_id {

	my $id = shift;
	$id =~ s/[^\w\-\:]/_/g;
	return $id;

}

################################################################################

sub _objectify_header {

	no warnings 'uninitialized';

	my $header = shift;

	my $oh = { };
	
	foreach my $line (@$header) {
		my %lo = ( );
		$line =~ s/^\#?\s+|\s+$//g;
		$line =~ s/=>/;/g;
		foreach (grep{/\=/} split(';', $line)) {
			my ($par, $val) = split('=', $_);
			$par =~ s/^\s+|\s+$//g;
			$val =~ s/^[\s\"\']+|[\s\"\']+$//g;
			if ($par eq 'biosample_id') {
				$val = _norm_sample_id($val) }
			$lo{$par} = $val;
		}
		# `biosample_id` containing lines are stored under this sample
		if (grep{ /^biosample_id$/} keys %lo) {
			foreach (keys %lo) {
				$oh->{samples}->{ $lo{biosample_id} }->{ $_ } = $lo{ $_ };
			}
		} elsif (grep{ /^group_id$/} keys %lo) {
			foreach (keys %lo) {
				$oh->{groups}->{ $lo{group_id} }->{ $_ } = $lo{ $_ };
			}
		} elsif	($line =~ /^plotpar/) {
			foreach (keys %lo) {
				$oh->{plotpars}->{ $_ } = $lo{ $_ };
			}	
		} else {
			foreach (keys %lo) {
				$oh->{meta}->{ $_ } = $lo{ $_ };
			}
		}
	}

	return $oh;

}

################################################################################

sub read_file_to_split_array {

	my $file = shift;

	if ($file =~ /\.(ods)|(xlsx?)$/i) {
		use	Spreadsheet::Read;
		use Spreadsheet::ReadSXC;
		my $book = ReadData($file);
		my $header = [ ];
		my $table = [ ];
		foreach my $currentRow (Spreadsheet::Read::rows($book->[1])) {
			if ($currentRow->[0] =~ /^\#/) {
				push( @$header, $currentRow) }
		else {
			push( @$table, $currentRow ) }
		}
		return ($header, $table);
	}
	
	return _l2t( _f2l( $file ) );

}

################################################################################

sub read_webfile_to_split_array {

	use LWP::UserAgent;
	use LWP::Simple;
	use	Spreadsheet::Read;
	use Spreadsheet::ReadSXC;

	my $web = shift;

	if ($web =~ /dropbox\.com/) {
	$web =~	s/(\?dl=\w)?$/?dl=1/ }

	$ENV{'PERL_LWP_SSL_VERIFY_HOSTNAME'} = 0;

	my $ua = new LWP::UserAgent;
	$ua->agent("Mozilla/8.0");

	my $req = new HTTP::Request 'GET' => $web;
	$req->header('Accept' => 'text/plain');

	my $res = $ua->request($req);
	my @content;

	if ($res =~ /\.(ods)|(xlsx?)$/i) {
		my $book = ReadData($res->{_content});
		foreach my $currentRow (Spreadsheet::Read::rows($book->[1])) {
			push( @content, join("\t", @{ $currentRow }) );
		}
	} else {
		@content = split("\n", $res->{_content});
		chomp @content;
	}

	if ($args{DELCOMMENT} =~ /^T/i) {
		@content = grep{ ! /^\#/ } @content;
		@content = grep{ /./ } @content;
	}
	
	return _l2t( \@content );

}



1;
