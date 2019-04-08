package PGX::IOUtilities::PGXfileReader;

use Data::Dumper;
use Math::Random qw(random_normal);
use PGX::GenomePlots::PlotParameters;
use PGX::Helpers::UtilityLibs;

require Exporter;
@ISA    =   qw(Exporter);
@EXPORT =   qw(
  read_probefile
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
        no              =>  __integer__,          # 1 -> n
        probe_id        =>  __string__,
        reference_name  =>  __string__,
        position        =>  __integer__,
        value           =>  __long__,
      },
      {
      ...
      },
    ]

=cut

########    ####    ####    ####    ####    ####    ####    ####    ####    ####

  my $pgx       =   shift;
  my $probeF    =   shift;
  my $probeT    =   shift;
  $probeT       ||= 'probedata';

  $pgx->{$probeT}    =   [];
  my @randomV;

  if (! -f $probeF) { return $pgx->{$probeT} }

  my $numfactor =   1;
  if (
    $pgx->{parameters}->{'reverse'} =~ /y/i
    &&
    $probeT !~ /frac/i
  ) { $numfactor = -1 }

  open  FILE, "$probeF" or die "No file $probeF $!";
  local   $/;                             # no input separator
  my $fContent  =   <FILE>;
  close FILE;
  my @probeData =   split(/\r\n?|\n/, $fContent);
  shift @probeData;

  my $i         =   0;
  foreach (@probeData) {

    $i++;
    my (
      $probe_id,
      $reference_name,
      $position,
      $value,
    )           =   split (/\s/, $_, 5);
    $probe_id   =~  s/[^\w\-\,]/_/g;
    $reference_name     =~ s/[^\dxXyY]//;
    $reference_name     =~ s/^23$/X/;
    $reference_name     =~ s/^24$/Y/;
    $position   =   sprintf "%.0f", $position;  # due to some erroneous .5 in-between pos.
    $value      =   sprintf "%.4f", $value;

    if ($reference_name !~ /^\w\d?$/)             { next }
    if ($position       !~ /^\d{1,9}$/)           { next }
    if ($value          !~ /^\-?\d+?(\.\d+?)?$/)  { next }

    push(
      @{ $pgx->{$probeT} },
      {
        no              =>  $i,
        probe_id        =>  $probe_id,
        reference_name  =>  $reference_name,
        position        =>  $position,
        value           =>  $numfactor * $value,
      }
    );
  }
  
  if ($probeT =~ /fracb/i) { return $pgx }

  # random values
  if ($pgx->{parameters}->{simulated_probes} =~ /y/i ) {
    my @randomV =   random_normal(scalar @{ $pgx->{$probeT} }, 0, 0.25);
    foreach my $n (0..$#{ $pgx->{$probeT} }) {
      $pgx->{$probeT}->[$n]->{value}  =   $randomV[$n];
    }
  }

  # baseline adjustment
  if ($pgx->{parameters}->{plot_adjust_baseline} =~ /[123456789]/) {
    foreach my $n (0..$#{ $pgx->{$probeT} }) {
      $pgx->{$probeT}->[$n]->{value}  +=   $pgx->{parameters}->{plot_adjust_baseline};
    }
  }

  return $pgx;

}

################################################################################

sub read_segmentfile {

=pod

Expects:
  - a standard tab-delimited Progenetix segments  file

  sample  chro  start stop  mean  probes
  GSM481286 1 742429  7883881 -0.1594 699
  GSM481286 1 115673158 115705254 -0.3829 8
  GSM481286 1 115722621 119771659 0.167 424
  GSM481286 1 119776776 162617092 0.4168  1587
  GSM481286 1 162621657 165278686 0.6508  350
  GSM481286 1 165280711 167221337 0.4056  241
  GSM481286 1 167248788 168289603 0.6784  130
  ...

Returns:
  - a list reference of genome CNV objects:
    [
      {
        no              =>  __integer__,    # 1 -> n
        callset_id      =>  __string__,
        reference_name  =>  __string__,
        start           =>  __integer__,
        end             =>  __integer__,
        variant_type    =>  __string__,     # DUP, DEL
        info            =>  {
          value           =>  __long__,
          svlen           =>  __integer__,
          probes          =>  __integer__,
          assembly_id     =>  __string__,     # GRCh36 ...
          experiment_type =>  __string__,     # aCGH ...
        },
      },
      {
      ...
      },
    ]

=cut

########    ####    ####    ####    ####    ####    ####    ####    ####    ####

  my $pgx       =   shift;
  my $segmentsF =   shift;
  my $segmentsT =   shift;
  $segmentsT    ||= 'segmentdata';
  $pgx->{$segmentsT}  =   [];

  if (! -f $segmentsF) { return $pgx }

  my $numfactor =   1;
  if (
    $pgx->{parameters}->{'reverse'} =~ /y/i
    &&
    $segmentsT !~ /frac/i
  ) { $numfactor = -1 }

  my %colOrder  =   (
    callset_id          =>  0,
    reference_name      =>  1,
    start               =>  2,
    end                 =>  3,
    value               =>  4,
    probes              =>  5,
  );

  if ($pgx->{parameters}->{format_inputfiles} =~ /tcga/i) {
    $colOrder{value}    =   5;
    $colOrder{probes}   =   4;
  };

  my $table     =   read_file_to_split_array($segmentsF);

  my $i         =   0;
  foreach my $segment (@$table) {
    $i++;
    if ($segment->[2] !~ /^\d+?$/) { next }
    my %segVals =   ();
    foreach (keys %colOrder) {
      $segVals{$_}  =   $segment->[$colOrder{$_}];
    };

    $segVals{callset_id}        =~  s/[^\w\-\,]/_/g;
    $segVals{reference_name}    =~ s/[^\dxXyY]//;
    $segVals{reference_name}    =~ s/^23$/X/;
    $segVals{reference_name}    =~ s/^24$/Y/;
    $segVals{start}     =   sprintf "%.0f", $segVals{start};
    $segVals{end}       =   sprintf "%.0f", $segVals{end};
    $segVals{probes}    =~  s/[^\d]//g;
    $segVals{value}     =   sprintf "%.4f", $segVals{value};

    if ($segVals{reference_name}!~ /^\w\d?$/)             { next }
    if ($segVals{start}         !~ /^\d{1,9}$/)           { next }
    if ($segVals{end}           !~ /^\d{1,9}$/)           { next }
    if ($segVals{value}         !~ /^\-?\d+?(\.\d+?)?$/)  { next }
    if (
      $segVals{probes} =~ /^\d\d?$/
      &&
      $segVals{probes} < $pgx->{parameters}->{segment_probecount_min}
    )                                                     { next }

    push(
      @{ $pgx->{$segmentsT} },
      {
        no              =>  $i,
        callset_id      =>  $segVals{callset_id},
        reference_name  =>  $segVals{reference_name},
        start           =>  [ 1 * $segVals{start} ],
        end             =>  [ 1 * $segVals{end} ],
        info            =>  {
          value         =>  $numfactor * $segVals{value},
          svlen         =>  1 * ($segVals{end} - $segVals{start}),
          probes        =>  1 * $segVals{probes},
        },
      }
    );

  }

  if ($segmentsT =~ /fracb/i) { return $pgx }

  # baseline adjustment
  if ($pgx->{parameters}->{plot_adjust_baseline} =~ /[123456789]/) {
    foreach my $i (0..$#{ $pgx->{$segmentsT} }) {
      $pgx->{$segmentsT}->[$i]->{info}->{value} +=   $pgx->{parameters}->{plot_adjust_baseline};
    }
  }

  for my $i (0..$#{ $pgx->{$segmentsT} }) {

    if ($pgx->{$segmentsT}->[$i]->{info}->{value} >= $pgx->{parameters}->{cna_gain_threshold}) {
      $pgx->{$segmentsT}->[$i]->{variant_type}  =   'DUP' }
    elsif ($pgx->{$segmentsT}->[$i]->{info}->{value} <= $pgx->{parameters}->{cna_loss_threshold}) {
      $pgx->{$segmentsT}->[$i]->{variant_type}  =   'DEL' }

  }

  return $pgx;

}

################################################################################

sub read_file_to_split_array {

	my $file      =   shift;
  my $table     =   [];

	if ($file =~ /\.(ods)|(xlsx?)$/i) {

    use	Spreadsheet::Read;
    use	Spreadsheet::XLSX;
    use Spreadsheet::ReadSXC;

    my $book			=	  ReadData($file);
    foreach my $currentRow (Spreadsheet::Read::rows($book->[1])) {
      push(
        @$table,
        $currentRow,
      );
    }

	} else {

    my $fContent  =   q{};
    open	FILE, "$file" or die "No file $file $!";
    local 	$/;															# no input separator
    $fContent   =	  <FILE>;
    close FILE;
    foreach my $line (split(/\r\n?|\n/, $fContent)) {
      push(
        @$table,
        [ split("\t", $line) ],
      );
    }
	}
	
	return	$table;

}

################################################################################

sub read_webfile_to_split_array {

	use LWP::UserAgent;
	use LWP::Simple;
	use	Spreadsheet::Read;
	use	Spreadsheet::XLSX;
	use Spreadsheet::ReadSXC;

	my $web      	=   shift;
  my $table     =   [];

  if ($web =~ /dropbox\.com/) {
  	$web 		=~	s/(\?dl=\w)?$/?dl=1/ }

	$ENV{'PERL_LWP_SSL_VERIFY_HOSTNAME'} 	= 	0;

	my $ua				=		new LWP::UserAgent;
  $ua->agent("Mozilla/8.0");

  my $req				=		new HTTP::Request 'GET' => $web;
  $req->header('Accept' => 'text/plain');

  my $res				=		$ua->request($req);
	my @filecontent;

	if ($res =~ /\.(ods)|(xlsx?)$/i) {
		my $book		=	ReadData($res->{_content});
		foreach my $currentRow (Spreadsheet::Read::rows($book->[1])) {
			push(
				@filecontent,
				join("\t", @{ $currentRow }),
			);
		}
	} else {
		@filecontent			  =		split("\n", $res->{_content});
		chomp	@filecontent;
	}

	if ($args{DELCOMMENT} =~ /^T/i) {
		@filecontent 			  = 	grep{ ! /^\#/ } @filecontent;
		@filecontent 			  = 	grep{ /./ } @filecontent;
	}
	
	foreach my $line (@filecontent) {
		push(
			@$table,
			[ split("\t", $line) ],
		);
	}
	

	return	$table;

}



1;
