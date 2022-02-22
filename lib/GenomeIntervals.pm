package lib::GenomeIntervals;

use Data::Dumper;

require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(
  make_genome_intervals
  get_reference_base_limits
  get_genome_basecount
);

################################################################################

sub make_genome_intervals {

=pod

Expects:
  - a list reference of genome interval objects, usually representing cytobands:
    [
      {
        no =>  __integer__,            # not used
        reference_name  =>  __string__,
        start =>  __integer__,
        end   =>  __integer__,
        stain =>  __string__,            # not used
        label =>  __string__,            # not used
      },
      {
      ...
      },
    ]
  - the binning size in bases (optional; defaults to 1000000)

Returns:
  - a list of genome intervals in a similar structure, based on the binning size
  - the intervals are per reference, i.e. each reference starts with a new interval,
    leading to the last interval < binning size
    [
      {
        no    =>  __integer__,          # 1 -> n
        reference_name  =>  __string__,
        start =>  __integer__,
        end   =>  __integer__,
        label =>  __string__,
      },
      {
      ...
      },
    ]

=cut

########    ####    ####    ####    ####    ####    ####    ####    ####    ####

	my $cytobands = shift;
	my $intSize = shift;
	
	my $terminal_intervals_soft_expansion = 100000;

	my $refLims = get_reference_base_limits($cytobands);
	my $gi = [];

	# references are sorted with numerical ones first, then others (e.g. 1 -> 22, X, Y)
	my @refNames = ((sort {$a <=> $b } grep{ /^\d\d?$/ } keys %$refLims), (sort grep{ ! /\d/ } keys %$refLims));

	my $i = 1;

	for my $chro (@refNames) {

		my $rbl = $refLims->{$chro};		
		my $p_max = $rbl->{"p"}->[1];
		my $q_max = $rbl->{"size"};
    	my $arm = "p";		
		my $start = 0;
		
    # calculate first interval to end p-arm with a full sized one
    my $p_first = $p_max;
    while ($p_first >= ($intSize + $terminal_intervals_soft_expansion)) {
        $p_first -= $intSize }
		
		my $end = $start + $p_first;

		while ($start < $q_max) {
		
			my $int_p = $intSize;

			# adjusting the end of the last interval
			if ($end > $q_max) {
				$end = $q_max }
			elsif ($q_max < ($end + $terminal_intervals_soft_expansion) ){
				$end = $q_max;
				$int_p += $terminal_intervals_soft_expansion;
			}

			if ($end >= $p_max) {
				$arm = "q" }

			my $size = $end - $start;
			
			push(
				@$gi,
				{
					no =>  $i,
					reference_name => $chro,
					arm => $arm,
					start => $start,
					end => $end,
					size => $size,
					label => $chro.$arm.':'.$start.'-'.$end,
				}
			);

			$start = $end;
			$end += $int_p;
			$i += 1;

		}
	}

	return $gi;

}

################################################################################

sub get_reference_base_limits {

=pod

Expects:
  - a list reference of genome interval objects:
    [
      {
        reference_name  =>  __string__,
        start           =>  __integer__,
        end             =>  __integer__,
        (others)...
      },
      {
        ...
      },
    ]

Returns:
  - a hash reference with keys for each chromosome and various [ start, end ]
  	arrays of integers representing the reference's (or arm's) min and max bases:
  	
    {
    	__chro__ : {
    		"size": __integer__,
    		"chro": [  __integer__,  __integer__ ],
    		"p": [  __integer__,  __integer__ ],
    		"q": [  __integer__,  __integer__ ],
    	}

      { __string__  =>  [  __integer__,  __integer__ ] }
     },
     { ... },
=cut

########    ####    ####    ####    ####    ####    ####    ####    ####    ####

	my $allRefs = shift;

	my $refLims = {};

	foreach my $ref (map{ $_->{"reference_name"} } @$allRefs) {

		my @refRefs = grep{ $_->{"reference_name"} =~ /^$ref$/i } @$allRefs;
		my @pRefs = grep{ $_->{"band"} =~ /^p/i } @refRefs;
		my @qRefs = grep{ $_->{"band"} =~ /^q/i } @refRefs;
		my @bases = sort { $a <=> $b } ((map{ $_->{start} } @refRefs), (map{ $_->{end} } @refRefs));
		$refLims->{$ref}->{"chro"} = [ $bases[0], $bases[-1] ];
		$refLims->{$ref}->{"size"} = $bases[-1] - $bases[0];
		@bases = sort { $a <=> $b } ((map{ $_->{start} } @pRefs), (map{ $_->{end} } @pRefs));
		$refLims->{$ref}->{"p"} = [ $bases[0], $bases[-1] ];
		@bases = sort { $a <=> $b } ((map{ $_->{start} } @qRefs), (map{ $_->{end} } @qRefs));
		$refLims->{$ref}->{"q"} = [ $bases[0], $bases[-1] ];

	}

	return $refLims;

}

################################################################################

sub get_genome_basecount {

	my $allRefs = shift;
	my $chr2plot = shift;
	
	my $genBases;

	my %refNames = map { $_ => 1 } @$chr2plot;

	foreach my $ref (keys %refNames) {
		my @refRefs = grep{ $_->{reference_name} =~ /^$ref$/i } @$allRefs;
		my @bases = sort { $a <=> $b } ((map{ $_->{start} } @refRefs), (map{ $_->{end} } @refRefs));
		$genBases += ($bases[-1] - $bases[0]);
	}

	return $genBases;

}

################################################################################
########    utility subs    ####    ####    ####    ####    ####    ####    ####
################################################################################

1;
