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
	
	my $last_interval_soft_expansion = 100000;

	my $refLims = get_reference_base_limits($cytobands);
	my $gi = [];

	# references are sorted with numerical ones first, then others (e.g. 1 -> 22, X, Y)
	my @refNames = ((sort {$a <=> $b } grep{ /^\d\d?$/ } keys %$refLims), (sort grep{ ! /\d/ } keys %$refLims));

	my $intI = 1;
	for my $i (0..$#refNames) {

		my $refName = $refNames[$i];
		my $rbl = $refLims->{$refName};
		
		my $p_max = $rbl->{"p"}->[1];
        my $arm = "p";
		
		my $start = $rbl->{"chro"}->[0];
		my $end = $start + $intSize;

		while ($start < $rbl->{"chro"}->[1]) {
				
			my $int_p = $intSize;

			# adjusting the end of the last interval
			if ($end > $rbl->{"size"}) {
				$end = $rbl->{"size"};
			} elsif ($rbl->{"size"} < ($end + $last_interval_soft_expansion) ){
				$end = $rbl->{"size"};
				$int_p += $last_interval_soft_expansion;
			} elsif ($p_max > 0) {
				if ($end >= $p_max) {
					$end = $p_max;
					$int_p = $end - $start;
					$p_max = 0;
					$arm = "q";
				}
			}

			my $thisSize = $end - $start;
			
			push(
				@$gi,
				{
					no =>  $intI,
					reference_name => $refName,
					arm => $arm,
					start => $start,
					end => $end,
					"length" => $thisSize,
					label => $refName.$arm.':'.$start.'-'.$end,
				}
			);

			$start += $thisSize;
			$end += $intSize;
			$intI++;

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
  - a hash reference with each value consisting of a list of 2 integers, representing
    the reference's min and max bases:
    [
      { __string__  =>  [  __integer__,  __integer__ ] },
      { ... },
    ]

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
