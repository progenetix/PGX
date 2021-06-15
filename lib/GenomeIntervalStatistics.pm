package lib::GenomeIntervalStatistics;

use lib::ClusterTree;

use Data::Dumper;

require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(
	segments_add_statusmaps
	interval_cnv_frequencies
	cluster_frequencymaps
	cluster_samples
	use_Eisen_tree
);

################################################################################

sub segments_add_statusmaps {
  
=podmd

### Sub "segments_add_statusmaps"

The subroutine returns an object containing maps for
* the (min, max) values of the overlapping variants
* the CNV coverage of the respective intervals

... foreach of the provided $pgx->{genomeintervals} genome intervals.

Structure:

maps:
  max:
    - 0 or pos. value x genomeintervals
  min:
    - 0 or neg. value x genomeintervals
  dup:
    - 0 or pos. value x genomeintervals
  del:
    - 0 or neg. value x genomeintervals

=cut

	no warnings 'uninitialized';

	my $pgx = shift;
	my $segments = shift;

	if (ref $segments ne 'ARRAY') {
		$segments = $pgx->{segmentdata} }
		

	my $maps = {
		intervals => scalar(@{ $pgx->{genomeintervals} }),
		binning => $pgx->{genomeintervals}->[0]->{end} - $pgx->{genomeintervals}->[0]->{start},
	};

	my %intCovLabs = ( DUP => 'dup', DEL => 'del' );
	my %intValLabs = ( DUP => 'max', DEL => 'min' );

	my $cnvStats = {
		cnvcoverage => 0,
		dupcoverage => 0,
		delcoverage => 0
	};

	foreach (values %intCovLabs) { $maps->{$_} = [ map{ 0 } 0..$#{ $pgx->{genomeintervals} } ] }
	foreach (values %intValLabs) { $maps->{$_} = [ map{ 0 } 0..$#{ $pgx->{genomeintervals} } ] }

	my $valueMap = [ map{[0]} 0..$#{ $pgx->{genomeintervals} } ];

	foreach my $csVar (@{ $segments }) {
	
		if (! grep{ $csVar->{variant_type} eq $_ } keys %intCovLabs) { next }

		# the index of intervals with a match to the current variant is created and used
		# to assign the status value and collect the segment value (since several variants
		# may overlap the same interval)
		foreach my $ind (grep{
			$csVar->{reference_name} eq	$pgx->{genomeintervals}->[ $_ ]->{reference_name}
			&&
			$csVar->{start} <= $pgx->{genomeintervals}->[ $_ ]->{end}
			&&
			$csVar->{end}	>= $pgx->{genomeintervals}->[ $_ ]->{start}
		} 0..$#{ $pgx->{genomeintervals} }) {

			my $ovEnd = (sort { $a <=> $b } ($pgx->{genomeintervals}->[ $ind ]->{end},  $csVar->{end} ) )[0];
			my $ovStart = (sort { $b <=> $a } ($pgx->{genomeintervals}->[ $ind ]->{start},  $csVar->{start} ) )[0];
			my $overlap = $ovEnd - $ovStart + 1;
			$maps->{ $intCovLabs{ $csVar->{variant_type } } }->[$ind] += $overlap;
			push( @{ $valueMap->[$ind] }, $csVar->{info}->{value} );
		  
		}
	}
  	
  	# statistics
	foreach my $cLab (values %intCovLabs) {
		foreach my $ind (grep{ $maps->{$cLab}->[$_] > 0 } 0..$#{ $pgx->{genomeintervals} }) {
			$cnvStats->{$cLab} += $maps->{$cLab}->[$ind];
			$cnvStats->{cnvcoverage} += $maps->{$cLab}->[$ind];
			$maps->{$cLab}->[$ind] = sprintf "%.5f", $maps->{$cLab}->[$ind] / ($pgx->{genomeintervals}->[$ind]->{end} - $pgx->{genomeintervals}->[$ind]->{start} );
		}
	}

	# the values for each interval are sorted, to allow extracting the min/max 
	# values by position
	$valueMap = [ map{ [ sort { $a <=> $b } @{ $valueMap->[$_] } ] } 0..$#{ $pgx->{genomeintervals} } ];

	# the last of the sorted values is assigned iF > 0
	foreach my $ind (grep{ $valueMap->[$_]->[-1] > 0 } 0..$#{ $pgx->{genomeintervals} }) {
		$maps->{max}->[$ind] = 1 * (sprintf "%.4f", $valueMap->[$ind]->[-1]) }

	# the first of the sorted values is assigned iF < 0
	foreach my $ind (grep{ $valueMap->[$_]->[0] < 0 } 0..$#{ $pgx->{genomeintervals} }) {
		$maps->{min}->[$ind] = 1 * (sprintf "%.4f", $valueMap->[$ind]->[0]) }

	$pgx->{statusmaps} = $maps;
	$pgx->{cnvstatistics} = $cnvStats;
  
	if ($pgx->{genomesize} > 1) {
		foreach my $covK (keys %$cnvStats) {
			my $fracK = $covK;
			$fracK =~ s/coverage/fraction/;
			$pgx->{cnvstatistics}->{$fracK}	= 1 * (sprintf "%.3f", $pgx->{cnvstatistics}->{$covK} / $pgx->{genomesize});
		} 
	}

	return $pgx;

}


################################################################################

sub interval_cnv_frequencies {

=pod

Expects:

Returns:

=cut

	no warnings 'uninitialized';

	my $pgx = shift;
	my $cnvmaps = shift;
	my $name = shift;
	my $labels = shift;
	
	my $sample_count = @{ $cnvmaps };

	my $maps = {
		interval_count => scalar(@{ $pgx->{genomeintervals} }),
		binning => $pgx->{genomeintervals}->[0]->{end} - $pgx->{genomeintervals}->[0]->{start},
		name => ($name =~ /\w/ ? $name : q{}),
		labels => (@$labels > 0 ? $labels : []),
		sample_count => $sample_count,
		intervals => [ ]
	};

	my %intLabs = ( DUP => 'dup', DEL => 'del' );
	my %freqLabs = ( DUP => 'gain_frequency', DEL => 'loss_frequency' );

	# avoiding division by 0 errors if improperly called
	my $fFactor = 100;
	if (@{ $cnvmaps } > 1) { $fFactor = 100 / @{ $cnvmaps } }
		
	$pgx->{parameters}->{bin_match_min} *= 1;
	for my $i (0..$#{ $pgx->{genomeintervals} }) {
		my $int = { "index" => $i };
		foreach my $type (keys %intLabs) {
			$int->{ $freqLabs{ $type } } = sprintf "%.3f", ( $fFactor * ( grep{ $_->{ $intLabs{ $type } }->[$i] >= $pgx->{parameters}->{bin_match_min} } @{ $cnvmaps } ));
		}
		push(@{$maps->{intervals}},  $int);
	}
	
	push(@{ $pgx->{frequencymaps} }, $maps);
	
	return $pgx;

}

################################################################################

sub use_Eisen_tree {

	use Algorithm::Cluster;
	no warnings;
	
	my $pgx = shift;
	my $matrix = shift;

	if ($pgx->{parameters}->{cluster_linkage_method} !~ /^[ascm]$/) {
		$pgx->{parameters}->{cluster_linkage_method} = 'm' }
	if ($pgx->{parameters}->{cluster_distance_metric} !~ /^[ecauxskb]$/) {
		$pgx->{parameters}->{cluster_distance_metric} = 'e' }

	my $EisenTree = Algorithm::Cluster::treecluster(
		transpose => 0,
		method => $pgx->{parameters}->{cluster_linkage_method},
		dist => $pgx->{parameters}->{cluster_distance_metric},
		data => $matrix,
	);

	return cluster_tree($EisenTree);

}

################################################################################

sub cluster_frequencymaps {

	my $pgx = shift;

	my @matrix = ();
	my $labels = [];
	my $order = [];

	foreach my $frequencymapsSet (@{ $pgx->{frequencymaps} }) {
		push(@{ $labels }, $frequencymapsSet->{name});
		push(
			@matrix,
			[
			  (map{ $frequencymapsSet->{intervals}->[$_]->{gain_frequency} + 0 } @{ $pgx->{matrixindex} }),
			  (map{ $frequencymapsSet->{intervals}->[$_]->{loss_frequency} + 0 } @{ $pgx->{matrixindex} })
			]
		);
	}
	
	
	if (scalar(@{ $labels }) < 3) { return $pgx }

	($pgx->{clustertree}, $order) = $pgx->use_Eisen_tree(\@matrix);
	
	$pgx->{frequencymaps} = [ map{ $pgx->{frequencymaps}->[$_] } reverse(@{ $order }) ];

	return $pgx;

}

################################################################################

sub cluster_samples {

	use Algorithm::Cluster;
	no warnings;

	my $pgx = shift;

	my @matrix = ();
	my $labels = [];
	my $order = [];

	my $i = 0;
	
	my $covThresh = $pgx->{parameters}->{bin_match_min} * 1;

	foreach my $sample (@{ $pgx->{samples} }) {
		$i++;
		my $label = 'sample_'.$i;
		if ($sample->{id} =~ /\w\w/) {
		  $label = $sample->{id} }
		elsif ($sample->{UID} =~ /\w\w/) {
		  $label = $sample->{UID} }
		$label =~ s/[^\w\,\-]/_/g;
		$label =~ s/_$//g;
		
		push(@{ $labels }, $label);
		push(
		  @matrix,
		  [
 			(map{ $sample->{statusmaps}->{dup}->[$_] >= $covThresh ? $sample->{statusmaps}->{dup}->[$_] : 0 } @{ $pgx->{matrixindex} }),
 			(map{ $sample->{statusmaps}->{del}->[$_] >= $covThresh ? $sample->{statusmaps}->{del}->[$_] : 0 } @{ $pgx->{matrixindex} })
		  ]
		);
	}

	if (scalar(@{ $labels }) < 3) { return $pgx }

	($pgx->{clustertree}, $order) = $pgx->use_Eisen_tree(\@matrix);

	$pgx->{samples} = [ map{ $pgx->{samples}->[$_] } reverse(@{ $order }) ];

	return $pgx;

}

1;
