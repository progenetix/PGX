package PGX::GenomeIntervals::IntervalStatistics;

use PGX::GenomeIntervals::ClusterTree;

use Data::Dumper;

require Exporter;
@ISA        =   qw(Exporter);
@EXPORT     =   qw(
  plot_segments_add_statusmap
  interval_cnv_frequencies
  cluster_frequencymaps
  cluster_samples
);

################################################################################

sub plot_segments_add_statusmap {

  no warnings 'uninitialized';
  
  # 
  my $pgx       =   shift;
  my $segments  =   shift;
  
  if (ref $segments ne 'ARRAY') {
    $segments   =   $pgx->{segmentdata} }
  
=pod

sub "plot_segments_add_statusmap"

The subroutine returns an object containing statusvalues (DUP, DEL) and the (min, max)
values of the overlapping variants, foreach of the provided $pgx->{genomeintervals} genome intervals.

Structure:

maps:
  dupmap:
    - null or DUP
    ...
  delmap:
    - null or DEL
    ...
  dupmax:
    - 0 or pos. value
    ...
  delmin:
    - 0 or neg. value
    ...

=cut

  my $maps      =   {
    intervals   =>  scalar(@{ $pgx->{genomeintervals} }),
    binning     =>  $pgx->{genomeintervals}->[0]->{end} - $pgx->{genomeintervals}->[0]->{start},
  };

  my %intStatLabs   =   (
    DUP         =>  'dupmap',
    DEL         =>  'delmap',
  );
  my %intCoverageLabs   =   (
    DUP         =>  'dupcoverage',
    DEL         =>  'delcoverage',
  );
  my %intValLabs   =   (
    DUP         =>  'dupmax',
    DEL         =>  'delmin',
  );

  foreach (values %intStatLabs) {
    $maps->{$_} =   [ map{''} 0..$#{ $pgx->{genomeintervals} } ] }
  foreach (values %intCoverageLabs) {
    $maps->{$_} =   [ map{ 0 } 0..$#{ $pgx->{genomeintervals} } ] }
  foreach (values %intValLabs) {
    $maps->{$_} =   [ map{ 0 } 0..$#{ $pgx->{genomeintervals} } ] }

  my $valueMap  =   [ map{[0]} 0..$#{ $pgx->{genomeintervals} } ];

	foreach my $csVar (@{ $segments }) {
	
	  if (! grep{ $csVar->{variant_type} eq $_ } keys %intStatLabs) { next }

    # the index of  intervals with a match to the current variant is created and used
    # to assign the status value and collect the segment value (since several variants
    # may overlap the same interval)
    foreach my $ind (grep{
      $csVar->{reference_name} eq	$pgx->{genomeintervals}->[ $_ ]->{reference_name}
      &&
      $csVar->{start}->[0] <=	$pgx->{genomeintervals}->[ $_ ]->{end}
      &&
      $csVar->{end}->[-1]	>=	$pgx->{genomeintervals}->[ $_ ]->{start}
    } 0..$#{ $pgx->{genomeintervals} }) {

      my $ovEnd     =   (sort { $a <=> $b } ($pgx->{genomeintervals}->[ $ind ]->{end},  $csVar->{end}->[-1]) )[0];
      my $ovStart   =   (sort { $b <=> $a } ($pgx->{genomeintervals}->[ $ind ]->{start},  $csVar->{start}->[0]) )[0];
      my $overlap   =   $ovEnd - $ovStart + 1;
    
      $maps->{ $intStatLabs{ $csVar->{variant_type } } }->[$ind] = $csVar->{variant_type};
      $maps->{ $intCoverageLabs{ $csVar->{variant_type } } }->[$ind]  += $overlap;
      push(
        @{ $valueMap->[$ind] },
        $csVar->{info}->{value},
      );

  }}
  
  foreach my $cLab (values %intCoverageLabs) {
    foreach my $ind (grep{ $maps->{$cLab}->[$_] > 0 } 0..$#{ $pgx->{genomeintervals} }) {
      $maps->{$cLab}->[$ind]  =  sprintf "%.3f", $maps->{$cLab}->[$ind] / ($pgx->{genomeintervals}->[$ind]->{end} - $pgx->{genomeintervals}->[$ind]->{start} + 1);
  }}

  # the values for each interval are sorted, to allow extracting the min/max 
  # values by position
  $valueMap     =   [ map{ [ sort { $a <=> $b } @{ $valueMap->[$_] } ] } 0..$#{ $pgx->{genomeintervals} } ];

  # the last of the sorted values is assigned iF > 0
  foreach my $ind (grep{ $valueMap->[$_]->[-1] > 0 } 0..$#{ $pgx->{genomeintervals} }) {
    $maps->{dupmax}->[$ind] =   1 * (sprintf "%.4f", $valueMap->[$ind]->[-1]) }

  # the first of the sorted values is assigned iF < 0
  foreach my $ind (grep{ $valueMap->[$_]->[0] < 0 } 0..$#{ $pgx->{genomeintervals} }) {
    $maps->{delmin}->[$ind] =   1 * (sprintf "%.4f", $valueMap->[$ind]->[0]) }

  $pgx->{statusmaps}    =   $maps;

  return $pgx;

}


################################################################################

sub interval_cnv_frequencies {

  no warnings 'uninitialized';
  
=pod

Expects:

Returns:

=cut

########    ####    ####    ####    ####    ####    ####    ####    ####    ####

  my $pgx       =   shift;
  my $cnvmaps   =   shift;
  my $name      =   shift;
  my $labels    =   shift;
  my $maps      =   {
    intervals   =>  scalar(@{ $pgx->{genomeintervals} }),
    binning     =>  $pgx->{genomeintervals}->[0]->{end} - $pgx->{genomeintervals}->[0]->{start} + 1,
    name        =>  ($name =~ /\w/ ? $name : q{}),
    labels      =>  (@$labels > 0 ? $labels : []),
  };
  my %intLabs   =   (
    DUP         =>  'dupmap',
    DEL         =>  'delmap',
  );
  my %freqLabs      =   (
    DUP         =>  'dupfrequencies',
    DEL         =>  'delfrequencies',
  );

  # avoiding division by 0 errors if improperly called
  my $fFactor   =   100;
  if (@{ $cnvmaps } > 1) { $fFactor = 100 / @{ $cnvmaps } }

  foreach my $type (keys %intLabs) {
    for my $i (0..$#{ $pgx->{genomeintervals} }) {
      $maps->{ $freqLabs{ $type } }->[$i]   =   sprintf "%.3f", ($fFactor * ( grep{ $_->{ $intLabs{ $type } }->[$i] eq $type } @{ $cnvmaps } ));
    }
  }

  push(
    @{ $pgx->{frequencymaps} },
    $maps,
  );
    
  return $pgx;

}

################################################################################

sub cluster_frequencymaps {

  use Algorithm::Cluster;
  no warnings 'uninitialized';

  my $pgx       =   shift;

  my @matrix    =   ();
  my $labels    =   [];
  my $order     =   [];
  
  if ($pgx->{parameters}->{cluster_linkage_method} !~ /^[ascm]$/) {
    $pgx->{parameters}->{cluster_linkage_method} = 'm' }
  if ($pgx->{parameters}->{cluster_distance_metric} !~ /^[ecauxskb]$/) {
    $pgx->{parameters}->{cluster_distance_metric} = 'e' }

  foreach my $frequencymapsSet (@{ $pgx->{frequencymaps} }) {
    push(@{ $labels }, $frequencymapsSet->{name});
    push(
        @matrix,
        [
          (map{ $frequencymapsSet->{dupfrequencies}->[$_] + 0 }  @{ $pgx->{matrixindex} }),
          (map{ $frequencymapsSet->{delfrequencies}->[$_] + 0 }  @{ $pgx->{matrixindex} })
        ]
      );
  }
  
  if (scalar(@{ $labels }) < 3) { return $pgx }

  my $EisenTree =   Algorithm::Cluster::treecluster(
                      transpose =>  0,
                      method    =>  $pgx->{parameters}->{cluster_linkage_method},
                      dist      =>  $pgx->{parameters}->{cluster_distance_metric},
                      data      =>  \@matrix,
                    );
  (
    $pgx->{clustertree},
    $order
  )             =   cluster_tree($EisenTree);

  $pgx->{frequencymaps}  =   [ map{ $pgx->{frequencymaps}->[$_] } reverse(@{ $order }) ];

  return $pgx;

}

################################################################################

sub cluster_samples {

  use Algorithm::Cluster;
  no warnings 'uninitialized';

  my $pgx       =   shift;

  my @matrix    =   ();
  my $labels    =   [];
  my $order     =   [];

  if ($pgx->{parameters}->{cluster_linkage_method} !~ /^[ascm]$/) {
    $pgx->{parameters}->{cluster_linkage_method} = 'm' }
  if ($pgx->{parameters}->{cluster_distance_metric} !~ /^[ecauxskb]$/) {
    $pgx->{parameters}->{cluster_distance_metric} = 'e' }

  my $i         =   0;
  foreach my $sample (@{ $pgx->{samples} }) {
    $i++;
    my $label   =   'sample_'.$i;
    if ($sample->{id} =~ /\w\w/) {
      $label    =   $sample->{id} }
    elsif ($sample->{UID} =~ /\w\w/) {
      $label    =   $sample->{UID} }
    $label      =~  s/[^\w\,\-]/_/g;
    $label      =~  s/_$//g;      

    push(@{ $labels }, $label);
    push(
      @matrix,
      [
        (map{ $sample->{statusmaps}->{dupmap}->[$_] eq 'DUP' ? 1 : 0 } @{ $pgx->{matrixindex} }),
        (map{ $sample->{statusmaps}->{delmap}->[$_] eq 'DEL' ? 1 : 0 } @{ $pgx->{matrixindex} })
      ]
    );
  }
  if (scalar(@{ $labels }) < 3) { return $pgx }

  my $EisenTree =   Algorithm::Cluster::treecluster(
                      transpose =>  0,
                      method    =>  $pgx->{parameters}->{cluster_linkage_method},
                      dist      =>  $pgx->{parameters}->{cluster_distance_metric},
                      data      =>  \@matrix,
                      mask      =>  [],
                    );
  (
    $pgx->{clustertree},
    $order
  )             =   cluster_tree($EisenTree);

  $pgx->{samples}  =   [ map{ $pgx->{samples}->[$_] } reverse(@{ $order }) ];

  return $pgx;

}

1;
