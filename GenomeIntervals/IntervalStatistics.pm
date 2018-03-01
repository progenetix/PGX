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

=pod

sub "callset_make_statusmap"

The subroutine returns an object containing statusvalues (DUP, DEL) and the (min, max)
values of the overlapping variants, foreach of the provided @{ $plotPars->{GENOINTVS} } genome intervals.

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

  my $plot      =   shift;
  my $maps      =   {
    intervals   =>  scalar(@{ $plot->{genomeintervals} }),
    binning     =>  $plot->{genomeintervals}->[0]->{end} - $plot->{genomeintervals}->[0]->{start},
  };

  my %intStatLabs   =   (
    DUP         =>  'dupmap',
    DEL         =>  'delmap',
  );
  my %intValLabs   =   (
    DUP         =>  'dupmax',
    DEL         =>  'delmin',
  );

  foreach (values %intStatLabs) {
    $maps->{$_} =   [ map{''} 0..$#{ $plot->{genomeintervals} } ] }
  foreach (values %intValLabs) {
    $maps->{$_} =   [ map{ 0 } 0..$#{ $plot->{genomeintervals} } ] }

  my $valueMap  =   [ map{[0]} 0..$#{ $plot->{genomeintervals} } ];

	foreach my $csVar (@{ $plot->{segmentdata} }) {
	
	  if (! grep{ $csVar->{variant_type} eq $_ } keys %intStatLabs) { next }

    # the index of  intervals with a match to the current variant is created and used
    # to assign the status value and collect the segment value (since several variants
    # may overlap the same interval)
    foreach my $ind (grep{
      $csVar->{reference_name} eq	$plot->{genomeintervals}->[ $_ ]->{reference_name}
      &&
      $csVar->{start} <=	$plot->{genomeintervals}->[ $_ ]->{end}
      &&
      $csVar->{end}	>=	$plot->{genomeintervals}->[ $_ ]->{start}
    } 0..$#{ $plot->{genomeintervals} }) {

      $maps->{ $intStatLabs{ $csVar->{variant_type } } }->[$ind] = $csVar->{variant_type};
      push(
        @{ $valueMap->[$ind] },
        $csVar->{info}->{value},
      );

  }}

  # the values for each interval are sorted, to allow extracting the min/max 
  # values by position
  $valueMap     =   [ map{ [ sort { $a <=> $b } @{ $valueMap->[$_] } ] } 0..$#{ $plot->{genomeintervals} } ];

  # the last of the sorted values is assigned iF > 0
  foreach my $ind (grep{ $valueMap->[$_]->[-1] > 0 } 0..$#{ $plot->{genomeintervals} }) {
    $maps->{dupmax}->[$ind] =   1 * (sprintf "%.4f", $valueMap->[$ind]->[-1]) }

  # the first of the sorted values is assigned iF < 0
  foreach my $ind (grep{ $valueMap->[$_]->[0] < 0 } 0..$#{ $plot->{genomeintervals} }) {
    $maps->{delmin}->[$ind] =   1 * (sprintf "%.4f", $valueMap->[$ind]->[0]) }

  $plot->{statusmaps}   =   $maps;

  return $plot;

}


################################################################################

sub interval_cnv_frequencies {

  no warnings 'uninitialized';
  
=pod

Expects:

Returns:

=cut

########    ####    ####    ####    ####    ####    ####    ####    ####    ####

  my $plot      =   shift;
  my $cnvmaps   =   shift;
  my $name      =   shift;
  my $labels    =   shift;
  my $maps      =   {
    intervals   =>  scalar(@{ $plot->{genomeintervals} }),
    binning     =>  $plot->{genomeintervals}->[0]->{end} - $plot->{genomeintervals}->[0]->{start} + 1,
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
    for my $i (0..$#{ $plot->{genomeintervals} }) {
      $maps->{ $freqLabs{ $type } }->[$i]   =   sprintf "%.3f", ($fFactor * ( grep{ $_->{ $intLabs{ $type } }->[$i] eq $type } @{ $cnvmaps } ));
    }
  }

  push(
    @{ $plot->{frequencymaps} },
    $maps,
  );
    
  return $plot;

}

################################################################################

sub cluster_frequencymaps {

  use Algorithm::Cluster;
  no warnings 'uninitialized';

  my $plot      =   shift;

  my @matrix    =   ();
  my $labels    =   [];
  my $order     =   [];

  foreach my $frequencymapsSet (@{ $plot->{frequencymaps} }) {
    push(@{ $labels }, $frequencymapsSet->{name});
    push(
        @matrix,
        [
          (map{ $frequencymapsSet->{dupfrequencies}->[$_] + 0 }  @{ $plot->{matrixindex} }),
          (map{ $frequencymapsSet->{delfrequencies}->[$_] + 0 }  @{ $plot->{matrixindex} })
        ]
      );
  }
  
  if (scalar(@{ $labels }) < 3) { return $plot }

  my $EisenTree =   Algorithm::Cluster::treecluster(
                      transpose =>  0,
                      method    =>  'm',
                      dist      =>  'e',
                      data      =>  \@matrix,
                    );
  (
    $plot->{clustertree},
    $order
  )             =   cluster_tree($EisenTree);

  $plot->{frequencymaps}  =   [ map{ $plot->{frequencymaps}->[$_] } reverse(@{ $order }) ];

  return $plot;

}

################################################################################

sub cluster_samples {

  use Algorithm::Cluster;
  no warnings 'uninitialized';

  my $plot      =   shift;

  my @matrix    =   ();
  my $labels    =   [];
  my $order     =   [];

  my $i         =   0;
  foreach my $sample (@{ $plot->{samples} }) {
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
        (map{ $sample->{statusmaps}->{dupmap}->[$_] eq 'DUP' ? 1 : 0 } @{ $plot->{matrixindex} }),
        (map{ $sample->{statusmaps}->{delmap}->[$_] eq 'DEL' ? 1 : 0 } @{ $plot->{matrixindex} })
      ]
    );
  }
  if (scalar(@{ $labels }) < 3) { return $plot }

  my $EisenTree =   Algorithm::Cluster::treecluster(
                      transpose =>  0,
                      method    =>  'm',
                      dist      =>  'e',
                      data      =>  \@matrix,
                      mask      =>  [],
                    );
  (
    $plot->{clustertree},
    $order
  )             =   cluster_tree($EisenTree);

  $plot->{samples}  =   [ map{ $plot->{samples}->[$_] } reverse(@{ $order }) ];

  return $plot;

}

1;
