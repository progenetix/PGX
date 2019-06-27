package PGX::GenomePlots::Genomeplot;

use Data::Dumper;
use PGX::GenomeIntervals::CytobandReader;
use PGX::GenomeIntervals::GenomeIntervals;
use PGX::GenomeIntervals::IntervalStatistics;
use PGX::GenomePlots::PlotParameters;
use PGX::GenomePlots::HistoPlotter;
use PGX::GenomePlots::ArrayPlotter;
use PGX::GenomePlots::StripPlotter;
use PGX::GenomePlots::CytobandsPlotter;
use PGX::IOUtilities::PGXfileReader;
use PGX::IOUtilities::PGXfileWriter;
use PGX::IOUtilities::PGXdataAggregation;

require Exporter;
@ISA    =   qw(Exporter);
@EXPORT =   qw(
  pgx_add_frequencymaps
  pgx_add_probes_from_file
  pgx_add_segments_from_file
  pgx_add_segmentsets_from_samples
  pgx_add_fracbprobes_from_file
  pgx_add_fracbsegments_from_file
  plot_adjust_random_probevalues
  pgx_get_genome_regions
);

########    ####    ####    ####    ####    ####    ####    ####    ####    ####
########    ####    ####    ####    ####    ####    ####    ####    ####    ####
########    ####    ####    ####    ####    ####    ####    ####    ####    ####

sub pgx_get_genome_regions {

  my $pgx       =   shift;
  my $regions   =   $pgx->{parameters}->{plotregions};
  my %chros     =   map{ $_->{reference_name} => 1 } @$regions;
  my @refNames  =   ((sort {$a <=> $b } grep{ /^\d\d?$/ } keys %chros), (sort grep{ ! /\d/ } keys %chros));

  if (! grep{ /^\d\w?$/ } @refNames) { return $pgx }
  if (! grep{ $_->{reference_name} =~ /^\d\w?$/ } @$regions) { return $pgx }

  my $refLims   =   {};
  my $baseCount =   0;

  foreach my $ref (@refNames) {
    my @allBounds       =   map{ $_->{start}, $_->{end} } (grep{ $_->{reference_name} eq $ref } @$regions);
    @allBounds          =   sort {$a <=> $b } @allBounds;
    $refLims->{$ref}    =   [ $allBounds[0], $allBounds[-1] ];
    $baseCount          +=  ($allBounds[-1] - $allBounds[0]);
  }

  $pgx->{parameters}->{do_chromosomes_proportional} = /n/;
  $pgx->{parameters}->{chr2plot}   =   [@refNames];
  $pgx->{referencebounds}  =   $refLims;
  $pgx->{genomesize}       =   $baseCount;

  $pgx->{matrixindex}  =   [];
  my @selIntI   =   ();  
  my $i         =   0;
  foreach my $int (@{$pgx->{genomeintervals}}) {
    if (
      $pgx->{referencebounds}->{ $int->{reference_name} }
      &&
      $int->{start} <= $pgx->{referencebounds}->{ $int->{reference_name} }->[1]
      &&
      $int->{end} >= $pgx->{referencebounds}->{ $int->{reference_name} }->[0]
    ) { 
      push(@{ $pgx->{matrixindex} }, $i);
    }
    $i++;

  }

  return $pgx;

}

########    ####    ####    ####    ####    ####    ####    ####    ####    ####
########    ####    ####    ####    ####    ####    ####    ####    ####    ####
########    ####    ####    ####    ####    ####    ####    ####    ####    ####

sub pgx_add_frequencymaps {

=pod
## pgx_add_frequencymaps

#### Expects:


#### Returns:


=cut

  my $pgx       =   shift;
  my $csColls   =   shift;

  $pgx->{frequencymaps}    =   [];

  foreach my $csColl (@$csColls) {
    $pgx->interval_cnv_frequencies(
      [ map{$_->{statusmaps}} @{ $csColl->{statusmapsets} } ],
      $csColl->{name},
      $csColl->{labels},
    );

  }

  return $pgx;

}

########    ####    ####    ####    ####    ####    ####    ####    ####    ####

sub pgx_add_probes_from_file {

=pod

This method adds probe values $pgx->{probedata}->[0] -> $pgx->{probedata}->[n]
to the $pgx object.
Probe data is handled as a single instance, assuming that the method won't be 
called for more than 1 array. However, it uses the standard "samples" construct
as input, which is an array of possibly multiple samples - the method should 
only be called when populated with one.


=cut

  my $pgx       =   shift;
  my $probefile =   shift;
  
  if (@{ $pgx->{samples} } != 1 ) {
  	push(
  		@{ $pgx->{errors} },
  		scalar(@{ $pgx->{samples} }).' samples in array probe reading attempt',
  	);
  }
  
  if ($probefile !~ /\w/) {
  	if ($pgx->{samples}->[0]->{paths}->{probefile} =~ /\w\.\w/) {
  		$probefile	=		$pgx->{samples}->[0]->{paths}->{probefile} }
  	else {
  		my $seriesId		=		(split('::', $pgx->{samples}->[0]->{id}))[1];
			my $arrayId			=		(split('::', $pgx->{samples}->[0]->{id}))[2];
			$probefile	=		$pgx->{parameters}->{genome}.'/'.$seriesId.'/'.$arrayId.'/probes,cn.tsv'
  	}
  	$probefile	=~ s/^\///;
  	$probefile	=			$pgx->{config}->{paths}->{dir_array_base_path}.'/'.$probefile;
  }
  
  $pgx->read_probefile($probefile);
  return $pgx;

}

########    ####    ####    ####    ####    ####    ####    ####    ####    ####

sub pgx_add_segments_from_file {

  my $pgx       =   shift;
  my $segfile   =   shift;

  $pgx->read_segmentfile($segfile);
  return $pgx;

}

########    ####    ####    ####    ####    ####    ####    ####    ####    ####

sub pgx_add_segmentsets_from_samples {

  my $pgx       =   shift;
  my $callsets  =   shift;
  my $idName    =   shift;
  if ($idName !~ /\w\w/) {
    $idName     =   'id'}
  
  $pgx->{segmentsets}  =   [];
       
  foreach my $cs (@$callsets) {
  
    if (! $cs->{name}) {
      $cs->{name}   =   $cs->{$idName} }
      
    push (
      @{ $pgx->{segmentsets} },
      {
        id          =>  $cs->{$idName},
        name        =>  $cs->{name},
        variants    =>  $cs->{variants},
        statusmaps  =>  $cs->{statusmaps},        
      }
    );

  }

  return $pgx;

}

########    ####    ####    ####    ####    ####    ####    ####    ####    ####

sub pgx_add_fracbprobes_from_file {

  my $pgx       =   shift;
  my $probefile =   shift;

  $pgx->read_probefile($probefile, 'probedata_fracb');
  return $pgx;

}

########    ####    ####    ####    ####    ####    ####    ####    ####    ####

sub pgx_add_fracbsegments_from_file {

  my $pgx       =   shift;
  my $segfile   =   shift;

  $pgx->read_segmentfile($segfile, 'segmentdata_fracb');
  return $pgx;

}

########    ####    ####    ####    ####    ####    ####    ####    ####    ####

sub plot_adjust_random_probevalues {

=pod

This method adjusts array probe values for the value of the segment they are
mapped to. The method is used for adjusting random probe values such as we are
using to simulate array data, in cases where only segments data is available.

The use of Term::ProgressBar here assumes that this function is only called in 
a local context (i.e. run in the terminal, not in web instances). 

=cut

  use Term::ProgressBar;

  my $pgx       =   shift;

  my $parent    = `ps -o ppid= -p $$ | xargs ps -o command= -p`;
  my $progBar;

  if ($pgx->{parameters}->{simulated_probes} =~ /y/i ) {

    my $i       =   0;
    
    if ($parent !~ /httpd/) {
      $progBar  =   Term::ProgressBar->new(
                      {
                        name  => 'Adjusting Simulated Values',
                        count => scalar @{ $pgx->{segmentdata} }
                      }
                    );
    }

    foreach my $seg (@{ $pgx->{segmentdata} }) {

      my @prI   =   map{ $_ } grep{
                      $pgx->{probedata}->[$_]->{reference_name} eq $seg->{reference_name}
                      &&
                      $pgx->{probedata}->[$_]->{position} >=  $seg->{start}
                      &&
                      $pgx->{probedata}->[$_]->{position} <=  $seg->{end}
                    } (0..$#{ $pgx->{probedata} });

      if ($parent !~ /httpd/) {
        $progBar->update($i++) }

      foreach (@prI) {
        $pgx->{probedata}->[$_]->{value}   +=  $seg->{info}->{value};
    }}

    if ($parent !~ /httpd/) {
      $progBar->update(scalar @{$pgx->{segmentdata}}) }

  }

  return $pgx;

}

1;
