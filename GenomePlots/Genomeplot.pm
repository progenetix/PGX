package PGX::GenomePlots::Genomeplot;

use Data::Dumper;
use PGX::GenomeIntervals::CytobandReader;
use PGX::GenomeIntervals::GenomeIntervals;
use PGX::GenomePlots::PlotParameters;
use PGX::FileUtilities::ArrayfileReader;

require Exporter;
@ISA    =   qw(Exporter);
@EXPORT =   qw(
  new
  plot_add_probedata
  plot_add_segmentdata
  plot_add_probedata_fracb
  plot_add_segmentdata_fracb
  plot_adjust_random_probevalues
  plot_get_plotregions
);

################################################################################

sub new {

  my $class     =   shift;
  my $args      =   shift;
  my $self      =   {
    parameters  =>  args_modify_plot_parameters(read_plot_defaults(), $args),
    cytobands   =>  read_cytobands($args->{'-genome'}),
    plotid      =>  ($args->{'-plotid'} !~ /^\w+?/ ? 'genomeplot' : $args->{'-plotid'}),
    svg         =>  q{},
    Y           =>  0,
  };
  bless $self, $class;
  $self->{genomeintervals}      =   make_genome_intervals(
                                      $self->{cytobands},
                                      $args->{'-binning'},
                                    );
  $self->{referencebounds}      =   get_reference_base_limits($self->{cytobands});
  $self->{genomesize}           =   get_genome_basecount(
                                      $self->{cytobands},
                                      $self->{parameters}->{chr2plot},
                                    );

  return $self;

}

################################################################################

sub plot_add_probedata {

  my $probefile =   shift;
  my $plot      =   shift;
  $plot->{probedata}    =   read_probefile($probefile, $plot);
  return $plot;

}

########    ####    ####    ####    ####    ####    ####    ####    ####    ####

sub plot_add_segmentdata {

  my $segfile   =   shift;
  my $plot      =   shift;
  $plot->{segmentdata}  =   read_segmentfile($segfile, $plot);
  return $plot;

}

########    ####    ####    ####    ####    ####    ####    ####    ####    ####
########    ####    ####    ####    ####    ####    ####    ####    ####    ####
########    ####    ####    ####    ####    ####    ####    ####    ####    ####

sub plot_add_probedata_fracb {

  my $probefile =   shift;
  my $plot      =   shift;
  $plot->{probedata_fracb}    =   read_probefile($probefile, $plot);
  return $plot;

}

########    ####    ####    ####    ####    ####    ####    ####    ####    ####
########    ####    ####    ####    ####    ####    ####    ####    ####    ####
########    ####    ####    ####    ####    ####    ####    ####    ####    ####

sub plot_add_segmentdata_fracb {

  my $segfile   =   shift;
  my $plot      =   shift;
  $plot->{segmentdata_fracb}  =   read_segmentfile($segfile, $plot);

  return $plot;

}

########    ####    ####    ####    ####    ####    ####    ####    ####    ####
########    ####    ####    ####    ####    ####    ####    ####    ####    ####
########    ####    ####    ####    ####    ####    ####    ####    ####    ####

sub plot_adjust_random_probevalues {

  use Term::ProgressBar;

  my $plot      =   shift;

	if ($plot->{parameters}->{simulated_probes} =~ /y/i ) {

    my $i				=		0;
    my $progBar =   Term::ProgressBar->new({name => 'Adjusting Simulated Values', count => scalar @{$plot->{segmentdata}}});

		foreach my $seg (@{ $plot->{segmentdata} }) {

			my @probeI				=		map{ $_ } grep{
															$plot->{probedata}->[$_]->{reference_name} eq $seg->{reference_name}
															&&
															$plot->{probedata}->[$_]->{position} >=  $seg->{start}
															&&
															$plot->{probedata}->[$_]->{position} <=  $seg->{end}
														} (0..$#{ $plot->{probedata} });

		  $progBar->update($i++);

			foreach (@probeI) {
				$plot->{probedata}->[$_]->{value}				+=	$seg->{info}->{value};
	  }}
	  $progBar->update(scalar @{$plot->{segmentdata}});
	}


  return $plot;

}

########    ####    ####    ####    ####    ####    ####    ####    ####    ####
########    ####    ####    ####    ####    ####    ####    ####    ####    ####
########    ####    ####    ####    ####    ####    ####    ####    ####    ####

sub plot_get_plotregions {

  my $regions   =   shift;
  my $plot      =   shift;

  my %chros     =   map{ $_->{reference_name} => 1 } @$regions;
  my @refNames  =   ((sort {$a <=> $b } grep{ /^\d\d?$/ } keys %chros), (sort grep{ ! /\d/ } keys %chros));

  if (! grep{ /^\d\w?$/ } @refNames) { return $plot }

  my $refLims   =   {};
  my $baseCount =   0;
  foreach my $ref (@refNames) {
    my @allBounds       =   map{ $_->{start}, $_->{end} } (grep{ $_->{reference_name} eq $ref } @$regions);
    @allBounds          =   sort {$a <=> $b } @allBounds;
    $refLims->{$ref}    =   [ $allBounds[0], $allBounds[-1] ];
    $baseCount          +=  ($allBounds[-1] - $allBounds[0]);
  }

  $plot->{parameters}->{do_chromosomes_proportional} = /n/;
  $plot->{parameters}->{chr2plot}   =   [@refNames];
  $plot->{referencebounds}      =   $refLims;
  $plot->{genomesize}           =   $baseCount;

  return $plot;

}

1;
