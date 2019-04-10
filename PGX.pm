package PGX;

use PGX::GenomeIntervals::CytobandReader;
use PGX::GenomeIntervals::GenomeIntervals;
use PGX::GenomeIntervals::IntervalStatistics;
use PGX::GenomePlots::Genomeplot;
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
@EXPORT =   qw(new);

########    ####    ####    ####    ####    ####    ####    ####    ####    ####
########    ####    ####    ####    ####    ####    ####    ####    ####    ####
########    ####    ####    ####    ####    ####    ####    ####    ####    ####

sub new {

  my $class     =   shift;
  my $args      =   shift;
  $args         =   args_modify_plot_parameters(read_plot_defaults(), $args);
  my $self      =   {
    parameters  =>  $args,
    cytobands   =>  read_cytobands($args->{genome}),
    plotid      =>  $args->{plotid},
    svg         =>  q{},
    Y           =>  0,
  };

  bless $self, $class;
  $self->{genomeintervals}  =   make_genome_intervals(
                                  $self->{cytobands},
                                  $self->{parameters}->{binning},
                                );
  $self->{referencebounds}  =   get_reference_base_limits($self->{cytobands});
  $self->{genomesize}       =   get_genome_basecount(
                                  $self->{cytobands},
                                  $self->{parameters}->{chr2plot},
                                );
  $self->{matrixindex}      =   [ 0..$#{ $self->{genomeintervals} } ];
  $self         =   pgx_get_genome_regions($self);
  return $self;

}

########    ####    ####    ####    ####    ####    ####    ####    ####    ####
########    ####    ####    ####    ####    ####    ####    ####    ####    ####
########    ####    ####    ####    ####    ####    ####    ####    ####    ####

1;
