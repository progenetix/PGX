package PGX;

use File::Basename;
use YAML::XS qw(LoadFile DumpFile);
use Data::Dumper;
use lib::CGItools;
use lib::CytobandReader;
use lib::GenomeIntervals;
use lib::GenomeIntervalStatistics;
use lib::PlotPrepareData;
use lib::PlotMakeParameters;
use lib::PlotHistoplots;
use lib::PlotArrays;
use lib::PlotErrorSVG;
use lib::PlotStripplots;
use lib::PlotCytobands;
use lib::Helpers;
use lib::ReadFiles;
use lib::WriteFiles;
use lib::AggregateData;

require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(
	new
	read_plot_defaults
);

########    ####    ####    ####    ####    ####    ####    ####    ####    ####

sub new {

	my $class = shift;
	my $args = shift;
	$args = args_modify_plot_parameters(read_plot_defaults(), $args);
	my $self = {
		parameters => $args,
		config => read_config(),
		cytobands => read_cytobands($args->{genome}),
		datasetid => q{},
		plotid => $args->{plotid},
		svg => q{},
		Y => 0,
		errors => [],
	};

	bless $self, $class;

	$self->{genomeintervals} = make_genome_intervals(
		$self->{cytobands},
		$self->{parameters}->{binning},
	);
	my $rbl = get_reference_base_limits($self->{cytobands});
	$self->{referencebounds} = _remap_referencebounds($rbl);
	$self->{genomesize} =  get_genome_basecount(
		$self->{cytobands},
		$self->{parameters}->{chr2plot},
	);
	$self->{matrixindex}  = [ 0..$#{ $self->{genomeintervals} } ];
	$self = pgx_get_genome_regions($self);

	return $self;

}

################################################################################

sub _remap_referencebounds {

	my $refLims = shift;
	my $rb = {};
	for my $rn (keys %$refLims) {
		$rb->{$rn} = $refLims->{$rn}->{"chro"};
	}
	
	return $rb;

}

################################################################################

sub read_config {

  my $path_of_this_module = File::Basename::dirname( eval { ( caller() )[1] } );
  my $config = LoadFile($path_of_this_module.'/config/config.yaml');
  return $config;

}

################################################################################

sub read_plot_defaults {

  my $path_of_this_module = File::Basename::dirname( eval { ( caller() )[1] } );
  my $plotPars = LoadFile($path_of_this_module.'/config/plotdefaults.yaml');
  return $plotPars;

}

################################################################################

1;
