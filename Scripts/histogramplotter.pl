#!/usr/bin/perl

# CPAN packages
use Data::Dumper;
use File::Basename;
use List::Util qw(shuffle);
use MongoDB;
use MongoDB::MongoClient;
use POSIX 'strftime';
use strict;
use Term::ProgressBar;

# local packages

my $path_of_this_module;
BEGIN {
  $path_of_this_module = File::Basename::dirname( eval { ( caller() )[1] } );
  push @INC, $path_of_this_module.'/../..';
}

use PGX::GenomeIntervals::CytobandReader;
use PGX::GenomeIntervals::GenomeIntervals;
use PGX::GenomeIntervals::IntervalStatistics;
use PGX::GenomePlots::CytobandsPlotter;
use PGX::GenomePlots::Genomeplot;
use PGX::GenomePlots::HistoPlotter;
use PGX::GenomePlots::PlotParameters;

# command line input
my %args        =   @ARGV;

$args{'-dataset'}       ||= 'arraymap_ga4gh';
$args{'-genome'}        ||= 'grch36';
$args{'-randno'}        ||= -1;
$args{'-binning'}       ||= 1000000;
$args{'-varquery'}      ||= {};
$args{'-csquery'}       ||= { biosample_id => { '$regex' => qr/GSM23/ } };
$args{'-biosquery'}     ||= { id => { '$regex' => qr/GSM23/ } };
$args{'-plottype'}      ||= 'histogram';
$args{'-plotregions'}   ||= q{};
$args{'-svg'}           ||= './histoplot.svg';       

_checkArgs();

# conventions
my $csColl      =   'callsets_cnv_'.genome_names_to_grch($args{'-genome'});
my $csQuery     =   $args{'-csquery'};

# predefined, for recycling
my $distincts;
my $cursor;
my $csIds       =   [];
my $callsets    =   [];

# preconfigured objects

########    ####    ####    ####    ####    ####    ####    ####    ####    ####
########    main    ####    ####    ####    ####    ####    ####    ####    ####
########    ####    ####    ####    ####    ####    ####    ####    ####    ####

my $dbconn  =   MongoDB::MongoClient->new()->get_database($args{'-dataset'});

########    ####    ####    ####    ####    ####    ####    ####    ####    ####

$cursor     =   $dbconn->get_collection($csColl)->find( $csQuery )->fields( { biosample_id => 1, info => 1 } );
$callsets   =   [ $cursor->all ];

if ($args{'-randno'} > 0) {
  $callsets =   [ (shuffle(@$callsets))[0..($args{'-randno'}-1)] ] }

my $csNo    =   scalar @$callsets;

my $timeLab =   strftime("%T", gmtime());

print <<END;

$csNo callsets will be processed for $csColl.

Parameters:

  -dataset  $args{'-dataset'}
  -genome   $args{'-genome'}

Start:      $timeLab

END

# biosamples query

$cursor     =   $dbconn->get_collection('biosamples')->find( $args{'-biosquery'} )->fields( { id => 1, 'bio_characteristics.ontology_terms' => 1, info => 1 } );
my $biosamples =   [ $cursor->all ];
if ($args{'-randno'} > 0) {
  $biosamples =   [ (shuffle(@$biosamples))[0..($args{'-randno'}-1)] ] }

for my $i (0..$#{ $biosamples }) {
  my @icdm      =  grep{ $_->{term_id} =~ /icdom/ } @{ $biosamples->[$i]->{bio_characteristics}->[0]->{ontology_terms} };
  $biosamples->[$i]->{term_id}     =   $icdm[0]->{term_id};
  $biosamples->[$i]->{term_label}  =   $icdm[0]->{term_label};
  delete $biosamples->[$i]->{bio_characteristics};
  print Dumper($biosamples->[$i]->{id}, $biosamples->[$i]->{term_id}, $biosamples->[$i]->{term_label})."\n";
}

#my $samleMap    =   { map{ $_}}

#exit;

$args{'-text_bottom_left'}      =   $csNo.' samples';

my $pgx        =   new PGX(\%args);
$pgx->pgx_add_frequencymaps( [ {
  statusmapsets =>  $callsets } ] );
$pgx->return_histoplot_svg();

open  (FILE, ">", $args{'-svg'}) || warn 'output file '.$args{'-svg'}.' could not be created.';
binmode(FILE, ":utf8");
print FILE  $pgx->{svg};
close FILE;

my $timeLab     =   strftime("%T", gmtime());

print <<END;

End:        $timeLab
SVG:        $args{'-svg'}

END

################################################################################
# subs #########################################################################
################################################################################

sub _checkArgs {

  # terminating if arraypath doesn't exist
  if (! $args{'-dataset'}) {

    print <<END;

Script parameters:

-dataset        Dataset; e.g. "arraymap_ga4gh".
                Required

-genome         Genome edition, for setting the correct coordinate space.
                Default: GRCh38
                Both "hg" and "GRCh" styles can be used.

-chr2plot       Chromosomes to be plotted (comma separated)
                Default: 1 => Y
                Example: "7,8,18"

-plotregions    For generating a single plot, with one or more plot regions
                Syntax: 9:20000000-22000000,17:0-25000000
                The filename will be adjusted accordingly (here
                  "arrayplot,chr9_20000000-22000000,17_0-25000000.svg").
                This overrides the "-chr2plot" parameter.

-value_plot_y_max
                Maximum Y value of the plot
                Default: 100
                Can be adjusted depending on probe dynamics.

-size_chromosome_w_px
                pixel width of the chromosomal bands
                Default: 12
                Cytoband plotting can be ski[pped by setting this to "0".

END

    exit;

  }


  # adjusting for special cases

  # / adjusting for special cases

}

################################################################################

sub _printFeedback {

  my @showArgs  =   qw(
title
text_bottom_left
genome
value_plot_y_max
colorschema
size_plotimage_w_px
size_chromosome_w_px
);
  print <<END;

################################################################################

You are processing the dataset

    $args{-dataset}

Some of the current parameters are shown below. Those may be changed by
providing them as command line arguments.

END

  foreach (@showArgs) {
    my $name    = '    -'.$_.(" " x 40);
    $name       =~  s/^(.{40}).*?$/$1/;
    print $name.$pgx->{parameters}->{$_}."\n";
  }
  print <<END;

################################################################################

END

}

1;