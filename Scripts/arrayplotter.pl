#!/usr/bin/perl

=pod

This is a utility script for plotting positional genome data:

- probes (log value and b-allele fraction)
- segments (segmented log value and b-allele fractions)
- labels

=cut

# CPAN packages
use Data::Dumper;
use File::Basename;
use POSIX 'strftime';
use strict;
use Term::ProgressBar;

# local packages

my $path_of_this_module;
BEGIN {
  $path_of_this_module = File::Basename::dirname( eval { ( caller() )[1] } );
  push @INC, $path_of_this_module.'/../..';
}

use PGX::GenomeIntervals::GenomeIntervals;
use PGX::GenomePlots::ArrayPlotter;
use PGX::GenomeIntervals::CytobandReader;
use PGX::GenomePlots::Genomeplot;

# command line input

our %args           =   @ARGV;
$args{'-plottype'}  =   'array';  # fixed
$args{'-genome'}    ||= 'GRCh38';

$args{'-do_allchros'}   ||= 'y';
$args{'-plotregions'}   ||= q{};

# (possibly) derived
if ($args{'-chr2plot'} =~ /\w/) { $args{'-do_allchros'} = 'n' }

$args{'-chr2plot'}      ||= join(',', 1..22, 'X');

# file names and paths
$args{'-probefilename'}     ||= $args{'-pfn'}   ||= 'probes,cn.tsv';
$args{'-segfilename'}       ||= $args{'-sfn'}   ||= 'segments,cn.tsv';
$args{'-fracbprobefilename'}||= $args{'-fbpfn'} ||= 'probes,fracb.tsv';
$args{'-fracbsegfilename'}  ||= $args{'-fbsfn'} ||= 'segments,fracb.tsv';
$args{'-defaultsfilename'}  ||= $args{'-dfn'}   ||= 'plotdefaults.yaml';

$args{'-arraypath'} ||= $args{'-in'}  ||= q{};
$args{'-arraypath'} =~  s/\/$//;

foreach (grep{ /filename$/} keys %args) {
  my $file      =   $args{'-arraypath'}.'/'.$args{$_};
  my $pathK     =   $_;
  $pathK        =~  s/name//;
  $args{$pathK} =   $file;
}

# the output files will be named later, if not given here as a single SVG
$args{'-out'}           ||= $args{'-arraypath'};
$args{'-svgfilename'}   ||= $args{'-svgfn'} ||= q{};

# check or feedback
_checkArgs();

# conventions


# predefined, for recycling


# preconfigured objects
our $plot       =   new PGX::GenomePlots::Genomeplot(\%args);

_printFeedback();

########    ####    ####    ####    ####    ####    ####    ####    ####    ####
########    main    ####    ####    ####    ####    ####    ####    ####    ####
########    ####    ####    ####    ####    ####    ####    ####    ####    ####

$plot->plot_add_fracbprobes_from_file($args{'-fracbprobefile'});
$plot->plot_add_fracbsegments_from_file($args{'-fracbsegfile'});
$plot->plot_add_probes_from_file($args{'-probefile'});
$plot->plot_add_segments_from_file($args{'-segfile'});
$plot->plot_adjust_random_probevalues();
$plot->return_arrayplot_svg();

my $plotfile;
if ($args{'-plotregions'} =~ /\w\:\d+?\-\d+?(?:\,|$)/) {
  $plotfile     =   'arrayplot,chr'.$args{'-plotregions'}.'.svg';
  $plotfile     =~  s/[\:]/_/g;
  $args{'-do_allchros'} =   'n';
}
elsif (scalar(@{$plot->{parameters}->{chr2plot}}) < 22) {
  $plotfile     =   'arrayplot,chr'.$args{'-chr2plot'}.'.svg' }
else {
  $plotfile     =   'arrayplot.svg' }

if ($args{'-svgfilename'} =~ /^[\w\,\-]+?\.svg/i) {
  $plotfile     =   $args{'-svgfilename'};
  $args{'-do_allchros'} =   'n';
}

open  (FILE, ">", $args{'-out'}.'/'.$plotfile);
binmode(FILE, ":utf8");
print FILE  $plot->{svg};
close FILE;

if ($args{'-do_allchros'} !~ /y/i) { exit }

my $i           =   0;
my $progBar     =   Term::ProgressBar->new(scalar @{$plot->{parameters}->{chr2plot}});
my $chr2plot    =   $plot->{parameters}->{chr2plot};

foreach my $chro (@$chr2plot) {

  # re-initializing some values for multiple plots ...
  $plot->{parameters}->{chr2plot} =   [$chro];
  $plot->{svg}                    =   q{};
  $plot->{genomesize}             =   get_genome_basecount(
                                        $plot->{cytobands},
                                        $plot->{parameters}->{chr2plot},
                                      );
  $plot         =   return_arrayplot_svg($plot);
  $plotfile     =   'arrayplot,chr'.$chro.'.svg';
  open  (FILE, ">", $args{'-out'}.'/'.$plotfile);
  binmode(FILE, ":utf8");
  print FILE  $plot->{svg};
  close FILE;

  $progBar->update($i++);
}

$progBar->update(scalar @$chr2plot);

################################################################################
# / main #######################################################################
################################################################################

################################################################################
# subs #########################################################################
################################################################################

sub _checkArgs {

  # terminating if arraypath doesn't exist
  if (
(! -d $args{'-arraypath'})
||
(! -d $args{'-out'})
) {

    print <<END;

Script parameters:

-arraypath      Path to the directory containging probe nad segment files.
(or -in)        Required

-out            Path to the output directory for the SVG files.
                Defaults to -arraypath value if not specified

-format_inputfiles
                Allows modification for different column order and file
                content. The only alternative option to the "progenetix"
                default is currently "TCGA", which will adjust the column order
                of probe number and segment value.

-genome         Genome edition, for setting the correct coordinate space.
                Default: GRCh38
                Both "hg" and "GRCh" styles can be used.

-cna_loss_threshold
                Theshold for calling loss segments.
                Default: -0.15
                Alternative values provided by local defaults file or through
                  command line parameters.

-cna_gain_threshold
                Theshold for calling gain segments.
                Default: 0.15
                ... as above ...
                
-plot_adjust_baseline
                Value is added to both probes and segment values, before 
                applying thresholds to the segments
                Default: 0

-plot_adjust_baseline
                Value is added to both probes and segment values, before
                applying thresholds to the segments
                Default: 0

-chr2plot       Chromosomes to be plotted (comma separated)
                Default: 1 => Y
                Example: "7,8,18"

-plotregions    For generating a single plot, with one or more plot regions
                Syntax: 9:20000000-22000000,17:0-25000000
                The filename will be adjusted accordingly (here
                  "arrayplot,chr9_20000000-22000000,17_0-25000000.svg").
                This overrides the "-chr2plot" parameter.

-markers        For adding colored overlay block on one or more specified
                regions, with an optional text label. Colors can be added, or
                will be ranomised.
                Example:
                  11:2000000-3000000:marker:#ffcc00,11:2900000-3400000:another

-value_plot_y_max
                Maximum Y value of the plot
                Default: 3.2
                Can be adjusted depending on probe dynamics.

-factor_probedots
                Relative size of the probe dots
                Default: 3
                Can be lowered or increased, depending on probe number

-size_chromosome_w_px
                pixel width of the chromosomal bands
                Default: 12
                Cytoband plotting can be ski[pped by setting this to "0".

-do_chromosomes_proportional
                When giving a single chromosome in "-chr2plot", the size of the
                  plot is adjusted relative to chromosome 1.
                Default: y
                Setting this to "n" will keep the full plot size.


END

    exit;

  }

  # title
  if ($args{'-title'} !~ /^[\w \,\-]+?$/) {
    if ($args{'-arraypath'} =~ /(?:^|\/)([\w\,\-]+?)$/) {
      $args{'-title'}   =   $1 }
  }

  # adjusting for special cases
  opendir DIR, $args{'-arraypath'};
  my @arrayDirF   =   grep{ /^\w[\w\.\,\-]+?\.\w\w\w\w?\w?\w?$/ } -f, readdir(DIR);
  close DIR;

  # TCGA
  if ($args{'-format_inputfiles'} =~ /tcga/i) {
    if (! -f $args{'-segfile'} ) {
      my @segGuess      =   grep{ /seg\.txt/ } @arrayDirF;
      if (@segGuess == 1) { $args{'-segfile'} = $args{'-arraypath'}.'/'.$segGuess[0] }
    }
    $args{'-probefile'}         =   $path_of_this_module.'/../rsrc/probemaps/'.genome_names_to_hg($args{'-genome'}).'/GPL6801,probes,cn.tsv';
    $args{'-simulated_probes'}  =   'y';
    $args{'-text_bottom_left'}  =   'TCGA (simulated probes)';

  }

  # / TCGA

  # / adjusting for special cases

}

################################################################################

sub _printFeedback {

  my @showArgs  =   qw(
title
text_bottom_left
genome
format_inputfiles
simulated_probes
cna_gain_threshold
cna_loss_threshold
segment_probecount_min
value_plot_y_max
factor_probedots
colorschema
size_plotimage_w_px
size_chromosome_w_px
size_segments_stroke_px
);
  print <<END;

################################################################################

You are processing the array files from

    $args{-arraypath}

Some of the current parameters are shown below. Those may be changed by
providing them as command line arguments.

END

  foreach (@showArgs) {
    my $name    = '    -'.$_.(" " x 40);
    $name       =~  s/^(.{40}).*?$/$1/;
    print $name.$plot->{parameters}->{$_}."\n";
  }
  print <<END;

################################################################################

END

}

1;
