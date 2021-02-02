#!/usr/bin/perl -w

$| = 1;

# CPAN packages
use Data::Dumper;
use File::Basename;
use POSIX 'strftime';
use strict;
use Term::ProgressBar;
use YAML::XS qw(LoadFile);

# use warnings;
no warnings 'recursion';
no warnings 'uninitialized';

=podmd

Input files can be either tab-separated plain text files (".tsv" or such), or
OpenOffice (".ods") spreadsheet files. In principle, Excel, too; but nobody 
should trust Excel.

Segment file name cannot contain ".", group_key has to be longer than 3 characters.

-f
  * a standard, single or multi-sample, tab-delimited Progenetix segments  file

  sample  chro  start stop  mean  probes
  GSM481286 1 742429  7883881 -0.1594 699
  GSM481286 1 115673158 115705254 -0.3829 8
  GSM481286 1 115722621 119771659 0.167 424
  GSM481286 1 119776776 162617092 0.4168  1587
  GSM481286 1 162621657 165278686 0.6508  350
  GSM481287 3 65280711 117221337 0.4056  20041
  ...

-sf
  * an optional tab-delimited file for labeling different sample groups
  
  sample  group_key group_label
  GSM481286 icdom-85003 breast cancer, NOS
  GSM481287 icdom-85003 breast cancer, NOS
  GSM533241 icdom-81403 gastric adenocarcinoma
  GSM533242 icdom-81403 gastric adenocarcinoma
  GSM533243 icdom-81403 gastric adenocarcinoma
  GSM126632 icdom-85003 breast cancer, NOS
  GSM126633 icdom-85003 breast cancer, NOS
  GSM126634 icdom-85003 breast cancer, NOS
  ...  
 
Examples:
  - perl pgx_plot_from_files.pl -f segments.tab -sf sortfile.ods

=cut

# local packages

my $path_of_this_module;
BEGIN {
  $path_of_this_module = File::Basename::dirname( eval { ( caller() )[1] } );
  push @INC, $path_of_this_module.'/..';
}

use PGX;

################################################################################
# data import & pre-parsing ####################################################
################################################################################

# predefined plot arguments
# those override the ones from insite PGX::rsrc::config
# but can be overridden again from command line input

my $plotargs = {
  -size_plotimage_w_px => 1024,
  -size_plotmargin_px => 25,
  -size_title_left_px => 200,
  -size_clustertree_w_px => 50,
  -size_plotarea_h_px => 40,
  -size_text_px => 12,
  -size_text_title_px => 48,
  -size_text_subtitle_px => 32,
  -size_text_title_left_px => 12,
  -label_y_m => '-50,0,50',
  -chr2plot => join(',', 1..22),
  -plotregions => q{},
  -genome => 'grch38',
  -plottype => 'histogram',
  -min_group_no => 3,
  -path_loc => $path_of_this_module.'/out',
};

# command line input
my %args = @ARGV;
$args{-sf} ||= q{};

if (-d $args{-outdir}) {
  $plotargs->{-path_loc} = $args{-outdir} }

foreach (keys %$plotargs) {
  if ($args{$_} =~ /\w/) {
    $plotargs->{$_} = $args{$_} }
}

mkdir $args{'-outdir'};

if (! -f $args{'-f'}) {
  print <<END;
No input file was specified:
  -f _path_to_my_file_/segments.tab

END
	exit;
}

# files and paths ##############################################################

my $outFileBase = $args{'-f'};
$outFileBase =~ s/^.*?\/([^\/]+?)\.\w{2,4}$/$1/;
my $sortFileBase = $args{'-sf'};
$sortFileBase =~ s/^.*?\/([^\/]+?)\.\w{2,4}$/,$1/;
$outFileBase = $outFileBase.$sortFileBase;
my $matrixplotF = $outFileBase.',samples_clustered.svg';
my $histoplotF = $outFileBase.',samples_histoplot.svg';
my $histoplotmultF = $outFileBase.',samples_histoplot_mult.svg';

# grouping of samples  #########################################################

my $pgx = new PGX($plotargs);
$pgx->{parameters}->{path_loc} = $args{'-outdir'};

# this uses the file reading routine; but multi-segment files have to be
# deconvoluted first ...
$pgx->pgx_add_segments_from_file($args{'-f'});
# print Dumper($pgx->{segmentdata});
# exit;
$pgx->pgx_create_samples_from_segments();
$pgx->pgx_callset_labels_from_file($args{'-sf'});
$pgx->pgx_create_sample_collections();

################################################################################

# new plot object, since 
# * $pgx->{samples}
# * $pgx->{samplecollections}
# are needed for the plots and would be deleted if recycling $pgx
my $plot;

$plot = new PGX($plotargs);
$plot->{parameters}->{plottype} = 'histogram';
$plot->{parameters}->{plotid} = 'histogram';
$plot->pgx_add_frequencymaps([ { statusmapsets => $pgx->{samples} } ]);
$plot->return_histoplot_svg();
$plot->write_svg($histoplotF);

print "Wrote histoplot to\n\t=> $histoplotF\n";

# sample matrix plot ##########################################################

$plotargs->{'-plottype'} = 'multistrip';
$plot = new PGX($plotargs);
$plot->{parameters}->{plotid} = 'multistripplot';
$plot->{samples} = $pgx->{samples};
$plot->cluster_samples();
$plot->return_stripplot_svg();
$plot->write_svg($matrixplotF);

print "Wrote sample plot to\n\t=> $matrixplotF\n";

# (clustered) CNA histograms

$plotargs->{'-plottype'} = 'multihistogram';
$plot = new PGX($plotargs);
$plot->pgx_add_frequencymaps($pgx->{samplecollections});
$plot->cluster_frequencymaps();
$plot->return_histoplot_svg();
$plot->write_svg($histoplotmultF);

print "Wrote grouped histograms to\n\t=> $histoplotmultF\n";

1;
