#!/usr/bin/perl -w

$| = 1;

# CPAN packages
use Data::Dumper;
use File::Basename;
use List::MoreUtils qw(any apply uniq first_index);
use MongoDB;
use MongoDB::MongoClient;
use POSIX 'strftime';
use strict;
use Term::ProgressBar;
use YAML::XS qw(LoadFile);

# use warnings;
no warnings 'recursion';
no warnings 'uninitialized';

=podmd

##### Sample numbers

Sample numbers are represented in 2 different ways:

* the overall number shown in the title is calculated as the sum of the
`code_matches` values of the subsets
* the numbers shown with the labels are from the `count` parameter, which for
hierarchical subsets may include the numbers of all child terms (which may show
up separately, too); so the sum of those may over-counting samples

Examples:
  - `perl subsetsPlotter.pl -size_title_left_px 800 -size_plotimage_w_px 2400 -mingroupno 50 -datasets arraymap,progenetix`
  - `perl subsetsPlotter.pl -size_title_left_px 800 -size_plotimage_w_px 2400 -mingroupno 50 -datasets arraymap,tcga,progenetix -inalldatasets y`
  - `perl PGX/Scripts/subsetsPlotter.pl -mingroupno 20 -datasets progenetix -outdir /Users/mbaudis/switchdrive/baudisgroup/artwork -size_plotarea_h_px -1 -query NCIT:C -excluded NCIT:C5226 -size_strip_h_px 15 -size_title_left_px 640 -color_plotarea_hex '#333333' -size_text_title_left_px 12`

=cut

# local packages

my $here_path;
BEGIN {
  $here_path = File::Basename::dirname( eval { ( caller() )[1] } );
  push @INC, $here_path.'/../..';
}

use PGX::PGX;
use PGX::GenomeIntervals::IntervalStatistics;
use PGX::GenomePlots::HistoPlotter;
use PGX::GenomePlots::StripPlotter;
use PGX::GenomePlots::Genomeplot;
use PGX::GenomePlots::PlotParameters;
use PGX::Helpers::UtilityLibs;

################################################################################
# data import & pre-parsing ####################################################
################################################################################

my $config     =   LoadFile($here_path.'/../config/config.yaml') or die print 'Content-type: text'."\n\n¡No config.yaml file in this path!";
bless $config;

# command line input
my %args        =   @ARGV;

if ($args{'-datasets'} =~ /^[\w\,]+?/) {
  $config->{ dataset_names }    =   [ split(',', $args{'-datasets'}) ] }

$args{'-datasets'}      ||= 'arraymap';
$args{'-split_sets'}    ||= 'y';
$args{'-genome'}        ||= 'grch38';
$args{'-randno'}        ||= -1;
$args{'-mingroupno'}    ||= 5;
$args{'-query'}         ||= 'icdom-8,icdom-9';
$config->{grouping}     =   $args{'-grouping'}  ||= 'icdom';
$args{'-plotregions'}   ||= q{};
$args{'-outdir'}        ||= $here_path.'/tmp';
$args{'-inalldatasets'} ||= '';
$args{'-excluded'}      ||= '';

# alternative plot defaults to the ones from PGX::rsrc::config::plotdefaults.yaml
$args{'-size_plotimage_w_px'}   ||=   1948;
$args{'-size_plotmargin_px'}    ||=   50;
$args{'-size_title_left_px'}    ||=   720;
$args{'-size_clustertree_w_px'} ||=   200;
$args{'-size_plotarea_h_px'}    ||=   40;
$args{'-size_text_px'}          ||=   12;
$args{'-size_text_title_px'}    ||=   48;
$args{'-size_text_subtitle_px'} ||=   32;
$args{'-size_text_title_left_px'}   ||=  16;
$args{'-label_y_m'}     ||=   '-50,0,50';
$args{'-chr2plot'}      ||=   '1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22';

if (@{ $config->{ dataset_names } } < 2) { $args{'-split_sets'} = 'n' }

my @plottypes   =   split(',', $args{'-plottypes'});
my @excluded    =   split(',', $args{'-excluded'});

mkdir $args{'-outdir'};
my $queryLabel  =   '';
my $subsetQ     =   {};
my $qList       =   [];

if (grep{ $args{'-query'} =~ /^($_[\:\-][\w\.]+?){1,99}/ } map{ $config->{biosubsets}->{$_}->{prefix} } keys %{ $config->{biosubsets} }) {
  $qList        =   [ @$qList, map{ { "id" => qr/$_/i } } split(',', $args{'-query'}) ];
  $queryLabel   =   ',select_'.$args{'-query'};
  $queryLabel   =~  s/[^\,\w]/_/g;
}

my %qBioPres    =   map{ $_ => 1 } ( grep{ $args{'-query'} =~ /$_/ } map{ $config->{biosubsets}->{$_}->{prefix} } keys %{ $config->{biosubsets} } );

my $sortValuesK =   "id";
my $sortValuesL =   "label";
my $sortValuesT =   $config->{biosubsets}->{ $config->{grouping} }->{ label }; #TODO: doesn't work for e.g. geogse

if (scalar @$qList == 1) {
  $subsetQ = $qList->[0] }
elsif (scalar @$qList > 1) {
  $subsetQ = { '$or' => $qList } }
else {
  print <<END;
No query in one of the formats was specified:
  -datasets arraymap,tcga -query icdom-817

END

  exit;
}

my @fileTags    = ();
for (qw(datasets grouping)) { push(@fileTags, $args{'-'.$_}) }
for (qw(query mingroupno excluded randno)) {
  if ($args{'-'.$_} =~ /^\w/) { push(@fileTags, $_.'__'.$args{'-'.$_}) }
}
if ($args{'-randno'} > 0) { push(@fileTags, 'randno__'.$args{ '-randno' }) }

my $fileTag     =   join(',', grep{ /\w/ } @fileTags);
$fileTag        =~  s/[^\w\-,\.]/_/g;

################################################################################

print <<END;

################################################################################
You are attempting to plot data from the $args{'-datasets'} dataset, with the
following parameters:

END

foreach (sort keys %args) {
  my $name      = '    '.$_.(" " x 40);
  $name         =~  s/^(.{40}).*?$/$1/;
  print $name.$args{$_}."\n";
}

print <<END;
################################################################################
END

my $subsets			=		[];
my %allCodes    =   ();
my $progBar;

foreach my $dataset (@{ $config->{ dataset_names } }) {

  my $dbconn    =   MongoDB::MongoClient->new()->get_database($dataset);
  my $cursor;

  $cursor       =   $dbconn->get_collection( 'biosubsets' )->find( $subsetQ );

  my $dssub     =   [ $cursor->all ];
  for my $i (0..$#{ $dssub }) {
    if ($dssub->[$i]->{count} >= $args{'-mingroupno'}) {
      if (any { $dssub->[$i]->{id} =~ /^$_/ } @excluded ) { next }
      $allCodes{ $dssub->[$i]->{id} } +=  1;
      $dssub->[$i]->{dataset} =   $dataset;
      $dssub->[$i]->{frequencymaps}->{name} =   $dssub->[$i]->{id};
      $dssub->[$i]->{frequencymaps}->{name} =~  s/^pgx\://;
      if (@{ $config->{ dataset_names } } > 1) {
        $dssub->[$i]->{frequencymaps}->{name} .=  ', '.$dataset }
      $dssub->[$i]->{frequencymaps}->{name} .=   ': '.$dssub->[$i]->{label}.' ('.$dssub->[$i]->{count}.')';
      push(@$subsets, $dssub->[$i]);
    }
  }
}

if ($args{'-inalldatasets'} =~ /y/i ) {
  my @subI      =   ();
  for my $i (0..$#{ $subsets }) {
    if (scalar(@{ $config->{ dataset_names } }) == $allCodes{$subsets->[$i]->{id}}) {
      push(@subI, $i) }
  }
  $subsets      =   [ (@$subsets)[@subI ] ];
}


if ($args{'-randno'} > 0) {
  $subsets      =   RandArr( $subsets, $args{'-randno'} ) }

print @$subsets." subsets will be processed\n";

if (@$subsets < 1) {

  print <<END;

No subsets could be identified.

END

  exit;

}

################################################################################
# histogram ####################################################################
################################################################################

my $usedSampleNo        =   0;
foreach (@$subsets) {
  if ( $_->{code_matches} ) {
    $usedSampleNo   +=    $_->{code_matches} }
  else {
    $usedSampleNo +=    $_->{count} }
}

print @$subsets." different subsets will be plotted\n";

$args{-title}   =   'CNA Frequencies in '.@$subsets.' Cancer Entities by '.$sortValuesT;
$args{-subtitle}=   'Analysis Based on '.$usedSampleNo.' Individual Genome Profiles from the '.join(', ', @{ $config->{ dataset_names } }).' Collection';
if (@{ $config->{ dataset_names } } > 1) { $args{-subtitle} .=  's' }

my $histoMultF;
my $fMatrixFile;
my $plotargs;
my $plot;

$histoMultF     =   $fileTag.',histoplot_mult.svg';
$fMatrixFile    =   $args{'-outdir'}.'/'.$fileTag.',frequencymatrix.tsv';

################################################################################
# (clustered) CNA histograms
################################################################################

$plotargs       =   \%args;
$plotargs->{'-plottype'}  =   'multihistogram';
$plotargs->{'-path_loc'}  =   $args{'-outdir'};

$plot           =   new PGX($plotargs);
$plot->{frequencymaps}  =   [ map{ $_->{frequencymaps} } @$subsets ];
$plot->cluster_frequencymaps();
$plot->return_histoplot_svg();
$plot->write_frequency_matrix($fMatrixFile);
$plot->write_svg($histoMultF);

print "The SVG was written to $histoMultF.\n";



1;
