#!/usr/bin/env perl -w

$| = 1;

# CPAN packages
use Data::Dumper;
use File::Basename;
use MongoDB;
use MongoDB::MongoClient;
use POSIX 'strftime';
use strict;
use Term::ProgressBar;
use YAML::XS qw(LoadFile);

=podmd

Examples:
  - perl CNVdbPlotter.pl -datasets arraymap -query icdom-817
  - perl CNVdbPlotter.pl -query icdom-817 -datasets arraymap -sortfile ./in/testsort817.tsv -plottypes multistrip -size_strip_h_px 30 -size_title_left_w_px 0
  - perl CNVdbPlotter.pl -datasets arraymap -query icdom-....3 -mingroupno 100

=cut

# local packages

my $path_of_this_module;
BEGIN {
  $path_of_this_module = File::Basename::dirname( eval { ( caller() )[1] } );
  push @INC, $path_of_this_module.'/..';
}

use PGX;
use lib::PlotMakeParameters;
use lib::Helpers;

################################################################################
# data import & pre-parsing ####################################################
################################################################################

my $config = LoadFile($path_of_this_module.'/../config/config.yaml') or die print 'Content-type: text'."\n\n¡No config.yaml file in this path!";
bless $config;

$config->{ERROR} = [];

# command line input
my %args = @ARGV;

$args{'-datasets'} ||= "progenetix";
$args{'-split_sets'} ||= 'y';
$args{'-stem_code'} ||= 'n';
$args{'-genome'} ||= 'grch38';
$args{'-randno'} ||= -1;
$args{'-mingroupno'} ||= 5;
$args{'-query'} ||= 'icdom-8,icdom-9';
$args{'-sortfile'} ||= '';
$args{'-plottype'} ||= 'histogram';
$args{'-grouping'} ||= 'icdom';
$args{'-plotregions'} ||= q{};
$args{'-outdir'} ||= $path_of_this_module.'/out';
$args{'-path_loc'} ||= $args{'-outdir'};
$args{'-plottypes'} ||= 'multistrip,multihistogram';
$args{'-size_plotimage_w_px'} ||= 1948;
$args{'-size_plotmargin_px'} ||= 50;
$args{'-size_title_left_w_px'} ||= 500;
$args{'-size_clustertree_w_px'} ||= 200;
$args{'-size_plotarea_h_px'} ||= 40;
$args{'-size_text_px'} ||= 12;
$args{'-size_text_title_px'} ||= 48;
$args{'-size_text_subtitle_px'} ||= 32;
$args{'-size_text_title_left_px'} ||= 16;
$args{'-label_y_m'} ||= '-50,0,50';
$args{'-chr2plot'} ||= join(',', 1..22);

$config->{ dataset_names } = [ split(',', $args{'-datasets'}) ];
if (@{ $config->{ dataset_names } } < 2) {
	$args{'-split_sets'} = 'n' }

$config->{grouping} = $args{'-grouping'};

my @plottypes = split(',', $args{'-plottypes'});

mkdir $args{'-outdir'};
my $queryLabel = '';
my $biosampleQ = {};
my $qList = [];

foreach my $qElement (split(',', $args{'-query'})) {
	$qElement =~ /^(\w+?)([\:\-](.*?))?$/;
	my ($pre, $other) = ($1, $2);
  if (grep{ $qElement =~ /$config->{datacollections}->{ $_ }->{greedymatch}/ } keys %{ $config->{datacollections} }) {
    push(@$qList, { $config->{datacollections}->{ $pre }->{ samplefield }.".id" => qr/$qElement/i });
    $queryLabel = ',select_'.$qElement;
  }
}

my $sortValuesK;
my $sortValuesL;
my $sortValuesT;

# external sortfile
if (-f $args{'-sortfile'}) {
  $config->{grouping} = 'custom';
  $qList = [];
}

# this assumes that the first column contains an entry for the selected id (or id)
# the second then the associated label
my %customSort;

if ($config->{grouping} =~ /custom/i) {
  open  FILE, "$args{'-sortfile'}" or die "No file $args{'-sortfile'} $!";
  local   $/;                             # no input separator
  my $fContent = <FILE>;
  close FILE;
  my @filedata = split(/\r\n?|\n/, $fContent);
  shift @filedata;

  foreach (split(/\r\n?|\n/, $fContent)) {
    my @fLine = split("\t", $_);
    print $fLine[0];
    print $fLine[1];
    if ($fLine[0] =~ /[\w\-]+/ && $fLine[1] =~ /\w+/) {
      $customSort{$fLine[0]} = $fLine[1];
      push(
        @$qList,
        { "id" => qr/(?:^|\-|\:)$fLine[0](?:$|\-|\:)/i },
      );
    }
  }
}

print "
$config->{grouping}

";

if (scalar @$qList == 1) {
  $biosampleQ = $qList->[0] }
elsif (scalar @$qList > 1) {
  $biosampleQ = { '$or' => $qList } }

################################################################################
################################################################################
################################################################################

$config->check_query($biosampleQ);
#$config->check_grouping();
if (grep{ /.../ } @{ $config->{ERROR} }) {
  print "\n".join("\n\n", @{ $config->{ERROR} })."\n";
  exit;
}

################################################################################
################################################################################
################################################################################

if ($config->{grouping} =~ /custom/i) {

  $sortValuesK = 'custom';
  $sortValuesL = 'custom';
  $sortValuesT = 'Custom Groups';

} elsif (grep{ $config->{grouping} eq $_ } keys %{ $config->{datacollections} }) {

  $sortValuesK = $config->{datacollections}->{ $config->{grouping} }->{ samplefield }.".id";
  $sortValuesL = $config->{datacollections}->{ $config->{grouping} }->{ samplefield }.".label";
  $sortValuesT = 'Cancer Subsets by '.$config->{datacollections}->{ $config->{grouping} }->{ label };

}

my $fileTag = $args{'-datasets'}.','.$config->{grouping}.$queryLabel.',mingroupno_'.$args{'-mingroupno'}.',stem_code_'.$args{'-stem_code'};
$fileTag =~  s/[^\w\,\-\.]/_/g;
if ($args{'-randno'} > 0) { $fileTag .= ',randno'.$args{'-randno'} }

################################################################################

print <<END;

################################################################################
You are attempting to plot data from the $args{'-datasets'} dataset, with the
following parameters:

END

foreach (sort keys %args) {
  my $name = '    '.$_.(" " x 40);
  $name =~  s/^(.{40}).*?$/$1/;
  print $name.$args{$_}."\n";
}

print <<END;
################################################################################
END

my $samples = [];
my $progBar;

foreach my $dataset (@{ $config->{ dataset_names } }) {

  my $dbconn = MongoDB::MongoClient->new()->get_database($dataset);
  my $cursor = $dbconn->get_collection( 'biosamples' )->find( $biosampleQ )->fields( { id => 1, histological_diagnosis => 1, icdo_morphology => 1, icdo_topography => 1, external_references => 1, "cnv_statistics" => 1 } );
  my $bioS = [ $cursor->all ];

  if ($args{'-randno'} > 0) {
    $bioS = RandArr( $bioS, $args{'-randno'} ) }

  my $bioByID = {};

  foreach my $thissample (@$bioS) {
    my $id = $thissample->{id};
    if ($sortValuesK =~ /custom/i) {
      if (grep{ $_ eq $id } keys %customSort) {
        $bioByID->{$id} = {
          id => $customSort{ $id },
          label => $customSort{ $id },
        };
        $bioByID->{ $id }->{id} =~  s/\s/_/g;
      }
    }
    elsif ($sortValuesK =~ /histological_diagnosis/) {
    	$bioByID->{$id} = $thissample->{histological_diagnosis} }
    elsif ($sortValuesK =~ /external_references/) {
      foreach my $thisbioc (@{ $thissample->{external_references}} ) {
        if ($thisbioc->{id} =~  /$config->{grouping}/) {
          $bioByID->{$id} = {
            id => $thisbioc->{id},
            label => $thisbioc->{id},
          };
        }
      }
    }
  }

  print @$bioS." biosamples have been found\n";

  $cursor = $dbconn->get_collection( 'callsets' )->find( { biosample_id => { '$in' => [ keys %$bioByID ] } } )->fields(  { id => 1, biosample_id => 1, info => 1 } );
  my $callS = [ $cursor->all ];

  print @$callS." callsets have been found\n";

  my $vars = [];
  if (grep{ /multistrip/ } @plottypes) {
    $cursor = $dbconn->get_collection( 'variants' )->find( { biosample_id => { '$in' => [ keys %$bioByID ] } } )->fields( { _id => 0, updated => 0, created => 0, digest => 0, reference_bases => 0, mate_name => 0} );
    $vars = [ $cursor->all ];

    print @$vars." variants will be represented\n";
  }

  my $csc = 0;

  $progBar = Term::ProgressBar->new(scalar @$callS);
  foreach my $cs (@$callS) {

    $csc++;
    $progBar->update($csc);

    my $sampleID = $cs->{biosample_id};
    my $label = $bioByID->{$sampleID}->{label};
    my $code = $bioByID->{$sampleID}->{id};

    $code =~ s/pgx\://;
    if ($args{'-stem_code'} =~ /y/i) {
      $code =~ s/[\.\/]\d?\d?$//;
      $label = $code;
    }
    if ($args{'-split_sets'} =~ /y/i) {
      $label .= ', '.$dataset;
      $code .= $dataset;
    }

    my $sample = {
      id => $sampleID,
      $sortValuesK => $code,
      $sortValuesL => $label,
      cnv_statusmaps => $cs->{cnv_statusmaps},
      variants => [ grep{ $_->{biosample_id} eq $sampleID } @$vars ],
    };

    push(@$samples, $sample);

  }

  $progBar->update(scalar @$callS);

}

$samples = [ grep{ defined($_->{ $sortValuesK }) } @{ $samples } ];

print @$samples." samples will be processed\n";

if (@$samples < 1) {

  print <<END;

No samples could be identified.

END

  exit;

}

################################################################################
# histogram ####################################################################
################################################################################

################################################################################
# individual interval frequencies are generated for all the sortgroup values
# above the sample_number_min member count
################################################################################

my $callsetCollections = [];
my $usedSampleNo = 0;

# if custom sort values have been defined, each sample is tested for a match of each value,
# in the sortcolumn; if one of the values matches, it is assigned to the new field "CUSTOMSORT"
# if more than one => "ambiguous" is used; if no match => "NA"
my @sortVals;

my %sortValues;

my %valBuffer = map{ lc($_->{ $sortValuesK }) => 1 } @{ $samples };
foreach (keys %valBuffer) {
  my $val = $_;
  $val =~ s/^ +?//;
  $val =~ s/ +?$//;
  $sortValues{$val} = 1;
}

foreach (sort grep{ /\w\w/ } keys %sortValues) {

  my $sortGroup = lc($_);
  $sortGroup =~  s/\[/\\\[/g;
  $sortGroup =~  s/\]/\\\]/g;
  $sortGroup =~  s/\(/\\\]/g;
  $sortGroup =~  s/\)/\\\]/g;
  $sortGroup =~  s/\//\\\//g;

  my $theseSamples;
  my @theirIndex;

  @theirIndex = grep{ $samples->[$_]->{ $sortValuesK } =~ /^ *?$sortGroup *?$/i } 0..$#{ $samples };
  my $label = {
    label_text => $sortGroup,
    label_link => q{},
    label_color => random_hexcolor(),
  };
  foreach (@theirIndex) {
    $samples->[$_]->{labels} = [ $label ];
    $samples->[$_]->{name} = $samples->[$_]->{id}.' ('.$sortGroup.')';
  }
  $theseSamples = [@{$samples}[@theirIndex]];

  if (scalar(@{ $theseSamples }) < $args{'-mingroupno'}) { next }

  $usedSampleNo += scalar @{ $theseSamples };

  my $groupTag = $sortGroup.':'.$theseSamples->[0]->{$sortValuesK}.': ';
  if ($groupTag =~ /custom/i) {
    $groupTag = q{} }
  my $name = $theseSamples->[0]->{$sortValuesL}.' ('.scalar(@$theseSamples).')';

  push(
    @$callsetCollections,
    {
      labels => [ $label ],
      name => $name,
      statusmapsets => [ map{  { cnv_statusmaps => $_->{cnv_statusmaps} } } @{ $theseSamples } ],
    },
  );

}

print @$callsetCollections." different subsets will be plotted\n";

$args{-title} = 'CNA Frequencies in '.@$callsetCollections.' '.$sortValuesT;
$args{-subtitle}= 'Analysis Based on '.$usedSampleNo.' Individual Genome Profiles from the '.join(', ', @{ $config->{ dataset_names } }).' Collection';
if (@{ $config->{ dataset_names } } > 1) { $args{-subtitle} .= 's' }

my $histoMultF;
my $stripF;
my $fMatrixFile;
my $matrixFile;
my $vMatrixFile;
my $plotargs;
my $plot;

$histoMultF = $fileTag.',histoplot_mult.svg';
$stripF = $fileTag.',stripplot_mult.svg';

$fMatrixFile = $args{'-path_loc'}.'/'.$fileTag.',frequencymatrix.tsv';
$matrixFile = $args{'-path_loc'}.'/'.$fileTag.',statusmatrix.tsv';
$vMatrixFile = $args{'-path_loc'}.'/'.$fileTag.',valuematrix.tsv';

################################################################################
# sample matrix ################################################################
################################################################################

$plotargs = { %args };
$plotargs->{-plottype} = 'multistrip';

if (grep{ $_ eq $plotargs->{-plottype} } @plottypes) {

  $plot = new PGX($plotargs);
  $plot->{parameters}->{plotid} = 'multistripplot';
  $plot->{parameters}->{subtitle} = 'Display of CNV Profiles from '.@$samples.' Individual Samples from the '.join(', ', @{ $config->{ dataset_names } }).' Collection';
  if ($plot->{parameters}->{size_title_left_w_px} > ($plot->{parameters}->{size_strip_h_px} - 2)) {
    $plot->{parameters}->{size_text_title_left_px} = $plot->{parameters}->{size_strip_h_px} - 2 }
  $plot->{parameters}->{text_bottom_left} = @{ $samples }.' samples';
  $plot->{samples} = $samples;
  $plot->cluster_samples();
  $plot->return_stripplot_svg();
  $plot->write_svg($stripF);

  print "The SVG was written to "."$args{'-path_loc'}"."/"."$stripF.\n";

}

################################################################################
# (clustered) CNA histograms
################################################################################

$plotargs = { %args };
$plotargs->{'-plottype'} = 'multihistogram';

if (grep{ $_ eq $plotargs->{-plottype} } @plottypes) {

#  $plotargs->{'-size_strip_h_px'} = 0;
  $plot = new PGX($plotargs);
  $plot->{samples} = $samples;
  $plot->pgx_add_frequencymaps($callsetCollections);
  $plot->cluster_frequencymaps();
  $plot->return_histoplot_svg();
  $plot->write_frequency_matrix($fMatrixFile);
  $plot->write_status_matrix($matrixFile);
  $plot->write_value_matrix($vMatrixFile);
  $plot->write_svg($histoMultF);

  print "The SVG was written to "."$args{'-path_loc'}"."/"."$histoMultF.\n";

}

exit;

################################################################################
################################################################################
#        #######################################################################
#  subs  #######################################################################
#        #######################################################################
################################################################################
################################################################################

sub check_query {

  my $config = shift;
  my $query = shift;

  if (! grep{ /.../ } keys %$query) {
    push(@{ $config->{ERROR} }, <<END
No query was specified:
  -datasets arraymap,tcga -query icdom-817 -grouping icdot

END
);
  }
  return $config;

}

sub check_grouping {

  my $config = shift;
  if (! grep{ $config->{grouping} eq $_ } (keys %{ $config->{biosubsets} }, keys %{ $config->{datacollections} })) {
    push(@{ $config->{ERROR} }, '
The "-grouping" value was not one of
'. join("\n", keys %{ $config->{biosubsets} })."\n");
  }
  return $config;
}


1;