#!/usr/bin/perl -w

$| = 1;

# CPAN packages
use Data::Dumper;
use File::Basename;
use POSIX 'strftime';
use strict;
#use Term::ProgressBar;
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

```
# biosample_id=GSM481286;group_id=NCIT:C4017;group_label=Ductal Breast Carcinoma
# biosample_id=GSM481418;group_id=NCIT:C3059;group_label=Glioma
biosample_id	chro	start	stop	mean	probes	variantType
GSM481286	1	742429	7883881	-0.1594	699	DEL
GSM481286	2	115673158	115705254	-0.3829	8	DEL
GSM481286	3	115722621	119771659	0.167	424	DUP
GSM481286	4	119776776	162617092	0.4168	1587	DUP
GSM481418	5	162621657	165278686	.	.	DUP
GSM481418	6	165280711	167221337	.	.	DUP
GSM481418	7	167248788	168289603	0.6784	.	DUP	
...
```

Examples:
  - perl CNVsegfilePlotter.pl -f segments.pgxseg

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
  -size_title_left_w_px => 200,
  -size_clustertree_w_px => 50,
  -size_text_px => 12,
  -size_text_title_px => 16,
  -size_text_subtitle_px => 14,
  -size_text_title_left_px => 10,
  -path_loc => $path_of_this_module.'/out',
};

# command line input
my %args = @ARGV;

if (-d $args{-o}) {
  $plotargs->{-path_loc} = $args{-o} }
else {
  print <<END;
No existing output path was specified:
  -o ___path_to_my_directory___

END
	exit;
}

foreach (keys %$plotargs) {
  if ($args{$_} =~ /\w/) {
    $plotargs->{$_} = $args{$_} }
}

if (! -f $args{'-f'}) {
  print <<END;
No input file was specified:
  -f ___path_to_my_file___

END
	exit;
}

# files and paths ##############################################################

my $outFileBase = $args{'-f'};
$outFileBase =~ s/^.*?\/([^\/]+?)\.\w{2,8}$/$1/;
my $matrixplotF = $outFileBase.',samples_clustered.svg';
my $histoplotF = $outFileBase.',samples_histoplot.svg';
my $histoplotmultF = $outFileBase.',samples_histoplot_mult.svg';

# grouping of samples  #########################################################

my $pgx = new PGX($plotargs);

# this uses the file reading routine; but multi-segment files have to be
# deconvoluted first ...
$pgx->pgx_add_segments_from_file($args{'-f'});

if (defined $pgx->{segfileheader}->{plotpars}) {
	foreach (keys %{ $pgx->{segfileheader}->{plotpars} }) {
		$plotargs->{$_} = $pgx->{segfileheader}->{plotpars}->{$_};
		$pgx->{parameters}->{$_} = $pgx->{segfileheader}->{plotpars}->{$_};
	}
}

$pgx->pgx_create_samples_from_segments();
$pgx->pgx_callset_labels_from_file($args{'-f'});
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
$plot->{parameters}->{size_title_left_w_px} = 0;
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
