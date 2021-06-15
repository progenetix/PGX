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

Examples:
  - perl CNVfrequencyfilePlotter.pl -f frequencies.pgxseg
  - perl CNVfrequencyfilePlotter.pl -f ../rsrc/test_files/testfile_interval_frequencies.pgxseg -o ~/Desktop

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
	-min_group_no => 5,
	-size_title_left_w_px => 250,
	-size_clustertree_w_px => 50,
	-size_text_px => 12,
	-size_text_title_px => 16,
	-size_text_subtitle_px => 14,
	-size_text_title_left_px => 14,
	-size_clustertree_w_px => 100,
	-path_loc => $path_of_this_module.'/out',
	-plottype => 'multihistogram'
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

foreach (keys %args) {
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
my $histoplotF = $outFileBase.',samples_histoplot.svg';
my $histoplotmultF = $outFileBase.',samples_histoplot_mult.svg';

my $pgx = new PGX($plotargs);

$pgx->pgx_add_frequencymaps_from_file($args{'-f'});

if (defined $pgx->{pgxfileheader}->{plotpars}) {
	foreach (keys %{ $pgx->{pgxfileheader}->{plotpars} }) {
		$plotargs->{$_} = $pgx->{pgxfileheader}->{plotpars}->{$_};
		$pgx->{parameters}->{$_} = $pgx->{pgxfileheader}->{plotpars}->{$_};
	}
}

$pgx->cluster_frequencymaps();
$pgx->return_histoplot_svg();
$pgx->write_svg($histoplotmultF);

print "Wrote grouped histograms to\n\t=> $histoplotmultF\n";

1;
