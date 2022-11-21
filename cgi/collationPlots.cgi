#!/usr/bin/env perl -w

$| = 1;

################################################################################
#                                                                              #
# Progenetix & arrayMap site scripts                                           #
#                                                                              #
# molecular cytogenetics, Comparative Genomic Hybridization, genomic arrays    #
# data analysis & visualization                                                #
#                                                                              #
# © 2000-2021 Michael Baudis: m@baud.is                                        #
#                                                                              #
################################################################################

=pod
https://progenetix.org/cgi/PGX/cgi/collationPlots.cgi?datasetIds=progenetix&id=NCIT:C3262&-size_plotimage_w_px=1084
=cut

use strict;
use CGI::Carp qw(fatalsToBrowser);
use CGI::Simple;
use Data::Dumper;
use MongoDB;

BEGIN { unshift @INC, ('..') };
use PGX;

my $config = PGX::read_config();
my $params = lib::CGItools::deparse_query_string();

my $debug_mode = 0;
if ($params->{debug}->[0] > 0) {
	$debug_mode = 1;
	print 'Content-type: text/plain'."\n\n";
}

################################################################################

if (scalar @{ $params->{id} } != 1) {
	print 'Content-type: text'."\n";
	print 'status: 422'."\n\n";
	print 'No single id provided...';
	exit;
}

if (! $params->{datasetIds}) {
	$params->{datasetIds} = [ 'progenetix' ] }
if (scalar @{ $params->{datasetIds} } != 1) {
	$params->{datasetIds} = [ 'progenetix' ] }
if (
	$params->{datasetIds}->[0] =~ /undefined/
	||
	$params->{datasetIds}->[0] =~ /null/
) {
	$params->{datasetIds} = [ 'progenetix' ] }

my $id = $params->{id}->[0];
my $dataset = $params->{datasetIds}->[0];
my $subset = MongoDB::MongoClient->new()->get_database( $dataset )->get_collection( "frequencymaps" )->find_one( { id => $id } );

my $plotargs = { map{ $_ => join(',', @{ $params->{$_} }) } (grep{ /^\-\w+?$/ } keys %{ $params }) };

if (! grep{ /id/ } keys %{ $subset }) {
	my $plot = new PGX($plotargs, $debug_mode);
	$plot->return_error_svg('No match for histogram with id "'.$id.'" in dataset "'.$dataset.'" ...');
	print 'Content-type: image/svg+xml'."\n\n";
	print $plot->{svg};
	exit;
}

if ($subset->{label} =~ /.../) {
	$plotargs->{-title} = $subset->{label}.' ('.$subset->{id}.')' }
else {
	$plotargs->{-title} = $subset->{id} }

my $plot = new PGX($plotargs, $debug_mode);
$plot->{rameters}->{plotid} = 'histoplot';
$plot->{parameters}->{text_bottom_left} = $subset->{counts}->{biosamples}.' samples';

$plot->{frequencymaps} = [ $subset->{ frequencymap } ];
$plot->return_histoplot_svg();

lib::CGItools::send_Google_tracking_no_log(
	$config->{cgi},
	"/cgi/PGX/cgi/collationPlots.cgi"
);

print 'Content-type: image/svg+xml'."\n\n";
print $plot->{svg};

exit;

################################################################################


1;
