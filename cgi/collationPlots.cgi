#!/usr/bin/perl -w

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
use lib::CGItools;

my $config = PGX::read_config();
my $params = lib::CGItools::deparse_query_string();

if ($params->{debug}->[0] > 0) {
	print 'Content-type: text/plain'."\n\n" }

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

my $id = $params->{id}->[0];
my $dataset = $params->{datasetIds}->[0];
my $subset = MongoDB::MongoClient->new()->get_database( $dataset )->get_collection( "frequencymaps" )->find_one( { id => $id } );

if (! grep{ /id/ } keys %{ $subset }) {
	print 'Content-type: text'."\n";
	print 'status: 422'."\n\n";
	print 'No match for histogram with id "'.$id.'" in dataset "'.$dataset.'" ...';
	exit;
}

my $plotargs = { map{ $_ => join(',', @{ $params->{$_} }) } (grep{ /^\-\w+?$/ } keys %{ $params }) };
if ($subset->{label} =~ /.../) {
	$plotargs->{-title} = $subset->{label}.' ('.$subset->{id}.')' }
else {
	$plotargs->{-title} = $subset->{id} }

my $plot = new PGX($plotargs);
$plot->{parameters}->{plotid} = 'histoplot';
$plot->{parameters}->{text_bottom_left} =   $subset->{counts}->{biosamples}.' samples';
$plot->{frequencymaps} = [ $subset->{ frequencymaps } ];
$plot->return_histoplot_svg();

send_Google_tracking_no_log(
	$config->{cgi},
	"/cgi/PGX/cgi/collationPlots.cgi"
);

print 'Content-type: image/svg+xml'."\n\n";
print $plot->{svg};

exit;

################################################################################


1;
