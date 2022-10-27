#!/usr/bin/perl -w

$| = 1;

################################################################################
#                                                                              #
# Progenetix site scripts                                                      #
#                                                                              #
# molecular cytogenetics, Comparative Genomic Hybridization, genomic arrays    #
# data analysis & visualization                                                #
#                                                                              #
# © 2000-2022 Michael Baudis: michael@baud.is                                  #
#                                                                              #
################################################################################

=podmd

=cut

use strict;
use CGI::Simple;
use CGI::Carp qw(fatalsToBrowser);
use Data::Dumper;

use File::Basename;
use JSON::XS;
use MongoDB;
$MongoDB::Cursor::timeout = 120000;

BEGIN { unshift @INC, ('..') };
use PGX;
use lib::CGItools;
use lib::Helpers;

my $config = PGX::read_config();
my $params = lib::CGItools::deparse_query_string();

my $debug_mode = 0;
if ($params->{debug}->[0] > 0) {
	$debug_mode = 1;
	print 'Content-type: text/plain'."\n\n";
}

my $accessid = $params->{accessid}->[0];

if (! $params->{datasetIds}) {
	$params->{datasetIds} = [ 'progenetix' ] }
if (
	$params->{datasetIds}->[0] =~ /undefined/
	||
	$params->{datasetIds}->[0] =~ /null/
) {
	$params->{datasetIds} = [ 'progenetix' ] }
	
my $api = {
	config => $config,
	debug_mode => $debug_mode,
	datasetid => '',
	coll => '',
	path_var => '/_process_'."$^T"."$$",
# 	plotargs =>	{ map{ $_ => join(',', @{ $params->{$_} }) } (grep{ /^\-\w+?$/ } keys %{ $params }) },
	plotargs =>	{ map{ $_ => join(',', @{ $params->{$_} }) } (grep{ /^\-?\w+?$/ } keys %{ $params }) },
	handover_db => $config->{handover_db},
	handover_coll => $config->{handover_coll},
	accessid => $accessid,
	segfile => $config->{paths}->{dir_tmp_base_path}.'/'.$accessid,
	technology_keys => $config->{technology_keys},
	server_link => ($ENV{HTTPS} ? 'https://' : 'http://').$ENV{HTTP_HOST},
	data => { },
	output => 'json',
	errors => [ ],
	warnings => [ ]
};
bless $api;

$api->{plotargs}->{-path_loc} = $api->{config}->{paths}->{dir_tmp_base_path}.$api->{path_var};
$api->{server_link} =~ s/^(https?\:\/\/)\w+?\.(\w+?\.\w+?(\/.*?)?)$/$1$2/;
$api->{plotargs}->{-path_web} = $api->{server_link}.$api->{config}->{paths}->{web_tmp_base_path}.$api->{path_var};

if ($accessid !~ /^[\w\-]+$/) {
	push(
  		@{ $api->{errors} },
		"Wrong or missing accessid parameter $api->{accessid}.",
		"Wrong or missing segments file."
	);
}

################################################################################

$api->_retrieve_samples();

=podmd

A single histogram SVG output can be forced with `&output=cnvhistogram`.

=cut

if (grep{ $params->{output}->[0] eq "cnvhistogram" } qw(output method)) {

	$api->{plotargs}->{ '-svg_embed' } = 1;
	$api->_return_histogram();
	print 'Content-type: image/svg+xml'."\n\n";
	print $api->{data}->{plots}->{ histogram }->{svg};
	exit;
	
}

################################################################################

mkdir $api->{plotargs}->{-path_loc};
$api->_return_histogram();
$api->_add_samplecollections();
$api->_return_multihistogram();
$api->_return_samplematrix();

send_Google_tracking_no_log(
	$config->{cgi},
	"/cgi/PGX/cgi/samplePlots.cgi"
);

$api->_return_json();

################################################################################
# subs #########################################################################
################################################################################

sub _retrieve_samples {

	my $api = shift;

	my $pgx = new PGX($api->{plotargs}, $api->{debug_mode});

	$pgx->pgx_add_segments_from_file($api->{segfile});

	if (defined $pgx->{pgxfileheader}->{plotpars}) {
		foreach (keys %{ $pgx->{pgxfileheader}->{plotpars} }) {
			$api->{plotargs}->{"-".$_} = $pgx->{pgxfileheader}->{plotpars}->{$_};
			$pgx->{parameters}->{$_} = $pgx->{pgxfileheader}->{plotpars}->{$_};
		}
	}
	$pgx->pgx_create_samples_from_segments();
	if (! -f $api->{segfile}) {
		$pgx->pgx_open_handover($api->{config}, $api->{accessid});
		$pgx->pgx_samples_from_handover();
	}
	if ($api->{plotargs}->{'-randno'} > 0) {
		$pgx->{samples} = RandArr($pgx->{samples}, $api->{plotargs}->{'-randno'}) }
	$pgx->pgx_callset_labels_from_biosamples($api->{config});
	$pgx->pgx_callset_labels_from_header();
	$pgx->pgx_add_variants_from_db();
	$pgx->pgx_remap_vrsified_segments();

	$api->{samples}	= $pgx->{samples};
	$api->{datasetid} = $pgx->{datasetid};
	
	return $api;

}

################################################################################

sub _return_histogram {

	my $api = shift;

	my $plotType = 'histogram';
	
	# modifying a copy of the standard plot arguments for the overview plot
	my $plotargs = bless { %{ $api->{plotargs} } }, ref $api->{plotargs};
	$plotargs->{'-plottype'} = $plotType;
	$plotargs->{-size_title_left_w_px} = 0;
	
	if (-f $api->{segfile}) {
		$plotargs->{-text_bottom_left} = 'Uploaded: '.scalar(@{ $api->{samples} }).' samples' }
	else {	
		$plotargs->{-text_bottom_left} = $api->{datasetid}.': '.scalar(@{ $api->{samples} }).' samples' }	

	my $pgx = new PGX($plotargs, $api->{debug_mode});

	$pgx->{datasetid} = $api->{datasetid};
	$pgx->{parameters}->{plotid} = 'histogram';
	$pgx->pgx_add_frequencymaps( [ { statusmapsets => $api->{samples} } ] );
	$pgx->return_histoplot_svg();
	$pgx->write_svg();
	
	if ($pgx->{parameters}->{svg_embed} > 0) {
		$api->{data}->{plots}->{ $plotType }->{svg} = $pgx->{svg} }
	$api->{data}->{plots}->{ $plotType }->{svg_link_tmp} = $pgx->{svg_path_web};
	
	return $api;
	
}

################################################################################

sub _add_samplecollections {

	my $api = shift;

	my $pgx = new PGX($api->{plotargs}, $api->{debug_mode});
	$pgx->{samples} = $api->{samples};

	$pgx->pgx_create_sample_collections();	
	$api->{samplecollections} = $pgx->{samplecollections};
	
	return $api;
	
}

################################################################################
# (clustered) CNA histograms
################################################################################

sub _return_multihistogram {

	my $api = shift;
	
	if (@{ $api->{samplecollections} } < 2) {
		return $api }
		
	my $plotType = 'multihistogram';
	
	my $plotargs = $api->{plotargs};
	$plotargs->{'-plottype'} = $plotType;	

	if (-f $api->{segfile}) {
		$plotargs->{-text_bottom_left} = 'Uploaded: '.scalar(@{ $api->{samples} }).' samples' }
	else {	
		$plotargs->{-text_bottom_left} = $api->{datasetid}.': '.scalar(@{ $api->{samples} }).' samples' }	
	
	my $pgx = new PGX($plotargs, $api->{debug_mode});
	$pgx->{samples} = $api->{samples};

	$pgx->{datasetid} = $api->{datasetid};
	$pgx->{parameters}->{plotid} = $plotType;
	$pgx->pgx_add_frequencymaps($api->{samplecollections});
	$pgx->cluster_frequencymaps();
	$pgx->return_histoplot_svg();
	$pgx->write_frequency_matrix($api->{plotargs}->{-path_loc}.'/frequencymatrix.tsv');
	$pgx->write_svg();

	if ($pgx->{parameters}->{svg_embed} > 0) {
		$api->{data}->{plots}->{ $plotType }->{svg} = $pgx->{svg} }

	$api->{data}->{plots}->{ $plotType }->{svg_link_tmp} = $pgx->{svg_path_web};
	$api->{data}->{data_files}->{frequencymatrix}->{"link"} = $api->{plotargs}->{-path_web}.'/frequencymatrix.tsv';

	return $api;
	
}

# ################################################################################
# sample matrix ##################################################################
# ################################################################################

sub _return_samplematrix {

	my $api = shift;

	my $plotType = 'multistrip';
	my $plotargs = $api->{plotargs};
	$plotargs->{'-plottype'} = $plotType;
	
	if (-f $api->{segfile}) {
		$plotargs->{-text_bottom_left} = 'Uploaded: '.scalar(@{ $api->{samples} }).' samples' }
	else {	
		$plotargs->{-text_bottom_left} = $api->{datasetid}.': '.scalar(@{ $api->{samples} }).' samples' }	
	
	my $pgx = new PGX($plotargs, $api->{debug_mode});
	$pgx->{samples} = $api->{samples};
	$pgx->{datasetid} = $api->{datasetid};
	$pgx->{parameters}->{plotid} = $plotType;
	$pgx->cluster_samples();
	$pgx->return_stripplot_svg();
	$pgx->write_svg();

	if ($pgx->{parameters}->{svg_embed} > 0) {
		$api->{data}->{plots}->{ $plotType }->{svg} = $pgx->{svg} }
	$api->{data}->{plots}->{ $plotType }->{svg_link_tmp} = $pgx->{svg_path_web};
	
	$pgx->write_status_matrix();
	$api->{data}->{samplematrix_link_tmp} = $pgx->{samplematrix_link_tmp};
	
	return $api;

}

################################################################################
################################################################################
################################################################################

sub _return_json {
	
	my $api = shift;
  
	if ($api->{output} =~ /json/i) {
	  print	'Content-type: application/json'."\n\n";
	  print JSON::XS->new->pretty( 1 )->allow_blessed()->convert_blessed()->encode( { data => $api->{data}, errrors => $api->{errors}, warnings => $api->{warnings} } )."\n";
	  exit;
	}
	
	return $api;

}


1;
