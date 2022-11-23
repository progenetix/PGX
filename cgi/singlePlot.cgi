#!/usr/bin/env perl -w

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

my $biosample_id = $params->{biosampleIds}->[0] || "__none__";
my $analysis_id = $params->{analysisIds}->[0] || "__none__";
my $probeF = $params->{probefile}->[0];


if (! $params->{datasetIds}) {
	$params->{datasetIds} = [ 'progenetix' ] }
if (
	$params->{datasetIds}->[0] =~ /undefined/
	||
	$params->{datasetIds}->[0] =~ /null/
) {
	$params->{datasetIds} = [ 'progenetix' ] }

my $db = $params->{datasetIds}->[0];

my $api = {
	config => $config,
	debug_mode => $debug_mode,
	datasetid => '',
	coll => '',
	path_var => '/_process_'."$^T"."$$",
	plotargs =>	{ map{ $_ => join(',', @{ $params->{$_} }) } (grep{ /^\-?\w+?$/ } keys %{ $params }) },
	handover_db => $config->{handover_db},
	handover_coll => $config->{handover_coll},
	technology_keys => $config->{technology_keys},
	server_link => ($ENV{HTTPS} ? 'https://' : 'http://').$ENV{HTTP_HOST},
	data => { },
	output => 'svg',
	errors => [ ],
	warnings => [ ]
};
bless $api;

$api->{plotargs}->{-path_loc} = $api->{config}->{paths}->{dir_tmp_base_path}.$api->{path_var};
$api->{plotargs}->{-plotid} = 'arrayplot';
$api->{plotargs}->{-plottype} = 'array';
$api->{server_link} =~ s/^(https?\:\/\/)\w+?\.(\w+?\.\w+?(\/.*?)?)$/$1$2/;
$api->{plotargs}->{-path_web} = $api->{server_link}.$api->{config}->{paths}->{web_tmp_base_path}.$api->{path_var};

mkdir $api->{plotargs}->{-path_loc};
	
################################################################################
################################################################################
################################################################################

if (! -f $probeF) {
	$probeF = "" }

my $query = {};
if ($biosample_id ne "__none__") {
	$query = { "biosample_id" => $biosample_id }
} elsif ($analysis_id ne "__none__") {
	$query = { "id" => $analysis_id }
} else {
	print 'Content-type: text/plain'."\n\n";
	print "no parameter for biosampleIds or callsetIds";
	exit;
}

my $pgx = new PGX($api->{plotargs}, $api->{debug_mode});
$pgx->{dataconn} = MongoDB::MongoClient->new()->get_database( $db );
$pgx->pgx_samples_from_callsets($query);
$pgx->pgx_add_variants_from_db($query);
$pgx->pgx_segmentdata_from_sample_0();
$pgx->pgx_remap_vrsified_segments( $pgx->{segmentdata} );
$pgx->pgx_add_probes_from_file($probeF);
$pgx->{parameters}->{text_bottom_left} = $pgx->{samples}->[0]->{id};
$pgx->return_arrayplot_svg();

if ($params->{print_file}->[0] !~ /^(?:y|1|t)/i) {
	print 'Content-type: image/svg+xml'."\n\n" ;
	print $pgx->{svg};
	exit;
}

mkdir $api->{plotargs}->{-path_loc};
$pgx->write_svg();
print <<END;
Content-type: text/html

<html>
	<img src="$api->{plotargs}->{-path_web}/$api->{plotargs}->{-plotid}.svg" />
</html>
END

exit;

1;
