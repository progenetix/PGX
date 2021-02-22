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
# uploader.cgi                                                                   #
#                                                                              #
################################################################################

use CGI::Carp qw(fatalsToBrowser);
use CGI::Upload;
use Data::Dumper;
use File::Basename ;
use JSON::XS;
$CGI::Simple::DISABLE_UPLOADS = 0;   # enable uploads
$CGI::Simple::POST_MAX = 1024 * 5000;

BEGIN { unshift @INC, ('..') };

use PGX;

my $upload = CGI::Upload->new;
my $config = PGX::read_config();

my $accessid = "$^T"."$$";

my $response = {
	error => q{},
	rel_path =>	$config->{paths}->{web_tmp_base_path}.'/'.$accessid,
	loc_path =>	$config->{paths}->{dir_tmp_base_path}.'/'.$accessid,
	accessid => $accessid,
	plot_link => $config->{url_base}.'/cgi/PGX/cgi/samplePlots.cgi?accessid='.$accessid,
	host => $config->{url_base},
	upload => $upload,
};

bless $response;

$response->_check_file_params();
$response->_return_on_error();
$response->_return_file_path();

exit;

################################################################################
################################################################################
################################################################################

sub _check_file_params {

	my $response = shift;	
	my $f_name = $response->{upload}->file_name('upload_file_name');

	if ($f_name =~ /[^\w\.\,\-]/) {
    $response->{error}->{error_message}	.=	'The file name may contain incorrect characters (e.g. spaces). Please rename your file using only word characters, dashes, commas and full stops.' }
	if ($f_name !~ /\w\.\w{2,7}$/) {
    $response->{error}->{error_message}	.=	"The file name $f_name seems incomplete, e.g. lacking an extension?" }

	return $response;
	
}

################################################################################

sub _return_on_error {

	my $response = shift;
	if ($response->{error} =~ /.../) {
		print	'Content-type: application/json'."\n\n";
		print JSON::XS->new->pretty( 0 )->allow_blessed->convert_blessed->encode( { error => $response->{error}, upload => $response->{upload}, param => $response->{param} } )."\n";
		exit;
	}	
	return $response;
}

################################################################################

sub _return_file_path {

	my $response = shift;

	my $f_handle = $response->{upload}->file_handle('upload_file_name');
	
	if ($f_handle) {

		my $count =  0;
		my %ids = ();
		open ( UPLOADFILE, ">$response->{loc_path}" ) or die "$!";
		while ( <$f_handle> ) { 
			my $in  = $_;
			chomp $in;
			if ($in =~ /^\#/) {
				print UPLOADFILE $in."\n";
				next;
			}
			my @line  =   split("\t", $in);
			$line[0]  =~  s/[^\w\,\.\-]//g;
			if ($line[1] =~ /^(chro?)?[\dXY]\d?$/) {
				$ids{ $line[0] } = 1;
				print UPLOADFILE $in."\n";
				$count++;
		}}
		close UPLOADFILE;
		print 'Content-type: application/json'."\n\n";
		print JSON::XS->new->pretty( 0 )->allow_blessed->convert_blessed->encode( { accessid => $response->{accessid}, file => $response->{host}.$response->{rel_path}, sample_count => scalar(keys %ids), seg_count => $count, actions => { plot_segment_file => $response->{plot_link} } } )."\n";
		exit;
	}	
	return $response;	
}

1;
