package lib::CGItools;

use CGI::Simple;
use Net::Google::Analytics::MeasurementProtocol;
use UUID::Tiny;

require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(
	deparse_query_string
	send_Google_tracking_no_log
);

################################################################################

sub deparse_query_string {

=podmd
#### Deparsing the query string

The query string is deparsed into a hash reference, in the "$query" object,
with the conventions of:

* each parameter is treated as containing a list of values
* values are __split into a list by the comma character__; so an example of
    `key=val1&key=val2,val3`
  would be deparsed to
    `key = [val1, val2, val3]`

The treatment of each attribute as pointing to an array leads to a consistent,
though sometimes awkward, access to the values; even consistently unary values
have to be addressed as (first) array element.

=cut

	my $cgi = CGI::Simple->new;
	my $params = {};
	foreach my $qkey ($cgi->param()) {
		my @qvalues = $cgi->param($qkey);
		$params->{ $qkey } = [];
		if ($qkey =~ /\w/ && grep{ /./ } @qvalues) {
			foreach my $val (grep{ /./} split(',', join(',', @qvalues))) {
				push(@{ $params->{ $qkey } }, $val);
	}}}

	return $params;

}

################################################################################

sub send_Google_tracking_no_log {

	my $conf_cgi = shift;
	my $web_root = shift;
	
	my $userAgent = $ENV{ HTTP_USER_AGENT };
	my $userIP = ( $ENV{ HTTP_CF_CONNECTING_IP } =~ /\d\d\d/ ? $ENV{ HTTP_CF_CONNECTING_IP } : $ENV{ REMOTE_ADDR } );
	my $lang = $ENV{ HTTP_ACCEPT_LANGUAGE };
	my $query = $ENV{ QUERY_STRING };

	my %googleParams = (
	tid => $conf_cgi->{google_tid},
	ua  => $userAgent,
	cid => create_UUID_as_string(UUID_V3, $userIP.'_'.$userAgent),
	uip => $userIP,
	);

	if ($lang =~ /^[^\,]*?(\w\w(\-\w\w)?)/) { $googleParams{ul} = lc($1) }

	my $ga =  Net::Google::Analytics::MeasurementProtocol->new( %googleParams );
	$ga->send(
		'pageview',
		{
			dr => $ENV{ HTTP_REFERER },
			dt => $conf_cgi->{google_dt},
			dp => $api_path.'?'.$query,
			dh => $web_root,
			ds => 'web',
		}
	);

}

1;
