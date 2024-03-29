package lib::CGItools;

use CGI::Simple;

require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(
	deparse_query_string
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

1;
