package PGX::WebUtilities::EuropmcRetriever;

use Data::Dumper;
use Encode qw(decode encode);
use LWP::Simple;


require Exporter;
@ISA        =   qw(Exporter);
@EXPORT     =   qw(retrieve_epmc_core);

sub retrieve_epmc_core {

=pod

expects:
  - one pmid

returns:
  - an object with article attributes extracted from the europepmc core data annotation

notes:
  - not using XML::Simple here, due to problems with the some types (e.g. "Abstract"
    sometimes contains the abstract, sometimes fields with abstract parts => funky ...;
    the specific solution for this type of file is also faster than the module call
  - still traditional Progenetix publication attributes ...

=cut

  my $pmid  =   $_[0];
  my $pub   =   {};

  my $epmc_web  =   'http://europepmc.org/abstract/MED/';
  my $epmc_q    =   'https://www.ebi.ac.uk/europepmc/webservices/rest/search?resulttype=core&query=ext_id:'.$pmid.'%20src:med';
  my $epmcxml   =   get($epmc_q);
  $epmcxml      =   $epmcxml;
  $epmcxml      =~  s/\n//gs;
  $epmcxml      =~  s/>\s+?</></gs;

  if ($epmcxml =~ /<doi>([^<]+?)<\/doi/) {
    $pub->{DOI} =   $1 }

  if ($epmcxml =~ /<pmcid>(PMC\d+?)<\/pmcid/) {
    $pub->{pmcid}   =   $1 }

  if ($epmcxml =~ /<pmid>(\d+?)<\/pmid/) {
    $pub->{pubmedid}    =   $1 }

  if ($epmcxml =~ /<affiliation>(.+?)<\/affiliation/) {
    $pub->{affiliation} =   $1 }

  if ($epmcxml =~ /<title>(.+?)<\/title/) {
    $pub->{TITLE}       =   $1 }

  if ($epmcxml =~ /<journalInfo>(.+?)<\/journalInfo>/) {

    my $journalData     =   $1;

    if ($journalData =~ /<title>(.+?)<\/title>/) {
      $pub->{JOURNALTITLE}      =   $1 }
    if ($journalData =~ /<ISOAbbreviation>([^>]+?)<\/ISOAbbreviation>/) {
      $pub->{journal}   = $1 }
    if ($journalData =~ /<volume>([^>]+?)<\/volume>/) {
      $pub->{journal}   .=  ' '.$1 }
    if ($journalData =~ /<issue>([^>]+?)<\/issue>/) {
      $pub->{journal}   .=  '('.$1.')' }
    if ($journalData =~ /<yearOfPublication>(\d\d\d\d)<\/yearOfPublication>/) {
      $pub->{journal}   .=  ', '.$1;
      $pub->{year}      =   $1;

  }}

  if ($epmcxml =~ /<abstractText>(.+?)<\/abstractText>/) {
    $pub->{abstract}    = $1;
    $pub->{abstract}    =~  s/<[^<]*?>/######/g;
    $pub->{abstract}    =~  s/^######//g;
    $pub->{abstract}    =~  s/######$//g;
    $pub->{abstract}    =~  s/######/ /g;
    $pub->{abstract}    =~  s/######/ /g;
  }

  if ($epmcxml =~ /<authorString>(.+?)<\/authorString>/) {
    $pub->{authors} =   $1 }

  $pub->{label}   =   $pub->{authors};
  $pub->{label}   =~  s/^(.+?\w\w\w) .*?$/$1 et al./;
  my $shortT    =   $pub->{TITLE};
  $shortT       =~    s/^(.{50,100} ).*?$/$1.../;;

  $pub->{label}   .=  ' ('.$pub->{year}.'): '.$shortT;

  return $pub;

}

1;
