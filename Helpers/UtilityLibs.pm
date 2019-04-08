package PGX::Helpers::UtilityLibs;

require Exporter;
@ISA    =   qw(Exporter);
@EXPORT =   qw(RandArr);

################################################################################

sub RandArr {

=pod
## RandArr

#### Expects:

* an array reference of arbitrary content
* the number of array elements to be returned (optional)

#### Return

* the re-shuffled array or a subset of its elements (as array reference)

This dependency-free array randomiser will return a re-shuffled array(ref) or
a slice of random $iL array elements.

The $overSamp factor is "empirical" and balances between oversampling with
out-of-range values + filtering, and cycling too many times to match all
index elements.

=cut

  my $arr       =   shift;
  my $iL        =   shift;
  my $overSamp  =   7;

  if (ref $arr ne 'ARRAY') { return \0 }

  # if no number of array elements => all
  if ($iL !~ /^\d+?$/) {
    $iL = scalar @$arr }
  # ... not more than all
  elsif ($iL > @$arr) {
    $iL = scalar @$arr }

  $overSamp *= $iL;

  # maximum index number, for filtering the oversampled values
  my $maxI      =   @$arr - 1;
  if ($maxI < 0) { return \0 }

  # use of a hash to have unique index numbers (keys of the hash)
  my %randNo    =   ();

  # adding to the hash keys until equal or more than needed
  while (keys %randNo < $iL) {
    %randNo     =   map{ $_ => 1 } (grep{ $_ <= $maxI } (keys %randNo, map{ int(rand($_)) } 0..$overSamp) );
  }

  return [ @$arr[ (keys %randNo)[0..($iL-1)] ] ];

}

1;
