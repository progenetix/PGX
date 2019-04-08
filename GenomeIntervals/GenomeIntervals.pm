package PGX::GenomeIntervals::GenomeIntervals;

use Data::Dumper;

require Exporter;
@ISA        =   qw(Exporter);
@EXPORT     =   qw(
  make_genome_intervals
  get_reference_base_limits
  get_genome_basecount
);

################################################################################

sub make_genome_intervals {

=pod

Expects:
  - a list reference of genome interval objects, usually representing cytobands:
    [
      {
        no    =>  __integer__,            # not used
        reference_name  =>  __string__,
        start =>  __integer__,
        end   =>  __integer__,
        stain =>  __string__,            # not used
        label =>  __string__,            # not used
      },
      {
      ...
      },
    ]
  - the binning size in bases (optional; defaults to 1000000)

Returns:
  - a list of genome intervals in a similar structure, based on the binning size
  - the intervals are per reference, i.e. each reference starts with a new interval,
    leading to the last interval < binning size
    [
      {
        no    =>  __integer__,          # 1 -> n
        reference_name  =>  __string__,
        start =>  __integer__,
        end   =>  __integer__,
        label =>  __string__,
      },
      {
      ...
      },
    ]

=cut

########    ####    ####    ####    ####    ####    ####    ####    ####    ####

  my (
    $cytobands,
    $intSize
  )             =   @_;

  if ($intSize !~ /^\d{3,9}$/) { $intSize = 1000000 }

  my $refLims   =   get_reference_base_limits($cytobands);
  my $gi        =   [];

  # references are sorted with numerical ones first, then others (e.g. 1 -> 22, X, Y)
  my @refNames  =   ((sort {$a <=> $b } grep{ /^\d\d?$/ } keys %$refLims), (sort grep{ ! /\d/ } keys %$refLims));

  my $intI      =   1;
  for my $i (0..$#refNames) {

    my $refName =   $refNames[$i];
    my $start   =   $refLims->{ $refName }->[0];
    my $end     =   $intSize - 1;

    while ($start < $refLims->{ $refName }->[1]) {

      # adjusting the end of the last interval
      if ($end > $refLims->{ $refName }->[1]) { $end = $refLims->{ $refName }->[1] };
      my $thisSize  =   $end - $start +1;
      push(
        @$gi,
        {
          no    =>  $intI,
          reference_name  =>  $refName,
          start =>  $start,
          end   =>  $end,
          length  =>  $thisSize,
          label =>  $refName.':'.$start.'-'.$end,
        }
      );

      $start    +=  $thisSize;
      $end      +=  $thisSize;
      $intI++;

  }}

  return $gi;

}

################################################################################

sub get_reference_base_limits {

=pod

Expects:
  - a list reference of genome interval objects:
    [
      {
        reference_name  =>  __string__,
        start           =>  __integer__,
        end             =>  __integer__,
        (others)...
      },
      {
        ...
      },
    ]

Returns:
  - a hash reference with each value consisting of a list of 2 integers, representing
    the reference's min and max bases:
    [
      { __string__  =>  [  __integer__,  __integer__ ] },
      { ... },
    ]

=cut

########    ####    ####    ####    ####    ####    ####    ####    ####    ####

  my $refLims   =   {};
  my $allRefs   =   $_[0];

  foreach my $ref (map{ $_->{reference_name} } @$allRefs) {
    my @refRefs =   grep{ $_->{reference_name} =~ /^$ref$/i } @$allRefs;
    my @bases   =   sort { $a <=> $b } ((map{ $_->{start} } @refRefs), (map{ $_->{end} } @refRefs));
    $refLims->{$ref}  =   [ $bases[0], $bases[-1] ];
  }

  return $refLims;

}

################################################################################

sub get_genome_basecount {

  my $genBases;
  my ($allRefs, $chr2plot)  =   @_;

  my %refNames  =   map { $_ => 1 } @$chr2plot;

  foreach my $ref (keys %refNames) {
    my @refRefs =   grep{ $_->{reference_name} =~ /^$ref$/i } @$allRefs;
    my @bases   =   sort { $a <=> $b } ((map{ $_->{start} } @refRefs), (map{ $_->{end} } @refRefs));
    $genBases   +=  ($bases[-1] - $bases[0]);
  }

  return $genBases;

}

################################################################################
########    utility subs    ####    ####    ####    ####    ####    ####    ####
################################################################################

1;
