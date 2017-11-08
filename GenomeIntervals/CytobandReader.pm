package PGX::GenomeIntervals::CytobandReader;

use Data::Dumper;

require Exporter;
@ISA    =   qw(Exporter);
@EXPORT =   qw(read_cytobands genome_names_to_hg genome_names_to_grch);

sub read_cytobands {

=pod

Expects:
  - as input a genome edition
  - a UCSC "cytoBandIdeo.txt" file of the given genome edition (hardcoded location):

chr1  0 2300000 p36.33  gneg
chr1  2300000 5300000 p36.32  gpos25
chr1  5300000 7100000 p36.31  gneg
chr1  7100000 9200000 p36.23  gpos25
chr1  9200000 12600000  p36.22  gneg
...

Returns:
  - a list reference of genome interval objects, usually representing cytobands:
    [
      {
        no    =>  __integer__,          # 1 -> n
        reference_name  =>  __string__,
        start =>  __integer__,
        end   =>  __integer__,
        stain =>  __string__,
        label =>  __string__,           # name of cytoband, e.g. "8q24.1"
      },
      {
      ...
      },
    ]

=cut

########    ####    ####    ####    ####    ####    ####    ####    ####    ####

  use File::Basename;

  my $path_of_this_module = File::Basename::dirname( eval { ( caller() )[1] } );
  my $cbF   =   $path_of_this_module.'/../rsrc/genomes/'.genome_names_to_hg($_[0]).'/cytoBandIdeo.txt';
  my @cb    =   ();

  open  FILE, "$cbF" or die "No file $cbF $!";
  my $i     =   0;
  foreach (grep{/^(chro?)?[\dXYZ]\d?\t/i} <FILE>) {
    $i++;
    chomp;
    my (
      $chro,
      $start,
      $end,
      $band,
      $staining,
    )       =   split (/\t/, $_, 6);

    $chro   =~ s/chr//;
    push(
      @cb,
      {
        no      =>  $i,
        reference_name  =>  $chro,
        start   =>  1 * $start,
        end     =>  1 * $end,
        stain   =>  $staining,
        label   =>  $chro.$band,
      }
    );
  }

  close FILE;

  return [@cb];

}

################################################################################
########    utility subs    ####    ####    ####    ####    ####    ####    ####
################################################################################

sub genome_names_to_grch {
  my %genome_names   =   (
    hg18    =>  'grch36',
    hg19    =>  'grch37',
    hg38    =>  'grch38',
  );
  if ($_[0] =~ /^grch\d\d$/i) {
    return lc($_[0]) }
  else {
    return $genome_names{lc($_[0])} }
}

################################################################################

sub genome_names_to_hg {
  my %genome_names   =   (
    grch36  =>  'hg18',
    grch37  =>  'hg19',
    grch38  =>  'hg38',
  );
  if ($_[0] =~ /^hg\d\d$/i) {
    return lc($_[0]) }
  else {
    return $genome_names{lc($_[0])} }
}

1;
