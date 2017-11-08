package PGX::GenomeIntervals::IntervalStatistics;

use Data::Dumper;

require Exporter;
@ISA        =   qw(Exporter);
@EXPORT     =   qw(interval_cnv_frequencies);

################################################################################

sub interval_cnv_frequencies {

=pod

Expects:

Returns:

=cut

########    ####    ####    ####    ####    ####    ####    ####    ####    ####

  my $cnvmaps   =   [];
  ($cnvmaps, $plotPars->{GENOINTVS})        =   @_;

  my $maps      =   {
    intervals   =>  scalar(@{ $plotPars->{GENOINTVS} }),
    binning     =>  $plotPars->{GENOINTVS}->[0]->{end} - $plotPars->{GENOINTVS}->[0]->{start} + 1,
  };
  my %intLabs   =   (
    DUP         =>  'dupmap',
    DEL         =>  'delmap',
  );
  my %freqLabs      =   (
    DUP         =>  'dupfrequencies',
    DEL         =>  'delfrequencies',
  );

  # avoiding division by 0 errors if improperly called
	my $fFactor		=   100;
	if (@{ $cnvmaps } > 1) { $fFactor = 100 / @{ $cnvmaps } }

  foreach my $type (keys %intLabs) {
  	for my $i (0..$#{ $plotPars->{GENOINTVS} }) {
      $maps->{ $freqLabs{ $type } }->[$i]	  =   sprintf "%.3f", ($fFactor * ( grep{ $_->{ $intLabs{ $type } }->[$i] eq $type } @{ $cnvmaps } ));
  	}
  }

  return $maps;

}

1;
