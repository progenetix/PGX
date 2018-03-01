package PGX::FileUtilities::PlotfileReader;

use Data::Dumper;
use Math::Random qw(random_normal);

require Exporter;
@ISA    =   qw(Exporter);
@EXPORT =   qw(read_probefile read_segmentfile);

################################################################################

sub read_probefile {

=pod

Expects:
  - a standard Progenetix style probe file

  ID  chro  pos log2
  cnvi0111187 17  35295593  0.0859121900
  cnvi0111188 8 65499402  -0.1438023000
  cnvi0111189 2 177061178 -0.0113166000
  cnvi0111190 5 70255894  0.0463862400
  ...

Returns:
  - a list reference of genome position / value objects:
    [
      {
        no              =>  __integer__,          # 1 -> n
        probe_id        =>  __string__,
        reference_name  =>  __string__,
        position        =>  __integer__,
        value           =>  __long__,
      },
      {
      ...
      },
    ]

=cut

########    ####    ####    ####    ####    ####    ####    ####    ####    ####

  my $plot      =   shift;
  my $probeF    =   shift;
  my $probeT    =   shift;
  $probeT       ||= 'probedata';

  $plot->{$probeT}    =   [];
  my @randomV;

  if (! -f $probeF) { return $plot->{$probeT} }

  my $numfactor =   1;
  if (
    $plot->{parameters}->{'reverse'} =~ /y/i
    &&
    $probeT !~ /frac/i
  ) { $numfactor = -1 }

  open  FILE, "$probeF" or die "No file $probeF $!";
  local   $/;                             # no input separator
  my $fContent  =   <FILE>;
  close FILE;
  my @probeData =   split(/\r\n?|\n/, $fContent);
  shift @probeData;

  my $i         =   0;
  foreach (@probeData) {

    $i++;
    my (
      $probe_id,
      $reference_name,
      $position,
      $value,
    )           =   split (/\s/, $_, 5);
    $probe_id   =~  s/[^\w\-\,]/_/g;
    $reference_name     =~ s/[^\dxXyY]//;
    $reference_name     =~ s/^23$/X/;
    $reference_name     =~ s/^24$/Y/;
    $position   =   sprintf "%.0f", $position;  # due to some erroneous .5 in-between pos.
    $value      =   sprintf "%.4f", $value;

    if ($reference_name !~ /^\w\d?$/)             { next }
    if ($position       !~ /^\d{1,9}$/)           { next }
    if ($value          !~ /^\-?\d+?(\.\d+?)?$/)  { next }

    push(
      @{ $plot->{$probeT} },
      {
        no              =>  $i,
        probe_id        =>  $probe_id,
        reference_name  =>  $reference_name,
        position        =>  $position,
        value           =>  $numfactor * $value,
      }
    );
  }
  
  if ($probeT =~ /fracb/i) { return $plot }

  # random values
  if ($plot->{parameters}->{simulated_probes} =~ /y/i ) {
    my @randomV =   random_normal(scalar @{ $plot->{$probeT} }, 0, 0.25);
    foreach my $n (0..$#{ $plot->{$probeT} }) {
      $plot->{$probeT}->[$n]->{value}  =   $randomV[$n];
    }
  }

  # baseline adjustment
  if ($plot->{parameters}->{plot_adjust_baseline} =~ /[123456789]/) {
    foreach my $n (0..$#{ $plot->{$probeT} }) {
      $plot->{$probeT}->[$n]->{value}  +=   $plot->{parameters}->{plot_adjust_baseline};
    }
  }

  return $plot;

}

################################################################################

sub read_segmentfile {

=pod

Expects:
  - a standard Progenetix segments  file

  sample  chro  start stop  mean  probes
  GSM481286 1 742429  7883881 -0.1594 699
  GSM481286 1 115673158 115705254 -0.3829 8
  GSM481286 1 115722621 119771659 0.167 424
  GSM481286 1 119776776 162617092 0.4168  1587
  GSM481286 1 162621657 165278686 0.6508  350
  GSM481286 1 165280711 167221337 0.4056  241
  GSM481286 1 167248788 168289603 0.6784  130
  ...

Returns:
  - a list reference of genome CNV objects:
    [
      {
        no              =>  __integer__,    # 1 -> n
        callset_id      =>  __string__,
        reference_name  =>  __string__,
        start           =>  __integer__,
        end             =>  __integer__,
        variant_type    =>  __string__,     # DUP, DEL
        info            =>  {
          value           =>  __long__,
          svlen           =>  __integer__,
          probes          =>  __integer__,
          assembly_id     =>  __string__,     # GRCh36 ...
          experiment_type =>  __string__,     # aCGH ...
        },
      },
      {
      ...
      },
    ]

=cut

########    ####    ####    ####    ####    ####    ####    ####    ####    ####

  my $plot      =   shift;
  my $segmentsF =   shift;
  my $segmentsT =   shift;
  $segmentsT    ||= 'segmentdata';
  $plot->{$segmentsT}  =   [];

  if (! -f $segmentsF) { return $plot->{$segmentsT} }

  my $numfactor =   1;
  if (
    $plot->{parameters}->{'reverse'} =~ /y/i
    &&
    $segmentsT !~ /frac/i
  ) { $numfactor = -1 }

  my %colOrder  =   (
    callset_id          =>  0,
    reference_name      =>  1,
    start               =>  2,
    end                 =>  3,
    value               =>  4,
    probes              =>  5,
  );

  if ($plot->{parameters}->{format_inputfiles} =~ /tcga/i) {
    $colOrder{value}    =   5;
    $colOrder{probes}   =   4;
  };

  open  FILE, "$segmentsF" or die "No file $segmentsF $!";
  my $i         =   0;
  foreach (<FILE>) {
    $i++;
    chomp;
    my @segment =   split (/\s/, $_, 6);
    my %segVals =   ();
    foreach (keys %colOrder) {
      $segVals{$_}  =   $segment[$colOrder{$_}];
    };

    $segVals{callset_id}        =~  s/[^\w\-\,]/_/g;
    $segVals{reference_name}    =~ s/[^\dxXyY]//;
    $segVals{reference_name}    =~ s/^23$/X/;
    $segVals{reference_name}    =~ s/^24$/Y/;
    $segVals{start}     =   sprintf "%.0f", $segVals{start};
    $segVals{end}       =   sprintf "%.0f", $segVals{end};
    $segVals{probes}    =~  s/[^\d]//g;
    $segVals{value}     =   sprintf "%.4f", $segVals{value};

    if ($segVals{reference_name}!~ /^\w\d?$/)             { next }
    if ($segVals{start}         !~ /^\d{1,9}$/)           { next }
    if ($segVals{end}           !~ /^\d{1,9}$/)           { next }
    if ($segVals{value}         !~ /^\-?\d+?(\.\d+?)?$/)  { next }
    if (
      $segVals{probes} =~ /^\d\d?$/
      &&
      $segVals{probes} < $plot->{parameters}->{segment_probecount_min}
    )                                                     { next }

    push(
      @{ $plot->{$segmentsT} },
      {
        no              =>  $i,
        callset_id      =>  $segVals{callset_id},
        reference_name  =>  $segVals{reference_name},
        start           =>  1 * $segVals{start},
        end             =>  1 * $segVals{end},
        info            =>  {
          value         =>  $numfactor * $segVals{value},
          svlen         =>  1 * ($segVals{end} - $segVals{start}),
          probes        =>  1 * $segVals{probes},
        },
      }
    );

  }

  if ($segmentsT =~ /fracb/i) { return $plot }

  # baseline adjustment
  if ($plot->{parameters}->{plot_adjust_baseline} =~ /[123456789]/) {
    foreach my $i (0..$#{ $plot->{$segmentsT} }) {
      $plot->{$segmentsT}->[$i]->{info}->{value}  +=   $plot->{parameters}->{plot_adjust_baseline};
    }
  }

  for my $i (0..$#{ $plot->{$segmentsT} }) {

    if ($plot->{$segmentsT}->[$i]->{info}->{value} >= $plot->{parameters}->{cna_gain_threshold}) {
      $plot->{$segmentsT}->[$i]->{variant_type}  =   'DUP' }
    elsif ($plot->{$segmentsT}->[$i]->{info}->{value} <= $plot->{parameters}->{cna_loss_threshold}) {
      $plot->{$segmentsT}->[$i]->{variant_type}  =   'DEL' }

  }

  return $plot;

}

1;
