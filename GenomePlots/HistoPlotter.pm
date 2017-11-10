package PGX::GenomePlots::HistoPlotter;

use Data::Dumper;
use PGX::GenomePlots::CytobandsPlotter;

require Exporter;
@ISA        =   qw(Exporter);
@EXPORT     =   qw(return_histoplot_svg get_histoplot_area);

################################################################################

sub return_histoplot_svg {

  our $plot     =   shift;

  $plot->{Y}    =   $plot->{parameters}->{size_plotmargin_top_px};
  my $plotW             =   $plot->{parameters}->{size_plotimage_w_px};
  $plot->{areawidth}    =   $plotW - 2 * $plot->{parameters}->{size_plotmargin_px};
  if (
    $plot->{parameters}->{do_chromosomes_proportional} =~ /y/i
    &&
    @{ $plot->{parameters}->{chr2plot} } == 1
  ) {
    $plot->{areawidth}  *=  ($plot->{referencebounds}->{ $plot->{parameters}->{chr2plot}->[0] }->[1] / $plot->{referencebounds}->{ '1' }->[1]);
    $plotW      =   $plot->{areawidth} + 2 * $plot->{parameters}->{size_plotmargin_px};
  }
  $plot->{basepixfrac}  =   ( $plot->{areawidth} - ($#{ $plot->{parameters}->{chr2plot} } * $plot->{parameters}->{size_chromosome_padding_px}) ) / $plot->{genomesize};
  $plot         =   get_title_svg($plot);
  $plot         =   get_cytobands_svg($plot);
  $plot         =   get_histoplot_area($plot);
  $plot         =   get_bottom_labels_svg($plot);
  $plot->{Y}    +=   $plot->{parameters}->{size_plotmargin_bottom_px};
  my $plotH     =   sprintf "%.0f", $plot->{Y};
  $plotW        =   sprintf "%.0f", $plotW;
  $plot->{svg}  =   '<svg
xmlns="http://www.w3.org/2000/svg"
xmlns:xlink="http://www.w3.org/1999/xlink"
version="1.1"
id="'.$plot->{plotid}.'"
width="'.$plotW.'px"
height="'.$plotH.'px"
style="margin: auto; font-family: Helvetica, sans-serif;">

<rect x="0" y="0" width="'.$plotW.'" height="'.$plotH.'" style="fill: '.$plot->{parameters}->{color_plotbackground_hex}.'; " />

'.$plot->{svg}.'
</svg>';

  return $plot;

}

################################################################################

sub get_histoplot_area {

=pod

Expects:
  - the current Y parameter for placing the plot elements on the SVG

Returns:
  - the increased end Y value as start for the next elements

=cut

  ######    ####    ####    ####    ####    ####    ####    ####    ####    ####

  get_labels_y_svg($plot);

  my $area_x0   =   $plot->{parameters}->{size_plotmargin_px};
  my $area_y0   =   $plot->{Y};
  my $area_yn   =   $plot->{Y} + $plot->{parameters}->{size_plotarea_h_px};
  my $area_ycen =   $plot->{Y} + $plot->{parameters}->{size_plotarea_h_px} / 2;

  foreach my $refName (@{ $plot->{parameters}->{chr2plot} }) {

    my $areaW   =  sprintf "%.1f", ($plot->{referencebounds}->{$refName}->[1] - $plot->{referencebounds}->{$refName}->[0]) * $plot->{basepixfrac};

      $plot->{svg}      .=  '
<rect x="'.$area_x0.'" y="'.$plot->{Y}.'" width="'.$areaW.'" height="'.$plot->{parameters}->{size_plotarea_h_px}.'" style="fill: '.$plot->{parameters}->{color_plotarea_hex}.'; fill-opacity: 0.8; " />';

     # intervals through index #    ####    ####    ####    ####    ####    ####
    my @ind     =  grep{
          $refName eq	$plot->{genomeintervals}->[$_]->{reference_name} } 0..$#{ $plot->{genomeintervals} };
    @ind        =  grep{
          $plot->{genomeintervals}->[$_]->{start} <=	$plot->{referencebounds}->{$refName}->[1]  } @ind;
    @ind        =  grep{
          $plot->{genomeintervals}->[$_]->{end} >=	$plot->{referencebounds}->{$refName}->[0]  } @ind;

    foreach my $GL (qw(dupfrequencies delfrequencies)) {
      $plot->{svg}	    .=	  '
<polygon points="'.$area_x0.' '.$area_ycen;

      foreach my $i (@ind) {

        my $start       =   $plot->{genomeintervals}->[$i]->{start};
        my $end         =   $plot->{genomeintervals}->[$i]->{end};
        if ($start < $plot->{referencebounds}->{$refName}->[0]) {
        $start  = $plot->{referencebounds}->{$refName}->[0] }
        if ($end > $plot->{referencebounds}->{$refName}->[1]) {
        $end   = $plot->{referencebounds}->{$refName}->[1] }

        my $X   =		sprintf "%.1f", $area_x0 + $plot->{basepixfrac} * ($end - ($end - $start) / 2);
        my $H	  =		sprintf "%.1f", $plot->{frequencymaps}->{$GL}->[$i] * $plot->{parameters}->{pixyfactor};
        $plot->{svg}	  .=	' '.$X.' '.(sprintf "%.1f", ( $GL eq 'delfrequencies' ? $area_ycen + $H : $area_ycen - $H) );
      }
      $plot->{svg}	    .=  ' '.(sprintf "%.1f", $area_x0 + $areaW ).' '.$area_ycen.'" fill="'.($GL =~ /del/i ? $plot->{parameters}->{color_var_del_hex} : $plot->{parameters}->{color_var_dup_hex}).'" stroke-width="0px" />';

    }
    $area_x0		+=	$areaW + $plot->{parameters}->{size_chromosome_padding_px};
  }

  # adding a baseline at 0
  $plot->{svg}  .=  '
<line x1="'.$plot->{parameters}->{size_plotmargin_px}.'"  y1="'.$area_ycen.'"  x2="'.($plot->{parameters}->{size_plotmargin_px} + $plot->{areawidth}).'"  y2="'.$area_ycen.'"  class="cen"  />';

  # TODO: grid, Y-labels

  $plot->{Y}    +=  $plot->{parameters}->{size_plotarea_h_px};

  return $plot;

}

1;
