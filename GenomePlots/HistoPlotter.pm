package PGX::GenomePlots::HistoPlotter;

use Data::Dumper;
use PGX::GenomePlots::CytobandsPlotter;

require Exporter;
@ISA        =   qw(Exporter);
@EXPORT     =   qw(return_histoplot_svg get_histoplot_area);

################################################################################

sub return_histoplot_svg {

  our ($plot, $fMaps)   =   @_;

  $plot->{Y}            =   $plot->{parameters}->{size_plotmargin_top_px};
  $plot->{areawidth}    =   $plot->{parameters}->{size_plotimage_w_px} - 2 * $plot->{parameters}->{size_plotmargin_px};
  $plot->{basepixfrac}  =   ( $plot->{areawidth} - ($#{ $plot->{parameters}->{chr2plot} } * $plot->{parameters}->{size_chromosome_padding_px}) ) / $plot->{genomesize};
  $plot         =   get_cytobands_svg($plot);
  $plot         =   get_histoplot_area($plot);
  $plot->{Y}    +=   $plot->{parameters}->{size_plotmargin_bottom_px};
  $plot->{svg}  =   '<svg
xmlns="http://www.w3.org/2000/svg"
xmlns:xlink="http://www.w3.org/1999/xlink"
version="1.1"
id="'.$plot->{plotid}.'"
onload="init(evt)"
width="'.$plot->{parameters}->{size_plotimage_w_px}.'px"
height="'.$plot->{Y}.'px"
style="margin: auto; font-family: Helvetica, sans-serif;">
'.$plot->{svg}.'
</svg>';

  return $plot;

}

sub get_histoplot_area {

=pod

Expects:
  - the current Y parameter for placing the plot elements on the SVG

Returns:
  - the increased end Y value as start for the next elements

=cut

########    ####    ####    ####    ####    ####    ####    ####    ####    ####

  my $areaX_0   =   $plot->{parameters}->{size_plotmargin_px};
  my $areaY_0   =   $plot->{Y} + $plot->{parameters}->{size_plotarea_h_px} / 2;

  foreach my $refName (@{ $plot->{parameters}->{chr2plot} }) {

    my $areaW   =  sprintf "%.1f", $plot->{referencebounds}->{$refName}->[1] * $plot->{basepixfrac};

      $plot->{svg}      .=  '
<rect x="'.$areaX_0.'" y="'.$plot->{Y}.'" width="'.$areaW.'" height="'.$plot->{parameters}->{size_plotarea_h_px}.'" style="fill: '.$plot->{parameters}->{color_plotarea_hex}.'; fill-opacity: 0.8; " />';

    my @ind     =  grep{
          $refName eq	$plot->{genomeintervals}->[ $_ ]->{reference_name} } 0..$#{ $plot->{genomeintervals} };
  #print Dumper($refName, join(',', @$plotInd));
    foreach my $GL (qw(dupfrequencies delfrequencies)) {
      $plot->{svg}	    .=	  '
<polygon points="'.$areaX_0.' '.$areaY_0;
      foreach my $i (@ind) {
        my $X   =		sprintf "%.1f", $areaX_0 + $plot->{basepixfrac} * ($plot->{genomeintervals}->[$i]->{end} - ($plot->{genomeintervals}->[$i]->{end} - $plot->{genomeintervals}->[$i]->{start}) / 2);
        my $H	  =		sprintf "%.1f", $fMaps->{$GL}->[$i] / 2 * $plot->{parameters}->{size_plotarea_h_px} / $plot->{parameters}->{value_plot_y_max};
        $plot->{svg}	  .=	' '.$X.' '.(sprintf "%.1f", ( $GL eq 'delfrequencies' ? $areaY_0 + $H : $areaY_0 - $H) );
      }
      $plot->{svg}	    .=  ' '.(sprintf "%.1f", $areaX_0 + $areaW ).' '.$areaY_0.'" fill="'.($GL =~ /del/i ? $plot->{parameters}->{color_var_del_hex} : $plot->{parameters}->{color_var_dup_hex}).'" stroke-width="0px" />';

    }
    $areaX_0		+=	$areaW + $plot->{parameters}->{size_chromosome_padding_px};
  }

  # TODO: grid, Y-labels

  $plot->{Y}    +=  $plot->{parameters}->{size_plotarea_h_px};

  return $plot;

}

1;
