package PGX::GenomePlots::HistoPlotter;

use Data::Dumper;
use PGX::GenomePlots::CytobandsPlotter;
use PGX::GenomePlots::StripPlotter;

require Exporter;
@ISA        =   qw(Exporter);
@EXPORT     =   qw(return_histoplot_svg get_histoplot_area);

################################################################################

sub return_histoplot_svg {

  my $plot      =   shift;

  $plot->{Y}    =   $plot->{parameters}->{size_plotmargin_top_px};
  my $plotW             =   $plot->{parameters}->{size_plotimage_w_px};
  $plot->{areastartx}   =   $plot->{parameters}->{size_plotmargin_px} + $plot->{parameters}->{size_title_left_px};
  $plot->{areawidth}    =   $plotW - ($plot->{areastartx} + $plot->{parameters}->{size_plotmargin_px});
  if (
    $plot->{parameters}->{do_chromosomes_proportional} =~ /y/i
    &&
    @{ $plot->{parameters}->{chr2plot} } == 1
  ) {
    $plot->{areawidth}  *=  ($plot->{referencebounds}->{ $plot->{parameters}->{chr2plot}->[0] }->[1] / $plot->{referencebounds}->{ '1' }->[1]);
    $plotW      =   $plot->{areawidth} + 2 * $plot->{parameters}->{size_plotmargin_px};
  }
  $plot->{basepixfrac}  =   ( $plot->{areawidth} - ($#{ $plot->{parameters}->{chr2plot} } * $plot->{parameters}->{size_chromosome_padding_px}) ) / $plot->{genomesize};
  $plotW        +=  $plot->{parameters}->{size_clustertree_w_px};
  $plotW        +=  $plot->{parameters}->{size_label_right_px};
  $plot->{areaendx}     =   $plot->{areastartx} + $plot->{areawidth};
  $plot->{areatreex}    =   $plot->{areaendx};
  if ($plot->{parameters}->{size_label_right_px} > 0) {
    $plot->{areatreex}  +=  $plot->{parameters}->{size_chromosome_padding_px} + $plot->{parameters}->{size_label_right_px} }
  $plot->svg_add_title();
  $plot->svg_add_cytobands();
  $plot->get_histoplot_area();
  $plot->svg_add_cluster_tree();
  if (
    @{ $plot->{frequencymaps} } > 1
    &&
    $plot->{parameters}->{size_strip_h_px} > 0
  ) {
    $plot->svg_add_cytobands();
    $plot->get_frequencystripplot_area();
    $plot->svg_add_cluster_tree();
  }
  $plot->svg_add_markers();
  $plot->svg_add_bottom_text();
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

  my $plot      =   shift;

  if ($plot->{parameters}->{size_plotarea_h_px} < 1) { return $plot }

  # preview of the first histogram area Y start (for use in cluster tree etc.)
  $plot->{areastarty}   =   $plot->{Y} + $plot->{parameters}->{size_plotarea_padding};

  my $defLabCol =   '#dddddd';
  my $altLabCol =   '#fefefe';
  my $labCol    =   '#dddddd';
  
  foreach my $frequencymapsSet (@{ $plot->{frequencymaps} }) {
    $plot->{Y}  +=  $plot->{parameters}->{size_plotarea_padding};

    if ($frequencymapsSet->{name} =~ /\w\w/) {
      $plot->{parameters}->{title_left} = $frequencymapsSet->{name} }

    my $area_x0     =   $plot->{areastartx};
    my $area_ycen   =   $plot->{Y} + $plot->{parameters}->{size_plotarea_h_px} / 2;

    $plot->svg_add_title_left($area_ycen);
    $plot->svg_add_labels_y();

    foreach my $refName (@{ $plot->{parameters}->{chr2plot} }) {

      my $areaW =  sprintf "%.1f", ($plot->{referencebounds}->{$refName}->[1] - $plot->{referencebounds}->{$refName}->[0]) * $plot->{basepixfrac};

      $plot->{svg}  .=  '
  <rect x="'.$area_x0.'" y="'.$plot->{Y}.'" width="'.$areaW.'" height="'.$plot->{parameters}->{size_plotarea_h_px}.'" style="fill: '.$plot->{parameters}->{color_plotarea_hex}.'; fill-opacity: 0.8; " />';

       # intervals through index #    ####    ####    ####    ####    ####    ####
      my @ind   =  grep{ $refName eq $plot->{genomeintervals}->[$_]->{reference_name} } 0..$#{ $plot->{genomeintervals} };
      @ind      =  grep{ $plot->{genomeintervals}->[$_]->{start} <=  $plot->{referencebounds}->{$refName}->[1]  } @ind;
      @ind      =  grep{ $plot->{genomeintervals}->[$_]->{end} >=  $plot->{referencebounds}->{$refName}->[0]  } @ind;

      foreach my $GL (qw(dupfrequencies delfrequencies)) {
        $plot->{svg}      .=    '
  <polygon points="'.$area_x0.' '.$area_ycen;

        foreach my $i (@ind) {

          my $start =   $plot->{genomeintervals}->[$i]->{start};
          my $end   =   $plot->{genomeintervals}->[$i]->{end};
          if ($start < $plot->{referencebounds}->{$refName}->[0]) {
            $start  =   $plot->{referencebounds}->{$refName}->[0] }
          if ($end > $plot->{referencebounds}->{$refName}->[1]) {
            $end    =   $plot->{referencebounds}->{$refName}->[1] }

          my $X =   sprintf "%.1f", $area_x0 + $plot->{basepixfrac} * ($end - ($end - $start) / 2);
          my $H =   sprintf "%.1f", $frequencymapsSet->{$GL}->[$i] * $plot->{parameters}->{pixyfactor};
          $plot->{svg}  .=  ' '.$X.' '.(sprintf "%.1f", ( $GL eq 'delfrequencies' ? $area_ycen + $H : $area_ycen - $H) );

        }

        $plot->{svg}    .=  ' '.(sprintf "%.1f", $area_x0 + $areaW ).' '.$area_ycen.'" fill="'.($GL =~ /del/i ? $plot->{parameters}->{color_var_del_hex} : $plot->{parameters}->{color_var_dup_hex}).'" stroke-width="0px" />';

      }

      $area_x0  +=  $areaW + $plot->{parameters}->{size_chromosome_padding_px};

    }

    # adding a baseline at 0
    $plot->{svg}    .=  '
  <line x1="'.$plot->{areastartx}.'"  y1="'.$area_ycen.'"  x2="'.($plot->{areastartx} + $plot->{areawidth}).'"  y2="'.$area_ycen.'"  class="cen"  />';


    my $labels_R    =   [];
    if ($frequencymapsSet->{labels}) { $labels_R = $frequencymapsSet->{labels} }
    # fallback color; defined here for alternation...
    if ($labCol eq $altLabCol) { 
      $labCol =   $defLabCol }
    else { 
      $labCol =   $altLabCol }
    $plot->svg_add_labels_right($labels_R, $plot->{parameters}->{size_plotarea_h_px}, $labCol);

    $plot->{Y}  +=  $plot->{parameters}->{size_plotarea_h_px};

  }

  return $plot;

}

1;
