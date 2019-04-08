package PGX::GenomePlots::HistoPlotter;

use Data::Dumper;
use PGX::GenomePlots::CytobandsPlotter;
use PGX::GenomePlots::StripPlotter;

require Exporter;
@ISA            =   qw(Exporter);
@EXPORT         =   qw(return_histoplot_svg get_histoplot_area);

################################################################################

sub return_histoplot_svg {

  my $pgx       =   shift;

  $pgx->{svg}   =   q{};

  $pgx->{Y}     =   $pgx->{parameters}->{size_plotmargin_top_px};
  my $plotW             =   $pgx->{parameters}->{size_plotimage_w_px};
  $pgx->{areastartx}    =   $pgx->{parameters}->{size_plotmargin_px} + $pgx->{parameters}->{size_title_left_px};
  if (
    $pgx->{parameters}->{size_clustertree_w_px} < 1
     ||
     (scalar @{ $pgx->{clustertree} } < 2)
  ) { $pgx->{parameters}->{size_clustertree_w_px} = 0 }
  $pgx->{areawidth} =   $plotW - ($pgx->{areastartx} + $pgx->{parameters}->{size_plotmargin_px} + $pgx->{parameters}->{size_label_right_px} + $pgx->{parameters}->{size_clustertree_w_px});
;
  if (
    $pgx->{parameters}->{do_chromosomes_proportional} =~ /y/i
    &&
    @{ $pgx->{parameters}->{chr2plot} } == 1
  ) {
    $pgx->{areawidth}   *=  ($pgx->{referencebounds}->{ $pgx->{parameters}->{chr2plot}->[0] }->[1] / $pgx->{referencebounds}->{ '1' }->[1]);
    $plotW      =   $pgx->{areawidth} + 2 * $pgx->{parameters}->{size_plotmargin_px};
  }
  $pgx->{basepixfrac}   =   ( $pgx->{areawidth} - ($#{ $pgx->{parameters}->{chr2plot} } * $pgx->{parameters}->{size_chromosome_padding_px}) ) / $pgx->{genomesize};
  $plotW        +=  $pgx->{parameters}->{size_label_right_px};
  $pgx->{areaendx}  =   $pgx->{areastartx} + $pgx->{areawidth};
  $pgx->{areatreex} =   $pgx->{areaendx};
  if ($pgx->{parameters}->{size_label_right_px} > 0) {
    $pgx->{areatreex}   +=  $pgx->{parameters}->{size_chromosome_padding_px} + $pgx->{parameters}->{size_label_right_px} }
  $pgx->svg_add_title();
  $pgx->svg_add_cytobands();
  $pgx->get_histoplot_area();
  $pgx->svg_add_cluster_tree();
  $pgx->{markerstarty}  =   $pgx->{areastarty}; # so that markers span the histograms
  if (
    @{ $pgx->{frequencymaps} } > 1
    &&
    $pgx->{parameters}->{size_strip_h_px} > 0
  ) {
    $pgx->svg_add_cytobands();
    $pgx->get_frequencystripplot_area_gd();
    $pgx->svg_add_cluster_tree();
  }
  $pgx->svg_add_markers();
  $pgx->svg_add_bottom_text();
  $pgx->{Y}     +=   $pgx->{parameters}->{size_plotmargin_bottom_px};
  my $plotH     =   sprintf "%.0f", $pgx->{Y};
  $plotW        =   sprintf "%.0f", $plotW;
  $pgx->{svg}   =   '<svg
xmlns="http://www.w3.org/2000/svg"
xmlns:xlink="http://www.w3.org/1999/xlink"
version="1.1"
id="'.$pgx->{plotid}.'"
width="'.$plotW.'px"
height="'.$plotH.'px"
style="margin: auto; font-family: Helvetica, sans-serif;">

<rect x="0" y="0" width="'.$plotW.'" height="'.$plotH.'" style="fill: '.$pgx->{parameters}->{color_plotbackground_hex}.'; " />

'.$pgx->{svg}.'
</svg>';

  return $pgx;

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

  my $pgx       =   shift;

  if ($pgx->{parameters}->{size_plotarea_h_px} < 1) { return $pgx }

  # preview of the first histogram area Y start (for use in cluster tree etc.)
  $pgx->{areastarty}   =   $pgx->{Y} + $pgx->{parameters}->{size_plotarea_padding};

  my $defLabCol =   '#dddddd';
  my $altLabCol =   '#fefefe';
  my $labCol    =   '#dddddd';
  
  foreach my $frequencymapsSet (@{ $pgx->{frequencymaps} }) {
    $pgx->{Y}  +=  $pgx->{parameters}->{size_plotarea_padding};

    if ($frequencymapsSet->{name} =~ /\w\w/) {
      $pgx->{parameters}->{title_left} = $frequencymapsSet->{name} }

    my $area_x0     =   $pgx->{areastartx};
    my $area_ycen   =   $pgx->{Y} + $pgx->{parameters}->{size_plotarea_h_px} / 2;

    $pgx->svg_add_title_left($area_ycen);
    $pgx->svg_add_labels_y();

    foreach my $refName (@{ $pgx->{parameters}->{chr2plot} }) {

      my $areaW =  sprintf "%.1f", ($pgx->{referencebounds}->{$refName}->[1] - $pgx->{referencebounds}->{$refName}->[0]) * $pgx->{basepixfrac};

      $pgx->{svg}  .=  '
  <rect x="'.$area_x0.'" y="'.$pgx->{Y}.'" width="'.$areaW.'" height="'.$pgx->{parameters}->{size_plotarea_h_px}.'" style="fill: '.$pgx->{parameters}->{color_plotarea_hex}.'; fill-opacity: 0.8; " />';

       # intervals through index #    ####    ####    ####    ####    ####    ####
      my @ind   =  grep{ $refName eq $pgx->{genomeintervals}->[$_]->{reference_name} } 0..$#{ $pgx->{genomeintervals} };
      @ind      =  grep{ $pgx->{genomeintervals}->[$_]->{start} <=  $pgx->{referencebounds}->{$refName}->[1]  } @ind;
      @ind      =  grep{ $pgx->{genomeintervals}->[$_]->{end} >=  $pgx->{referencebounds}->{$refName}->[0]  } @ind;

      foreach my $GL (qw(dupfrequencies delfrequencies)) {
        $pgx->{svg}      .=    '
  <polygon points="'.$area_x0.' '.$area_ycen;

        foreach my $i (@ind) {
        
          my $x_corr    =   $pgx->{referencebounds}->{$refName}->[0] * $pgx->{basepixfrac};

          my $start =   $pgx->{genomeintervals}->[$i]->{start};
          my $end   =   $pgx->{genomeintervals}->[$i]->{end};
          if ($start < $pgx->{referencebounds}->{$refName}->[0]) {
            $start  =   $pgx->{referencebounds}->{$refName}->[0] }
          if ($end > $pgx->{referencebounds}->{$refName}->[1]) {
            $end    =   $pgx->{referencebounds}->{$refName}->[1] }

          my $X =   sprintf "%.1f", $area_x0 - $x_corr + $pgx->{basepixfrac} * ($end - ($end - $start) / 2);
          my $H =   sprintf "%.1f", $frequencymapsSet->{$GL}->[$i] * $pgx->{parameters}->{pixyfactor};
          $pgx->{svg}  .=  ' '.$X.' '.(sprintf "%.1f", ( $GL eq 'delfrequencies' ? $area_ycen + $H : $area_ycen - $H) );

        }

        $pgx->{svg}    .=  ' '.(sprintf "%.1f", $area_x0 + $areaW ).' '.$area_ycen.'" fill="'.($GL =~ /del/i ? $pgx->{parameters}->{color_var_del_hex} : $pgx->{parameters}->{color_var_dup_hex}).'" stroke-width="0px" />';

      }

      $area_x0  +=  $areaW + $pgx->{parameters}->{size_chromosome_padding_px};

    }

    # adding a baseline at 0
    $pgx->{svg}    .=  '
  <line x1="'.$pgx->{areastartx}.'"  y1="'.$area_ycen.'"  x2="'.($pgx->{areastartx} + $pgx->{areawidth}).'"  y2="'.$area_ycen.'"  class="cen"  />';

    my $labels_R    =   [];
    if ($frequencymapsSet->{labels}) { $labels_R = $frequencymapsSet->{labels} }
    # fallback color; defined here for alternation...
    if ($labCol eq $altLabCol) { 
      $labCol   =   $defLabCol }
    else { 
      $labCol   =   $altLabCol }
    $pgx->svg_add_labels_right($labels_R, $pgx->{parameters}->{size_plotarea_h_px}, $labCol);

    $pgx->{Y}   +=  $pgx->{parameters}->{size_plotarea_h_px};

  }

  return $pgx;

}

1;
