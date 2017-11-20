package PGX::GenomePlots::CytobandsPlotter;

use Data::Dumper;

require Exporter;
@ISA        =   qw(Exporter);
@EXPORT     =   qw(
  svg_add_title
  svg_add_cytobands
  svg_add_labels_y
  svg_add_markers
  svg_add_bottom_text
);

################################################################################

################################################################################

sub svg_add_title {

=pod

Expects:

Returns:

=cut

########    ####    ####    ####    ####    ####    ####    ####    ####    ####

  my $plot      =   shift;

  if ($plot->{parameters}->{title} !~ /\w+?/) { return $plot }

  my $titleSize =   sprintf "%.0f", $plot->{parameters}->{size_text_px} * 1.3;
  $plot->{Y}    +=   $titleSize;
  $plot->{svg}  .=  '
<text x="'.($plot->{areawidth} / 2 + $plot->{parameters}->{size_plotmargin_px}).'" y="'.$plot->{Y}.'" style="text-anchor: middle; font-size: '.$titleSize.'px">'.$plot->{parameters}->{title}.'</text>';

  $plot->{Y}    +=  $plot->{parameters}->{size_chromosome_padding_px};
  return $plot;


}

################################################################################

sub svg_add_labels_y {

=pod

Expects:

Returns:

=cut

########    ####    ####    ####    ####    ####    ####    ####    ####    ####

  my $plot      =   shift;

  if (@{ $plot->{parameters}->{label_y_m} } < 1) { return $plot }

  $plot->{svg}  .=  '
<style type="text/css"><![CDATA[
  .cen {stroke-width: '.$plot->{parameters}->{size_centerline_stroke_px}.'px; stroke: '.$plot->{parameters}->{color_plotgrid_hex}.'; opacity: 0.8 ; }
  .tick {stroke-width: 1px; stroke: '.$plot->{parameters}->{color_label_y_hex}.'; opacity: 0.8 ; }
  .ylabs {text-anchor: end; font-size: '.$plot->{parameters}->{size_text_lab_px}.'px; fill: '.$plot->{parameters}->{color_label_y_hex}.';}
  .ylabe {text-anchor: start; font-size: '.$plot->{parameters}->{size_text_lab_px}.'px; fill: '.$plot->{parameters}->{color_label_y_hex}.';}
]]></style>';

  my $area_x0   =   $plot->{parameters}->{size_plotmargin_px};
  my $area_y0   =   $plot->{Y};
  my $area_ycen =   $plot->{Y} + $plot->{parameters}->{size_plotarea_h_px} / 2;
  my $area_yn   =   $plot->{Y} + $plot->{parameters}->{size_plotarea_h_px};

  foreach my $lab (@{ $plot->{parameters}->{label_y_m} }) {

    my $lab_y   =   sprintf "%.1f", $area_ycen - $lab * $plot->{parameters}->{pixyfactor};

    # checking area boundaries
    if ($lab_y < $area_y0 || $lab_y > $area_yn) { next }

    $plot->{svg}        .=  '
<line x1="'.($area_x0 - 1).'"  y1="'.$lab_y.'"  x2="'.$area_x0.'"  y2="'.$lab_y.'"  class="tick"  />
<line x1="'.($area_x0 + $plot->{areawidth}).'"  y1="'.$lab_y.'"  x2="'.($area_x0 + $plot->{areawidth} + 1).'"  y2="'.$lab_y.'"  class="tick"  />
<line x1="'.$area_x0.'"  y1="'.$lab_y.'"  x2="'.($area_x0 + $plot->{areawidth}).'"  y2="'.$lab_y.'"  class="cen"  />';

    # avoiding too dense labels
    if (@{ $plot->{parameters}->{label_y_m} } > 9 && $lab !~ /^\-?\d\d?\d?\%?$/ ) { next }
    # positioning the label text
    $lab_y              +=  ($plot->{parameters}->{size_text_lab_px} / 2) - 1;
    $plot->{svg}        .=  '
<text x="'.($area_x0 - 2).'" y="'.$lab_y.'" class="ylabs">'.$lab.$plot->{parameters}->{plot_unit_y}.'</text>
<text x="'.($area_x0 + $plot->{areawidth} + 2).'" y="'.$lab_y.'" class="ylabe">'.$lab.$plot->{parameters}->{plot_unit_y}.'</text>';

  }

  return $plot;


}

################################################################################

sub svg_add_markers {

=pod

Expects:

Returns:

=cut

########    ####    ####    ####    ####    ####    ####    ####    ####    ####

  my $plot      =   shift;

  if (@{ $plot->{parameters}->{markers} } < 1) { return $plot }

  $plot->{Y}      +=  $plot->{parameters}->{size_chromosome_padding_px};

  $plot->{svg}  .=  '
<style type="text/css">
<![CDATA[
.marker { text-anchor: middle; font-size: '.$plot->{parameters}->{size_text_marker_px}.'px }
]]>
</style>';

  my $areaX_0   =   $plot->{areastartx};

  my $markerLineHeight  =   $plot->{parameters}->{size_text_marker_px} + 4;
  my %markerLineNo  =   (0 => $areaX_0);

  foreach my $refName (@{ $plot->{parameters}->{chr2plot} }) {

    my $areaW   =  sprintf "%.1f", ($plot->{referencebounds}->{$refName}->[1] - $plot->{referencebounds}->{$refName}->[0]) * $plot->{basepixfrac};

    my $areamarkers     =   [ grep{ $_->{reference_name} eq $refName } @{ $plot->{parameters}->{markers} } ];
    $areamarkers        =   [ grep{ $_->{start} <= $plot->{referencebounds}->{$refName}->[1] } @$areamarkers];
    $areamarkers        =   [ grep{ $_->{end} >= $plot->{referencebounds}->{$refName}->[0] } @$areamarkers ];
    %markerLineNo       =   (1 => $areaX_0);

    foreach my $marker (@$areamarkers) {
      if ($marker->{start} < $plot->{referencebounds}->{$refName}->[0]) {
        $marker->{start}    =   $plot->{referencebounds}->{$refName}->[0] }
      if ($marker->{end} > $plot->{referencebounds}->{$refName}->[1]) {
        $marker->{end}      =   $plot->{referencebounds}->{$refName}->[1] }

      my $mark_X0     =   sprintf "%.1f", $areaX_0 + $marker->{start} * $plot->{basepixfrac};
      my $mark_W      =   sprintf "%.2f", ($marker->{end} - $marker->{start}) * $plot->{basepixfrac};
      if ($mark_W < 0.5) {$mark_W = 0.5}
      my $mark_H      =   sprintf "%.0f", ($plot->{areaendy} - $plot->{areastarty});
      $plot->{svg}    .=  '
<rect x="'.$mark_X0.'" y="'.$plot->{areastarty}.'" width="'.$mark_W.'" height="'.$mark_H.'" style="fill: '.$marker->{color}.'; fill-opacity: 0.3; " />';


      if ($marker->{label} =~ /\w/) {
      
        my $marklablen  =   sprintf "%.0f", length($marker->{label}) * $plot->{parameters}->{size_text_marker_px} * 0.7 + 4;
        my $marklab_Xcen =   sprintf "%.1f", $mark_X0 + $mark_W / 2;
        my $marklab_X0  =   $marklab_Xcen - $marklablen / 2;
        my $marklab_Xn  =   $marklab_X0 + $marklablen;
        my $markbox_Y   =   $plot->{Y};
        my $marklab_Y   =   $markbox_Y + $markerLineHeight - 3;
        foreach my $i (1..100) {
          if ($markerLineNo{$i} > $marklab_X0) {
            $marklab_Y  +=  $markerLineHeight + 1;
            $markbox_Y  +=  $markerLineHeight + 1;
          } else {
            $markerLineNo{$i} =   $marklab_Xn;
            last;
          }
        }
        $plot->{svg}    .=  '
<rect x="'.$marklab_X0.'" y="'.$markbox_Y.'" width="'.$marklablen.'" height="'.$markerLineHeight.'" style="fill: '.$marker->{color}.'; fill-opacity: 0.3; " />
<text x="'.$marklab_Xcen.'" y="'.$marklab_Y.'" class="marker">'.$marker->{label}.'</text>';

      }

    }
    $areaX_0    +=  $areaW + $plot->{parameters}->{size_chromosome_padding_px};
  }
  
  my $maxline   =   (sort { $a <=> $b } keys %markerLineNo)[-1];
  
  $plot->{Y}    +=  $maxline * ($markerLineHeight + 1) - 1;

  return $plot;

}

################################################################################

sub svg_add_bottom_text {

  use Time::Piece;

=pod

Expects:

Returns:

=cut

########    ####    ####    ####    ####    ####    ####    ####    ####    ####

  my $plot      =   shift;

  $plot{Y}      +=  $plot->{parameters}->{size_chromosome_padding_px};

  $plot->{parameters}->{text_bottom_right} ||= '&#x24B8; '.(localtime->strftime('%Y')).' progenetix.org';
  if (
    $plot->{parameters}->{text_bottom_left} !~ /\w+?/
    &&
    $plot->{parameters}->{text_bottom_right} !~ /\w+?/
  ) { return $plot }

  $plot->{Y}    +=   $plot->{parameters}->{size_text_px};
  $plot->{svg}  .=  '
<text x="'.$plot->{parameters}->{size_plotmargin_px}.'" y="'.$plot->{Y}.'" style="text-anchor: start; font-size: '.$plot->{parameters}->{size_text_px}.'px; fill: '.$plot->{parameters}->{color_label_bottom_hex}.'">'.$plot->{parameters}->{text_bottom_left}.'</text>
<text x="'.($plot->{parameters}->{size_plotmargin_px} + $plot->{areawidth}).'" y="'.$plot->{Y}.'" style="text-anchor: end; font-size: '.$plot->{parameters}->{size_text_px}.'px; fill: '.$plot->{parameters}->{color_label_bottom_hex}.'">'.$plot->{parameters}->{text_bottom_right}.'</text>';

  $plot->{Y}    +=  $plot->{parameters}->{size_chromosome_padding_px};
  return $plot;

}

################################################################################

sub svg_add_cytobands {

=pod

Expects:

Returns:

=cut

########    ####    ####    ####    ####    ####    ####    ####    ####    ####

  my $plot      =   shift;

  # one can skip cytoband plotting by setting the "size_chromosome_w_px"
  # parameter < 1
  if ($plot->{parameters}->{size_chromosome_w_px} < 1) { return $plot }

  my $areaX_0   =   $plot->{parameters}->{size_plotmargin_px};

  $plot->{svg}  .=  '
<defs>';
  $plot->{svg}  .=  _cytoband_svg_gradients($plot->{plotid});
  $plot->{svg}  .=  '
<style type="text/css">
<![CDATA[
.chrlab { text-anchor: middle; font-size: '.$plot->{parameters}->{size_text_px}.'px }
]]>
</style>';
  $plot->{svg}  .=  '
</defs>';

  my $chrolabY  =   $plot->{Y} + $plot->{parameters}->{size_text_px} + $plot->{parameters}->{size_chromosome_padding_px};
  my $chroBandY =   $chrolabY + $plot->{parameters}->{size_chromosome_padding_px};

  foreach my $refName (@{ $plot->{parameters}->{chr2plot} }) {

    my $areaW   =  sprintf "%.1f", ($plot->{referencebounds}->{$refName}->[1] - $plot->{referencebounds}->{$refName}->[0]) * $plot->{basepixfrac};
    my $chroX   =  sprintf "%.1f", $areaX_0 + $areaW / 2;

    my $areabands   =   [ grep{ $_->{reference_name} eq $refName } @{ $plot->{cytobands} } ];

    my $chrolabel   =   $refName;

    # adding the base boundaries if incomplete chromosome (from first and last
    # cytoband )
    if (
      $plot->{referencebounds}->{$refName}->[0] > $areabands->[0]->{start}
      ||
      $plot->{referencebounds}->{$refName}->[1] < $areabands->[-1]->{end}
    ) {
      $chrolabel        .=  ' ('.$plot->{referencebounds}->{$refName}->[0].'-'.$plot->{referencebounds}->{$refName}->[1].')';
    }
    $plot->{svg}        .=  '
<text x="'.$chroX.'" y="'.$chrolabY.'" class="chrlab">'.$chrolabel.'</text>';

    $areabands  =   [ grep{ $_->{start} <= $plot->{referencebounds}->{$refName}->[1] } @$areabands];
    $areabands  =   [ grep{ $_->{end} >= $plot->{referencebounds}->{$refName}->[0] } @$areabands ];

    foreach my $cb (@$areabands) {

      if ($cb->{start} < $plot->{referencebounds}->{$refName}->[0]) {
        $cb->{start}    =   $plot->{referencebounds}->{$refName}->[0] }
      if ($cb->{end} > $plot->{referencebounds}->{$refName}->[1]) {
        $cb->{end}      =   $plot->{referencebounds}->{$refName}->[1] }

      my $cbX   =   sprintf "%.1f", $areaX_0 + $plot->{basepixfrac} * ($cb->{start} - $plot->{referencebounds}->{$refName}->[0]);
      my $cbW   =   sprintf "%.1f", $plot->{basepixfrac} * ($cb->{end} - $cb->{start});
      my $cbY   =   $chroBandY;
      my $cbH   =   $plot->{parameters}->{size_chromosome_w_px};

      # terminal and centromere bands are more narrow
      if (
        ($cb->{stain} =~ /cen/i)
        ||
        $cb >= ($plot->{referencebounds}->{$refName}->[1] - 500)
        ||
        $cb->{start} <= ($plot->{referencebounds}->{$refName}->[0] + 500) ) {
        $cbY    +=  (0.1*$plot->{parameters}->{size_chromosome_w_px});
        $cbH    -=  (0.2*$plot->{parameters}->{size_chromosome_w_px});
      }

      # acrocentric "stalk" bands are slim
      if ($cb =~ /stalk/i) {
        $cbY    +=  (0.3*$plot->{parameters}->{size_chromosome_w_px});
        $cbH    -=  (0.6*$plot->{parameters}->{size_chromosome_w_px});
      }

      $plot->{svg}      .=  '
<rect x="'.$cbX.'" y="'.$cbY.'" width="'.$cbW.'" height="'.$cbH.'" style="fill: url(#'.$plot->{plotid}.$cb->{stain}.'); " />';

    }
    $areaX_0    +=  $areaW + $plot->{parameters}->{size_chromosome_padding_px};
  }

  $plot->{Y}    =  $chroBandY + $plot->{parameters}->{size_chromosome_w_px};

  return $plot;

}

################################################################################

sub _cytoband_svg_gradients {

  my $id        =   $_[0];

  return '
<linearGradient id="'.$id.'gpos100" x1="0%" x2="0%" y1="0%" y2="80%" spreadMethod="reflect">
    <stop offset="0%" stop-color="rgb(39,39,39)" />
    <stop offset="100%" stop-color="rgb(0,0,0)" />
</linearGradient>
<linearGradient id="'.$id.'gpos75" x1="0%" x2="0%" y1="0%" y2="80%" spreadMethod="reflect">
    <stop offset="0%" stop-color="rgb(87,87,87)" />
    <stop offset="100%" stop-color="rgb(39,39,39)" />
</linearGradient>
<linearGradient id="'.$id.'gpos50" x1="0%" x2="0%" y1="0%" y2="80%" spreadMethod="reflect">
    <stop offset="0%" stop-color="rgb(196,196,196)" />
    <stop offset="100%" stop-color="rgb(111,111,111)" />
</linearGradient>
<linearGradient id="'.$id.'gpos25" x1="0%" x2="0%" y1="0%" y2="80%" spreadMethod="reflect">
    <stop offset="0%" stop-color="rgb(223,223,223)" />
    <stop offset="100%" stop-color="rgb(196,196,196)" />
</linearGradient>
<linearGradient id="'.$id.'gneg" x1="0%" x2="0%" y1="0%" y2="80%" spreadMethod="reflect">
    <stop offset="0%" stop-color="white" />
    <stop offset="100%" stop-color="rgb(223,223,223)" />
</linearGradient>
<linearGradient id="'.$id.'gvar" x1="0%" x2="0%" y1="0%" y2="80%" spreadMethod="reflect">
    <stop offset="0%" stop-color="rgb(196,196,196)" />
    <stop offset="100%" stop-color="rgb(111,111,111)" />
</linearGradient>
<linearGradient id="'.$id.'stalk" x1="0%" x2="0%" y1="0%" y2="80%" spreadMethod="reflect">
    <stop offset="0%" stop-color="rgb(39,39,39)" />
    <stop offset="100%" stop-color="rgb(0,0,0)" />
</linearGradient>
<linearGradient id="'.$id.'acen" x1="0%" x2="0%" y1="0%" y2="80%" spreadMethod="reflect">
    <stop offset="0%" stop-color="rgb(163,55,247)" />
    <stop offset="100%" stop-color="rgb(138,43,226)" />
</linearGradient>
';

}


1;
