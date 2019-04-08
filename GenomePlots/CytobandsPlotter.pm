package PGX::GenomePlots::CytobandsPlotter;

use Data::Dumper;

require Exporter;
@ISA        =   qw(Exporter);
@EXPORT     =   qw(
  svg_add_title
  svg_add_cytobands
  svg_add_title_left
  svg_add_labels_y
  svg_add_markers
  svg_add_bottom_text
  svg_add_labels_right
  svg_add_cluster_tree
);

################################################################################

################################################################################

sub svg_add_title {

=pod

Expects:

Returns:

=cut

########    ####    ####    ####    ####    ####    ####    ####    ####    ####

  my $pgx       =   shift;

  if ($pgx->{parameters}->{title} !~ /\w+?/) { return $pgx }

  $pgx->{Y}    +=   $pgx->{parameters}->{size_text_title_px};
  $pgx->{svg}  .=  '
<text x="'.($pgx->{areastartx} + $pgx->{areawidth} / 2).'" y="'.$pgx->{Y}.'" style="text-anchor: middle; font-size: '.$pgx->{parameters}->{size_text_title_px}.'px">'.$pgx->{parameters}->{title}.'</text>';

  if ($pgx->{parameters}->{subtitle} !~ /\w+?/) { 
    $pgx->{Y}  +=  $pgx->{parameters}->{size_text_title_px}; 
    return $pgx;
  }

  $pgx->{Y}    +=   sprintf "%.0f", $pgx->{parameters}->{size_text_subtitle_px} * 1.5;
  $pgx->{svg}  .=  '
<text x="'.($pgx->{areastartx} + $pgx->{areawidth} / 2).'" y="'.$pgx->{Y}.'" style="text-anchor: middle; font-size: '.$pgx->{parameters}->{size_text_subtitle_px}.'px">'.$pgx->{parameters}->{subtitle}.'</text>';

  $pgx->{Y}    +=  $pgx->{parameters}->{size_text_title_px};
  
  return $pgx;

}

################################################################################
################################################################################
################################################################################

sub svg_add_cytobands {

=pod

Expects:

Returns:

=cut

########    ####    ####    ####    ####    ####    ####    ####    ####    ####

  my $pgx       =   shift;

  # one can skip cytoband plotting by setting the "size_chromosome_w_px"
  # parameter < 1
  if ($pgx->{parameters}->{size_chromosome_w_px} < 1) { return $pgx }

  my $areaX_0   =   $pgx->{areastartx};

  $pgx->{svg}  .=  '
<defs>';
  $pgx->{svg}  .=  _cytoband_svg_gradients($pgx->{plotid});
  $pgx->{svg}  .=  '
<style type="text/css">
<![CDATA[
.chrlab { text-anchor: middle; font-size: '.$pgx->{parameters}->{size_text_px}.'px; fill: '.$pgx->{parameters}->{color_text_hex}.' }
]]>
</style>';
  $pgx->{svg}  .=  '
</defs>';

  my $chrolabY  =   $pgx->{Y} + $pgx->{parameters}->{size_text_px} + $pgx->{parameters}->{size_plotarea_padding};
  my $chroBandY =   $chrolabY + $pgx->{parameters}->{size_plotarea_padding};

  foreach my $refName (@{ $pgx->{parameters}->{chr2plot} }) {

    my $areaW   =  sprintf "%.1f", ($pgx->{referencebounds}->{$refName}->[1] - $pgx->{referencebounds}->{$refName}->[0]) * $pgx->{basepixfrac};
    my $chroX   =  sprintf "%.1f", $areaX_0 + $areaW / 2;

    my $areabands   =   [ grep{ $_->{reference_name} eq $refName } @{ $pgx->{cytobands} } ];
    my $chrolabel   =   $refName;

    # adding the base boundaries if incomplete chromosome (from first and last
    # cytoband )
    if (
      $pgx->{referencebounds}->{$refName}->[0] > $areabands->[0]->{start}
      ||
      $pgx->{referencebounds}->{$refName}->[1] < $areabands->[-1]->{end}
    ) {
      $chrolabel        .=  ' ('.$pgx->{referencebounds}->{$refName}->[0].'-'.$pgx->{referencebounds}->{$refName}->[1].')';
    }

    $pgx->{svg}        .=  '
<text x="'.$chroX.'" y="'.$chrolabY.'" class="chrlab">'.$chrolabel.'</text>';

    $areabands  =   [ grep{ $_->{start} <= $pgx->{referencebounds}->{$refName}->[1] } @$areabands];
    $areabands  =   [ grep{ $_->{end} >= $pgx->{referencebounds}->{$refName}->[0] } @$areabands ];

    foreach my $cb (@$areabands) {

      if ($cb->{start} < $pgx->{referencebounds}->{$refName}->[0]) {
        $cb->{start}    =   $pgx->{referencebounds}->{$refName}->[0] }
      if ($cb->{end} > $pgx->{referencebounds}->{$refName}->[1]) {
        $cb->{end}      =   $pgx->{referencebounds}->{$refName}->[1] }

      my $cbX   =   sprintf "%.1f", $areaX_0 + $pgx->{basepixfrac} * ($cb->{start} - $pgx->{referencebounds}->{$refName}->[0]);
      my $cbW   =   sprintf "%.1f", $pgx->{basepixfrac} * ($cb->{end} - $cb->{start});
      my $cbY   =   $chroBandY;
      my $cbH   =   $pgx->{parameters}->{size_chromosome_w_px};

      # terminal and centromere bands are more narrow
      if (
        ($cb->{stain} =~ /cen/i)
        ||
        $cb->{end} >= ($pgx->{referencebounds}->{$refName}->[1] - 500)
        ||
        $cb->{start} <= ($pgx->{referencebounds}->{$refName}->[0] + 500) ) {
        $cbY    +=  (0.1*$pgx->{parameters}->{size_chromosome_w_px});
        $cbH    -=  (0.2*$pgx->{parameters}->{size_chromosome_w_px});
      }

      # acrocentric "stalk" bands are slim
      if ($cb =~ /stalk/i) {
        $cbY    +=  (0.3*$pgx->{parameters}->{size_chromosome_w_px});
        $cbH    -=  (0.6*$pgx->{parameters}->{size_chromosome_w_px});
      }

      $pgx->{svg}      .=  '
<rect x="'.$cbX.'" y="'.$cbY.'" width="'.$cbW.'" height="'.$cbH.'" style="fill: url(#'.$pgx->{plotid}.$cb->{stain}.'); " />';

    }
    $areaX_0    +=  $areaW + $pgx->{parameters}->{size_chromosome_padding_px};
  }

  $pgx->{Y}    =  $chroBandY + $pgx->{parameters}->{size_chromosome_w_px};

  return $pgx;

}

################################################################################
################################################################################
################################################################################

sub svg_add_title_left {

  my $pgx       =   shift;
  my $titleY    =   shift;

  if ($pgx->{parameters}->{title_left} !~ /\w+?/) { return $pgx }
  if ($pgx->{parameters}->{size_title_left_px} < 1) { return $pgx }

  my $plotareaHeight  =   $pgx->{parameters}->{size_plotarea_h_px};
  if ($plotareaHeight < $pgx->{parameters}->{size_strip_h_px}) {
    $plotareaHeight  =   $pgx->{parameters}->{size_strip_h_px} }

  my $text_x    =   $pgx->{areastartx} /2;
  my $text_y    =   $titleY + $pgx->{parameters}->{size_text_title_left_px} / 2;

  $pgx->{svg}  .=  '
<text x="'.$text_x.'" y="'.$text_y.'"  transform="rotate('.$pgx->{parameters}->{title_left_rotation}.','.$text_x.','.$text_y.')" style="text-anchor: middle; fill: '.$pgx->{parameters}->{color_text_hex}.'; font-size: '.$pgx->{parameters}->{size_text_title_left_px}.'px; ">'.$pgx->{parameters}->{title_left}.'</text>';

  return $pgx;

}

################################################################################
################################################################################
################################################################################################################################################################

sub svg_add_labels_y {

=pod

Expects:

Returns:

=cut

########    ####    ####    ####    ####    ####    ####    ####    ####    ####

  my $pgx       =   shift;

  if (@{ $pgx->{parameters}->{label_y_m} } < 1) { return $pgx }

  $pgx->{svg}  .=  '
<style type="text/css"><![CDATA[
  .cen {stroke-width: '.$pgx->{parameters}->{size_centerline_stroke_px}.'px; stroke: '.$pgx->{parameters}->{color_plotgrid_hex}.'; opacity: 0.8 ; }
  .tick {stroke-width: 1px; stroke: '.$pgx->{parameters}->{color_label_y_hex}.'; opacity: 0.8 ; }
  .ylabs {text-anchor: end; font-size: '.$pgx->{parameters}->{size_text_lab_px}.'px; fill: '.$pgx->{parameters}->{color_label_y_hex}.';}
  .ylabe {text-anchor: start; font-size: '.$pgx->{parameters}->{size_text_lab_px}.'px; fill: '.$pgx->{parameters}->{color_label_y_hex}.';}
]]></style>';

  my $area_x0   =   $pgx->{areastartx};
  my $area_y0   =   $pgx->{Y};
  my $area_ycen =   $pgx->{Y} + $pgx->{parameters}->{size_plotarea_h_px} / 2;
  my $area_yn   =   $pgx->{Y} + $pgx->{parameters}->{size_plotarea_h_px};

  foreach my $lab (@{ $pgx->{parameters}->{label_y_m} }) {

    my $lab_y   =   sprintf "%.1f", $area_ycen - $lab * $pgx->{parameters}->{pixyfactor};

    # checking area boundaries
    if ($lab_y < $area_y0 || $lab_y > $area_yn) { next }

    $pgx->{svg}        .=  '
<line x1="'.($area_x0 - 1).'"  y1="'.$lab_y.'"  x2="'.$area_x0.'"  y2="'.$lab_y.'"  class="tick"  />
<line x1="'.($area_x0 + $pgx->{areawidth}).'"  y1="'.$lab_y.'"  x2="'.($area_x0 + $pgx->{areawidth} + 1).'"  y2="'.$lab_y.'"  class="tick"  />
<line x1="'.$area_x0.'"  y1="'.$lab_y.'"  x2="'.($area_x0 + $pgx->{areawidth}).'"  y2="'.$lab_y.'"  class="cen"  />';

    # avoiding too dense labels
    if (@{ $pgx->{parameters}->{label_y_m} } > 9 && $lab !~ /^\-?\d\d?\d?\%?$/ ) { next }
    # positioning the label text
    $lab_y              +=  ($pgx->{parameters}->{size_text_lab_px} / 2) - 1;
    $pgx->{svg}        .=  '
<text x="'.($area_x0 - 2).'" y="'.$lab_y.'" class="ylabs">'.$lab.$pgx->{parameters}->{plot_unit_y}.'</text>';
    if ($pgx->{parameters}->{size_clustertree_w_px} < 1) {
      $pgx->{svg}      .=  '
<text x="'.($area_x0 + $pgx->{areawidth} + 2).'" y="'.$lab_y.'" class="ylabe">'.$lab.$pgx->{parameters}->{plot_unit_y}.'</text>';
    }
  }

  return $pgx;

}

################################################################################
################################################################################
################################################################################

sub svg_add_markers {

=pod

Expects:

Returns:

=cut

########    ####    ####    ####    ####    ####    ####    ####    ####    ####

  my $pgx       =   shift;

  if (@{ $pgx->{parameters}->{markers} } < 1) { return $pgx }

  $pgx->{Y}    +=  $pgx->{parameters}->{size_chromosome_padding_px};
  $pgx->{svg}  .=  '
<style type="text/css">
<![CDATA[
.marker { text-anchor: middle; font-size: '.$pgx->{parameters}->{marker_text_px}.'px }
]]>
</style>';

  my $areaX_0   =   $pgx->{areastartx};
  my $markerY_0 =   $pgx->{areastarty};
  if ($pgx->{markerstarty} > 0) {
    $markerY_0  =   $pgx->{markerstarty} }

  # stores the marker line index and the last X value there
  my %markerLineNo      =   (1 => $pgx->{areastartx});
  my $markerLineHeight  =   $pgx->{parameters}->{marker_text_px} + 4;

  foreach my $refName (@{ $pgx->{parameters}->{chr2plot} }) {

    my $areaW   =  sprintf "%.1f", ($pgx->{referencebounds}->{$refName}->[1] - $pgx->{referencebounds}->{$refName}->[0]) * $pgx->{basepixfrac};

    my $areamarkers =   [ grep{ $_->{reference_name} eq $refName } @{ $pgx->{parameters}->{markers} } ];
    $areamarkers    =   [ grep{ $_->{start} <= $pgx->{referencebounds}->{$refName}->[1] } @$areamarkers];
    $areamarkers    =   [ grep{ $_->{end} >= $pgx->{referencebounds}->{$refName}->[0] } @$areamarkers ];

    foreach my $marker (@$areamarkers) {

      if ($marker->{start} < $pgx->{referencebounds}->{$refName}->[0]) {
        $marker->{start}    =   $pgx->{referencebounds}->{$refName}->[0] }
      if ($marker->{end} > $pgx->{referencebounds}->{$refName}->[1]) {
        $marker->{end}      =   $pgx->{referencebounds}->{$refName}->[1] }

      my $mark_X0   =   sprintf "%.1f", $areaX_0 + ($marker->{start} - $pgx->{referencebounds}->{$refName}->[0]) * $pgx->{basepixfrac};
      my $mark_W    =   sprintf "%.2f", ($marker->{end} - $marker->{start}) * $pgx->{basepixfrac};
      if ($mark_W < 0.5) {$mark_W = 0.5}
      my $mark_H    =   sprintf "%.0f", ($pgx->{Y} - $pgx->{markerstarty});
      $pgx->{svg}  .=  '
<rect x="'.$mark_X0.'" y="'.$markerY_0.'" width="'.$mark_W.'" height="'.$mark_H.'" style="fill: '.$marker->{color}.'; fill-opacity: '.$pgx->{parameters}->{marker_opacity}.'; " />';

      # adding a text label ####################################################

      if ($marker->{label} =~ /\w/) {

        my $mLablen     =   sprintf "%.0f", length($marker->{label}) * $pgx->{parameters}->{marker_text_px} * 0.7 + 4;
        my $mLab_Xcen   =   sprintf "%.1f", $mark_X0 + $mark_W / 2;
        my $mLab_X0     =   $mLab_Xcen - $mLablen / 2;
        my $mLab_Xn     =   $mLab_X0 + $mLablen;
        my $markbox_Y   =   $pgx->{Y};
        my $mLab_Y      =   $markbox_Y + $markerLineHeight - 3;
        foreach my $i (1..100) {
          if ($markerLineNo{$i} > $mLab_X0) {
            $mLab_Y  +=  $markerLineHeight + 1;
            $markbox_Y  +=  $markerLineHeight + 1;
          } else {
            $markerLineNo{$i} =   $mLab_Xn;
            last;
          }
        }
        $pgx->{svg}    .=  '
<rect x="'.$mLab_X0.'" y="'.$markbox_Y.'" width="'.$mLablen.'" height="'.$markerLineHeight.'" style="fill: '.$marker->{color}.'; fill-opacity: '.$pgx->{parameters}->{marker_label_opacity}.'; " />
<text x="'.$mLab_Xcen.'" y="'.$mLab_Y.'" class="marker">'.$marker->{label}.'</text>';

      }

      # / text label ###########################################################
      
    }
    $areaX_0    +=  $areaW + $pgx->{parameters}->{size_chromosome_padding_px};
  }

  my $maxline   =   (sort { $a <=> $b } keys %markerLineNo)[-1];
  $pgx->{Y}     +=  $maxline * ($markerLineHeight + 1) - 1;

  return $pgx;

}

################################################################################
################################################################################
################################################################################

sub svg_add_bottom_text {

  use Time::Piece;

=pod

Expects:

Returns:

=cut

########    ####    ####    ####    ####    ####    ####    ####    ####    ####

  my $pgx       =   shift;

  if (
    $pgx->{parameters}->{text_bottom_left} !~ /\w+?/
    &&
    $pgx->{parameters}->{text_bottom_right} !~ /\w+?/
  ) { return $pgx }

  $pgx->{Y}    +=  $pgx->{parameters}->{size_chromosome_padding_px};

  if ($pgx->{parameters}->{text_bottom_right} =~ /copy/i) {
    $pgx->{parameters}->{text_bottom_right} = '&#x24B8; '.(localtime->strftime('%Y')).' progenetix.org'
}
  $pgx->{Y}    +=   $pgx->{parameters}->{size_text_px};
  $pgx->{svg}  .=  '
<text x="'.$pgx->{parameters}->{size_plotmargin_px}.'" y="'.$pgx->{Y}.'" style="text-anchor: start; font-size: '.$pgx->{parameters}->{size_text_px}.'px; fill: '.$pgx->{parameters}->{color_label_bottom_hex}.'">'.$pgx->{parameters}->{text_bottom_left}.'</text>
<text x="'.($pgx->{areastartx} + $pgx->{areawidth}).'" y="'.$pgx->{Y}.'" style="text-anchor: end; font-size: '.$pgx->{parameters}->{size_text_px}.'px; fill: '.$pgx->{parameters}->{color_label_bottom_hex}.'">'.$pgx->{parameters}->{text_bottom_right}.'</text>';

  $pgx->{Y}    +=  $pgx->{parameters}->{size_chromosome_padding_px};
  return $pgx;

}

################################################################################
################################################################################
################################################################################

sub svg_add_labels_right {

  my $pgx       =   shift;
  my $labels    =   shift;
  my $thisH     =   shift;
  my $labCol    =   shift;
  
  if ($pgx->{parameters}->{size_clustertree_w_px} < 1) { return $pgx }
  # adding labels to the right side of the plot
  if ($pgx->{parameters}->{size_label_right_px} > 0) {
    if (@$labels < 1) {
      push (
        @$labels,
        {
          label_text  =>  q{},
          label_link  =>  q{},
          label_color =>  $labCol,     
        }
      );
    }
    my $labW    =   sprintf "%.1f", ( $pgx->{parameters}->{size_label_right_px} / @$labels);
    my $labX    =   $pgx->{areaendx} + $pgx->{parameters}->{size_chromosome_padding_px};
    foreach my $lab (@$labels) {
      $pgx->{svg}  .=  '
<rect x="'.$labX.'" y="'.$pgx->{Y}.'" width="'.$labW.'" height="'.$thisH.'" style="fill: '.$lab->{label_color}.'; " />';
      $labX     +=  $labW; 
    }
  }
  # / labels

	return	$pgx;

}

################################################################################
################################################################################
################################################################################

sub svg_add_cluster_tree {

  my $pgx       =   shift;

  if (! defined($pgx->{clustertree})) { return $pgx }

  if (
    $pgx->{parameters}->{size_clustertree_w_px} < 1
     ||
     (scalar @{$pgx->{clustertree}} < 2)
  ) { return $pgx }

  my $nodeN     =   scalar @{ $pgx->{clustertree} } + 1;
  my $treePixH  =   $pgx->{Y} - $pgx->{areastarty};
	my $yPixF			=		$treePixH	/ $nodeN;
	my $yCorr			=		$pgx->{areastarty} + ($treePixH / $nodeN / 2);
	my $maxTreeX	=		(sort {$a <=> $b } ( map{ $_->{NODEX} } @{ $pgx->{clustertree} } ) )[-1];
	my $xPixF			=		$pgx->{parameters}->{size_clustertree_w_px} / ($maxTreeX == 0 ? 1 : $maxTreeX);

	foreach (@{ $pgx->{clustertree} }) {

		my $xStartLpix	=		sprintf "%.1f", $pgx->{areatreex} + $_->{LX} * $xPixF;
		my $xStartRpix	=		sprintf "%.1f", $pgx->{areatreex} + $_->{RX} * $xPixF;
		my $xNodePix	  =		sprintf "%.1f", $pgx->{areatreex} + $_->{NODEX} * $xPixF;
		my $yLpix		    =		sprintf "%.1f", $yCorr + $_->{LY} * $yPixF;
		my $yRpix		    =		sprintf "%.1f", $yCorr + $_->{RY} * $yPixF;

		$pgx->{svg} .=	'
<line x1="'.$xStartLpix.'" y1="'.$yLpix.'" x2="'.$xNodePix.'" y2="'.$yLpix.'" stroke="#666666" />
<line x1="'.$xStartRpix.'" y1="'.$yRpix.'" x2="'.$xNodePix.'" y2="'.$yRpix.'" stroke="#666666" />
<line x1="'.$xNodePix.'" y1="'.$yLpix.'" x2="'.$xNodePix.'" y2="'.$yRpix.'" stroke="#666666" />';

	}

	return	$pgx;

}

################################################################################
################################################################################
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
