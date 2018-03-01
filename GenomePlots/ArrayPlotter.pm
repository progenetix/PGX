package PGX::GenomePlots::ArrayPlotter;

use Data::Dumper;
use GD::Simple;
use MIME::Base64 qw(encode_base64);
use PGX::GenomePlots::CytobandsPlotter;
use PGX::GenomePlots::PlotParameters;

require Exporter;
@ISA            =   qw(Exporter);
@EXPORT         =   qw(return_arrayplot_svg get_arrayplot_area);

################################################################################

sub return_arrayplot_svg {

  our $plot     =   shift;

  $plot->{Y}            =   $plot->{parameters}->{size_plotmargin_top_px};
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
  $plot         =   svg_add_title($plot);
  $plot         =   svg_add_cytobands($plot);
  $plot->{areastarty}   =   $plot->{Y};
  $plot         =   get_arrayplot_area($plot);
  $plot         =   get_fracb_area($plot);
  $plot         =   svg_add_markers($plot);
  $plot         =   svg_add_bottom_text($plot);
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

sub get_arrayplot_area {

=pod

Expects:
  - the plot object with current Y parameter for placing the plot elements on
    the SVG

Returns:
  - the  extended SVG and the increased end Y value as start for
    the next elements

=cut

########    ####    ####    ####    ####    ####    ####    ####    ####    ####

  if ($plot->{parameters}->{size_plotarea_h_px} < 1) { return $plot }

  $plot->{svg}    .=  '
<style type="text/css"><![CDATA[
  .DUP {stroke-width: '.$plot->{parameters}->{size_segments_stroke_px}.'px; stroke: '.$plot->{parameters}->{color_var_dup_hex}.'; opacity: 0.8  }
  .DEL {stroke-width: '.$plot->{parameters}->{size_segments_stroke_px}.'px; stroke: '.$plot->{parameters}->{color_var_del_hex}.'; opacity: 0.8  }
  .cen {stroke-width: '.$plot->{parameters}->{size_centerline_stroke_px}.'px; stroke: '.$plot->{parameters}->{color_plotgrid_hex}.'; opacity: 0.8 ; }
]]></style>';

  $plot->{Y}      +=  $plot->{parameters}->{size_plotarea_padding};

  svg_add_labels_y($plot);

  my $area_x0   =   $plot->{areastartx};
  my $area_y0   =   $plot->{Y};
  my $area_ycen =   $plot->{Y} + $plot->{parameters}->{size_plotarea_h_px} / 2;
  my $area_yn   =   $plot->{Y} + $plot->{parameters}->{size_plotarea_h_px};
  my $lowSegY   =   $area_yn + $plot->{parameters}->{size_segments_stroke_px};

  $plot->{areaendy} =   $area_yn;

  my $probeSize =   $plot->{parameters}->{factor_probedots};
  if (scalar @{ $plot->{probedata} } < 200000) { $probeSize *= 2 }
  if (scalar @{ $plot->{probedata} } < 8000)  { $probeSize *= 2 }

  # probe area, probes & segments   ####    ####    ####    ####    ####    ####

  foreach my $refName (@{ $plot->{parameters}->{chr2plot} }) {

    my $areaW   =  sprintf "%.1f", ($plot->{referencebounds}->{$refName}->[1] - $plot->{referencebounds}->{$refName}->[0]) * $plot->{basepixfrac};

    $plot->{svg}        .=  '
<rect x="'.$area_x0.'" y="'.$plot->{Y}.'" width="'.$areaW.'" height="'.$plot->{parameters}->{size_plotarea_h_px}.'" style="fill: '.$plot->{parameters}->{color_plotarea_hex}.'; fill-opacity: 0.8; " />';

    # probes ###    ####    ####    ####    ####    ####    ####    ####    ####
    my $areaProbes      =  [ grep{ $_->{reference_name} eq $refName } @{ $plot->{probedata} } ];
    $areaProbes         =  [ grep{ $_->{position} <= $plot->{referencebounds}->{$refName}->[1] } @$areaProbes];
    $areaProbes         =  [ grep{ $_->{position} >= $plot->{referencebounds}->{$refName}->[0] } @$areaProbes ];

    # probes are plotted using GD
    my $probeArea       =   GD::Image->new($areaW, $plot->{parameters}->{size_plotarea_h_px}, 1);
    my $gdDotS          =   1 * $probeSize;
    my $gdAreaCol       =   $probeArea->colorAllocate( @{ hex2rgb($plot->{parameters}->{color_plotarea_hex}) } );
    $probeArea->transparent($gdAreaCol);
    my $gdDotcol        = $probeArea->colorAllocateAlpha(32,32,32,63);
    $probeArea->filledRectangle(0, 0, $areaW, $plot->{parameters}->{size_plotarea_h_px}, $gdAreaCol);

    foreach (@$areaProbes) {
      my $dotX  =   sprintf "%.2f", $plot->{basepixfrac} * ($_->{position} - $plot->{referencebounds}->{$refName}->[0]);
      my $dotY  =   sprintf "%.2f", ($plot->{parameters}->{size_plotarea_h_px} / 2 - $_->{value} * $plot->{parameters}->{pixyfactor});
      $probeArea->filledEllipse($dotX, $dotY, $gdDotS, $gdDotS, $gdDotcol);
    }

    # the GD object is encoded as base64 and embedded into the svg
    $plot->{svg}        .=  '
<image
  x="'.$area_x0.'"
  y="'.$plot->{Y}.'"
  width="'.$areaW.'"
  height="'.$plot->{parameters}->{size_plotarea_h_px}.'"
  xlink:href="data:image/png;base64,'.encode_base64($probeArea->png).'"
/>';

    # / probes #    ####    ####    ####    ####    ####    ####    ####    ####

    # segments #    ####    ####    ####    ####    ####    ####    ####    ####

    my $areaSegments    =   [ grep{ $_->{reference_name} eq $refName } @{ $plot->{segmentdata} } ];
    $areaSegments       =   [ grep{ $_->{variant_type} =~ /\w/ } @$areaSegments];
    $areaSegments       =   [ grep{ $_->{start} <= $plot->{referencebounds}->{$refName}->[1] } @$areaSegments];
    $areaSegments       =   [ grep{ $_->{end} >= $plot->{referencebounds}->{$refName}->[0] } @$areaSegments ];

    foreach my $seg (@$areaSegments) {

      if ($seg->{start} < $plot->{referencebounds}->{$refName}->[0]) {
        $seg->{start}   = $plot->{referencebounds}->{$refName}->[0] }
      if ($seg->{end} > $plot->{referencebounds}->{$refName}->[1]) {
        $seg->{end}     = $plot->{referencebounds}->{$refName}->[1] }

      # providing a minimum sub-pixel segment plot length
      my $segPixLen     =   sprintf "%.1f", ($plot->{basepixfrac} * ($seg->{end} - $seg->{start}));
      if ($segPixLen < 0.2) { $segPixLen = 0.2 }

      my $seg_x0        =   sprintf "%.1f", $area_x0 + $plot->{basepixfrac} * ($seg->{start} - $plot->{referencebounds}->{$refName}->[0]);
      my $seg_xn        =   $seg_x0 + $segPixLen;
      my $seg_y         =   sprintf "%.1f", $area_ycen - $seg->{info}->{value} * $plot->{parameters}->{pixyfactor};

      $plot->{svg}      .=  '
<line x1="'.$seg_x0.'"  y1="'.$seg_y.'"  x2="'.$seg_xn.'"  y2="'.$seg_y.'"  class="'.$seg->{variant_type}.'"  />';
      if (scalar @{$plot->{parameters}->{chr2plot}} == 1) {
        $plot->{svg}    .=  '
<a
xlink:href="http://genome.ucsc.edu/cgi-bin/hgTracks?db='.$plot->{parameters}->{genome}.'&amp;position=chr'.$refName.'%3A'.$seg->{start}.'-'.$seg->{end}.'"
xlink:show="new"
xlink:title="'.$seg->{info}->{value}.' at '.$refName.':'.$seg->{start}.'-'.$seg->{end}.'">';
      }
      $plot->{svg}      .=  '
<line x1="'.$seg_x0.'"  y1="'.$lowSegY.'"  x2="'.$seg_xn.'"  y2="'.$lowSegY.'"  class="'.$seg->{variant_type}.'"  />';
      if (scalar @{$plot->{parameters}->{chr2plot}} == 1) {
        $plot->{svg}    .=  '</a>' }

    }

    # / segments    ####    ####    ####    ####    ####    ####    ####    ####

    # moving x to the next chromosome area
    $area_x0    +=  $areaW + $plot->{parameters}->{size_chromosome_padding_px};

 }

  # adding a baseline at 0
  $plot->{svg}  .=  '
<line x1="'.$plot->{parameters}->{size_plotmargin_px}.'"  y1="'.$area_ycen.'"  x2="'.($plot->{parameters}->{size_plotmargin_px} + $plot->{areawidth}).'"  y2="'.$area_ycen.'"  class="cen"  />';

  $plot->{Y}    =  $lowSegY + $plot->{parameters}->{size_segments_stroke_px} / 2;

  return $plot;

}

################################################################################

sub get_fracb_area {

=pod

Expects:
  - the plot object with current Y parameter for placing the plot elements on
    the SVG

Returns:
  - the  extended SVG and the increased end Y value as start for
    the next elements

=cut

########    ####    ####    ####    ####    ####    ####    ####    ####    ####

  if (
    @{ $plot->{probedata_fracb} } < 1
    &&
    @{ $plot->{segmentdata_fracb} } < 1
  ) { return $plot }

  if ($plot->{parameters}->{size_fracbarea_h_px} < 1) { return $plot }

########    ####    ####    ####    ####    ####    ####    ####    ####    ####

  $plot->{Y}      +=  $plot->{parameters}->{size_chromosome_padding_px};

  $plot->{svg}    .=  '
<style type="text/css"><![CDATA[
  .fb {stroke-width: '.$plot->{parameters}->{size_segments_stroke_px}.'px; stroke: #FF3333; opacity: 0.8  }
]]></style>';

  # area init; here, $area_y0 is the lower margin, corresponding then to
  # a 0-value
  my $area_x0   =   $plot->{areastartx};
  my $area_yn   =   $plot->{Y} + $plot->{parameters}->{size_fracbarea_h_px};
  my $area_ycen =   $plot->{Y} + $plot->{parameters}->{size_fracbarea_h_px} / 2;
  my $fbPixYfac =   1 / $plot->{parameters}->{size_fracbarea_h_px};

  $plot->{areaendy} =   $area_yn;

  my $probeSize =   $plot->{parameters}->{factor_probedots};
  if (scalar @{ $plot->{probedata_fracb} } < 200000) { $probeSize *= 2 }
  if (scalar @{ $plot->{probedata_fracb} } < 8000)  { $probeSize *= 2 }

  foreach my $refName (@{ $plot->{parameters}->{chr2plot} }) {

    my $areaW   =  sprintf "%.1f", ($plot->{referencebounds}->{$refName}->[1] - $plot->{referencebounds}->{$refName}->[0]) * $plot->{basepixfrac};

    $plot->{svg}    .=  '
<rect x="'.$area_x0.'" y="'.$plot->{Y}.'" width="'.$areaW.'" height="'.$plot->{parameters}->{size_fracbarea_h_px}.'" style="fill: '.$plot->{parameters}->{color_plotarea_hex}.'; fill-opacity: 0.8; " />';

    # probes ###    ####    ####    ####    ####    ####    ####    ####    ####

    my $areaProbes  =  [ grep{ $_->{reference_name} eq $refName } @{ $plot->{probedata_fracb} } ];
    $areaProbes     =  [ grep{ $_->{position} <= $plot->{referencebounds}->{$refName}->[1] } @$areaProbes];
    $areaProbes     =  [ grep{ $_->{position} >= $plot->{referencebounds}->{$refName}->[0] } @$areaProbes ];

    # probes are plotted using GD

    my $probeArea   =   GD::Image->new($areaW, $plot->{parameters}->{size_fracbarea_h_px}, 1);
    my $gdDotS      =   1 * $probeSize;
    my $gdAreaCol   =   $probeArea->colorAllocate( @{ hex2rgb($plot->{parameters}->{color_plotarea_hex}) } );
    $probeArea->transparent($gdAreaCol);
    my $gdDotcol    = $probeArea->colorAllocateAlpha(0,0,0,111);
    $probeArea->filledRectangle(0, 0, $areaW, $plot->{parameters}->{size_fracbarea_h_px}, $gdAreaCol);

    foreach (@$areaProbes) {
      my $dotX  =   sprintf "%.2f", $plot->{basepixfrac} * ($_->{position} - $plot->{referencebounds}->{$refName}->[0]);
      my $dotY  =   sprintf "%.2f", ($plot->{parameters}->{size_fracbarea_h_px} - $_->{value} / $fbPixYfac);
      $probeArea->filledEllipse($dotX, $dotY, $gdDotS, $gdDotS, $gdDotcol);
    }

    # the GD object is encoded as base64 and embedded into the svg
    $plot->{svg}    .=  '
<image
  x="'.$area_x0.'"
  y="'.$plot->{Y}.'"
  width="'.$areaW.'"
  height="'.$plot->{parameters}->{size_fracbarea_h_px}.'"
  xlink:href="data:image/png;base64,'.encode_base64($probeArea->png).'"
/>';

    # / probes #    ####    ####    ####    ####    ####    ####    ####    ####

    # fracb segments s###    ####    ####    ####    ####    ####    ####    ####
    my $areaSegments  =  [ grep{ $_->{reference_name} eq $refName } @{ $plot->{segmentdata_fracb} } ];
    $areaSegments     =  [ grep{ $_->{start} <= $plot->{referencebounds}->{$refName}->[1] } @$areaSegments];
    $areaSegments     =  [ grep{ $_->{end} >= $plot->{referencebounds}->{$refName}->[0] } @$areaSegments ];

    foreach my $seg (@$areaSegments) {

      if ($seg->{start} < $plot->{referencebounds}->{$refName}->[0]) {
        $seg->{start} = $plot->{referencebounds}->{$refName}->[0] }
      if ($seg->{end} > $plot->{referencebounds}->{$refName}->[1]) {
        $seg->{end} = $plot->{referencebounds}->{$refName}->[1] }

      # providing a minimum sub-pixel segment plot length
      my $segPixLen     =   sprintf "%.1f", ($plot->{basepixfrac} * ($seg->{end} - $seg->{start}));
      if ($segPixLen < 0.2) { $segPixLen = 0.2 }

      my $seg_x0        =   sprintf "%.1f", $area_x0 + $plot->{basepixfrac} * ($seg->{start} - $plot->{referencebounds}->{$refName}->[0]);
      my $seg_xn        =   $seg_x0 + $segPixLen;

      my $val     =   $seg->{info}->{value} / $fbPixYfac;
      my $seg_y   =   sprintf "%.2f", ($area_yn - $val);

      $plot->{svg}      .=  '
<line x1="'.$seg_x0.'" y1="'.$seg_y.'" x2="'.$seg_xn.'" y2="'.$seg_y.'" class="fb" />';

    }

    # / segments    ####    ####    ####    ####    ####    ####    ####    ####

    $area_x0    +=  $areaW + $plot->{parameters}->{size_chromosome_padding_px};
 
  }

  # adding a baseline at 0.5
  $plot->{svg}  .=  '
<line x1="'.$plot->{parameters}->{size_plotmargin_px}.'"  y1="'.$area_ycen.'"  x2="'.($plot->{parameters}->{size_plotmargin_px} + $plot->{areawidth}).'"  y2="'.$area_ycen.'"  class="cen"  />';

  $plot->{Y}    +=  $plot->{parameters}->{size_fracbarea_h_px};

  return $plot;

}

1;
