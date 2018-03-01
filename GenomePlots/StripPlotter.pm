package PGX::GenomePlots::StripPlotter;

use Data::Dumper;
use PGX::GenomePlots::CytobandsPlotter;
use PGX::GenomePlots::PlotParameters;

require Exporter;
@ISA            =   qw(Exporter);
@EXPORT         =   qw(
  return_stripplot_svg
  get_stripplot_area
  get_stripplot_area_gd
  get_frequencystripplot_area
);

################################################################################

sub return_stripplot_svg {

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
  $plotW        +=  $plot->{parameters}->{size_clustertree_w_px};
  $plotW        +=  $plot->{parameters}->{size_label_right_px};
  $plot->{areaendx}   =   $plot->{areastartx} + $plot->{areawidth};
  $plot->{areatreex}  =   $plot->{areaendx};
  if ($plot->{parameters}->{size_label_right_px} > 0) {
    $plot->{areatreex}  +=  $plot->{parameters}->{size_chromosome_padding_px} + $plot->{parameters}->{size_label_right_px} }
  $plot->{basepixfrac}  =   ( $plot->{areawidth} - ($#{ $plot->{parameters}->{chr2plot} } * $plot->{parameters}->{size_chromosome_padding_px}) ) / $plot->{genomesize};
  $plot->svg_add_title();
  $plot->svg_add_cytobands();
  if (! $plot->{samples}) { $plot->{samples} = [ $plot->{segmentdata} ] }
  $plot->get_stripplot_area_gd();
  $plot->get_frequencystripplot_area();
  $plot->svg_add_cluster_tree();
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

sub get_stripplot_area_gd {

  use GD::Simple;
  use MIME::Base64 qw(encode_base64);
  no warnings 'uninitialized';

=pod

Expects:
  - the current Y parameter for placing the plot elements on the SVG

Returns:
  - the increased end Y value as start for the next elements

=cut

  ######    ####    ####    ####    ####    ####    ####    ####    ####    ####

  my $plot      =   shift;

  if ($plot->{parameters}->{size_strip_h_px} < 1) { return $plot }
  if (! $plot->{samples}) { return $plot }

  $plot->{Y}    +=  $plot->{parameters}->{size_plotarea_padding};
  $plot->{areastarty}   =   $plot->{Y};

  my $defLabCol =   '#dddddd';
  my $altLabCol =   '#fefefe';
  my $labCol    =   '#dddddd';
  
  my $areaH     =   $plot->{parameters}->{size_strip_h_px} * @{$plot->{samples}};
  my $stripArea =   GD::Image->new($plot->{areawidth}, $areaH, 1);
  my $gdBgCol   =   $stripArea->colorAllocate( @{ hex2rgb($plot->{parameters}->{color_plotbackground_hex}) } );
  $stripArea->filledRectangle(0, 0, $plot->{areawidth}, $areaH, $gdBgCol);
  my $gdAreaCol =   $stripArea->colorAllocate( @{ hex2rgb($plot->{parameters}->{color_plotarea_hex}) } );
  $stripArea->transparent($gdBgCol);
  my $gdDupCol  = $stripArea->colorAllocate( @{ hex2rgb($plot->{parameters}->{color_var_dup_hex}) } );
  my $gdDelCol  = $stripArea->colorAllocate( @{ hex2rgb($plot->{parameters}->{color_var_del_hex}) } );

  my $gd_y0     =   0;
  my $gd_yn;
  my $area_x0;
  
  foreach my $sample (@{$plot->{samples}}) {

    $gd_yn      =   $gd_y0 + $plot->{parameters}->{size_strip_h_px};
    $area_x0    =   0;

    my $segSet  =   $sample->{variants};
 
   if ($sample->{name} =~ /\w\w/) {
      $plot->{parameters}->{title_left} = $sample->{name} }
    $plot->svg_add_title_left(($plot->{Y} + $plot->{parameters}->{size_strip_h_px} / 2));

    foreach my $refName (@{ $plot->{parameters}->{chr2plot} }) {
      my $areaW =  sprintf "%.1f", ($plot->{referencebounds}->{$refName}->[1] - $plot->{referencebounds}->{$refName}->[0]) * $plot->{basepixfrac};
      $stripArea->filledRectangle($area_x0, $gd_y0, ($area_x0 + $areaW), $gd_yn, $gdAreaCol);
#print Dumper(join(' :: ', $sample->{name}, $refName, $area_x0, $gd_y0, $gd_yn));

      my $areaSegs  =   [ grep{ $_->{reference_name} eq $refName } @{ $segSet } ];
      $areaSegs     =   [ grep{ $_->{variant_type} =~ /\w/ } @$areaSegs];
      $areaSegs     =   [ grep{ $_->{start} <= $plot->{referencebounds}->{$refName}->[1] } @$areaSegs];
      $areaSegs     =   [ grep{ $_->{end} >= $plot->{referencebounds}->{$refName}->[0] } @$areaSegs ];

      foreach my $seg (@$areaSegs) {

        if ($seg->{start} < $plot->{referencebounds}->{$refName}->[0]) {
          $seg->{start} = $plot->{referencebounds}->{$refName}->[0] }
        if ($seg->{end} > $plot->{referencebounds}->{$refName}->[1]) {
          $seg->{end}   = $plot->{referencebounds}->{$refName}->[1] }

        my $seg_x0  =   sprintf "%.1f", $area_x0 + $plot->{basepixfrac} * ($seg->{start} - $plot->{referencebounds}->{$refName}->[0]);
        my $segPixEnd   =   sprintf "%.2f", ($seg_x0 + $plot->{basepixfrac} * ($seg->{end} - $seg->{start}));
        # providing a minimum sub-pixel segment plot length
        if ($segPixEnd < 0.02) { $segPixLen = 0.02 }

        if ($seg->{variant_type} =~ /DEL/i) {
          $stripArea->filledRectangle($seg_x0, $gd_y0, $segPixEnd, $gd_yn, $gdDelCol) }
        elsif ($seg->{variant_type} =~ /DUP/i) {
          $stripArea->filledRectangle($seg_x0, $gd_y0, $segPixEnd, $gd_yn, $gdDupCol) }

      }

      $area_x0  +=  $areaW + $plot->{parameters}->{size_chromosome_padding_px};

    }

    my $labels_R    =   [];
    if ($sample->{labels}) { $labels_R = $sample->{labels} }
    # fallback color; defined here for alternation...
    if ($labCol eq $altLabCol) { 
      $labCol =   $defLabCol }
    else { 
      $labCol =   $altLabCol }
    $plot->svg_add_labels_right($labels_R, $plot->{parameters}->{size_strip_h_px}, $labCol);

    $plot->{Y}  +=  $plot->{parameters}->{size_strip_h_px};
    $gd_y0      +=  $plot->{parameters}->{size_strip_h_px};

  }
  
  $plot->{svg}    .=  '
<image
  x="'.$plot->{areastartx}.'"
  y="'.$plot->{areastarty}.'"
  width="'.$plot->{areawidth}.'"
  height="'.$areaH.'"
  xlink:href="data:image/png;base64,'.encode_base64($stripArea->png).'"
/>';


  return $plot;

}
################################################################################

sub get_stripplot_area {

  no warnings 'uninitialized';

=pod

Expects:
  - the current Y parameter for placing the plot elements on the SVG

Returns:
  - the increased end Y value as start for the next elements

=cut

  ######    ####    ####    ####    ####    ####    ####    ####    ####    ####

  my $plot      =   shift;

  if ($plot->{parameters}->{size_strip_h_px} < 1) { return $plot }
  if (! $plot->{samples}) { return $plot }

  $plot->{Y}    +=  $plot->{parameters}->{size_plotarea_padding};
  $plot->{areastarty}   =   $plot->{Y};

  my $defLabCol =   '#dddddd';
  my $altLabCol =   '#fefefe';
  my $labCol    =   '#dddddd';
  
  foreach my $sample (@{$plot->{samples}}) {

    my $area_x0 =   $plot->{areastartx};
    my $segSet  =   $sample->{variants};
 
   if ($sample->{name} =~ /\w\w/) {
      $plot->{parameters}->{title_left} = $sample->{name} }
    $plot->svg_add_title_left(($plot->{Y} + $plot->{parameters}->{size_strip_h_px} / 2));
    foreach my $refName (@{ $plot->{parameters}->{chr2plot} }) {

      my $areaW =  sprintf "%.1f", ($plot->{referencebounds}->{$refName}->[1] - $plot->{referencebounds}->{$refName}->[0]) * $plot->{basepixfrac};

      $plot->{svg}  .=  '
  <rect x="'.$area_x0.'" y="'.$plot->{Y}.'" width="'.$areaW.'" height="'.$plot->{parameters}->{size_strip_h_px}.'" style="fill: '.$plot->{parameters}->{color_plotarea_hex}.'; fill-opacity: 0.8; " />';

      my $areaSegs  =   [ grep{ $_->{reference_name} eq $refName } @{ $segSet } ];
      $areaSegs     =   [ grep{ $_->{variant_type} =~ /\w/ } @$areaSegs];
      $areaSegs     =   [ grep{ $_->{start} <= $plot->{referencebounds}->{$refName}->[1] } @$areaSegs];
      $areaSegs     =   [ grep{ $_->{end} >= $plot->{referencebounds}->{$refName}->[0] } @$areaSegs ];

      foreach my $seg (@$areaSegs) {

        if ($seg->{start} < $plot->{referencebounds}->{$refName}->[0]) {
          $seg->{start} = $plot->{referencebounds}->{$refName}->[0] }
        if ($seg->{end} > $plot->{referencebounds}->{$refName}->[1]) {
          $seg->{end}   = $plot->{referencebounds}->{$refName}->[1] }

        # providing a minimum sub-pixel segment plot length
        my $segPixLen   =   sprintf "%.1f", ($plot->{basepixfrac} * ($seg->{end} - $seg->{start}));
        if ($segPixLen < 0.2) { $segPixLen = 0.2 }

        my $seg_x0  =   sprintf "%.1f", $area_x0 + $plot->{basepixfrac} * ($seg->{start} - $plot->{referencebounds}->{$refName}->[0]);

          $plot->{svg}  .=  '
<a xlink:href="http://genome.ucsc.edu/cgi-bin/hgTracks?db='.$plot->{parameters}->{genome}.'&amp;position=chr'.$refName.'%3A'.$seg->{start}.'-'.$seg->{end}.'" xlink:show="new" xlink:title="'.$seg->{info}->{value}.' at  '.$refName.':'.$seg->{start}.'-'.$seg->{end}.'">';

        $plot->{svg}    .=  '
  <rect x="'.$seg_x0.'" y="'.$plot->{Y}.'" width="'.$segPixLen.'" height="'.$plot->{parameters}->{size_strip_h_px}.'" fill="'.($seg->{variant_type} =~ /DEL/i ? $plot->{parameters}->{color_var_del_hex} : $plot->{parameters}->{color_var_dup_hex}).'" stroke-width="0px" />';
          $plot->{svg}  .=  '</a>';

      }

      $area_x0  +=  $areaW + $plot->{parameters}->{size_chromosome_padding_px};

    }

    my $labels_R    =   [];
    if ($sample->{labels}) { $labels_R = $sample->{labels} }
    # fallback color; defined here for alternation...
    if ($labCol eq $altLabCol) { 
      $labCol =   $defLabCol }
    else { 
      $labCol =   $altLabCol }
    $plot->svg_add_labels_right($labels_R, $plot->{parameters}->{size_strip_h_px}, $labCol);

    $plot->{Y}  +=  $plot->{parameters}->{size_strip_h_px};

  }

  return $plot;

}

################################################################################

sub get_frequencystripplot_area {

=pod

Expects:
  - the current Y parameter for placing the plot elements on the SVG

Returns:
  - the increased end Y value as start for the next elements

=cut

  ######    ####    ####    ####    ####    ####    ####    ####    ####    ####

  my $plot      =   shift;

  if ($plot->{parameters}->{size_strip_h_px} < 1) { return $plot }
  if (! $plot->{frequencymaps}->[0]) { return $plot }

  $plot->{Y}    +=  $plot->{parameters}->{size_plotarea_padding};
  $plot->{areastarty}   =   $plot->{Y};

  my $fMapIndex =  0;
  my $defLabCol =   '#dddddd';
  my $altLabCol =   '#fefefe';
  my $labCol    =   '#dddddd';
  
  foreach my $frequencymapsSet (@{ $plot->{frequencymaps} }) {

    $fMapIndex++;

    if ($frequencymapsSet->{name} =~ /\w\w/) {
      $plot->{parameters}->{title_left}   =   $frequencymapsSet->{name} }

    my $area_x0 =   $plot->{areastartx};

    $plot->svg_add_title_left(($plot->{Y} + $plot->{parameters}->{size_strip_h_px} / 2));

    my $gradientDefs;

    my $maxF		=		(sort {$a <=> $b} (@{ $frequencymapsSet->{dupfrequencies} }, @{ $frequencymapsSet->{delfrequencies} }) )[-1];

    foreach my $refName (@{ $plot->{parameters}->{chr2plot} }) {

      $gradientDefs .=   '
<linearGradient id="frequencygradient_'.$fMapIndex.'_'.$refName.'"  x1="0" x2="1" y1="0" y2="0">';

      my $areaBases =   $plot->{referencebounds}->{$refName}->[1] - $plot->{referencebounds}->{$refName}->[0];
      my $areaW     =   sprintf "%.1f", $areaBases * $plot->{basepixfrac};

      $plot->{svg}  .=  '
<rect x="'.$area_x0.'" y="'.$plot->{Y}.'" width="'.$areaW.'" height="'.$plot->{parameters}->{size_strip_h_px}.'" style="fill: url(#frequencygradient_'.$fMapIndex.'_'.$refName.');" />';

      # intervals through index #    ####    ####    ####    ####    ####    ####
      my @ind   =  grep{ $refName eq $plot->{genomeintervals}->[$_]->{reference_name} } 0..$#{ $plot->{genomeintervals} };
      @ind      =  grep{ $plot->{genomeintervals}->[$_]->{start} <=  $plot->{referencebounds}->{$refName}->[1]  } @ind;
      @ind      =  grep{ $plot->{genomeintervals}->[$_]->{end} >=  $plot->{referencebounds}->{$refName}->[0]  } @ind;

      foreach my $i (@ind) {

        my $start =   $plot->{genomeintervals}->[$i]->{start};
        my $end   =   $plot->{genomeintervals}->[$i]->{end};
        if ($start < $plot->{referencebounds}->{$refName}->[0]) {
          $start  = $plot->{referencebounds}->{$refName}->[0] }
        if ($end > $plot->{referencebounds}->{$refName}->[1]) {
          $end    = $plot->{referencebounds}->{$refName}->[1] }

        my $currentFraction		=	  sprintf "%.3f", 1 - ($plot->{referencebounds}->{$refName}->[1] - ($start + $end) / 2 ) / $areaBases;

        my $fill  =   frequencies2rgb(
                        $plot->{parameters},
                        $frequencymapsSet->{dupfrequencies}->[$i],
                        $frequencymapsSet->{delfrequencies}->[$i],
                        $maxF,
                      );

        $gradientDefs .=  '
<stop offset="'.$currentFraction.'" stop-color="rgb('.$fill.')"/>';

      }

      $gradientDefs   .=   '
</linearGradient>';
      $area_x0    +=  $areaW + $plot->{parameters}->{size_chromosome_padding_px};

    }
    
    my $labels_R    =   [];
    if ($frequencymapsSet->{labels}) { $labels_R = $frequencymapsSet->{labels} }
    # fallback color; defined here for alternation...
    if ($labCol eq $altLabCol) { 
      $labCol =   $defLabCol }
    else { 
      $labCol =   $altLabCol }
    $plot->svg_add_labels_right($labels_R, $plot->{parameters}->{size_strip_h_px}, $labCol);

    $plot->{svg}  .=  $gradientDefs;
    $plot->{Y}    +=  $plot->{parameters}->{size_strip_h_px};

  }

  return $plot;

}


1;
