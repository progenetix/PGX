package PGX::GenomePlots::StripPlotter;

use GD::Simple;
use Data::Dumper;
use PGX::GenomePlots::CytobandsPlotter;
use PGX::GenomePlots::PlotParameters;

require Exporter;
@ISA            =   qw(Exporter);
@EXPORT         =   qw(
  return_stripplot_svg
  get_stripplot_area_gd
  get_frequencystripplot_area_gd
);

################################################################################

sub return_stripplot_svg {

  my $pgx       =   shift;

  $pgx->{svg}   =   q{};

  $pgx->{Y}     =   $pgx->{parameters}->{size_plotmargin_top_px};
  my $plotW     =   $pgx->{parameters}->{size_plotimage_w_px};
  $pgx->{areastartx}   =   $pgx->{parameters}->{size_plotmargin_px} + $pgx->{parameters}->{size_title_left_px};
  if (
    $pgx->{parameters}->{size_clustertree_w_px} < 1
     ||
     (scalar @{ $pgx->{clustertree} } < 2)
  ) { $pgx->{parameters}->{size_clustertree_w_px} = 0 }
  $pgx->{areawidth}    =   $plotW - ($pgx->{areastartx} + $pgx->{parameters}->{size_plotmargin_px} + $pgx->{parameters}->{size_label_right_px} + $pgx->{parameters}->{size_clustertree_w_px});
;
  if (
    $pgx->{parameters}->{do_chromosomes_proportional} =~ /y/i
    &&
    @{ $pgx->{parameters}->{chr2plot} } == 1
  ) {
    $pgx->{areawidth}  *=  ($pgx->{referencebounds}->{ $pgx->{parameters}->{chr2plot}->[0] }->[1] / $pgx->{referencebounds}->{ '1' }->[1]);
    $plotW      =   $pgx->{areawidth} + 2 * $pgx->{parameters}->{size_plotmargin_px};
  }
  $pgx->{areaendx}   =   $pgx->{areastartx} + $pgx->{areawidth};
  $pgx->{areatreex}  =   $pgx->{areaendx};
  if ($pgx->{parameters}->{size_label_right_px} > 0) {
    $pgx->{areatreex}  +=  $pgx->{parameters}->{size_chromosome_padding_px} + $pgx->{parameters}->{size_label_right_px} }
  $pgx->{basepixfrac}  =   ( $pgx->{areawidth} - ($#{ $pgx->{parameters}->{chr2plot} } * $pgx->{parameters}->{size_chromosome_padding_px}) ) / $pgx->{genomesize};
  $pgx->svg_add_title();
  $pgx->svg_add_cytobands();
  if (! $pgx->{samples}) { $pgx->{samples} = [ $pgx->{segmentdata} ] }
  $pgx->get_stripplot_area_gd();
  $pgx->get_frequencystripplot_area_gd();
  $pgx->{markerstarty}  =   $pgx->{areastarty}; # so that markers span the histograms
  $pgx->svg_add_cluster_tree();
  $pgx->svg_add_markers();
  $pgx->svg_add_bottom_text();
  $pgx->{Y}     +=  $pgx->{parameters}->{size_plotmargin_bottom_px};
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

sub get_stripplot_area_gd {

  use GD::Simple;
  use MIME::Base64 qw(encode_base64);
  no warnings 'uninitialized';

=pod

Expects:
  - the current Y parameter for placing the plot elements on the SVG
  - a "sample" object, which has a "name" and "variants":
  
    name : "this-sample"
    variants :
      - start :
          - 50601893
            51088076
        end :
          - 61109900
        reference_name : 2,
        variant_type : "DUP",
      - start : 
      (...)

Returns:
  - the increased end Y value as start for the next elements

=cut

  ######    ####    ####    ####    ####    ####    ####    ####    ####    ####

  my $pgx       =   shift;

  if ($pgx->{parameters}->{size_strip_h_px} < 1) { return $pgx }
  if (! $pgx->{samples}) { return $pgx }

  $pgx->{Y}     +=  $pgx->{parameters}->{size_plotarea_padding};
  $pgx->{areastarty}   =   $pgx->{Y};

  my $defLabCol =   '#dddddd';
  my $altLabCol =   '#fefefe';
  my $labCol    =   '#dddddd';
  
  my $areaH     =   $pgx->{parameters}->{size_strip_h_px} * @{ $pgx->{samples} };

  my $stripArea =   GD::Image->new($pgx->{areawidth}, $areaH, 1);
  my $gdBgCol   =   $stripArea->colorAllocate( @{ hex2rgb($pgx->{parameters}->{color_plotbackground_hex}) } );
  $stripArea->filledRectangle(0, 0, $pgx->{areawidth}, $areaH, $gdBgCol);
  my $gdAreaCol =   $stripArea->colorAllocate( @{ hex2rgb($pgx->{parameters}->{color_plotarea_hex}) } );
  $stripArea->transparent($gdBgCol);
  my $gdDupCol  = $stripArea->colorAllocate( @{ hex2rgb($pgx->{parameters}->{color_var_dup_hex}) } );
  my $gdDelCol  = $stripArea->colorAllocate( @{ hex2rgb($pgx->{parameters}->{color_var_del_hex}) } );

  my $gd_y0     =   0;
  my $gd_yn;
  my $area_x0;
  
  foreach my $sample (@{ $pgx->{samples} }) {

    $gd_yn      =   $gd_y0 + $pgx->{parameters}->{size_strip_h_px};
    $area_x0    =   0;

    my $segSet  =   $sample->{variants};

   if ($sample->{name} =~ /\w\w/) {
      $pgx->{parameters}->{title_left} = $sample->{name} }
    $pgx->svg_add_title_left(($pgx->{Y} + $pgx->{parameters}->{size_strip_h_px} / 2));

    foreach my $refName (@{ $pgx->{parameters}->{chr2plot} }) {

      my $areaW =  sprintf "%.1f", ($pgx->{referencebounds}->{$refName}->[1] - $pgx->{referencebounds}->{$refName}->[0]) * $pgx->{basepixfrac};
      $stripArea->filledRectangle($area_x0, $gd_y0, ($area_x0 + $areaW), $gd_yn, $gdAreaCol);

      my $areaSegs  =   [ grep{ $_->{reference_name} eq $refName } @{ $segSet } ];
      $areaSegs     =   [ grep{ $_->{variant_type} =~ /\w/ } @$areaSegs];
      $areaSegs     =   [ grep{ $_->{start}->[0] <= $pgx->{referencebounds}->{$refName}->[1] } @$areaSegs];
      $areaSegs     =   [ grep{ $_->{end}->[-1] >= $pgx->{referencebounds}->{$refName}->[0] } @$areaSegs ];
      foreach my $seg (@$areaSegs) {
        if ($seg->{start}->[0] < $pgx->{referencebounds}->{$refName}->[0]) {
          $seg->{start}->[0] = $pgx->{referencebounds}->{$refName}->[0] }
        if ($seg->{end}->[-1] > $pgx->{referencebounds}->{$refName}->[1]) {
          $seg->{end}->[-1]   = $pgx->{referencebounds}->{$refName}->[1] }

        my $seg_x0  =   sprintf "%.1f", $area_x0 + $pgx->{basepixfrac} * ($seg->{start}->[0] - $pgx->{referencebounds}->{$refName}->[0]);
        my $segPixEnd   =   sprintf "%.2f", ($seg_x0 + $pgx->{basepixfrac} * ($seg->{end}->[-1] - $seg->{start}->[0]));
        # providing a minimum sub-pixel segment plot length
        if ($segPixEnd < 0.02) { $segPixLen = 0.02 }

        if ($seg->{variant_type} =~ /DEL/i) {
          $stripArea->filledRectangle($seg_x0, $gd_y0, $segPixEnd, $gd_yn, $gdDelCol) }
        elsif ($seg->{variant_type} =~ /DUP/i) {
          $stripArea->filledRectangle($seg_x0, $gd_y0, $segPixEnd, $gd_yn, $gdDupCol) }

      }

      $area_x0  +=  $areaW + $pgx->{parameters}->{size_chromosome_padding_px};

    }

    my $labels_R    =   [];
    if ($sample->{labels}) { $labels_R = $sample->{labels} }
    # fallback color; defined here for alternation...
    if ($labCol eq $altLabCol) { 
      $labCol   =   $defLabCol }
    else { 
      $labCol   =   $altLabCol }
    $pgx->svg_add_labels_right($labels_R, $pgx->{parameters}->{size_strip_h_px}, $labCol);

    $pgx->{Y}   +=  $pgx->{parameters}->{size_strip_h_px};
    $gd_y0      +=  $pgx->{parameters}->{size_strip_h_px};

  }
  
  $pgx->{svg}   .=  '
<image
  x="'.$pgx->{areastartx}.'"
  y="'.$pgx->{areastarty}.'"
  width="'.$pgx->{areawidth}.'"
  height="'.$areaH.'"
  xlink:href="data:image/png;base64,'.encode_base64($stripArea->png).'"
/>';


  return $pgx;

}
################################################################################

sub get_frequencystripplot_area_gd {

  use GD::Simple;
  use MIME::Base64 qw(encode_base64);

=pod

Expects:
  - the current Y parameter for placing the plot elements on the SVG

Returns:
  - the increased end Y value as start for the next elements

=cut

  ######    ####    ####    ####    ####    ####    ####    ####    ####    ####

  my $pgx       =   shift;

  if ($pgx->{parameters}->{size_strip_h_px} < 1) { return $pgx }
  if (! $pgx->{frequencymaps}->[0]) { return $pgx }

  $pgx->{Y}    +=  $pgx->{parameters}->{size_plotarea_padding};
  $pgx->{areastarty}   =   $pgx->{Y};

  my $defLabCol =   '#dddddd';
  my $altLabCol =   '#fefefe';
  my $labCol    =   '#dddddd';

  my $areaH     =   $pgx->{parameters}->{size_strip_h_px} * @{ $pgx->{frequencymaps} };
  my $stripArea =   GD::Image->new($pgx->{areawidth}, $areaH, 1);
  my $gdBgCol   =   $stripArea->colorAllocate( @{ hex2rgb($pgx->{parameters}->{color_plotbackground_hex}) } );
  $stripArea->filledRectangle(0, 0, $pgx->{areawidth}, $areaH, $gdBgCol);
  my $gdAreaCol =   $stripArea->colorAllocate( @{ hex2rgb($pgx->{parameters}->{color_plotarea_hex}) } );
  $stripArea->transparent($gdBgCol);

  my $gd_y0     =   0;
  my $gd_yn;
  my $area_x0;
    
  foreach my $frequencymapsSet (@{ $pgx->{frequencymaps} }) {

    $fMapIndex++;

    if ($frequencymapsSet->{name} =~ /\w\w/) {
      $pgx->{parameters}->{title_left}   =   $frequencymapsSet->{name} }

    $gd_yn      =   $gd_y0 + $pgx->{parameters}->{size_strip_h_px};
    $area_x0    =   0;

    $pgx->svg_add_title_left(($pgx->{Y} + $pgx->{parameters}->{size_strip_h_px} / 2));

    my $maxF		=		(sort {$a <=> $b} (@{ $frequencymapsSet->{dupfrequencies} }, @{ $frequencymapsSet->{delfrequencies} }) )[-1];

    foreach my $refName (@{ $pgx->{parameters}->{chr2plot} }) {

      my $areaBases =   $pgx->{referencebounds}->{$refName}->[1] - $pgx->{referencebounds}->{$refName}->[0];
      my $areaW     =   sprintf "%.1f", $areaBases * $pgx->{basepixfrac};

      # intervals through index #    ####    ####    ####    ####    ####    ####
      my @ind   =  grep{ $refName eq $pgx->{genomeintervals}->[$_]->{reference_name} } 0..$#{ $pgx->{genomeintervals} };
      @ind      =  grep{ $pgx->{genomeintervals}->[$_]->{start} <=  $pgx->{referencebounds}->{$refName}->[1]  } @ind;
      @ind      =  grep{ $pgx->{genomeintervals}->[$_]->{end} >=  $pgx->{referencebounds}->{$refName}->[0]  } @ind;

      foreach my $i (@ind) {

        my $start =   $pgx->{genomeintervals}->[$i]->{start};
        my $end   =   $pgx->{genomeintervals}->[$i]->{end};
        if ($start < $pgx->{referencebounds}->{$refName}->[0]) {
          $start  = $pgx->{referencebounds}->{$refName}->[0] }
        if ($end > $pgx->{referencebounds}->{$refName}->[1]) {
          $end    = $pgx->{referencebounds}->{$refName}->[1] }

        my $seg_x0  =   sprintf "%.1f", $area_x0 + $pgx->{basepixfrac} * ($start - $pgx->{referencebounds}->{$refName}->[0]);
        my $segPixEnd   =   sprintf "%.2f", ($seg_x0 + $pgx->{basepixfrac} * ($end - $start));
        # providing a minimum sub-pixel segment plot length

        my $fill  =   frequencies2rgb(
                        $pgx->{parameters},
                        $frequencymapsSet->{dupfrequencies}->[$i],
                        $frequencymapsSet->{delfrequencies}->[$i],
                        $maxF,
                      );
        my $gdCol =   $stripArea->colorAllocate( split(',', $fill) );

        $stripArea->filledRectangle($seg_x0, $gd_y0, $segPixEnd, $gd_yn, $gdCol);

      }
      
      $area_x0    +=  $areaW + $pgx->{parameters}->{size_chromosome_padding_px};

    }
    
    my $labels_R  =   [];
    if ($frequencymapsSet->{labels}) { $labels_R = $frequencymapsSet->{labels} }
    # fallback color; defined here for alternation...
    if ($labCol eq $altLabCol) { 
      $labCol   =   $defLabCol }
    else { 
      $labCol   =   $altLabCol }
    $pgx->svg_add_labels_right($labels_R, $pgx->{parameters}->{size_strip_h_px}, $labCol);

    $pgx->{Y}   +=  $pgx->{parameters}->{size_strip_h_px};
    $gd_y0      +=   $pgx->{parameters}->{size_strip_h_px};

  }

  $pgx->{svg}   .=  '
<image
  x="'.$pgx->{areastartx}.'"
  y="'.$pgx->{areastarty}.'"
  width="'.$pgx->{areawidth}.'"
  height="'.$areaH.'"
  xlink:href="data:image/png;base64,'.encode_base64($stripArea->png).'"
/>';

  return $pgx;

}

################################################################################


1;
