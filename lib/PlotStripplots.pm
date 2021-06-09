package lib::PlotStripplots;

use GD::Simple;
use Data::Dumper;
use lib::PlotCytobands;
use lib::PlotMakeParameters;

require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(
  return_stripplot_svg
  get_stripplot_area_gd
  get_frequencystripplot_area_gd
);

################################################################################

sub return_stripplot_svg {

	my $pgx = shift;
	
	my $pp = $pgx->{parameters};

	$pgx->{svg} = q{};

	$pgx->{Y} = $pp->{size_plotmargin_top_px};

	my $plotW = $pp->{size_plotimage_w_px};
	$pgx->{areastartx} = $pp->{size_plotmargin_px} + $pp->{size_title_left_w_px};

	if (
		$pp->{size_clustertree_w_px} < 1
		 ||
		 (scalar @{ $pgx->{clustertree} } < 2)
	) { $pp->{size_clustertree_w_px} = 0 }
	$pgx->{areawidth} = $plotW - ($pgx->{areastartx} + $pp->{size_plotmargin_px} + $pp->{size_label_right_px} + $pp->{size_clustertree_w_px});
	;
	
	if (
		$pp->{do_chromosomes_proportional} =~ /y/i
		&&
		@{ $pp->{chr2plot} } == 1
	) {
		$pgx->{areawidth}  *= ($pgx->{referencebounds}->{ $pp->{chr2plot}->[0] }->[1] / $pgx->{referencebounds}->{ '1' }->[1]);
		$plotW = $pgx->{areawidth} + $pgx->{areastartx} + $pp->{size_plotmargin_px} + $pp->{size_label_right_px} + $pp->{size_clustertree_w_px};
	}
	
	$pgx->{areaendx} = $pgx->{areastartx} + $pgx->{areawidth};
	$pgx->{areatreex} = $pgx->{areaendx};
	
	if ($pp->{size_label_right_px} > 0) {
		$pgx->{areatreex}  += $pp->{size_chromosome_padding_px} + $pp->{size_label_right_px} }

	$pgx->{basepixfrac} = ( $pgx->{areawidth} - ($#{ $pp->{chr2plot} } * $pp->{size_chromosome_padding_px}) ) / $pgx->{genomesize};
	
	$pgx->svg_add_title();
	$pgx->svg_add_cytobands();
	
	if (! $pgx->{samples}) {
		$pgx->{samples} = [ $pgx->{segmentdata} ] }
		
	$pgx->get_stripplot_area_gd();
	$pgx->get_frequencystripplot_area_gd();
	
	$pgx->{markerstarty} = $pgx->{areastarty}; # so that markers span the histograms
	
	$pgx->svg_add_cluster_tree();
	$pgx->svg_add_markers();
	$pgx->svg_add_bottom_text();
	
	$pgx->{Y} += $pp->{size_plotmargin_bottom_px};

	my $plotH = sprintf "%.0f", $pgx->{Y};
	$plotW = sprintf "%.0f", $plotW;
	$pgx->{svg} = '<svg
xmlns="http://www.w3.org/2000/svg"
xmlns:xlink="http://www.w3.org/1999/xlink"
version="1.1"
id="'.$pgx->{plotid}.'"
width="'.$plotW.'px"
height="'.$plotH.'px"
style="margin: auto; font-family: Helvetica, sans-serif;">
<style type="text/css"><![CDATA[
  .title-left {text-anchor: middle; fill: '.$pp->{color_text_hex}.'; font-size: '.$pp->{size_text_title_left_px}.'px;}
]]></style>

<rect x="0" y="0" width="'.$plotW.'" height="'.$plotH.'" style="fill: '.$pp->{color_plotbackground_hex}.'; " />

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

	my $pgx = shift;
	my $pp = $pgx->{parameters};

	if ($pp->{size_strip_h_px} < 1) { return $pgx }
	if (! $pgx->{samples}) { return $pgx }

	$pgx->{Y} += $pp->{size_plotarea_padding};
	$pgx->{areastarty} = $pgx->{Y};

	my $defLabCol = '#dddddd';
	my $altLabCol = '#fefefe';
	my $labCol = '#dddddd';

	my $areaH = $pp->{size_strip_h_px} * @{ $pgx->{samples} };

	my $gdArea = GD::Image->new($pgx->{areawidth}, $areaH, 1);

	my $gdBgCol = $gdArea->colorAllocate( @{ hex2rgb($pp->{color_plotbackground_hex}) } );
	$gdArea->filledRectangle(0, 0, $pgx->{areawidth}, $areaH, $gdBgCol);
	my $gdAreaCol = $gdArea->colorAllocate( @{ hex2rgb($pp->{color_plotarea_hex}) } );
	$gdArea->transparent($gdBgCol);
	my $gdDupCol = $gdArea->colorAllocate( @{ hex2rgb($pp->{color_var_dup_hex}) } );
	my $gdDelCol = $gdArea->colorAllocate( @{ hex2rgb($pp->{color_var_del_hex}) } );

	my $gd_y0 = 0;
	my $gd_yn;
	my $area_x0;
	
	# sample label size adjustment	
	for my $i (0..$#{ $pgx->{samples} }) {
		if ($pgx->{samples}->[$i]->{name} =~ /\w\w/) {
			my $t_l =  length($pgx->{samples}->[$i]->{name}) * 0.5 * $pgx->{parameters}->{size_text_title_left_px};
			if ($t_l > $pgx->{parameters}->{size_title_left_w_px}) {
				my ($t_text, $t_code) = ( $pgx->{samples}->[$i]->{name} =~ /^(.*?)( \([^\()]+?\))?$/ );
				if (length($t_text) > $pgx->{parameters}->{text_labels_left_max_letters}) {		
					$t_maxl = $pgx->{parameters}->{text_labels_left_max_letters} - length($t_code);
					$t_text =~ s/^(.{6,$t_maxl}[\w\.]) .*?$/$1.../;
					$pgx->{samples}->[$i]->{name} = $t_text.$t_code;
				}	
			}
		}	
	}
	
	# sample label text size adjustment	
	$pgx->{parameters}->{size_text_title_left_px} = text_width_from_array(
		[ map{$_->{name}} @{ $pgx->{samples} } ],
		$pgx->{parameters}->{size_text_title_left_px},
		$pgx->{parameters}->{size_title_left_w_px}
	);

	foreach my $sample (@{ $pgx->{samples} }) {

		$gd_yn = $gd_y0 + $pp->{size_strip_h_px};
		$area_x0 = 0;

		my $segSet = $sample->{variants};

		if ($sample->{name} =~ /\w\w/) {
			my $titeL = {
				text => $sample->{name},
				pos_y => $pgx->{Y} + $pp->{size_strip_h_px} * 0.5,
				linkout => q{},
			};
			if ($pgx->{datasetid} =~ /.../ && $sample->{id} =~ /.../) {
				$titeL->{linkout} = '/callsets/details?id='.$sample->{id} }
			$pgx->svg_add_title_left($titeL);
		}

		foreach my $chro (@{ $pp->{chr2plot} }) {

		  my $areaW = sprintf "%.1f", ($pgx->{referencebounds}->{$chro}->[1] - $pgx->{referencebounds}->{$chro}->[0]) * $pgx->{basepixfrac};
		  $gdArea->filledRectangle($area_x0, $gd_y0, ($area_x0 + $areaW), $gd_yn, $gdAreaCol);

		  my $areaSegs = [ grep{ $_->{reference_name} eq $chro } @{ $segSet } ];
		  $areaSegs = [ grep{ $_->{variant_type} =~ /\w/ } @$areaSegs];
		  $areaSegs = [ grep{ $_->{start} <= $pgx->{referencebounds}->{$chro}->[1] } @$areaSegs];
		  $areaSegs = [ grep{ $_->{end} >= $pgx->{referencebounds}->{$chro}->[0] } @$areaSegs ];
		  foreach my $seg (@$areaSegs) {
			if ($seg->{start} < $pgx->{referencebounds}->{$chro}->[0]) {
			  $seg->{start} = $pgx->{referencebounds}->{$chro}->[0] }
			if ($seg->{end} > $pgx->{referencebounds}->{$chro}->[1]) {
			  $seg->{end} = $pgx->{referencebounds}->{$chro}->[1] }

			my $seg_x0 = sprintf "%.1f", $area_x0 + $pgx->{basepixfrac} * ($seg->{start} - $pgx->{referencebounds}->{$chro}->[0]);
			my $segPixEnd = sprintf "%.2f", ($seg_x0 + $pgx->{basepixfrac} * ($seg->{end} - $seg->{start} ));
			# providing a minimum sub-pixel segment plot length
			if ($segPixEnd < 0.02) { $segPixLen = 0.02 }

			if ($seg->{variant_type} =~ /DEL/i) {
			  $gdArea->filledRectangle($seg_x0, $gd_y0, $segPixEnd, $gd_yn, $gdDelCol) }
			elsif ($seg->{variant_type} =~ /DUP/i) {
			  $gdArea->filledRectangle($seg_x0, $gd_y0, $segPixEnd, $gd_yn, $gdDupCol) }

		  }

		  $area_x0 += $areaW + $pp->{size_chromosome_padding_px};

		}

		my $labels_R = [];
		if ($sample->{labels}) { $labels_R = $sample->{labels} }
		# fallback color; defined here for alternation...
		if ($labCol eq $altLabCol) { 
		  $labCol = $defLabCol }
		else { 
		  $labCol = $altLabCol }
		$pgx->svg_add_labels_right($labels_R, $pp->{size_strip_h_px}, $labCol);

		$pgx->{Y} += $pp->{size_strip_h_px};
		$gd_y0 += $pp->{size_strip_h_px};

	}
  
	$pgx->{svg} .= '
<image
  x="'.$pgx->{areastartx}.'"
  y="'.$pgx->{areastarty}.'"
  width="'.$pgx->{areawidth}.'"
  height="'.$areaH.'"
  xlink:href="data:image/png;base64,'.encode_base64($gdArea->png).'"
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

  my $pgx = shift;
  my $pp = $pgx->{parameters};

  if ($pp->{size_strip_h_px} < 1) { return $pgx }
  if (! $pgx->{frequencymaps}->[0]) { return $pgx }

  $pgx->{Y} += $pp->{size_plotarea_padding};
  $pgx->{areastarty} = $pgx->{Y};

  my $defLabCol = '#dddddd';
  my $altLabCol = '#fefefe';
  my $labCol = '#dddddd';

  my $areaH = $pp->{size_strip_h_px} * @{ $pgx->{frequencymaps} };
  my $gdArea = GD::Image->new($pgx->{areawidth}, $areaH, 1);
  my $gdBgCol = $gdArea->colorAllocate( @{ hex2rgb($pp->{color_plotarea_hex}) } );
  $gdArea->filledRectangle(0, 0, $pgx->{areawidth}, $areaH, $gdBgCol);

  my $gd_y0 = 0;
  my $gd_yn;
  my $area_x0;
    
  foreach my $frequencymapsSet (@{ $pgx->{frequencymaps} }) {

    $fMapIndex++;

    if ($frequencymapsSet->{name} =~ /\w\w/) {
      my $titeL = {
        text => $frequencymapsSet->{name},
        pos_y => $pgx->{Y} + $pp->{size_strip_h_px} * 0.5,
        linkout => q{},
      };
      $pgx->svg_add_title_left($titeL);
    }

    $gd_yn = $gd_y0 + $pp->{size_strip_h_px};
    $area_x0 = 0;

    my $maxF = (sort {$a <=> $b} (map{ $_->{gain_frequency}, $_->{loss_frequency} } @{ $frequencymapsSet->{intervals} }) )[-1];

    foreach my $chro (@{ $pp->{chr2plot} }) {

      my $areaBases = $pgx->{referencebounds}->{$chro}->[1] - $pgx->{referencebounds}->{$chro}->[0];
      my $areaW = sprintf "%.1f", $areaBases * $pgx->{basepixfrac};

      # intervals through index #    ####    ####    ####    ####    ####    ####
      my @ind = grep{ $chro eq $pgx->{genomeintervals}->[$_]->{reference_name} } 0..$#{ $pgx->{genomeintervals} };
      @ind = grep{ $pgx->{genomeintervals}->[$_]->{start} <= $pgx->{referencebounds}->{$chro}->[1]  } @ind;
      @ind = grep{ $pgx->{genomeintervals}->[$_]->{end} >= $pgx->{referencebounds}->{$chro}->[0]  } @ind;

      foreach my $i (@ind) {

        my $start = $pgx->{genomeintervals}->[$i]->{start};
        my $end = $pgx->{genomeintervals}->[$i]->{end};

        if ($start < $pgx->{referencebounds}->{$chro}->[0]) {
          $start = $pgx->{referencebounds}->{$chro}->[0] }
        if ($end > $pgx->{referencebounds}->{$chro}->[1]) {
          $end = $pgx->{referencebounds}->{$chro}->[1] }

        my $seg_x0 = sprintf "%.1f", $area_x0 + $pgx->{basepixfrac} * ($start - $pgx->{referencebounds}->{$chro}->[0]);
        my $segPixEnd = sprintf "%.1f", ($seg_x0 + $pgx->{basepixfrac} * ($end - $start));
        # providing a minimum sub-pixel segment plot length
        if ($segPixEnd - $seg_x0 < 0.2) {
        	$segPixEnd = $seg_x0 + 0.2 }

        my $fill = frequencies2rgb(
			$pp,
			$frequencymapsSet->{intervals}->[$i]->{gain_frequency},
			$frequencymapsSet->{intervals}->[$i]->{loss_frequency},
			$maxF,
		);

        my $gdCol = $gdArea->colorAllocate( split(',', $fill) );
        $gdArea->filledRectangle($seg_x0, $gd_y0, $segPixEnd, $gd_yn, $gdCol);

      }
      
      $area_x0 += $areaW + $pp->{size_chromosome_padding_px};

    }
    
    my $labels_R = [];
    if ($frequencymapsSet->{labels}) { $labels_R = $frequencymapsSet->{labels} }

    # fallback color; defined here for alternation...
    $labCol = _alt_col($labCol, $defLabCol, $altLabCol);
    $pgx->svg_add_labels_right($labels_R, $pp->{size_strip_h_px}, $labCol);

    $pgx->{Y} += $pp->{size_strip_h_px};
    $gd_y0 += $pp->{size_strip_h_px};

  }

  $pgx->{svg} .= '
<image
  x="'.$pgx->{areastartx}.'"
  y="'.$pgx->{areastarty}.'"
  width="'.$pgx->{areawidth}.'"
  height="'.$areaH.'"
  xlink:href="data:image/png;base64,'.encode_base64($gdArea->png).'"
/>';

  return $pgx;

}

################################################################################

sub _alt_col {

	my $labCol = shift;
	my $defLabCol = shift;
	my $altLabCol = shift;
	
    if ($labCol eq $altLabCol) { 
      return $defLabCol }
    else { 
      return $altLabCol }
		
}

################################################################################


1;
