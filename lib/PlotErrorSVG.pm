package lib::PlotErrorSVG;

use Data::Dumper;
require Exporter;

@ISA = qw(Exporter);
@EXPORT = qw(return_error_svg);

################################################################################

sub return_error_svg {

	no warnings 'uninitialized';

	my $pgx = shift;
	my $error = shift;
	
	$pgx->{svg} = q{};

	$pgx->{Y} = $pgx->{parameters}->{size_plotmargin_top_px};
	my $plotW = $pgx->{parameters}->{size_plotimage_w_px};

	$pgx->{areastartx} = $pgx->{parameters}->{size_plotmargin_px};
	$pgx->{areawidth} = $plotW - $pgx->{areastartx} - $pgx->{parameters}->{size_plotmargin_px};
	$pgx->{basepixfrac} = ( $pgx->{areawidth} - ($#{ $pgx->{parameters}->{chr2plot} } * $pgx->{parameters}->{size_chromosome_padding_px}) ) / $pgx->{genomesize};

	my $e_px = 16;
	my $t_w = $pgx->{areawidth} * 0.9;

	$pgx->svg_add_cytobands();
	
	my $e_l = length($error) * $e_px * 0.5;
	
	while ($e_l > $t_w) {		
		$e_px--;
		$e_l = length($error) * $e_px * 0.5;
	}
	
	$pgx->{Y} += $pgx->{parameters}->{size_plotarea_padding};
	my $p_a_h = $e_px + $pgx->{parameters}->{size_plotarea_padding} * 3;

	$pgx->{svg}  .= '
<rect x="'.$pgx->{areastartx}.'" y="'.$pgx->{Y}.'" width="'.$pgx->{areawidth}.'" height="'.$p_a_h.'" style="fill: #fec; fill-opacity: 0.6; " />';

	$pgx->{Y} += $p_a_h - $pgx->{parameters}->{size_plotarea_padding} * 2;

	$pgx->{svg} .= '
<text x="'.($plotW / 2).'" y="'.$pgx->{Y}.'" style="text-anchor: middle; font-size: '.$e_px.'px; fill: #dd3333">'.$error.'</text>';

	$pgx->{Y} += $pgx->{parameters}->{size_plotarea_padding};
	$pgx->{markerstarty} = $pgx->{Y};
	$pgx->svg_add_markers();
	$pgx->{Y} += $pgx->{parameters}->{size_plotmargin_top_px};
	
	my $plotH = sprintf "%.0f", $pgx->{Y};
	$plotW = sprintf "%.0f", $plotW;
	
	$pgx->{svg} = '<svg
xmlns="http://www.w3.org/2000/svg"
xmlns:xlink="http://www.w3.org/1999/xlink"
version="1.1"
id="error"
width="'.$plotW.'px"
height="'.$plotH.'px"
style="margin: auto; font-family: Helvetica, sans-serif;">
<rect x="0" y="0" width="'.$plotW.'" height="'.$plotH.'" style="fill: '.$pgx->{parameters}->{color_plotbackground_hex}.'; " />
'.$pgx->{svg}.'
</svg>';

	return $pgx;

}

1;
