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
	my $e_px = 16;
	my $t_w = $plotW * 0.8;
	
	my $e_l = length($error) * $e_px * 0.5;
	
	while ($e_l > $t_w) {		
		$e_px--;
		$e_l = length($error) * $e_px * 0.5;
	}
	
	$pgx->{svg} .= '
<text x="'.($plotW / 2).'" y="'.$pgx->{Y}.'" style="text-anchor: middle; font-size: '.$e_px.'px; fill: #dd3333">'.$error.'</text>';

	$pgx->{Y} += $e_px;
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
