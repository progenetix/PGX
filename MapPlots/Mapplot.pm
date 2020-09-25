package PGX::MapPlots::Mapplot;

use Data::Dumper;
use PGX;
use PGX::GenomePlots::PlotParameters;
use PGX::IOUtilities::PGXfileReader;
use PGX::IOUtilities::PGXfileWriter;
use PGX::IOUtilities::PGXdataAggregation;
require Exporter;
@ISA    =   qw(Exporter);
@EXPORT =   qw(
  new
  pgx_get_web_geomap
);

########    ####    ####    ####    ####    ####    ####    ####    ####    ####
########    ####    ####    ####    ####    ####    ####    ####    ####    ####
########    ####    ####    ####    ####    ####    ####    ####    ####    ####


sub new {

  my $class     =   shift;
  my $args      =   shift;
  $args         =   args_modify_plot_parameters(PGX::read_plot_defaults(), $args);
  my $self      =   {
    parameters  =>  $args,
    plotid      =>  $args->{plotid},
    map         =>  q{},
  };

  bless $self, $class;
  return $self;

}

########    ####    ####    ####    ####    ####    ####    ####    ####    ####
########    ####    ####    ####    ####    ####    ####    ####    ####    ####
########    ####    ####    ####    ####    ####    ####    ####    ####    ####

sub pgx_get_web_geomap {
    
	my $pgx				=		shift;

	# stop processing now if nothing to do
	if ( scalar @{$pgx->{geomarkers}} < 1 ) {
			return $pgx }

	# get max marker size
	my @markerS		=		map{ $_->[3] } @{$pgx->{geomarkers}};
	my $markerMax	=		(sort {$b <=> $a} @markerS)[0];
	if ($markerMax < 1) {
		$markerMax 	= 	1 }
  my $locsizeF 	= 	( 50000000000 * $pgx->{parameters}->{map_marker_scale} / $markerMax );
  if (@{$pgx->{geomarkers}} < 2) {
  	$pgx->{parameters}->{map_marker_type}		=		'marker' }
  	
  my @markersJs;
  
  
  foreach my $marker (@{$pgx->{geomarkers}}) {

    my ($title, $lat, $lon, $size, $type) 	= 	@{$marker};

    if (! grep{ $type eq $_ } qw(circle marker)) {
    	$type			=		$pgx->{parameters}->{map_marker_type} }
    
  	if ($type eq 'marker') {
    	push @markersJs, qq!
L.$type([$lat, $lon]).bindPopup('$title').addTo(map)
        ! }
    else {
    	my $radius 	= 	sprintf "%.0f", sqrt($size / 3.14 * $locsizeF);
#print Dumper($title)."\n";
    	push @markersJs, qq!
L.$type([$lat, $lon], {
    stroke: true,
    color: '$pgx->{parameters}->{map_bubble_stroke_color}',
    weight: $pgx->{parameters}->{map_bubble_stroke_weight},
    fillColor: '$pgx->{parameters}->{map_bubble_fill_color}',
    fillOpacity: $pgx->{parameters}->{map_bubble_opacity},
    radius: $radius,
    count: $size
}).bindPopup('$title').addTo(map)
        !
      }  
  }
  
  my $_markersJs 	= 	'' . join(';', @markersJs) . '';

  $pgx->{map} = 	<< "__HTML__";
$pgx->{parameters}->{map_head}

<!-- map needs to exist before we load leaflet -->
<div id="map-canvas" style="width: $pgx->{parameters}->{size_plotimage_w_px}px; height: $pgx->{parameters}->{size_plotimage_h_px}px;"></div>

<!-- Make sure you put this AFTER Leaflet's CSS -->
<script src="https://unpkg.com/leaflet\@1.3.4/dist/leaflet.js"
      integrity="sha512-nMMmRyTVoLYqjP9hrbed9S+FzjZHW5gY1TWCHA5ckwXZBadntCNs8kEqAWdrb9O7rxbCaA4lKTIWjDXZxflOcA=="
      crossorigin=""></script>
<script>
  // Create the map.
  var map = L.map('map-canvas', { renderer: L.svg() } ).setView([$pgx->{parameters}->{map_latitude}, $pgx->{parameters}->{map_longitude}], $pgx->{parameters}->{map_zoom});

  L.tileLayer('$pgx->{parameters}->{map_tiles}', {
      minZoom: $pgx->{parameters}->{map_zoom_min},
      maxZoom: $pgx->{parameters}->{map_zoom_max},
      $pgx->{parameters}->{map_extra_JS}
      attribution: '$pgx->{parameters}->{map_attribution}'
  }).addTo(map);

  $_markersJs;
</script>
__HTML__

  return $pgx;

}


1;
