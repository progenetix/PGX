package PGX::IOUtilities::PGXdataAggregation;

use Data::Dumper;
use Math::Random qw(random_normal);
use PGX::GenomePlots::PlotParameters;
use PGX::Helpers::UtilityLibs;

require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(
  pgx_open_handover
  pgx_samples_from_handover
  pgx_segmentdata_from_sample_0
  pgx_add_variants_from_db
  pgx_create_samples_from_segments
  pgx_callset_labels_from_biosamples
  pgx_callset_labels_from_file
  pgx_create_sample_collections
  pgx_geo_markers_from_provenance
);

################################################################################

sub pgx_open_handover {

	my $pgx = shift;
	my $config = shift;
	my $accessid = shift;

	if ($accessid !~ /..../) { return $pgx }

	$pgx->{handover} = MongoDB::MongoClient->new()->get_database( $config->{handover_db} )->get_collection( $config->{handover_coll} )->find_one( { id => $accessid } );
	$pgx->{dataconn} = MongoDB::MongoClient->new()->get_database( $pgx->{handover}->{source_db} );

	return $pgx;

}

################################################################################

sub pgx_samples_from_handover {

	my $pgx = shift;

	if (! $pgx->{handover}) { return $pgx }
	if ($pgx->{handover}->{target_collection} ne 'callsets') { return $pgx }
	my $cscoll = $pgx->{dataconn}->get_collection( $pgx->{handover}->{target_collection} );
	my $dataQuery = { $pgx->{handover}->{target_key} => { '$in' => $pgx->{handover}->{target_values} } };
	my $cursor = $cscoll->find( $dataQuery )->fields( { _id => 1, id => 1, biosample_id => 1, info => 1 } );
	my $callsets = [ $cursor->all ];
	$callsets = [ grep{ exists $_->{info}->{statusmaps} } @$callsets ];

	$pgx->{datasetid} = $pgx->{handover}->{source_db};

	$pgx->{samples} = [
		map{
			{
				id => $_->{id},
				biosample_id => $_->{biosample_id},
				statusmaps => $_->{info}->{statusmaps},
				info => { cnvstatistics => $_->{info}->{cnvstatistics} },
				paths => $_->{info}->{paths},
			}
		} @$callsets
	];

	return $pgx;

}

################################################################################

sub pgx_add_variants_from_db {

	my $pgx = shift;
	my $query = shift;

	if (! $pgx->{handover}) { return $pgx }
	if (! $pgx->{dataconn}) { return $pgx }

	my $vcoll = $pgx->{dataconn}->get_collection('variants');

	for my $i (0..$#{ $pgx->{samples} }) {
		my $dataQuery = { 'callset_id' => $pgx->{samples}->[$i]->{id} };
		my $cursor = $vcoll->find( $dataQuery )->fields( { _id => 0, updated => 0, created => 0 } );
		$pgx->{samples}->[$i]->{variants} = [ $cursor->all ];
	}

	return $pgx;

}

################################################################################

sub pgx_segmentdata_from_sample_0 {

	my $pgx = shift;
	
	if (! $pgx->{samples}) { return $pgx }
	if (! $pgx->{samples}->[0]->{variants}) { return $pgx }
	if ($pgx->{segmentdata}) { return $pgx }
	
	$pgx->{segmentdata} = $pgx->{samples}->[0]->{variants};
	
	return $pgx;

}

################################################################################

sub pgx_create_samples_from_segments {

	my $pgx = shift;

	if (! $pgx->{segmentdata}) { return $pgx }

	my %csIds = map{ $_->{callset_id} => 1 } @{ $pgx->{segmentdata} };

	$pgx->{samples} ||= [];

	foreach my $csId (keys %csIds) {

		# a bit awkward re-assignment of segmentsdata
		my $segments = [ grep{ $_->{callset_id} eq $csId } @{ $pgx->{segmentdata} } ];	
		$pgx->segments_add_statusmaps($segments);
		push(
		  @{ $pgx->{samples} },
		  {
			id => $csId,
			statusmaps => $pgx->{statusmaps},
			variants => $segments,
			name => $csId,
		  }
		);
	}

	return $pgx;

}

################################################################################

sub pgx_create_sample_collections {

	no warnings 'uninitialized';

	my $pgx = shift;

	$pgx->{samplecollections} = [];

	my %sortKeys = map{ $_->{sortkey} => 1 } @{ $pgx->{samples} };

	# creation of the groups
	foreach my $sortKey (keys %sortKeys) {

	my @theirIndex = grep{ $pgx->{samples}->[$_]->{sortkey} eq $sortKey } 0..$#{ $pgx->{samples} };

	my $labelColor = random_hexcolor();
	if ( @theirIndex < $pgx->{parameters}->{min_group_no} ) {
	  $labelColor = '#cccccc' }

	my $label = {
	  label_text => $sortKey,
	  label_link => q{},
	  label_color => $labelColor,
	};

	$sortKey =~ s/^.+?\:\:?//;
	foreach (@theirIndex) {
	  $pgx->{samples}->[$_]->{labels} = [ $label ];
	  $pgx->{samples}->[$_]->{name} = $pgx->{samples}->[$_]->{id}.($sortKey =~ /.../ ? ' ('.$sortKey.')' : q{});
	  $pgx->{samples}->[$_]->{name} =~ s/^.+?\:\://g;
	  $pgx->{samples}->[$_]->{name} =~ s/^.+?\:\://g;
	}

	if ( @theirIndex < $pgx->{parameters}->{min_group_no} ) { next }

	my $theseSamples = [@{ $pgx->{samples} }[@theirIndex]];
	my $name = $theseSamples->[0]->{sortlabel}.' ('.scalar(@$theseSamples).')';

	push(
	  @{$pgx->{samplecollections}},
	  {
		labels => [ $label ],
		name => $name,
		statusmapsets => [ map{  { statusmaps => $_->{statusmaps} } } @{ $theseSamples } ],
	  },
	);

	}

	return $pgx;

}

################################################################################

sub pgx_callset_labels_from_biosamples {

	my $pgx = shift;
	my $config = shift;

	if (! $pgx->{dataconn}) { return $pgx }

	my ($groupAttr, $groupType);

	# CAVE: both "$config" (environment) and "$pgx->{config}" used here ...
	if ($config->{param}->{group_by}->[0] =~ /.../) {
		for my $gra (qw(biosubsets datacollections)) {
			for my $grt (keys %{ $pgx->{config}->{$gra} }) {
				if ($grt eq $config->{param}->{group_by}->[0]) {
					$groupType = $config->{param}->{group_by}->[0];
					$groupAttr = $pgx->{config}->{$gra}->{$grt}->{samplefield};			
	} } } }

	if ($groupAttr !~ /.../)  { $groupAttr = 'biocharacteristics' }
	if ($groupType !~ /.../)  { $groupType = 'xxxx' }

	my %biosIds = map{ $_->{biosample_id} => 1 } @{ $pgx->{samples} };
	my $bscoll = $pgx->{dataconn}->get_collection( 'biosamples' );
	my $cursor = $bscoll->find( { id => { '$in' => [ keys %biosIds ] } } )->fields( { id => 1, $groupAttr => 1 } );
	my $bioS = [ $cursor->all ];

	for my $i (0..$#{ $pgx->{samples} }) {

		my $csId = $pgx->{samples}->[$i]->{id};
		my $bsId = $pgx->{samples}->[$i]->{biosample_id};

		$pgx->{samples}->[$i]->{sortkey} = 'NA';
		$pgx->{samples}->[$i]->{sortlabel} = 'not specified';

		if ($bsId !~ /.../) { next }
		my ($thisBios) = grep{ $_->{id} eq $bsId } @$bioS;
		my ($thisbioc) = grep{ $_->{id} =~ /$groupType/ } @{ $thisBios->{$groupAttr} };
		if ($thisbioc->{id} !~ /.../) { next }

		$pgx->{samples}->[$i]->{sortkey} = $thisbioc->{id};
		$pgx->{samples}->[$i]->{sortlabel} = ( $thisbioc->{label} =~ /.../ ? $thisbioc->{label} : $thisbioc->{id} );
		$pgx->{samples}->[$i]->{sortlabel} =~ s/^.+?\:\:?//;

	}

	return $pgx;

}

################################################################################

sub pgx_callset_labels_from_file {

	my $pgx = shift;
	my $sortFile = shift;

	my $fallbackK = 'NA';
	my $fallbackL = 'not specified';

	for my $i (0..$#{ $pgx->{samples} }) {
		$pgx->{samples}->[$i]->{sortkey} = $fallbackK;
		$pgx->{samples}->[$i]->{sortlabel} = $fallbackL;
	}

	if (! -f $sortFile)  { return $pgx }

	# this assumes that the first column contains an entry for the selected id (or id)
	# the second then the associated label
	my $customSort = {};

	my $table = read_file_to_split_array($sortFile);
	foreach (@$table) {
	  if ($_->[0] =~ /\w\w/ && $_->[1] =~ /\w/) {
		$customSort->{$_->[0]} = {
		  sortkey => $_->[1] =~ /\w\w/ ? $_->[1] : $fallbackK,
		  sortlabel => $_->[2] =~ /\w\w/ ? $_->[2] : $fallbackL,
		};
	}}

	for my $i (0..$#{ $pgx->{samples} }) {

		my $csId = $pgx->{samples}->[$i]->{id};

		if ($customSort->{$csId}->{sortkey} !~ /.../) { next }

		$pgx->{samples}->[$i]->{sortkey} = $customSort->{$csId}->{sortkey};
		$pgx->{samples}->[$i]->{sortlabel} = $customSort->{$csId}->{sortlabel};

	}

	return $pgx;

}

################################################################################

sub pgx_geo_markers_from_provenance {

=pod

input format (other attributes optional):

[
	"provenance" : {
		"geo" : {
			"latitude" : 38.98,
			"longitude" : -77.1,
			"city" : "Bethesda",
			"country" : "United States",
			"continent" : "North America",
			"precision" : "city",
			"label" : "Bethesda, United States, North America"
		}
	},
	"counts" : {
		"ccgh" : 9,
		"acgh" : 0,
		"wes" : 0,
		"wgs" : 0
	},
	"id" : "10027410",
]

=cut

	my $pgx = shift;
	my $data = shift;
	my %locData = ();

	# the following re-assigns city etc. each time, assuming that there are
	# always values
	# the single key assignment is necessary so as not to lose the counts
	foreach my $this (@$data) {
    my $locKey = join(
			'',
			'la',
			$this->{provenance}->{geo}->{latitude},
			'lo',
			$this->{provenance}->{geo}->{longitude},
		);
		$locData{$locKey}->{label} = $this->{provenance}->{geo}->{label};
		$locData{$locKey}->{latitude} = $this->{provenance}->{geo}->{latitude};
		$locData{$locKey}->{longitude} = $this->{provenance}->{geo}->{longitude};
		$locData{$locKey}->{city} = $this->{provenance}->{geo}->{city};
		$locData{$locKey}->{country} = $this->{provenance}->{geo}->{country};
		$locData{$locKey}->{counts}->{items} += 1;
		push(@{$locData{$locKey}->{ids}}, $this->{id});
		foreach my $countKey (grep{ /.../ } keys %{ $this->{counts} }) {
			$locData{$locKey}->{counts}->{$countKey} += $this->{counts}->{$countKey};
		}
	}

	# compute marker data
	my $maxsamples = 1;

	foreach my $locKey ( keys %locData ) {

		if (
			$locData{$locKey}->{latitude} !~ /\d/
			||
			$locData{$locKey}->{longitude} !~ /\d/
			||
			$locData{$locKey}->{counts}->{items} < 1
		) { next }

		my $itemNo = $locData{$locKey}->{counts}->{items};
    my $circleSize = $itemNo;    # * $scalingF;
		if ($circleSize > $maxsamples) {
			$maxsamples = $circleSize }

    my @techNos;

		foreach ( grep{ ! /items/ } keys %{ $locData{$locKey}->{counts} } ) {
			if ( $locData{$locKey}->{counts}->{$_} > 0 ) {
				push(
					@techNos,
					$locData{$locKey}->{counts}->{$_} . ' ' . $_
				);
			}
		}

		my $infoText = $locData{$locKey}->{city} . ', '
			. $locData{$locKey}->{country}
			. ':<br />'
			. $locData{$locKey}->{counts}->{items}
			. ' item'
			. ( $locData{$locKey}->{counts}->{items} > 1 ? 's' : q{} );

		if (grep{ /\d/ } @techNos) {
			$infoText .= ' with '
			. join( ' and ', @techNos )
			. ' samples<ul>' }

		$infoText =~ s/\'/\\\'/g;
		if (grep{ /PMID/ } @{$locData{$locKey}->{ids}}) {
			$infoText .= '<ul>';
      foreach my $pmid (grep{ /PMID/ } @{$locData{$locKey}->{ids}}) {
				$pmid =~ s/pubmed\://;
				 $infoText .= '<li><a href="/cgi-bin/publications.cgi?id=PMID:'
					. $pmid
					. '" target="_blank">pubmedid '
					. $pmid
					. '</a></li>';
			}
			$infoText .= '</ul>';
		}

    push(
    	@{ $pgx->{geomarkers} },
    	[$infoText, $locData{$locKey}->{latitude}, $locData{$locKey}->{longitude}, $circleSize]
    );

  }

	return	$pgx;

}

################################################################################



1;
