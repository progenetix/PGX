package lib::AggregateData;

use Data::Dumper;
use lib::PlotMakeParameters;
use lib::Helpers;

use lib::ReadFiles;

require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(
  pgx_open_handover
  pgx_samples_from_handover
  pgx_segmentdata_from_sample_0
  pgx_add_variants_from_db
  pgx_create_samples_from_segments
  pgx_handover_pagination_range
  pgx_remap_vrsified_segments
  pgx_callset_labels_from_header
  pgx_callset_labels_from_biosamples
  pgx_callset_labels_from_file
  pgx_create_sample_collections
  pgx_geo_markers_from_provenance
  pgx_paginate_this
  pgx_samples_from_callsets
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

	$pgx->{datasetid} = $pgx->{handover}->{source_db};

	$pgx->pgx_handover_pagination_range();
	my $t_vs = $pgx->pgx_paginate_this($pgx->{handover}->{target_values});

	if ($pgx->{debug_mode} > 0) {		
		print Dumper("IDs from HANDOVER",$t_vs);
	}

	my $t_k = "id";

	if ($pgx->{handover}->{target_collection} eq 'biosamples') {
		# get the callset ids
		# $csIds = ...
		if ($pgx->{handover}->{target_key} eq "id") {
			$t_k = "biosample_id";
		} else {
			return $pgx;
		}
	} elsif ($pgx->{handover}->{target_collection} eq 'callsets') {
		$t_vs = $pgx->{handover}->{target_values};
		$t_k = $pgx->{handover}->{target_key};
	} else {
		return $pgx }

	my $query = { $t_k => { '$in' => $t_vs } };

	$pgx->pgx_samples_from_callsets($query);

	return $pgx;

}

################################################################################

sub pgx_samples_from_callsets {

	my $pgx = shift;
	my $query = shift;

	my $cscoll = $pgx->{dataconn}->get_collection( 'callsets' );
	my $cursor = $cscoll->find( $query );
	my $callsets = [ $cursor->all ];

	$callsets = [ grep{ exists $_->{cnv_statusmaps} } @$callsets ];

	if ($pgx->{debug_mode} > 0) {		
		print Dumper("QUERY", $query);
	}
# print Dumper($callsets);	

	$pgx->{samples} = [
		map{
			{
				id => $_->{id},
				biosample_id => $_->{biosample_id},
				cnv_statusmaps => $_->{cnv_statusmaps},
				cnv_statistics => $_->{cnv_statistics},
				paths => $_->{info}->{paths},
			}
		} @$callsets
	];

	if ($pgx->{debug_mode} > 0) {		
		print Dumper("FIRST SAMPLE", $pgx->{samples}->[0]);
	}

	return $pgx;

}

################################################################################

sub pgx_handover_pagination_range {

	my $pgx = shift;
	
	if (! $pgx->{handover}) { return $pgx }

	$pgx->{datasetid} = $pgx->{handover}->{source_db};

	my $skip_docs = $pgx->{parameters}->{skip} || 0;
	my $limit_docs = $pgx->{parameters}->{limit} || 0;

	$skip_docs *= $limit_docs;

	my $doc_no = $pgx->{handover}->{target_count};
	my $left_no = $doc_no - $skip_docs;

	if ($left_no < 1) {
		$pgx->{samples} = [];
		return $pgx;
	} elsif ($limit_docs > $left_no) {
		$limit_docs = $left_no;
	}

	my $range_i = $skip_docs + $limit_docs - 1;

	$pgx->{handover_pagination_range} = ([$skip_docs, $range_i]);
		
	return $pgx;

}

################################################################################

sub pgx_paginate_this {

	my $pgx = shift;
	my $this = shift || [];
	my $range = $pgx->{handover_pagination_range};

	if ($range->[0] > scalar @$this) {
		return [] }

	if ($range->[1] > scalar @$this) {
		$range->[1] = scalar @$this }

	return [ @$this[ $range->[0]..$range->[1] ] ]

}

################################################################################

sub pgx_add_variants_from_db {

	my $pgx = shift;
	my $query = shift;

	if (! $pgx->{dataconn}) { return $pgx }

	my $vcoll = $pgx->{dataconn}->get_collection('variants');

	for my $i (0..$#{ $pgx->{samples} }) {
		my $dataQuery = { 'callset_id' => $pgx->{samples}->[$i]->{id} };
		my $cursor = $vcoll->find( $dataQuery )->fields( { _id => 0, updated => 0, created => 0 } );
		$pgx->{samples}->[$i]->{variants} = [ $cursor->all ];
		$pgx->{samples}->[$i]->{variants} = $pgx->pgx_remap_vrsified_segments($pgx->{samples}->[$i]->{variants});
	}

	return $pgx;

}

################################################################################

sub pgx_segmentdata_from_sample_0 {

	my $pgx = shift;
	
	if (! $pgx->{samples}) { return $pgx }
	if (! $pgx->{samples}->[0]->{variants}) { return $pgx }
	if ($pgx->{segmentdata}) { return $pgx }
	
	$pgx->{segmentdata} = $pgx->pgx_remap_vrsified_segments($pgx->{samples}->[0]->{variants});
	
	return $pgx;
}

################################################################################

sub pgx_remap_vrsified_segments {

	my $pgx = shift;
	my $segmentdata = shift || [];
	
	my %refseqs = (
		"refseq:NC_000001.11" => "1",
		"refseq:NC_000002.12" => "2",
		"refseq:NC_000003.12" => "3",
		"refseq:NC_000004.12" => "4",
		"refseq:NC_000005.10" => "5",
		"refseq:NC_000006.12" => "6",
		"refseq:NC_000007.14" => "7",
		"refseq:NC_000008.11" => "8",
		"refseq:NC_000009.12" => "9",
		"refseq:NC_000010.11" => "10",
		"refseq:NC_000011.10" => "11",
		"refseq:NC_000012.12" => "12",
		"refseq:NC_000013.11" => "13",
		"refseq:NC_000014.9" => "14",
		"refseq:NC_000015.10" => "15",
		"refseq:NC_000016.10" => "16",
		"refseq:NC_000017.11" => "17",
		"refseq:NC_000018.10" => "18",
		"refseq:NC_000019.10" => "19",
		"refseq:NC_000020.11" => "20",
		"refseq:NC_000021.9" => "21",
		"refseq:NC_000022.11" => "22",
		"refseq:NC_000023.11" => "X",
		"refseq:NC_000024.10" => "Y"
	);
	
	my %varTypes = (
		"EFO:0030070" => "DUP",
		"EFO:0030067" => "DEL",
		"EFO:0030073" => "DUP", #AMP
		"EFO:0030069" => "DEL" #HOMODEL
	);
	
	for my $i (0..$#{ $segmentdata }) {	
		if ($segmentdata->[$i]->{location}) {
			my $loc = $segmentdata->[$i]->{location};
			my $efo = $segmentdata->[$i]->{variant_state};
			$segmentdata->[$i]->{variant_type} = $varTypes{ $efo->{id} };
			$segmentdata->[$i]->{reference_name} = $refseqs{ $loc->{sequence_id} };
			$segmentdata->[$i]->{start} = $loc->{interval}->{start}->{value};
			$segmentdata->[$i]->{end} = $loc->{interval}->{end}->{value};
			delete $segmentdata->[$i]->{location};
		}
	}
	
	return $segmentdata;

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
			cnv_statusmaps => $pgx->{cnv_statusmaps},
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
			statusmapsets => [ map{  { cnv_statusmaps => $_->{cnv_statusmaps} } } @{ $theseSamples } ],
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
	if ($pgx->{parameters}->{group_by} =~ /.../) {
		for my $grt (keys %{ $pgx->{config}->{datacollections} }) {
			if ($grt eq $pgx->{parameters}->{group_by}) {
				$groupType = $pgx->{parameters}->{group_by};
				$groupAttr = $pgx->{config}->{datacollections}->{$grt}->{samplefield};			
	}}}

	if ($groupAttr !~ /.../)  { $groupAttr = 'histological_diagnosis' }
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
		my $thisbioc = {};
		if (ref $thisBios->{$groupAttr} eq 'ARRAY') {  
			($thisbioc) = grep{ $_->{id} =~ /$groupType/ } @{ $thisBios->{$groupAttr} } }
		else {
			$thisbioc = $thisBios->{$groupAttr} }
		if ($thisbioc->{id} !~ /.../) { next }

		$pgx->{samples}->[$i]->{sortkey} = $thisbioc->{id};
		$pgx->{samples}->[$i]->{sortlabel} = ( $thisbioc->{label} =~ /.../ ? $thisbioc->{label} : $thisbioc->{id} );
		$pgx->{samples}->[$i]->{sortlabel} =~ s/^.+?\:\:?//;

	}

	return $pgx;

}

################################################################################

sub pgx_callset_labels_from_header {

	my $pgx = shift;

	if (! defined $pgx->{pgxfileheader}) { 
		return $pgx }

	my $fallbackK = 'NA';
	my $fallbackL = 'not specified';

	my $hv = $pgx->{pgxfileheader};
	for my $i (0..$#{ $pgx->{samples} }) {
		my $csId = $pgx->{samples}->[$i]->{id};
		if (grep{ /^group_id$/} keys %{ $hv->{samples}->{ $csId } } ) {
			$pgx->{samples}->[$i]->{sortkey} = $hv->{samples}->{ $csId }->{group_id} }
		if (grep{ /^group_label$/} keys %{ $hv->{samples}->{ $csId } } ) {
			$pgx->{samples}->[$i]->{sortlabel} = $hv->{samples}->{ $csId }->{group_label} }
		if (! defined $pgx->{samples}->[$i]->{sortkey}) {
			$pgx->{samples}->[$i]->{sortlabel} = $fallbackK }
		if (! defined $pgx->{samples}->[$i]->{sortlabel}) {
			$pgx->{samples}->[$i]->{sortlabel} = $fallbackL }	
	}

	return $pgx;

}

################################################################################

sub pgx_callset_labels_from_file {

	my $pgx = shift;
	my $sortFile = shift;

	if (! -f $sortFile)  { return $pgx }

	my $fallbackK = 'NA';
	my $fallbackL = 'not specified';

	for my $i (0..$#{ $pgx->{samples} }) {
		$pgx->{samples}->[$i]->{sortkey} = $fallbackK;
		$pgx->{samples}->[$i]->{sortlabel} = $fallbackL;
	}
	
	if ( defined $pgx->{pgxfileheader}) {
		my $hv = $pgx->{pgxfileheader};
		for my $i (0..$#{ $pgx->{samples} }) {
			my $csId = $pgx->{samples}->[$i]->{id};
			if (grep{ /^group_id$/} keys %{ $hv->{samples}->{ $csId } } ) {
				$pgx->{samples}->[$i]->{sortkey} = $hv->{samples}->{ $csId }->{group_id} }
			if (grep{ /^group_label$/} keys %{ $hv->{samples}->{ $csId } } ) {
				$pgx->{samples}->[$i]->{sortlabel} = $hv->{samples}->{ $csId }->{group_label} }
	}}


	# this assumes that the first column contains an entry for the selected id (or id)
	# the second then the associated label
	my $customSort = {};

	my ($header, $table) = read_file_to_split_array($sortFile);
	foreach (@$table) {
	  if ($_->[0] =~ /\w/ && $_->[1] =~ /\w/) {
		$customSort->{$_->[0]} = {
		  sortkey => $_->[1] =~ /\w/ ? $_->[1] : $fallbackK,
		  sortlabel => $_->[2] =~ /\w/ ? $_->[2] : $fallbackL,
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
