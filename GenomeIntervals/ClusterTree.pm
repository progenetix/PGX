package PGX::GenomeIntervals::ClusterTree;

use Data::Dumper;

require Exporter;
@ISA        =   qw(Exporter);
@EXPORT     =   qw(cluster_tree);

=pod

cluster_tree implements a way to retrieve an ordered tree structure including branch lengths
from a tree structure produced by Algorithm::Cluster::treecluster

The function will return 2 objects:

-	a list reference, where each item is a hash reference containing the X and Y coordinates for
	each node, which can then be used for drawing the dendrogram by assigning some pixel conversion
	factors:

	$_->{LX}, $_->{LY}	=>		|-----------|
												|
	$_->{NODEY}			=>					|---------------|
												|				|
	$_->{RX}, $_->{RY}	=>		|-----------|				|
																|-------- ...
												|			   ...
										$_->{NODEX}

	The corresponding coordinates e.g. in an SVG plot are then calculated like this (example assumes
	dendrogram on the right, with $dendroStart being the start X coordinate and $plotYstart being
	the, well, Y pixel coordinate of the top of the corresponding image/heatmap, and $plotHeightPix
	being the, well, plot area pixel height):

	my $maxTreeX	=	  max( map{ $_->{NODEX} } @{ $clusterTree } );
	my $xPixF		  =	  $totalTreeWidthPix / ($maxTreeX == 0 ? 1 : $maxTreeX);
	my $yPixF			=	 $plotHeightPix	/ (scalar(@{ $clusterTree }) + 1);

	foreach (@{ $clusterTree }) {

		my $xLpix0	=   sprintf "%.1f", $dendroStart + $_->{LX} * $xPixF;
		my $xRpix0	=   sprintf "%.1f", $dendroStart + $_->{RX} * $xPixF;
		my $xNdePix	=	  sprintf "%.1f", $dendroStart + $_->{NODEX} * $xPixF;
		my $yLpix		=   sprintf "%.1f", $yCorr + $_->{LY} * $yPixF;
		my $yRpix		=   sprintf "%.1f", $yCorr + $_->{RY} * $yPixF;

		$SVG				.=	'
	<line x1="'.$xStartLpix.'" y1="'.$yLpix.'" x2="'.$xNodePix.'" y2="'.$yLpix.'" stroke="#666666" />
	<line x1="'.$xStartRpix.'" y1="'.$yRpix.'" x2="'.$xNodePix.'" y2="'.$yRpix.'" stroke="#666666" />
	<line x1="'.$xNodePix.'" y1="'.$yLpix.'" x2="'.$xNodePix.'" y2="'.$yRpix.'" stroke="#666666" />';

	}

-	an index corresponding to the cluster matrix line numbers, which can be used to order the
	items used for clustering

=cut

sub cluster_tree {

	my $EisenTree	=	  shift;
	my $tree			=	  {};

	for (my $i = 0; $i < $EisenTree->length; $i++) {

		my $node		=   $EisenTree->get($i);
		my $nodeID	=   1 * (-1-$i);

		$tree->{ $nodeID }->{ID}		=	1 * $nodeID;
		$tree->{ $nodeID }->{LEFT}	= 1 * 	$node->left;
		$tree->{ $nodeID }->{RIGHT}	=	1 * $node->right;
		$tree->{ $nodeID }->{DIST}	=	1 * $node->distance;

		# making sure that nodes start with a leaf if there is one

		if (
			$tree->{ $nodeID }->{LEFT} =~ /^\-/
			&&
			$tree->{ $nodeID }->{RIGHT} !~ /^\-/
		) {
			$tree->{ $nodeID }->{LEFT}	  =	  $tree->{ $nodeID }->{RIGHT};
			$tree->{ $nodeID }->{RIGHT}	  =	  $node->left;
		}

	}

	my $nodes			=   [ sort keys %{ $tree } ];

	# declaring this as an empty list ref to have something blessed to submit ...
	my $sortNodes	=   [];
	$sortNodes		=	  _checkNode(
    									$nodes->[0],
    									$tree,
    									$nodes,
    									$sortNodes
    								);

	return	_drawNodes(
    				$tree,
    				$sortNodes,
    			);

}

################################################################################
# local subs ###################################################################
################################################################################

sub	_getParent {

	my ($nodeID, $cTree, $nodes)	=	@_;

	my $pNode		  =	(
		grep{
			$cTree->{ $_ }->{LEFT} eq $nodeID
			||
			$cTree->{ $_ }->{RIGHT} eq $nodeID
		} @$nodes
	)[0];
	return $pNode;

}

################################################################################

sub _addNode {

	my (
    $nodeID,
    $sortNodes
  )	            =   @_;

	unless ( grep{ $_ eq $nodeID } @$sortNodes ) {
		push(@$sortNodes, $nodeID) }

	return	$sortNodes;

}

################################################################################

sub _checkNode {

	my (
		$nodeID,
		$cTree,
		$nodes,
		$sortNodes,
	)					    =   @_;

	if (! $nodeID) { return	$sortNodes }
	if (! grep{ $_ eq $nodeID } @$nodes) { return	$sortNodes }

	if (
    (
		  $cTree->{ $nodeID }->{LEFT} !~ /^\-/
			||
			(grep{ $_ eq $cTree->{ $nodeID }->{LEFT} } @$sortNodes )
		)
		&&
		(
			$cTree->{ $nodeID }->{RIGHT} !~ /^\-/
			||
			(grep{ $_ eq $cTree->{ $nodeID }->{RIGHT} } @$sortNodes )
		)
	) {
		$sortNodes  =   _addNode($nodeID, $sortNodes);
		$sortNodes  =	  _checkNode(
        							_getParent(
                        $nodeID,
                        $cTree,
                        $nodes
                      ),
        							$cTree,
        							$nodes,
        							$sortNodes
        						);
	} else {

		if (
			$cTree->{ $nodeID }->{LEFT} =~ /^\-/
			&&
			(! grep{ $_ eq $cTree->{ $nodeID }->{LEFT} } @$sortNodes)
		) {
			$sortNodes    =	  _checkNode(
                          $cTree->{ $nodeID }->{LEFT},
                          $cTree,
                          $nodes,
                          $sortNodes
                        );
		}

		if (
			$cTree->{ $nodeID }->{RIGHT} =~ /^\-/
			&&
			(! grep{ $_ eq $cTree->{ $nodeID }->{RIGHT} } @$sortNodes)
		) {
			$sortNodes    =	  _checkNode(
                          $cTree->{ $nodeID }->{RIGHT},
                          $cTree,
                          $nodes,
                          $sortNodes
                        );
	}}

	return	$sortNodes;

}

################################################################################

sub _drawNodes {

  $DB::deep = 1000;
  
	my (
		$tree,
		$sortNodes,
	)					    =   @_;

	my $order	    =	  [];
	my $sampleI		=	  scalar(@$sortNodes) + 1;

	# not using pixels here; easy to multiply later by sample width (Y) or distance factor (X)

	foreach (@$sortNodes) {

		if ($tree->{ $_ }->{LEFT} =~ /^\-/) {

			$tree->{ $_ }->{LY}	=	1 * $tree->{ $tree->{ $_ }->{LEFT} }->{NODEY};
			$tree->{ $_ }->{LX}	=	1 * $tree->{ $tree->{ $_ }->{LEFT} }->{NODEX};

		} else {

			$sampleI	-=	1;
			$tree->{ $_ }->{LY}	=	1 * $sampleI;
			$tree->{ $_ }->{LX}	=	0;

			push(@$order, $tree->{ $_ }->{LEFT});

		}

		 if ($tree->{ $_ }->{RIGHT} =~ /^\-/) {

			$tree->{ $_ }->{RY}	=	1 * $tree->{ $tree->{ $_ }->{RIGHT} }{ NODEY };
			$tree->{ $_ }->{RX}	=	1 * $tree->{ $tree->{ $_ }->{RIGHT} }{ NODEX };

		} else {

			$sampleI	-=	1;
			$tree->{ $_ }->{RY}	=	1 * $sampleI;
			$tree->{ $_ }->{RX}	=	0;

			push(@$order, $tree->{ $_ }->{RIGHT});

		}

		$tree->{ $_ }->{NODEX}	=	1 * (sort {$a <=> $b } ($tree->{ $_ }->{LX}, $tree->{ $_ }->{RX}))[-1] + $tree->{ $_ }->{DIST};
		$tree->{ $_ }->{NODEY}	=	($tree->{ $_ }->{LY} + $tree->{ $_ }->{RY}) / 2;

	}

	return	(
    [ map{ $tree->{ $_ } } sort { $a <=> $b } keys %{ $tree } ],
    $order
  );

}

################################################################################

1;
