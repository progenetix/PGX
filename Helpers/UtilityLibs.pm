package PGX::Helpers::UtilityLibs;

require Exporter;
@ISA    =   qw(Exporter);
@EXPORT =   qw(
	MaxTextWidthPix
);

################################################################################

sub MaxTextWidthPix {

	my $texts			=		shift;
	my $fontSize	=		shift;

	return	( sort {$a <=> $b} (map{ length($_) } @$texts) )[-1] * $fontSize * 0.5;

}

1;
