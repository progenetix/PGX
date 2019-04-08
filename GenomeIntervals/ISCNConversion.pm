package PGX::GenomeIntervals::ISCNConversion;

use PGX::GenomeIntervals::CytobandReader;
use Data::Dumper;

require Exporter;
@ISA    =   qw(Exporter);
@EXPORT =   qw(
  deparse_karyotype
  deparse_cgh
);

sub deparse_karyotype {

=pod

Expects:


Returns:


=cut

########    ####    ####    ####    ####    ####    ####    ####    ####    ####

  my $pgx       =   shift;
	my $karyo     =   shift();

	my @deparsed_banding;

	# catching of odd annotations

	$karyo .= ',';
	$karyo =~ s/\.$//g;
	$karyo =~ s/[\–\—]/-/g;					# converts different dashes
	$karyo =~ s/[\s\"]//g;					# kills whitespace etc.
	$karyo =~ s/\{/\[/g;					  # sometimes a typing error (same key on American keyboard)
	$karyo =~ s/\}/\]/g;					  # sometimes a typing error (same key on American keyboard)
	$karyo =~ s/tt\(/t\(/g;					# typing error of double t was found in Mitelman database
	$karyo =~ s/(\([\w\.\;\-\~\?]*?)\,/$1\),/g;	# closing if forgotten
	$karyo =~ s/(\w)\:(\w)/$1;$2/g;
	$karyo =~ s/\[.*?\]/,/g;

	# topics one may treat differently
	$karyo = ','.$karyo;
	$karyo =~ s/,c(\[|,)/,$1/g;				# deletes the "c" of constitutional changes
	$karyo =~ s/\)c/\)/g;					    # deletes the "c" of constitutional changes
	$karyo =~ s/or/,/g;						    # indecision is converted to both
	$karyo =~ s/(\w)(\+\w)/$1,$2/g;		# lost comma is inserted; does not work for other lost commas (only the "+" context is fairly specific)
	$karyo =~ s/trpdup/dup/g;				  # trpdup found once => either or...
	$karyo =~ s/derdic/dic/g;				  # derdic found once => either or...
	$karyo =~ s/(X|Y)\/(X|Y)/$1-$2/g;	# XY number indecision is changed, due to "/" used as clone separator
	$karyo =~ s/^,[\d\-]+?(<[\dn]+?>)?,/,/;	# removes the chromosome number
	$karyo =~ s/(\<d.*?\>)/,/g;				# killing ploidy
	$karyo =~ s/^,[XY\-\/].*?,/,/;		# removes the XY annotation
	$karyo =~ s/(\w)\?/$1/g;				  # removes question marks from bands
	$karyo =~ s/\?([\w\;])/$1/g;			# removes question marks from bands
	$karyo =~ s/\;\?\??\)/;Z9)/g;			# replaces question marks which were
													# in place of chromosomes/bands to keep the pattern alive
#print qq($karyo\n);
	$karyo =~ s/,.*?mar.*?,/,/g;			# removes markers
	$karyo =~ s/,.*?mar.*?,/,/g;			# removes markers
	$karyo =~ s/,.*?dmin.*?,/,/g;			# removes dmins
	$karyo =~ s/inc//g;						    # incomplete label is removed

	$karyo =~ s/,trp(.*?),/,dup$1,dup$1,/ig;	# splits multiple occurrences

	$karyo =~ s/,([\w\.\;\-\(\)\+]*?)x(1\-)?2/,$1,$1,/ig;	      # splits multiple occurrences
	$karyo =~ s/,([\w\.\;\-\(\)\+]*?)x(2\-)?3/,$1,$1,$1,/ig;    # splits multiple occurrences
	$karyo =~ s/,([\w\.\;\-\(\)\+]*?)x(3\-)?4/,$1,$1,$1,$1,/ig;	# splits multiple occurrences

	$karyo =~ s/,([\w\.\;\-\(\)\+]*?)x(3\-)?4/,$1,$1,$1,$1,/ig;	# splits multiple occurrences

	$karyo =~ s/,der\(\w\d?\)(((add)|(inv)|t)[\w\.\;\-\(\)]*?)(((add)|(inv)|t)[\w\.\;\-\(\)]*?),/,$1,$5,/g;	# non-standard triplication of a chromosome

	$karyo =~ s/,hsr[^,]*?,/,/g;	# hsr cannot be interpreted

	$karyo .= ',';
# this is contentious: here, the minor karyotypes are killed
	$karyo =~ s/^(.*?)\/[\d\-\~]*?,idem(x\d)?,/$1,/g;
	$karyo =~ s/^(.*?)\/(.*?)$/$1/g;
	$karyo =~ s/^(.*?)\|\|(.*?)/$1/;

# conversion of special annotations ("Philadelphia chromosome" etc. ############################

	$karyo .= ',';

	$karyo =~ s/ph,/der(22)t(9;22)(q11;q34),/ig;
	$karyo =~ s/t\(9;22\),/t\(9;22\)\(q34;q11\),/g;		# assumption, because usually this is meant...
	$karyo =~ s/t\(14;18\),/t\(14;18\)\(q32;q21\),/g;
	$karyo =~ s/(t\(11;14\)),/$1\(q13;q32\),/g;
# c-myc/Ig
	$karyo =~ s/t\(8;14\),/t\(8;14\)\(q24;q32\),/g;
	$karyo =~ s/t\(8;22\),/t\(8;22\)\(q24;q11\),/g;
	$karyo =~ s/t\(2;8\),/t\(2;8\)\(p12;q24\),/g;

	# marker with 2 translocations, e.g. der(14)t(1;14)(q11;p11)t(8;14)(q24;q32), is split
	$karyo =~ s/,(der\(\w\d?\))(t\(\w*?\;\w*?\)\([\w\.]*?\;[\w\.]*?\))(t\(\w*?\;\w*?\)\([\w\.]*?\;[\w\.]*?\)),/,$1$2,$1$3,/g;

	$karyo =~ s/,(der\(\w\d?\))(t\(\w*?\;\w*?\)\(\w*?\;\w*?\))(del\(\w*?\)\([\w\.]+?\)),/,$1$2,$3,/g;

#	$karyo =~ s/,(\+?der\(\w\d?\))(t\([\w\;]*?\)\([\w\;]*?\))(t\([\w\;]*?\)\([\w\;]*?\)),/,$1$2,$1$3,/g;


	################################################################################################
	# end of pre-treatment #########################################################################
	################################################################################################

	my @remaining;

	foreach my $currentannotation (grep { /\w/ } split(/\,/, $karyo)) {

		my $foundmarker = 0;

		# single losses
		if ($currentannotation =~ /^\-(\w\d?)$/) {
			push(@deparsed_banding, $1.':-1');
			$foundmarker++;
		}

		# single gains
		if ($currentannotation =~ /^\+(\w\d?)$/) {
			push(@deparsed_banding, $1.':1');
			$foundmarker++;
		}

		# ring chromosomes
		if ($currentannotation =~ /^\+?r\((\w\d?)\)$/) {
			my $chro = $1;
			push(@deparsed_banding, '='.$chro.'pter=='.$chro.'qter=');
			if ($currentannotation =~ /\+/) {
				push(@deparsed_banding, $chro.':1');
			}
			$foundmarker++;
		}

		# simple translocations and derivatives
		if ($currentannotation =~ /^(\+?)(der\((\w\d?)\))?t\((\w\d?)\;(\w\d?)\)(\(([\w\.\-]*?)\;([\w\.\-]*?)\))?$/) {
			my $plus = $1;
			my $derchro = $3;
			my $chro1 = $4;
			my $chro2 = $5;
			my $band1 = $7;
			my $band2 = $8;
			my ($derband, $addchro, $addband);
			push(@deparsed_banding, '='.$chro1.$band1.'=='.$chro2.$band2.'=');
			if ($derchro =~ /\w/) {
				if ($chro1 eq $derchro) {
					$derband = $band1;
					$addchro = $chro2;
					$addband = $band2;
				} else {
					$derband = $band2;
					$addchro = $chro1;
					$addband = $band1;
				}
				if ($derband =~ /p\d/) {
					push(@deparsed_banding, $derchro.'pter-'.$derband.':-1');
				} elsif ($derband =~ /q\d/) {
					push(@deparsed_banding, $derchro.$derband.'-qter:-1');
				}
				if ($addband =~ /p\d/) {
					push(@deparsed_banding, $addchro.'pter-'.$addband.':1');
				} elsif ($addband =~ /q\d/) {
					push(@deparsed_banding, $addchro.$addband.'-qter:1');
				}

				if ($plus eq '+') {
					push(@deparsed_banding, $derchro.':1');
				}
			}
			$foundmarker++;
		}


		# three-way translocations
		if ($currentannotation =~ /^t\((\w\d?)\;(\w\d?)\;(\w\d?)\)(\(([\w\.\-]*?)\;([\w\.\-]*?)\;([\w\.\-]*?)\))?$/) {
			push(@deparsed_banding, '='.$1.$5.'=='.$2.$6.'=::='.$2.$6.'=='.$3.$7.'=');
			$foundmarker++;
		}

		# deletion of bands
		if ($currentannotation =~ /^(\+)?del\((\w\d?)\)\(([\w\.\-]*?)\)$/) {
			my $plus = $1;
			my $chro = $2;
			my $bands = $3;
			if ($bands =~ /([pq])([\d\.]+?)\-[pq]?([\d\.]+?)/) {
				push(@deparsed_banding, '='.$chro.$1.$2.'=='.$chro.$1.$3.'=');
			}
			push(@deparsed_banding, $chro.$bands.':DEL');
			if ($plus eq '+') {
				push(@deparsed_banding, $chro.':DUP');
			}
			$foundmarker++;
		}

		# duplication of bands
		if ($currentannotation =~ /^(\+)?dup\((\w\d?)\)\(([\w\.\-]*?)\)$/) {
			my $plus = $1;
			my $chro = $2;
			my $bands = $3;
			if ($bands =~ /([pq])([\d\.]+?)\-[pq]?([\d\.]+?)/) {
				push(@deparsed_banding, '='.$chro.$1.$3.'=='.$chro.$1.$2.'=');
			}
			push(@deparsed_banding, $chro.$bands.':1');
			if ($plus eq '+') {
				push(@deparsed_banding, $chro.':1');
			}
			$foundmarker++;
		}

		# "add" is treated like "del" to terminal band, since it results in deletion of a part of the
		# target chromosome, with addition of unknown origin
		if ($currentannotation =~ /^(\+)?add\((\w\d?)\)\(([\w\.\-]*?)\)$/) {
			my $plus = $1;
			my $chro = $2;
			my $bands = $3;
			push(@deparsed_banding, $chro.$bands.':-1');
			if ($plus eq '+') {
				push(@deparsed_banding, $chro.':1');
			}
			$foundmarker++;
		}

		# simple inversion
		if ($currentannotation =~ /^inv\((\w\d?)\)\(([pq][\d\.\-]*?)([pq][\d\.\-]*?)\)?$/) {
			push(@deparsed_banding, '='.$1.$2.'=='.$1.$3.'=::='.$1.$3.'=='.$1.$2.'=');
			$foundmarker++;
		}

		# isochromosome
		if ($currentannotation =~ /^(\+?)i\((\w\d?)(\)\()?([pq]).*?\)$/) {
			my $plus = $1;
			my $chro = $2;
			my $arm = $4;
			if ($plus eq '+') {
				if ($arm eq 'p') {
					push(@deparsed_banding, ($chro.'p:1', $chro.'p:1'));
				} elsif  ($arm eq 'q') {
					push(@deparsed_banding, ($chro.'q:1', $chro.'q:1'));
				}
			} else {
				if ($arm eq 'p') {
					push(@deparsed_banding, ($chro.'p:1', $chro.'q:-1'));
				} elsif  ($arm eq 'q') {
					push(@deparsed_banding, ($chro.'p:-1', $chro.'q:1'));
				}
			}
			$foundmarker++;
		}

		# centromeric fusion
		if ($currentannotation =~ /^(\+?)der\((\w\d?)\;(\w\d?)\)(\(([pq])(1\d?)?\;([pq])(1\d?)?\))?$/) {
			my $plus = $1;
			my $chro1 = $2;
			my $chro2 = $3;
			my $arm1 = $4;
			my $arm2 = $5;
			push(@deparsed_banding,  '='.$chro1.$arm1.'11=='.$chro2.$arm2.'11=');
			if ($plus eq '+') {
				if ($arm1 eq 'p') {
					push(@deparsed_banding, $chro1.'p:1');
				} elsif ($arm1 eq 'q') {
					push(@deparsed_banding, $chro1.'q:1');
				}
				if ($arm2 eq 'p') {
					push(@deparsed_banding, $chro2.'p:1');
				} elsif ($arm2 eq 'q') {
					push(@deparsed_banding, $chro2.'q:1');
				}
			} else {
				if ($arm1 eq 'p') {
					push(@deparsed_banding, $chro1.'q:-1');
				} elsif  ($arm1 eq 'q') {
					push(@deparsed_banding, $chro1.'p:-1');
				}
				if ($arm2 eq 'p') {
					push(@deparsed_banding, $chro2.'q:-1');
				} elsif  ($arm2 eq 'q') {
					push(@deparsed_banding, $chro2.'p:-1');
				}
			}
			$foundmarker++;
		}

		# dicentric chromosome
		if ($currentannotation =~ /^(\+?)dic\((\w\d?)\;(\w\d?)\)\(([pq][\d\.\-]*?)\;([pq][\d\.\-]*?)\)$/) {
			my $plus = $1;
			my $chro1 = $2;
			my $chro2 = $3;
			my $band1 = $4;
			my $band2 = $5;
			push(@deparsed_banding, '='.$chro1.$band1.'=='.$chro2.$band2.'=');
			if ($plus eq '+') {
				if ($arm1 eq 'p') {
					push(@deparsed_banding, $chro1.$band1.'qter:1');
				} elsif ($arm1 eq 'q') {
					push(@deparsed_banding, $chro1.'pter'.$band1.':1');
				}
				if ($arm2 eq 'p') {
					push(@deparsed_banding, $chro2.$band2.'qter:1');
				} elsif ($arm2 eq 'q') {
					push(@deparsed_banding, $chro2.'pter'.$band2.':1');
				}
			} else {
				if ($arm1 eq 'p') {
					push(@deparsed_banding, $chro1.'pter'.$band1.':-1');
				} elsif  ($arm1 eq 'q') {
					push(@deparsed_banding, $chro1.$band1.'qter:-1');
				}
				if ($arm2 eq 'p') {
					push(@deparsed_banding, $chro2.'pter'.$band2.':-1');
				} elsif  ($arm2 eq 'q') {
					push(@deparsed_banding, $chro2.$band2.'qter:-1');
				}
			}
			$foundmarker++;
		}

		# for logging
		if ($foundmarker == 0) {
			push (@remaining, "not found: $currentannotation");
		} elsif ($foundmarker > 1) {
			push (@remaining, "multiple ($foundmarker) matches");
		}
	}

	my $cna       =   {};
	my $index     =   0;

#	foreach my $CNAtype (qw(DUP DEL)) {
#		foreach my $change (grep{ /\:$CNAtype$/ } @deparsed_banding) {
#			$change   =~ s/\:\w\w\w$//;
#			$index++;
#			$cna->{ $index }  = {
#				CYTOREGION  => $change,
#				CNATYPE     => $CNAtype,
#			};
#			(
#				$cna->{ $index }->{ CHRO },
#				$cna->{ $index }->{ BASESTART },
#				$cna->{ $index }->{ BASESTOP },
#			) = split(/[\-\:]/, getGP4interval($change));
#		}
#	}
#
#	foreach my $change (grep{ /=/ }  @deparsed_banding) {
#		$index++;
#		$cna->{ $index } = {
#			CYTOREGION => $change,
#			CNATYPE => 'STRUCTURAL',
#			CNSTATUS => 0,
#		};
#	}

	return $pgx;

}

1;
