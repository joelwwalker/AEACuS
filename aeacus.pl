#!/usr/bin/perl

#*******************************#
# aeacus.pl Version 3.034	#
# March '11 - July '21		#
# Joel W. Walker		#
# Sam Houston State University	#
# jwalker@shsu.edu		#
# arXiv:1207.3383		#
# Copy: GNU Public License V3	#
#*******************************#

=pod
* OUT issues + sorting, especially MET -> MPX MPY MAP
* external / multiple + interleave ...
tree implementation & jets
issue of undef propagation
MT2 validation
weights, mixing, uncertainties, correlation matrix
=cut

# Require minimal perl version and specify AEACuS package version
{; package Local::AEACuS; require 5.008_009; our ($VERSION) = 3.033; }

# Apply a strict coding pragma and define numerical constants
use strict; use sort q(stable); use constant +{
	NIL => 10**-8, EPS => 0.0005, ONE => 1.0, BIG => 500, HUG => 10**+8, TRY => 25, PRT => 12,
	TUP => 64, RPT => 12, LPR => 100, BLK => 6, BMX => 52, SET => 10**+4, CAP => 250_000,
	EXP => qr/[-+]?(?:\d+\.\d*|\d*\.\d+|\d+)(?:E[-+]?\d+)?/i, IMX => ( ~0 >> 1 ), PIE => 4*( atan2 (1,1)),
	do { my ($inf) = 9**9**9; ( INF => ($inf), NIN => ( -1*($inf)), NAN => (($inf)/($inf))) }};

# Register named and ordered parameter options from the command line input
our ($OPT) = ( &LOAD_OPTS(@ARGV));

# Begin local enclosing scope for main program routine; bypass if SUBROUTINE library is included by external program
if ( $0 =~ /\/aeacus\.pl$/ ) {

# Read event filtering, object reconstruction and event selection specifications from cardfile
my ($crd) = map { (/^(.*\/)?([^\/]*?)(?:\.dat)?$/); my ($crd,$err,$fil) =
	( &LOAD_CARD([ (($1) or ( q(./Cards/))), [ ((defined) ? ($2) : ( q(cut_card))), q(), q(.dat) ]]));
	($crd) or ((defined) ? ( die 'Cannot open card file for read' ) : (exit));
	( die ( join "\n", ( 'Malformed instruction(s) in card '.($$fil[0].$$fil[1]),
		( map {( "\t".'* Line '.$$_[0].':'."\t".'>> '.$$_[1].' <<' )} (@$err)), q()))) if (@$err);
	($crd) } ( &$OPT( q(crd)));
# If the CRD parameter includes the '/' character then the path is treated literally; it otherwise defaults to './Cards/'
# If the CRD parameter was undefined then the default cut card file name is 'cut_card'; otherwise it must be given explicitly
# The cut card file must use the extension '.dat', but it is not necessary to write the extention explicitly in the CRD parameter
# Program terminates if the card file cannot be found; an exception is raised only if the CRD parameter was defined
# Program terminates with exception if lexical errors exist in the card file

# Establish the event statistic and cut flow sequences and Curry the primary analysis function
my ($flw,$stc,$paf) = ((%$crd) ? ( &ANALYSIS_CODE($crd)) : ());

# Establish primary file name base and cross section specifications
my ($fil,$xsc,$err,$pre,$out,$ipb,$cap,$lhc,$r2l) = map {(
	( do { (/^(.*\/)?([^\/*]*)$/) or ( die 'Invalid file specification' ); [[ (($1) or q(./Events/)), ($2) ]] } ),
	( map {( &$OPT($_))} ( qw( xsc err pre out ipb )))
	)} grep {(defined)} ( &$OPT( q(fil)));
# If the FIL parameter includes the '/' character then the path is treated literally; it otherwise defaults to './Events/' if FIL is defined
# Wildcards are not permissible in an explicit singular filename at this level
# If FIL was defined then the cross section and its error [pb], the tag prefix, the output path, and the target luminisity [/pb] are read

# Establish global event channel management specifications and locate (possibly compressed) .lhco or .root event file(s)
my (@fil); my (@chn) = ( map { my ($i,$chn) = (($_), (($$crd{chn}[$_]) or (($_ > 0) ? (undef) : +{} ))); !($chn) ? () : do {
	# The zeroth channel, wherein statistics are computed from .lhco events and output to .cut files, runs automatically
	# Positively numbered channels, wherein sorting is performed relative to previously computed statistics, short circuit if not elaborated in the card
	my ($dir,$mrg) = q().( &DEFINED(($$chn{dir}[0]), (($i > 0) ? ($out) : ( q(./Events/)))));
	# The active directory DIR is read from the card file; if undefined it defaults to './Events/' for channel zero, or otherwise to the output directory
	my ($fil) = ((($i == 0) && ($fil)) or ( grep {(@$_)}
		map {[ map {[ $dir, $_ ]} (($fil) ? () : ((@$_) or ($i > 0)) ? (@$_) : (($mrg=1),(undef))[1] ) ]}
		[ map { (defined) ? qq($_) : () } (@{$$chn{fil}||[]}) ] ));
	# The active file is the command line file, including its bundled directory, if it exists and the zeroth channel is processing
	# The active file is undefined for numbered channels if the command line file exists or if the array of defined FIL card entries is empty
	# If that array is empty for the zeroth channel then "merge" mode is triggered with an undefined file name tied to the active directory
	# Otherwise, the non-empty array of defined FIL card entries is bundled with the active directory
	# FOR THE ZEROTH CHANNEL:
	(($i == 0) ? ( do { ($mrg) && do {
			$cap = ( map {((defined) ? 0+(($_ >= 0) && ((int) or (CAP))) : (undef))} ($$chn{cap}[0]))[0];
			# If merging, the file size event CAP is read (UNDEF for 1-to-1 file mapping,
			# or negative for unlimited single-file merging; 0 defaults to 250,000 events per file)
			$lhc = [ (( map {((defined) ? ( qq($_)) : ())} ($$chn{lhc}[0])), ($dir)) ];
			# Also, the LHC merged Olympics directory is read (defaults to active)
			do { $r2l = +{ wgt => (1,1,!1)[$$_[0]], gen => (!1,1,!1)[$$_[1]], fat => (!1,1,!1)[$$_[2]] }}
				for ([ map {( $$_[0] <=> 0 )} @$chn{( qw( wgt gen fat ))} ]);
			# Also, the root2lhco conversion parameters WGT, GEN, FAT, and AUX (defaulting to True, False, False, and Empty) are read if merging
			$$r2l{aux} = [ grep { s/\s//g; 1 }
				map {( m/^((?!\.)(?:(?:^|\.)\s*[A-Za-z][A-Za-z0-9]*(?:\s*(?:\[\s*\d+\s*]|\(\s*\)))?\s*)+)$/ )}
				map {( qq($_))} @{$$chn{aux}||[]} ]; };
		$pre = ( map { (/^([A-Za-z][\w-]*?)_*$/) ? ( qq($1).q(_)) :
			do { (defined) && ( print STDERR 'INVALID SELECTION PREFIX '.$_."\n" ); (undef) }} ( &DEFINED( $pre, $$chn{pre}[0] )))[0];
		# The prefix PRE is read; it is leading alpha plus word characters and dashes
		# A (consolidated) trailing underscore is appended if defined (defaults to empty)
		$out = (( &Local::FILE::PATH(( $out = q().( &DEFINED( $out, $$chn{out}[0], q(./Cuts/)))), 2 )) or ( die 'Cannot write to directory '.$out ));
		# The OUT directory is read, defaulting to './Cuts/' if undefined, and created (recursively) if it does not exist
		$ipb = ( &DEFINED( $ipb, ( do { my ($lum,$ipb) = [ @$chn{( qw( ipb ifb iab izb iyb ))} ]; for my $i (0..4) {
			($ipb) = (${(($$lum[$i]) or (next))}[0])*((10)**(3*$i)); (last) } (( &MAX(0,$ipb)) or (undef)) } )));
		# The positive luminosity scale is read in a known unit ( IPB, IFB, IAB, IZB, IYB ) or left undefined
		my (%fil); do { push @{ $fil{( shift @$_ )}[( shift @$_ )] }, $_ } for grep { pop @$_; pop @$_; 1 } values %{{
			map {(($$_[2].$$_[4]) => ($_))} sort { our ($a,$b); ( $$a[5] <=> $$b[5] ) } (
			( map {(( $$_[1] !~ /_uncleaned_events\.lhco(?:\.gz)?/ ) && ( $$_[1] =~ /^((.*?)(?:_000)*((?:_\d{3})*))\.lhco(\.gz)?$/ ) ?
				[ qq($2), (0+ !(length $3)), @$_[0,1], qq($1), (!$4)*(1<<0) ] : ())}
				map {( &Local::FILE::LIST($_))} (@{$fil||[]})),
			( map {(( $$_[1] =~ /^(([A-Za-z][\w-]*?)?(?:_+\d+)*)_+(?:delphes|(pgs))_events\.(?:lhco(\.gz)|(root))?$/ ) ? (($5) && ($3)) ? () :
				[ ((length $2) ? ( qq($2)) : ( q(TAG))), (2), @$_[0,1], qq($1), ((!$3)*(1<<2) + (!$5)*(1<<1) + (!$4)*(1<<0)) ] : ())}
				(($mrg) ? ( &Local::FILE::LIST($$fil[0],1)) : ()))) }};
		# Any .lhco[.gz] files (except '_uncleaned_events') in the active dir consistent with fil are matched for a key and a strict index
		# If merging, the fil matching criterion is guaranteed to accept all files consistent with basic naming requirements of the FILE class
		# Also, _events.lhco[.gz] and delphes_events.root files one dir below are matched for a key (defaults to TAG) and a non-strict index
		# If a key+index string repeats in a directory, Delphes files are valued over PGS files, .lhco over .root, and unzipped over zipped
		# Located files retain dir and fil info, along with the key and a flag for lower-level (2) or upper-level w/wout (0,1) an explicit index
		# Only the zon-zero set of files having the lowest numbered flag is stored for each key
		(@fil) = ( map { my ($k) = $_; ( grep {($$_[2])} map {[ $k, $_, $fil{$k}[$_]]} (0..2))[0] }
			( sort { our ($a,$b); ( $a cmp $b ) } keys %fil )); () } ) :
	# FOR NON-ZERO CHANNELS:
	( map { +{ cid => ($i), fil => ($fil), esc => [ map {[ ($_ < 0), ( abs ), ( @{ (($_) && ( $$crd{esc}[( abs )] )) or
		do { print STDERR 'INVALID EVENT SELECTION CUT SPECIFICATION IN CHANNEL '.$i."\n"; +{}}}{( qw( key cut ))} ) ]} (@$_) ] }}
		grep {(@$_)} [ map {((defined) ? ( int ) : ())} (@{$$chn{esc}||[]}) ] )) }} ( 0, (1..(@{$$crd{chn}||[]}-1))));
		# Processing proceeds if there are non-empty ESC specifications
		# A hash is pushed onto the channels list with the channel index, the active file, and a list of processed event selection cuts
		# Stored here is the sign (negative to reject or positive to accept), the cut index, and the referenced key/cut pair from the cardfile

# LOOP over LINEs from .lhco input FILEs and perform sequential mode channel filtering on processed .cut files
LOOP: while ( my ($k,$m,$f) = @{(( shift @fil ) || [] )} ) {
	use Fcntl qw(:seek); my (@FHT,@xsc); (($m == 2) or ($paf) or (next LOOP));
	# Processing starts with any fil sets associated with the zeroth channel
	# The loop advances unless merging ( flag 2 ) or the primary analysis function exists
	# An inner loop processes for each file associated with the given key
	FILE: for my $fil ( sort SORT_LIST_ALPHA (@$f)) { my ($lns,$cut);

		# A new output group is initialzed if none exists or if merging with 1-to-1 file mapping
#		if ( !(@FHT) or (($m == 2) && ( defined $cap ))) { push @FHT, [ 0, []]; }

		# A conversion from .root files to the .lhco format is attempted, if applicable
		if (($m == 2) && ( $$fil[1] =~ /\.root$/ )) { ( &ROOT_2_LHCO($fil,$r2l)) or
			do { print STDERR 'CANNOT CONVERT TO LHCO FROM FILE '.$$fil[0].$$fil[1]."\n"; ( next FILE ) };
			( $$fil[1] =~ s/\.root$/.lhco.gz/ ); }
		my ($FHI) = ( &Local::FILE::HANDLE($fil)) or
			do { print STDERR 'CANNOT READ FROM FILE '.$$fil[0].$$fil[1]."\n"; ( next FILE ) };
		# The cross-section record associated with the key is initialized with the command line value or the Merged_XSEC.TXT value or UNDEF
		push @xsc, [ my ($xsc) = (( defined $xsc ) ? ($xsc) : ( &MERGED_XSEC($fil))) ];
		# The cross-section record associated with the key is initialized with the command line value or the SUMMARY.TXT value or UNDEF
		push @xsc, [ my ($xsc) = (( defined $xsc ) ? ($xsc) : ( &SUMMARY_XSEC($fil))) ];
		# The cross-section record associated with the key is initialized with the command line value or the Pythia.LOG value or UNDEF
		push @xsc, [ my ($xsc) = (( defined $xsc ) ? ($xsc) : ( &PYTHIA_XSEC($fil))) ]; local ($_);
		# An inner loop processes each line in the current file
		LINE: while (<$FHI>) {
			# Any proprietary or MadGraph cross-section line is appended to the list unless the initial file-level cross section line was defined
			( m/^\s*#*\s*(${\EXP})\s+PB\s+CROSS\s+SECTION/i ) ? do { ( push @xsc, [ (0+ $1) ] ) unless ( defined $xsc ) } :
			( m/^\s*#*\s*<MGGenerationInfo>/i ) ? do { ( push @xsc, [ ( &MADGRAPH_XSEC($FHI))[0,1]] ) unless ( defined $xsc ) } :
			# The LHCO formatted events are read
			( m/^\s*(\d+)\s+(\d+)/ ) && do {
				($1 eq q(0)) ? ( $lns = [ $_ ] ) :
				($2 ne q(6)) ? ( push @{$lns||[]}, $_ ) :
				($lns) && do {
					# An event number, trigger bit mask, and partial weight are extracted from the event header
					my ($nbr,$trg,$wgt) = ( map {((length) ? (0+ $_) : (undef))}
						( $$lns[0] =~ m/^\s*0\s+(\d+)(?:(?:\s+(\d+))(?:\s+(${\EXP}))?)?/ ));
					if ( !(@FHT) or (($m == 2) && ($cap) && ( $FHT[-1][1] == $cap ))) { push @FHT, [
						grep {((defined) or ( die 'Cannot open temporary file for read/write' ))}
							( &Local::FILE::HANDLE()) ]; }
					# If merging, events are copied to a set of sized temporary file handles
					if ($m == 2) { print "\n", @$lns, $_; }
					# Otherwise, if processing, the MET and event objects are run through the primary analysis function
					else {
						(( my ($xit,$vls) = (($paf)->(
							[ $nbr, $trg, $wgt ],
							( scalar &LORENTZ(( shift @{( &EVENT_OBJECTS($_) || [] )} ),1,1,!1)),
							( scalar &EVENT_OBJECTS($lns))))) or ( next LINE ));
						# A counter corresponding to the level of exit from the cut flow (or -1 for survival) is incremented
						${ $cut ||= [ map {[0,!1]} (0..(0+@{$flw||[]})) ] }[$xit][0]++;
						# Computed statistics for surviving events are printed to the temporary filehandle
						if ($xit < 0) { print q().( join q( ), map {( sprintf (( shift @$_), ( shift @$_ )))}
							( [ q(%07.7i), 1+$FHT[-1][1]], @$vls ))."\n"; }
						# A correlation matrix is computed for failing events with multiple statistics at a given flow
						elsif (@$vls > 1) { $$cut[$xit][1] += ( &Local::VECTOR::OUTER_PRODUCT($vls,$vls)); }}
					# The event count and event count applicable to the active cross-section are incremented
					$FHT[-1][1]++; $xsc[-1][2]++; ( undef $lns ); }}}
		if ($m == 2) {}
		# After completion of the loop over file lines, processing continues if not in merging mode
		else {
			# Processed .cut files are generated with a header, prepending any prefix string, in the OUT directory
			my ($FHT,$hdr,$xsc) = map {(( grep {(( seek $_, 0, SEEK_SET ) or ( die 'Cannot rewind temporary file' ))} ( shift @$_ )),
				( &FORMAT_HEADER( 1, [ ( shift @$_ ), \@xsc, $ipb ], ($cut,$flw,$stc))))} (( shift @FHT ) || ()); ( undef @FHT ); ( undef @xsc );
			# Cut files are run through higher order channel filters for channels with no files specified in the card (forced for command-line files) 
			do { my ($f) = ($_); for my $c (@chn) { for ( @{ ( &AUXILIARY_CHANNEL($c,$f)) || [] } ) {( print STDERR $_."\n" )}}} for
				map { my ($h,$f) = (( shift @$_ ), ( shift @$_ )); print $hdr; local ($_); while (<$FHT>) {( print )}; print "\n"; ( close $h );
					(( defined $xsc ) or ( print STDERR 'CANNOT ESTABLISH CROSS SECTION FOR FILE '.$$f[0].$$f[1]."\n" )); ($f) }
				grep {((@$_) or ( die 'Cannot open file in directory '.$out.' for write' ))}
				# Non-indexed files ( flag 1 ) are incremented relative to existing files for non-destructive output
				map {[ ($m == 1) ? ( &Local::FILE::NEXT($_)) : ( &Local::FILE::HANDLE( $_, [1,1,1,0] )) ]}
				map {[ ($out), [ ($pre.$$_[0]), ($$_[1]), ( q(.cut)) ]]}
				(((defined $hdr) && ( &Local::FILE::KEYS($$fil[1]))) or ()); }}
	# After completion of the loop over key files, processing continues if in merging mode
	if ($m == 2) {
		# Each of the temporary file handles, along with a new header, is dumped to a suitably incremented file name in the LHC directory
		my ($lhc,$lnk) = map {(( &Local::FILE::PATH($_,2)) or ( die 'Cannot write to directory '.$_ ))} (@$lhc);
		# The corresponding file is appended to the list for processing with flag 0, as a numbered lower-level file
 		my ($n) = [ $lhc, [ $k, 0, q(.lhco) ]]; push @fil, [ $k, 0, [ map {
			my ($FHT,$hdr,$xsc) = (( grep {(( seek $_, 0, SEEK_SET ) or ( die 'Cannot rewind temporary file' ))} ( shift @$_ )),
				( &FORMAT_HEADER( 1, [ ( shift @$_ ), \@xsc ] ))); (
			# A linked copy of this file is stored in ./Events to document the assigned file indexing and block regeneration
			map { my ($h,$f) = (( shift @$_ ), ( shift @$_ )); print $hdr; local ($_); while (<$FHT>) {( print )}; print "\n"; ( close $h );
				$n = [ @$f ]; if ($lnk) { ( link ( $lhc.$$f[1], $lnk.$$f[1] )); ( $$f[0] = $lnk ) }
				(( defined $xsc ) or ( print STDERR 'CANNOT ESTABLISH CROSS SECTION FOR FILE '.$$f[0].$$f[1]."\n" )); ($f) }
			grep {(( defined $$_[0] ) or ( die 'Cannot open file in directory '.$lhc.' for write' ))}
				(( defined $hdr ) ? [ &Local::FILE::NEXT($n) ] : ())) } (@FHT) ]]; }}

# Perform batch mode channel filtering on processed .cut files
for my $c (@chn) { for ( @{ ( &AUXILIARY_CHANNEL($c)) || [] } ) {( print STDERR $_."\n" )}}
# Files specified explicitly in the card are processed at this stage (unless in command-line file mode)

}	# End local enclosing scope for main program routine

#*************#
# SUBROUTINES #
#*************#

# Returns a subroutine closure encapsulating access by index or hashed key to runtime command line parameters
sub LOAD_OPTS { my (@opts); my (%opts) = map { ( /^-(-)?([A-Za-z]\w*)(?:=(.*))?$/ ) ? do { ($1) ? ((lc $2) => ((defined $3) ? qq($3) : q(1))) :
	( map {((lc) => q(1))} ( $2 =~ /([A-Za-z])/g )) } : do { push @opts, $_; () }} (@_); sub { my ($key) = ( lc shift );
		(return $_) for ( map { (m/^(DEFAULT)|(TRUE)|(FALSE)|(UNDEF)$/i) ? (($1)?(0):($2)?(+1):($3)?(-1):(undef)) : ($_) }
		(( exists $opts{$key} ) ? ($opts{$key}) : ( shift @opts ))); }}
# --KEY=VAL : KEY begins alpha and is subsequently word characters; VAL is arbitrary; without the '=', VAL is '1'
# -ABC : each of the flags keyed by 'A', 'B', and 'C' is set to '1'
# Other: string is pushed onto the ordered queue
# On query, existing (case insensitive) KEYs return VAL and missing keys shift values from the front of the queue
# Special strings DEFAULT, TRUE, FALSE, and UNDEF are mapped internally to 0, +1, -1, and undef

# Returns a data structure encoding user specified AEACuS meta language instructions
{; my ($abc,$idx,$val); sub LOAD_CARD { my ($crd,$err,$fil) = (+{},[]);
	$abc ||= qr'(?i:[A-Z][A-Z\d]{2})'; $idx ||= qr'\d{1,3}'; $val ||= do { my ($key) =
		[ qr"(?:(${idx})|(${abc})(?:_(${idx}))?)", sub { (shift); ( map { (length) ? (0+ $_) : ( +{ ( lc (shift)) => (0+ (shift)) } ) } (shift))[0] } ]; [
		[ qr'(?i:DEFAULT|TRUE|FALSE|UNDEF)', sub { ${{ default => (0), true => (+1), false => (-1), undef => (undef) }}{ lc (shift) }} ],
		[ qr"${\EXP}", sub { 0+(shift) } ],
		[ qr'"([^"]*)"', sub {( &UNESCAPE_STRING((shift,shift)[1] ))} ], ($key),
		[ qr'\{([^}]*)}', sub { my ($q,@t) = ((shift,shift)[1] );
			while ( my $ref = ( &REX_LIST( 2, $q, [ qr",\s*$$key[0]\s*$", $$key[1]] ))) {( unshift @t, $$ref )}
			( map { (defined) ? [$_,@t] : (undef) } ( &STRING_FUNCTIONAL($q)))[0] } ]] }; (($crd) -> {$$_[0]} -> [$$_[1]] = ($$_[2])) for
	map { my ($l) = $_; ( $$l[1] =~ m/^(?:\*|\s*$)/ ) ? () : ( $$l[1] =~ m/^(?:${abc}_)?(${abc})(?:_(${idx}))?\s*=\s*(.*?)$/ ) ? [ (lc $1), (0+ $2),
		do { my ($h,$q,$e) = (+{},qq($3)); while ( $q =~ m/\G(?(?!^),\s*)(${abc})\s*:\s*/gc ) { $$h{lc $1} = do { my ($a,$b) = ([],!!( $q =~ m/\G\[\s*/gc ));
		while ((($b) || !(@$a)) and (!(@$a) or ( $q =~ m/\G,\s*/gc )) and ( my $ref = ( &REX_LIST( 1, $q, @$val )))) { ( push @$a, $$ref ) && ( $q =~ m/\G\s+/gc ) }
		($e) ||= (($b) and ( $q !~ m/\G]\s*/gc )); ($a) }} ( push @$err, $l ) if (($e) || ( $q !~ m/\G$/gc )); ($h) } ] : do { push @$err, $l; () }}
	do { my (@t,$i,$c) = ([1]); ((( my ($FHI)),($fil)) = ( &Local::FILE::HANDLE(shift))) or (return); local ($_); while (<$FHI>) { $i++; chomp; ($c = 1) if (s/^!.*//); s/#.*//;
		(( push @t, [$i] ) && ( undef $c )) if (m/^\S/); ( $t[-1][1] .= $_ ) unless ($c); } (@t) }; ((wantarray) ? ($crd,$err,$fil) : ($crd)) }}
# The card data structure is an array (list) of hashes (dictionaries) of hashes, corresponding to lines of the form "ABC_IDX = ... "
# ABC is a letter followed by precisedly two letters or digits (case-insensitive); IDX is a number with one to three digits (defaults to "_000")
# The list elements are derived from comma-separated ABC:VAL or ABC:[VAL,VAL,...] expressions after the equality assignment; indent for line continuation
# VAL may be one of DEFAULT, TRUE, FALSE, UNDEF, which map internally to 0, +1, -1, and undef
# VAL may be a signed number, possibly in scientific notation, which maps to the corresponding numerical value
# VAL may be a string inside double quotes which goes through an unescape and substitution cycle prior to storage
# VAL may be a KEY of the form IDX or ABC or ABC_IDX, which is converted internally to a number or to the hash reference { ABC => undef } or { ABC => IDX }
# VAL may be a function { f($1,$2,...), VAL, VAL, ... }, which is converted internally to an array reference holding a closure and keys for value substition

# Return list of event objects encapsulating kinematic and tagging variables converted from list of LHCO formatted strings
sub EVENT_OBJECTS { map {((wantarray) ? (@$_) : (return $_))} [ map { my ($o,$t,$s) = ( $_, (( map {(int)} ( split /\./, ( shift @$_ ))), 0 ));
	((( $t < 0 ) or ( $t == 5 ) or ( $t > 6 )) ? () : ( grep { (( @$o[1,2,3] ) = (( &PRINCIPAL_RAD($$o[1])), ( map {( &MAX(0,$_))} ( @$o[2,3] ))));
	( @$_{ qw( typ sub eta phi ptm mas trk muo sgn hft fem ptc etr ep0 ep1 ep2 ep3 dm1 dm2 aux )} ) = ( $t, ( scalar ( &INT_2_FLAGS($s))), @$o[0..3],
	( map {(($t == 4) ? (((( map {(int)} ( split /\./ )), 0 )[0,1] ), 0 ) : (( int abs ), 0, (((undef),1,1,1)[$t] ? ( $_ <=> 0 ) : (0))))} ($$o[4])), (($t == 4) ? ( abs $$o[5] ) : (0)),
	( map {(($t == 2) ? ((undef), (int), (( &RATIO((int), $$o[2] )) + ((abs) - ( int abs )))) : ((1/( 1 + (abs))), (undef,undef)))} ($$o[6])),
	( map {(( sqrt (($$o[3])**2 + ($$_[0])**2 + ($$_[1])**2 + ($$_[2])**2)), @$_ )} [ map {($_*$$o[2])} (( cos $$o[1] ), ( sin $$o[1] ),
	( map {( &RATIO(( cos $_ ), ( sin $_ )))} 2*( atan2 (( exp ( -1 * $$o[0] )), 1 )))[0] )] ), (@$o[7,8]), [( splice @$o, 9 )] ); 1 } +{} )) }
	grep {(( shift @$_ ) > 0 )} map {[ map {(( m/^NAN$/i ) ? (undef) : (0+ $_))} ( split ) ]}
	grep {( m/^\s*(\d+)\s+(\d+)/ )} map {( qq($_))} map {(( ref eq q(ARRAY)) ? (@$_) : ($_))} (@_) ] }

# Returns code reference for generation of primary object reconstruction and event selection statistics
sub ANALYSIS_CODE { my ($crd,@k,@i,@c) = (shift); do { my ($t) = [ grep {( defined $$_[0] )} (@$_) ];
	push @k, [ map {($$_[0])} (@$t) ]; push @i, [ grep {($$t[$_][1])} (0..(@$t-1)) ]; push @c, [ map {($$_[2])} (@$_) ]; } for
	map { my (%t); ( push @{ $t{(shift @$_)} ||= [] }, $_ ) for (@$_); map {( $t{$_} )} sort { our ($a,$b); ( $a <=> $b ) } ( keys %t ) }
	map { my (@t); ( push @{ $t[0+(( shift @$_ ) > 0 )] ||= [] }, $_ ) for ( sort { our ($a,$b); ( $$a[0] <=> $$b[0] ) } (@$_)); ( grep {($_)} (@t)) }
	map { my ($h,$l,$p) = @$_; [ map { my ($k,$c,$f,$a,$z) = @$_; map { my ($i) = $_; map { my ($e,$m,@o) = ( $_,
		( !( &MATCH_VALUE( $$_{cut}, (undef))) or ((0) <=> (($p) && ( &OUTPUT_VALUE( $$_{out} ))))), (($p) ? () : ( &OUTPUT_OBJECT( $$_{out} ))));
		map { my ($r,$s,$n,$q,$j,$w,$x) = (!!($_),@{$_||[]}); my ($w,$f) = map {((((!1)x(4)),1)[$_],$_)} (($r) ? (( ${{ eta => 3, phi => 3 }}{$q} ) or (($w)?(-1):(2))) : ( 0+(0..4)[$f] ));
			[ 0+($s), (($r) ? ( $j, [(undef),$q,$j,$w], !1 ) : ( $i, (($m) ? ( [$h,$k,$i,$w], ($m > 0)) : ((undef), !1 )))), ( sub { my ($o,$v) = (shift,shift);
			( return ((defined) ? [ q().( q(       ), q(%7.1i), q(%7.1f), q(%7.3f), q(%+10.3E))[ (($f) or (((($_) == (int)) and (($_) < (10**(+7)))) ? (1) :
				((abs) >= (10**(+3))) ? (2) : (3))) ], 0+($_) ] : [ (($w) ? q(   ) : q()).q(  UNDEF), (undef) ] )) for
			(($r) ? (($$v{$q}[$j]) = ( map {(($x) ? do { my ($k,$v) = (%$x); ($$_{$k}[$v]) } : ($$_{$q}))} ( $$o{$k}[$i][$n] || +{} ))) :
			( map { ($m) ? (($m < 0) or ( &MATCH_VALUE( $$e{cut}, $_ )) or (return undef)) : (return); ($$_) } \( $$v{$k}[$i] = (($c)->($e,$o,$v,$i))))); } ) ] } ((undef),@o) }
		(( &CLONE($$crd{$k}[$i])) || (($i > 0) ? () : +{} )) } (( &DEFINED($a,1))..( &DEFINED($z,( &MAX(0,(@{$$crd{$k}||[]}-1)))))) } (@$l) ] }
				( [ obj => [
		# Group photons and individual lepton flavors, and filter on kinematics, sign, and isolation
	( map { my ($k,$src) = (($_), ((undef), qw( ele muo tau all ))[$_]); [ ($src) => sub {
		my ($t) = [ grep { ($k < 0) || ($$_{typ} == $k) } (@{ $_[1]{(($k < 0)?q(obj):q(all))}[0] || [] }) ];
		( map {( 0+@$_ )} grep { !($k < 0) && ( 0+(0,1,0)[($_[0]{jet}[0])] ) && do { do { $$_{typ} = 4 } for ( &EXCLUDE_OBJECTS($_,(@$t))) }; 1 }
			grep { $_[1]{$src}[$_[3]] = $_; 1 } ( &SELECT_OBJECTS( $_[0], $t, (undef), $_[3], ((undef),0,0,0,-1)[$k], $k )))[0] }, 1, 0, 0 ] } (-1,(1..3))),
		# Reconstruct photon, lepton, and jet subsets consistent with various kinematic, tagging, and back-reference criteria
# HERE ... allow for selection of comparison group ... include explicit, implicit, allowed keys globally
	( map { my ($j,$k,$src,$cmp) = (($_), ((undef),4,0)[$_], ( qw( lep jet pho ))[(($_),(1,0,1)[$_])] ); [ ($src) => sub {
		my ($t) = [ ($_[3] > 0) ? (( &INCLUDE_OBJECTS( $_[0]{src}, $_[3], $_[1]{$src} )), (($j == 1) ? (( &IPHO ),( &ILEP )) : ())) :
			(($j == 0) ? ( map {(@{ $$_[0] || [] })} (@{ $_[1] || +{}}{ qw( ele muo tau ) })) : ( grep {($$_{typ} == $k)} (@{ $_[1]{all}[0] || [] }))) ];
		( map {( 0+@$_ )} grep { ($_[3] == 0) && ($j != 1) && ( 0+(0,1,0)[($_[0]{jet}[0])] ) && do { do { $$_{typ} = 4 } for ( &EXCLUDE_OBJECTS($_,(@$t))) }; 1 }
			grep { $_[1]{$src}[$_[3]] = $_; 1 } ( &SELECT_OBJECTS( $_[0], $t, [ &INCLUDE_OBJECTS( $_[0]{cmp}, $_[3], $_[1]{$cmp} ) ], $_[3], $j )))[0] }, 1, 0 ] } (-1..1)),
				], !1 ], [ evt => [
		# Filter on the missing transverse momentum, scalar sum on transverse energy, and effective mass of various lepton and jet reconstructions
	[ cal => sub {( ${( $_[1]{met} = [ $_[1]{cal}[0] || []] )}[0][0] )}, 2, 0, 0 ],
	[ met => sub {(	${(( grep { my ($a,$i) = (($_[1]{met}||[]), $_[3] );
		(($i < (@$a - 1)) ? ( splice @$a, $i, 1, $_ ) : ( splice @$a, (-1), 0, (((undef)x($i-@$a+1)),$_))); 1 }
		(( scalar &MET(($_[3] == 0) ? (@{ $_[1]{all}[0] || [] }) : ( &IOBJ ))) || [] ))[0] )}[0] )}, 2, 0 ],
	[ mht => sub {(	&MHT( $_[0]{msv}, (($_[3] == 0) ? (@{ $_[1]{all}[0] || [] }) : ( &IOBJ ))))}, 2, 0 ],
	[ mef => sub {(	&SUM(( scalar &INDEXED_VALUES((($_[3] > 0) ? $_[0]{met} : [0]), $_[2]{met} )),
		( scalar &INDEXED_VALUES((($_[3] > 0) ? $_[0]{mht} : [0]), $_[2]{mht} ))))}, 2, 0 ],
		# Filter on various compound ratios and differences composed out of the previously computed mass dimensioned statistics
	[ ret => sub {( &RATIO( &INDEXED_VALUES( $_[0]{num}, $_[2]{met}, $_[0]{den}, $_[2]{met} )))}, 3 ],
	[ rhr => sub {( &RATIO(( scalar &INDEXED_VALUES( $_[0]{num}, $_[2]{met} )),
		( map { (defined) ? (sqrt $_) : (undef) } ( scalar &INDEXED_VALUES( $_[0]{den}, $_[2]{mht} )))[0] ))}, 3 ],
	[ ref => sub {( &RATIO(( scalar &INDEXED_VALUES( $_[0]{num}, $_[2]{met} )), ( scalar &INDEXED_VALUES( $_[0]{den}, $_[2]{mef} ))))}, 3 ],
	[ rhh => sub {( &RATIO( &INDEXED_VALUES( $_[0]{num}, $_[2]{mht}, $_[0]{den}, $_[2]{mht} )))}, 3 ],
	[ det => sub {( ${ &LORENTZ_DIFFERENCE( &INDEXED_VALUES( $_[0]{one}, $_[1]{met}, $_[0]{two}, $_[1]{met} )) || []}[0] )}, 2 ],
		# Filter on various specialized discovery statistics employed by the LHC experimental collaborations
	[ ote => sub {( &TRANSVERSE_ENERGY((( &IMET ) || ()), ( &IOBJ )))}, 2 ],
	[ oim => sub {( &INVARIANT_MASS((( &IMET ) || ()), ( &IOBJ )))}, 2 ],
	[ otm => sub {( &TRANSVERSE_MASS(( &IMET ), ( scalar &LORENTZ_SUM( undef, ( &IOBJ )))))}, 2 ],
	[ stm => sub {( &S_TRANSVERSE_MASS( &HEMISPHERES( 0, q(LND), ( &IOBJ ))))}, 2 ],
	[ atm => sub {( &MIN( grep {(defined)} map {( &A_TRANSVERSE_MASS(( &IMET ), (@$_)))}
		( &AMT2_ROLES( $_[0]{mod}, $_[0]{lep}, $_[1]{lep}, $_[0]{jet}, $_[1]{jet} ))))}, 2 ],
	[ tjm => sub {( &TRI_JET_MASS( $_[0]{lim}, ( &IJET )))}, 2 ],
	[ ttm => sub {( &TAU_TAU_MASS( [ &ILEP ], [ &IJET ] ))}, 2 ],
	[ jzb => sub {( &JET_Z_BALANCE( [ (( &ILEP ), undef )[0]], [ &IJET ]))}, 2 ],
	[ jrm => sub {( &ALPHA_R(( &IMET ), ( &HEMISPHERES( 0, q(MIM), ( &IOBJ )))))[1] }, 2 ],
	[ odr => sub {( &DELTA_RPA( &LORENTZ_MERGE([ 2, 1, +1 ], ( &IOBJ ))))}, 3 ],
	[ oda => sub {( &DELTA_RSA( &LORENTZ_MERGE([ 2, 2, +1 ], ( &IOBJ ))))}, 3 ],
	[ odp => sub {( &DELTA_PHI( &LORENTZ_MERGE([ 2, 3, +1 ], ( &IOBJ ))))}, 3 ],
	[ ode => sub {( &DELTA_ETA( &LORENTZ_MERGE([ 2, 4, +1 ], ( &IOBJ ))))}, 3 ],
	[ alr => sub {( scalar &ALPHA_R(( &IMET ), ( &HEMISPHERES( 0, q(MIM), ( &IOBJ )))))}, 3 ],
	[ alt => sub {( &ALPHA_T(( &IMET ), ( scalar &INDEXED_VALUES( $_[0]{mht}, $_[2]{mht} )), ( &HEMISPHERES( 0, q(MDH), ( &IOBJ )))))}, 3 ],
	[ mdp => sub {( &MET_DELTA_PHI(( &IMET ), ( &IOBJ )))}, 3 ],
	[ bdp => sub {( &BIASED_DELTA_PHI(( &IMET ), ( &IOBJ )))}, 3 ],
	[ cts => sub {( &COSINE_THETA_STAR( &ILEP ))}, 3 ],
	[ lwp => sub {( &LEP_W_PROJECTION(( &IMET ), (( &ILEP ), undef )[0] ))}, 3 ],
	[ tts => sub {( scalar &THRUST_SHAPE(( scalar &INDEXED_VALUES( $_[0]{mht}, $_[2]{mht} )), ( &IOBJ )))}, 3 ],
	[ tms => sub {( &THRUST_SHAPE(( scalar &INDEXED_VALUES( $_[0]{mht}, $_[2]{mht} )), ( &IOBJ )))[1] }, 3 ],
	[ tos => sub {( &SPHEROCITY_SHAPE(( scalar &INDEXED_VALUES( $_[0]{mht}, $_[2]{mht} )), ( &IOBJ )))}, 3 ],
	[ tss => sub {( &SPHERICITY_SHAPE( &IOBJ ))}, 3 ],
	[ tfs => sub {( &F_MATRIX_SHAPE( &IOBJ ))}, 3 ],
				], 1 ], [ usr => [
		# Filter on user-defined composite event statistics
	[ var => sub { my ($sub,@val) = map { ( ref eq q(ARRAY)) ? (@$_) : ( sub {(shift)} , $_ ) } (${$_[0]{val}||$_[0]{key}||[]}[0]);
		(($sub)->( map { ( ref eq q(HASH)) ? do { my ($k,$v) = %$_; ${$_[2]{$k}||[]}[$v] } : (undef) } (@val)))}, -1 ],
		# Filter on externally defined event statistics
# HERE ... itemize die, use system/IPC ... safety wrapper
	[ ext => sub { (( my ($exe) = ( grep {((length) && ( -x ))} (( ${$_[0]{exe}||[]}[0] =~ /^(.*\/)?([^\/]+)$/ ) ?
		( join q(), @{ ( &Local::FILE::HANDLE([ (($1) or ( q(./External/))), ($2) ]))[1]||[] } ) : ()))) or
		( die 'Cannot execute external routine' ));
		( map {(((length) && !(m/UNDEF/i)) ? ( sprintf q(%+10.3E), (0+ $_)) : (undef))}
		map {`($_)`} ( join q( ), (($exe), ( map {((defined) ? ( sprintf q(%+10.3E), $_ ) : q(UNDEF))}
		map { my ($sub,@val) = ( ref eq q(ARRAY)) ? (@$_) : ( sub {(shift)} , $_ );
			(($sub)->( map { ( ref eq q(HASH)) ? do { my ($k,$v) = %$_; ${$_[2]{$k}||[]}[$v] } : (undef) } (@val))) } (@{$_[0]{key}||[]})))))[0]
		}, -1 ],
				], 1 ], [ wgt => [
		# Filter on user-defined event weight statistics
	[ wgt => sub {( $_[1]{wgt} )}, -1, 0 ],
				], 1 ] );
		# Construct event analysis closure and list of reported statistics
	( [ map {[ @{ $k[$_] }[ @{ $i[$_] } ]]} (0..(@i-1)) ], [ map {(@$_)} (@k) ], sub { my ($o,$v) = (+{},+{});
		@$o{( qw( nbr trg wgt cal obj ))} = ( @{(shift)||[]}[0..2], [((shift) or (return))], [[ @{((shift) || (return))} ]] ); ( -1,
			[ map { my ($i) = $_; map {(@$_)} grep { my ($t); do { ( return ( $i, $_ )) if ($t) } for
			[ map {( !(defined) && ($t = 1))} ( @$_[ @{ $i[$i] } ] ) ]; 1 }
			[ map {(($_)->($o,$v))} (@{ $c[$i] }) ] } (0..(@i-1)) ] ) } ) }

# Returns a formatted header summarizing processed and surviving event counts, and selection cut activity rates
sub FORMAT_HEADER { my ($prm) = !(!(shift)); my ($evt,$xsc,$lum,$srv,$res) = ( map {(( ref eq q(ARRAY)) ? (( &MAX(0,( int $$_[0] ))),
	( map {(( ref eq q(ARRAY)) ? ( &MERGE_XSEC(@$_)) : (defined) ? (0+$_) : (undef))} ($$_[1])), (0+ $$_[2]),
	(( &BOUNDED([0,( int $$_[0])],( int $$_[(($prm) ? (0) : (3))] ))) or (return)), (0+ $$_[4])) : (return))} (shift));
	my ($cut,$flw,$stc) = ( map {(( ref eq q(ARRAY)) ? ($_) : (undef))} (shift,shift,shift));
	my ($pas) = ( &MAX(0,( int ${${$cut||[]}[-1]||[]}[0] ))); ( map { (wantarray) ? ( $_, $xsc ) : (return $_) } (
	(($prm) ? "\n".(($cut) ? q() : q(# )).($evt).' EVENT SAMPLE'.(($evt==1) ? q() : q(S)).' PROCESSED IN TOTAL'."\n" .
		(( defined $xsc ) ? "\n".(($cut) ? q() : q(# )).( sprintf q(%+10.3E), $xsc ).' PB CROSS SECTION IMPLIES ' .
		( sprintf q(%+10.3E), ($evt/$xsc)).' PER PB LUMINOSITY'."\n".(($lum && $cut) ? "\n" .
		( sprintf q(%+10.3E), $lum ).' PER PB SCALED LUMINOSITY YIELDS ' .
		( sprintf q(%+10.3E), ($xsc*$lum)).' EFFECTIVE EVENTS'."\n" : q()) : q()) : q()) .
	(($cut) ? do { $lum *= ($xsc/$evt); my ($s); my ($t) = ( join q(), map { my ($t) = $_; ( join q( ), (
		( map {(( uc $$_[1] ).( q(_)).( sprintf q(%3.3i), $$_[2] ))} map {((@$_==1) ? ($$_[0]) : [(undef),@{$$_[0]}[0,2]] )} ($$t[0])),
		( map {((($lum) ? ( sprintf q(%7.1E), $_*$lum ) : ( sprintf q(% 7i), $_ )), (( defined $xsc ) ? ( sprintf q(%7.1E), $xsc*($_/$evt)) : ()),
		( map {((defined) ? ( sprintf q(%07.3f), 100*$_ ) : q(  UNDEF))} (( &RATIO($$t[1],($_+$$t[1]))), (1-($_/$srv)))))}
		($s+$pas-$$t[2]))))."\n" } map { (($srv == ($s + $pas)) or (return)); (@$_) } [ map {((@{$$flw[$_]||[]}) ?
		[ $$flw[$_], ( map {( $_, $s += $_ )} ( &MAX(0,( int $$cut[$_][0] )))) ] : ())} (0..(@{$flw||[]}-1)) ] );
		"\n".($pas).' SAMPLE'.(($pas==1) ? q( SURVIVES) : q(S SURVIVE)).' '.(($prm) ? q(PRIMARY) : q(AUXILIARY)).' SELECTION WITH ' .
		(( defined $xsc ) ? ( sprintf q(%+10.3E), ($xsc*($pas/$evt))).' PB RESIDUAL CROSS SECTION AND ' : q()) .
		( sprintf q(%07.3f), ( 100*(1-($pas/$evt)))).' PERCENT CUMULATIVE EFFICIENCY'."\n" .
		(($s > 0) ? "\n".'CUT_KEY SURVIVE'.(( defined $xsc ) ? ' PB_AREA' : q()).' %_LOCAL %_TALLY'."\n".($t) : q()) } : q()) .
	((($prm) && ($cut) && ($pas)) ? ( "\n".( join q( ), ( q(EVENT_N),
		( map {((($$_[3]) ? q(   ) : q()).( uc $$_[1] ).( q(_)).( sprintf q(%3.3i), $$_[2] ))} ( @{$stc||[]} ))))."\n" ) : q()))) }

# Applies auxiliary channel filtering specification to an existing .cut event selection file
sub AUXILIARY_CHANNEL { my (@err); my ($chn) = grep {(( ref eq q(HASH)) or (return undef))} (shift); my ($fil) = grep {( ref eq q(ARRAY))} (shift);
	my ($cid,$did) = map {( $_, (( sprintf q(CHN_%3.3i), $_ ).q(/)))} (( int abs ($$chn{cid})) or (return undef));
	CHN: for my $fil ( sort SORT_LIST_ALPHA ( values %{ +{ map {(($$_[0].$$_[1]) => ($_))} grep {( $$_[1] =~ /^[\w-]+\.cut$/ )}
		(($fil) ? ($$chn{fil}) ? () : ($fil) : ( map {( &Local::FILE::LIST($_))} (@{$$chn{fil}||[]}))) }} )) {
	( my ($FHI) = ( &Local::FILE::HANDLE($fil))) or do { push @err, 'CANNOT READ FROM FILE '.$$fil[0].$$fil[1]; (next CHN) };
	( my ($xsc,$idx,$txt) = ( &AEACUS_XSEC($FHI))) or do { push @err, 'CANNOT CUE HANDLE OF FILE '.$$fil[0].$$fil[1]; (next CHN) };
	(($$xsc[3] > 0) && ($$xsc[0] >= $$xsc[3])) or do { push @err, 'NO SURVIVING EVENTS IN FILE '.$$fil[0].$$fil[1]; (next CHN) };
	( defined $idx ) or do { push @err, 'CANNOT INTERPRET HEADER OF FILE '.$$fil[0].$$fil[1]; (next CHN) }; ( my (@cut) = map {[$_]}
		sort { our ($a,$b); ( $$a[1] <=> $$b[1] ) } grep { ( $$_[2] = ( &HASHED_FUNCTIONAL( $idx, ( map { ( ref eq q(ARRAY)) ? (@$_) : (undef,$_) }
			($$_[2][0]))))) or do { print STDERR 'INVALID KEY IN SELECTION '.$$_[1].' FOR CHANNEL '.$cid.' ON FILE '.$$fil[0].$$fil[1]."\n"; !1 }}
		grep { !( &MATCH_VALUE( $$_[3], (undef))) } map {[(@$_)]} (@{$$chn{esc}||[]})) or
		do { push @err, 'NO ACTIVE SELECTION FOR FILE'.$$fil[0].$$fil[1].' IN CHANNEL '.$cid; (next CHN) };
	( my ($FHT) = ( &Local::FILE::HANDLE())) or do { push @err, 'CANNOT OPEN TEMPORARY FILE FOR READ/WRITE'; (next CHN) };
	my ($cut,$flw) = ( [ map {[0,!1]} (0..(0+@cut)) ], [ map {[ map {[ q(cut), (($$_[0]) ? q(inv) : q(esc)), 0+$$_[1], 0 ]} (@$_) ]} (@cut) ] );
	local ($_); while (<$FHI>) { my ($l,$v) = ((/^\s*(\d+)/) ? ( $_, [ map { (/^UNDEF$/) ? (undef) : (0+ $_) } ( split ) ] ) :
		do { print $FHT $_; (last) } ); for my $i ((0..(@cut-1)),-1) { my ($t); if ($i < 0) { print $FHT $l } else {
		my ($v) = [ map { my ($x,(undef),$k,$c) = @$_; !(($x) xor ( &MATCH_VALUE( $c, (($k)->($v))))) && ($t = 1) } (@{$cut[$i]}) ];
		(next) unless ($t); ( $$cut[$i][1] += ( &Local::VECTOR::OUTER_PRODUCT($v,$v))) if (@$v > 1) }; $$cut[$i][0]++; (last) }}
	( my ($hdr) = ( &FORMAT_HEADER(!1,$xsc,$cut,$flw))) or do { push @err, 'CANNOT FORMAT HEADER FOR CHANNEL '.$cid.' OF '.$$fil[0].$$fil[1]; (next CHN) };
	( seek $FHT, 0, SEEK_SET ) or do { push @err, 'CANNOT REWIND TEMPORARY FILE WHILE FILTERING '.$$fil[0].$$fil[1].' IN CHANNEL '.$cid; (next CHN) };
	( my ($FHO) = ( &Local::FILE::HANDLE([[$$fil[0],$did],$$fil[1]],1))) or do { push @err, 'CANNOT WRITE TO FILE '.$$fil[0].$did.$$fil[1]; (next CHN) };
	print $FHO q(), ( splice @$txt, 0, -2 ), ($hdr), (@$txt); local ($_); while (
		defined ( $_ = ((<$FHT>) or (<$FHI>)))) {( print )}} ((@err) ? (\@err) : (undef)) }

# HERE ... unify the below.  Include pythia8 / other sources?
# Returns the same-directory Merged_XSEC.TXT event production cross section (mean) and error (RMS) associated with a standardized .lhco file
sub MERGED_XSEC { my ($FHI,@xsc,@dev) = ( &Local::FILE::HANDLE(
	grep {(( $$_[-1] =~ s/(?:^|\/)([\w-]+?)_(?:delphes|pgs)_events\.lhco(?:\.gz)?$/${1}_merged_xsecs.txt/ ) or (return))}
	map {[ ( ref eq q(ARRAY)) ? (@$_) ? (@$_) : (return) : (($_),((@_)?(shift):())) ]} (shift)) or (return));
	local ($_); while (<$FHI>) {(( m/^\s*${\EXP}\s+(${\EXP})\s+(${\EXP})\s*$/ ) && ( push @xsc, (0+ $1)) && ( push @dev, (0+ $2)))}
	map { (wantarray) ? ( $_, ( &RMS( @dev ))) : ( return $_ ) } ( &ARITHMETIC( @xsc )) }

# Returns the same-directory SUMMARY.TXT event production cross section and error associated with a standardized .lhco file
sub SUMMARY_XSEC { my ($FHI,$xsc,$dev) = ( &Local::FILE::HANDLE(
	grep {(( $$_[-1] =~ s/(?:^|\/)([\w-]+?)_(?:delphes|pgs)_events\.lhco(?:\.gz)?$/summary.txt/ ) or (return))}
	map {[ ( ref eq q(ARRAY)) ? (@$_) ? (@$_) : (return) : (($_),((@_)?(shift):())) ]} (shift)) or (return));
	local ($_); while (<$FHI>) {((($xsc,$dev) = ( m/^\s*Total\s+cross\s+section:\s+(${\EXP})\s+\+-\s+(${\EXP})\s+pb/i )) && (last))}
	map { (wantarray) ? ( $_, ( 0+ $dev )) : ( return $_ ) } ( 0+ $xsc ) }

# Returns the same-directory Pythia.LOG event production cross section associated with a standardized .lhco file (fragmentation corrected, old file format)
sub PYTHIA_XSEC { use Fcntl qw(:seek); my ($FHI) = grep {( seek ($_,-80,SEEK_END))} ( &Local::FILE::HANDLE(
	grep {(( $$_[-1] =~ s/(?:^|\/)([\w-]+?)_(?:delphes|pgs)_events\.lhco(?:\.gz)?$/${1}_pythia.log/ ) or (return undef))}
	map {[ ( ref eq q(ARRAY)) ? (@$_) ? (@$_) : (return undef) : (($_),((@_)?(shift):())) ]} (shift)) or (return undef));
	local ($_); while (<$FHI>) {(( m/^\s*Cross\s+section\s+\(pb\)\s*:\s*(${\EXP})/i ) && ( return (0+ $1)))} (undef) }

# Returns the <MGGenerationInfo> tag event cross section associated with an open, cued filehandle
sub MADGRAPH_XSEC { my ($FHI,$evt,$mch,$wgt) = grep {((( ref eq q(GLOB)) or
	( &ISA( 1, $_, q(Local::FILE)))) && ((tell $_) >= 0) or (return))} (shift);
	local ($_); while (<$FHI>) {(
		( m/^\s*#*\s*<\/MGGenerationInfo>/i ) ? (last) :
		( m/^\s*#*\s*Number\s+of\s+Events\s*:\s*(\d+)/i ) ? ( $evt = (0+ $1)) :
		( m/^\s*#*\s*(Matched\s+)?Integrated\s+weight\s+\(pb\)\s*:\s*(${\EXP})/i ) &&
			(($mch,$wgt) = (((($mch) && !($1)) or (($2) eq q(-1.0))) ? (next) : (!!($1),(0+ $2)))))}
	map {((wantarray) ? (@$_) : (return $$_[0]))} ((( defined $wgt ) && (($mch) or ($evt > 0))) ? [ ($wgt), (($mch) ? () : ($evt)) ] : ()) }

# Returns the AEACuS header information, event count and cross section associated with an open, cued filehandle
sub AEACUS_XSEC { my ($FHI,$e,$x,$l,$s,$r,$i,@t) =
	grep { (( ref eq q(GLOB)) or ( &ISA( 1, $_, q(Local::FILE)))) && ((tell $_) >= 0) or (return) } (shift);
	local ($_); while (<$FHI>) { push @t, $_;
		( m/^\s*(\d+)\s+EVENT\s+SAMPLES?/i ) ? do { $e = (0+ $1) } :
		( m/^\s*(${\EXP})\s+PB\s+CROSS/i ) ? do { $x = (0+ $1) } :
		( m/^\s*(${\EXP})\s+PER\s+PB\s+SCALED/i ) ? do { $l = (0+ $1) } :
		( m/^\s*(\d+)\s+SAMPLES?\s+SURVIVES?(?:\s+(?:PRIMARY|AUXILIARY)\s+SELECTION\s+WITH\s+(${\EXP})\s+PB\s+RESIDUAL)?/i ) ?
			do { $s = (0+ $1); ( $r = (0+ $2)) if ( length $2 ) } :
		( m/^\s*EVENT_(?:N|#)/i ) ? do { my ($j); $i = +{ map {((lc) => ($j++))} ( split ) }; (last) } : (next) }
	((wantarray) ? ([$e,$x,$l,$s,$r],$i,\@t) : ($x)) }

# Returns the unified cross section for a merged data set
sub MERGE_XSEC { my (@x) = grep {($$_[2] > 0)} grep {(( ref eq q(ARRAY)) or (return undef))} (@_); ( &RATIO(
	( &SUM( map { ( defined $$_[0] ) ? (($$_[0])*(($$_[1] > $$_[2]) ? ($$_[2]/$$_[1]) : (1))*($$_[2])) : (undef) } (@x))), ( &SUM( map {($$_[2])} (@x))))) }

# Returns the VECTOR of cut counts, and MATRIX cut correlation metric for an input weighted, normalized Boolean MATRIX of event-by-event cuts
sub CUT_METRIC { my ($e); my ($m) = (
	map {( scalar &Local::TENSOR::MAP( $_, ( sub {( &FLUSH((shift), NIL ))} )))}
	map { my ($t) = $_; map {(( $t * ( $_ x $_ )) or (return))}
		( &Local::VECTOR::OBJECT([ map {(( &RATIO(1,( sqrt $$t[$_][$_] ))) or (1))} (0..(@$t-1)) ]) or (return)) }
	map { my ($t) = $_; map {(( $t - ( $_ x $_ )) or (return))}
		( $e = &Local::VECTOR::OBJECT([ map {($$t[$_][$_])} (0..(@$t-1)) ]) or (return)) }
	( &Local::MATRIX::OBJECT(shift) or (return))); ( $$m[$_][$_] = 1 ) for (0..(@$m-1));
	((wantarray) or ( return $m )); ($e,$m) }

# Attempts to generates an LHCO event file from a DELPHES event file in the fashion of the standard MadGraph script "run_delphes3"
sub ROOT_2_LHCO { my ($fil,$r2l,%cfg,$xpb) = (( grep {(( ref eq q(ARRAY)) or (return undef))} (shift)), (shift));
	do { my ($FHI) = (( &Local::FILE::HANDLE([[ $$fil[0], q(../../Cards/) ], q(me5_configuration.txt) ])) or
		( &Local::FILE::HANDLE([[ $$fil[0], q(../../Cards/) ], q(amcatnlo_configuration.txt) ])) or (return undef));
		local ($_); while (<$FHI>) { ( $cfg{$1} = $2 ) if ( m/^\s*([\w-]+)\s*=\s*(\S*)\s*$/ ) }};
	my ($dph) = (( &Local::FILE::PATH([ $cfg{mg5_path}, ( &DEFINED( $cfg{delphes_path}, q(./Delphes))) ])) or (return undef));
	my ($pth,$run) = ( $$fil[0] =~ /^(.*\/)([^\/]+)\/$/ ) or (return undef);
	my ($tag) = ( $$fil[1] =~ /^([^\/]+)_delphes_events.root$/ ) or (return undef);
	my ($wgt,$gen,$fat) = ( map {(($_) ? q(True) : q(False))} @$r2l{( qw( wgt gen fat ))} );
	my ($aux) = ( join q( ), ( q([), ( grep {(length)} ( join q(, ), ( map {( qq('${_}'))} @{$$r2l{aux}||[]} ))), q(])));
	local ($?); do { open my $FHO, q(|-), q(python >/dev/null 2>&1); print $FHO <<DONE;
# Adapted for Python from "root2lhco.cpp", Delphes, (C) 2012-20XX UCL Belgium, GNU GPL V3ff
import sys
if (( sys.version_info[0] < 2 ) or (( sys.version_info[0] == 2 ) and ( sys.version_info[1] < 7 ))) : sys.exit( 1 )
import ROOT
ROOT.gROOT.ProcessLine( '.include ${dph}external' )
ROOT.gSystem.Load( '${dph}libDelphes' )
wgt = ${wgt}; gen = ${gen}; fat = ${fat} 
aux = ${aux}; amx = min( 100, len( aux))
def auxKey( idx=0 ) : return (( 'X{:02d}'.format( int( idx ))) if (( idx >= 0 ) and ( idx <= 99 )) else None )
def printObject( obj=None, num=0, typ=0, **kwargs ) :
  vls = { 'num':num, 'typ':typ, 'eta':0.0, 'phi':0.0, 'pt':0.0, 'jmas':0.0, 'ntrk':0.0, 'btag':0.0, 'hadem':0.0, 'dum1':0.0, 'dum2':0.0 }
  for k in kwargs.keys() :
    if ( k in vls ) : vls[k] = kwargs[k]
  for i in range( amx ) :
    key = auxKey( i ); atr = str( aux[i] ); vls[key] = float( 'NAN' )
    if ( len( atr ) == 0 ) : continue
    try : vls[key] = float( eval(( 'obj.' + atr ), {'__builtins__':{}}, {'obj':obj} ))
    except : pass
  handle.write(( '{num:4d} ' + ( '{typ:4d} {eta:8.3F}' if isinstance( typ, int ) else '{typ:6.1F} {eta:6.3F}' ) +
    ' {phi:8.3F} {pt:7.2F} {jmas:7.2F} {ntrk:6.1F} {btag:6.1F} {hadem:7.2F} {dum1:6.1F} {dum2:6.1F}' +
    ( ''.join(( ' {' + auxKey( i ) + ':10.3E}' ) for i in range( amx ))) + '\\n' ).format( **vls ))
  return num + 1
tree = ROOT.TFile.Open( '${pth}${run}/${tag}_delphes_events.root', 'READ' )
handle = open( '${pth}${run}/delphes_events.lhco', 'w' )
nE = 1; handle.write( '   #  typ      eta      phi      pt    jmas   ntrk   btag  had/em   dum1   dum2' + ( ''.join(( '        ' + auxKey( i )) for i in range( amx ))) + '\\n' )
for e in tree.Delphes :
  nO = 1; handle.write(( '{:4d} {:13d} {:8d}' + ( '      {:10.3E}' if wgt else '' ) + '\\n' ).format( 0, nE, 0, ( e.Weight[ e.Weight_size - 1 ].Weight if wgt else None )))
  for o in e.Photon :
    nO = printObject( o, nO, 0, eta=o.Eta, phi=o.Phi, pt=o.PT, hadem=o.EhadOverEem )
  for o in e.Electron :
    nO = printObject( o, nO, 1, eta=o.Eta, phi=o.Phi, pt=o.PT, ntrk=float( o.Charge ), hadem=o.EhadOverEem )
  muons = [ o for o in e.Muon ]; taus = [ ]; jets = [ ]
  for o in e.Jet :
    jets.append( o ) if ( o.TauTag == 0 ) else taus.append( o )
  iJet0 = nO + len( muons ) + len( taus )
  for o in muons :
    iJet = 0; drMin = None; ptSum = 0.0; etSum = 0.0
    for i in range( len( jets )) :
      drJet = o.P4().DeltaR( jets[i].P4())
      if (( drMin is None ) or ( drJet < drMin )) : iJet = i + iJet0; drMin = drJet
    for t in e.Track :
      if ( o.P4().DeltaR( t.P4()) < 0.5 ) : ptSum += t.PT
    for t in e.Tower :
      if ( o.P4().DeltaR( t.P4()) < 0.5 ) : etSum += t.ET
    nO = printObject( o, nO, 2, eta=o.Eta, phi=o.Phi, pt=o.PT, jmas=0.11, ntrk=float( o.Charge ),
      btag=float( iJet ), hadem=( round( ptSum ) + min( 0.99, etSum/o.PT )))
  for o in taus :
    nO = printObject( o, nO, 3, eta=o.Eta, phi=o.Phi, pt=o.PT, jmas=o.Mass, ntrk=float( o.Charge ), hadem=o.EhadOverEem )
  for typ, trk, lst in (( 4, True, jets ), ( 4.1, False, e.GenJet if gen else [ ] ), ( 4.2, False, e.FatJet if fat else [ ] )) :
    for o in lst :
      nO = printObject( o, nO, typ, eta=o.Eta, phi=o.Phi, pt=o.PT, jmas=o.Mass,
        ntrk=( float( sum( 1 for t in e.Track if ( o.P4().DeltaR( t.P4()) < 0.5 ))) if trk else o.NCharged ),
        btag=o.BTag, hadem=o.EhadOverEem )
  for o in e.MissingET :
    nO = printObject( o, nO, 6, phi=o.Phi, pt=o.MET )
  nE += 1
handle.close()
sys.exit(0)
DONE
	(( close $FHO ) && (($? >> 8) == 0 )) } and do {
	open my $FHO, q(|-), q(bash >/dev/null 2>&1);
	print $FHO <<DONE;
cd ${pth}
if [ -e ./${run}/delphes_events.lhco ]; then
  sed -e "s/^/#/g" ./${run}/${run}_${tag}_banner.txt > ./${run}/${tag}_delphes_events.lhco
  # echo "##  Integrated weight (pb)  : ${xpb}" >> ./${run}/${tag}_delphes_events.lhco
  cat ./${run}/delphes_events.lhco >> ./${run}/${tag}_delphes_events.lhco
  gzip ./${run}/${tag}_delphes_events.lhco
  rm -f ./${run}/delphes_events.lhco
fi
DONE
	(( close $FHO ) && (($? >> 8) == 0 )) }}

# Attempts to generates an LHCO event file from a parton-level LHE event file
sub LHE_2_LHCO { my ($fin,$out,$pre) = (@_);
	$fin =~ /(?:^|\/)([^\/]+)\.lhe$/ or die q(bad file name for input);
	open my $FHI, q(<), $fin or die q(cannot open file ${fin} for read);
	( my ($h,$f) = ( &Local::FILE::NEXT([ ( &DEFINED( $out, q(./LHCO))), [ ( &DEFINED( $pre, $1 )), 0, q(lhco) ]]))) or die q(cannot open file for output);
	my ($xsc,@evt); while (<$FHI>) {
		next unless /^\s*<event>/;
		my ($lin) = ( scalar <$FHI> );
		next unless $lin =~ /^\s*\d+/;
		my ($obj) = [ map {(0+ $_ )} ( split q( ), $lin ) ];
		next unless ( @$obj == 6 );
		my ($wgt) = 0+$$obj[2];
		my ($met,$str,$i,@obj) = [0,0];
		while (<$FHI>) {
			last if /^\s*<\/event>/;
			next unless /^\s*(?:-|\+)?\d+/;
			my ($obj) = [ map {(0+ $_ )} ( split ) ];
			next unless (( @$obj == 13 ) && ( $$obj[1] == 1 ));
			my ($pid,$trk,$hft) = (5,0,0);
			LOOP: for my $pth ( [[22],[0,0,0]], [[11],[1,-1,0]], [[13],[2,-1,0]], [[15],[3,-1,0]],
					[[1,2,3,4],[4,1,0]], [[5,6],[4,1,1]], [[21],[4,0,0]] ) {
				for (@{$$pth[0]}) {
					next unless ( abs $$obj[0] == $_ );
					($pid,$trk,$hft) = (@{$$pth[1]});
					if ( $trk < 0 ) { $trk = ( $$obj[0] <=> 0 ); }
					last LOOP; }}
			push @obj, [ $pid, @{(( &ETA_PHI_PTM_MAS([ @$obj[(9,6..8)]] )) || (next))}[0..3], $trk, $hft, (( $pid == 4 ) ? (999.9) : (0)), 0, 0 ];
			if ( $pid == 5 ) { $$met[0] += $$obj[6]; $$met[1] += $$obj[7]; }}
		next unless @obj;
		$xsc += $wgt;
		push @obj, [ 6, 0, @{(( &ETA_PHI_PTM_MAS([ 0, @$met[0,1], 0 ])) || (next))}[1..2], 0, 0, 0, 0, 0, 0 ];
		for ( sort { $$a[0] <=> $$b[0] } @obj ) {
			$str .= ( sprintf qq(%4i %4i %8.3f %8.3f %7.2f %7.2f %6.1f %6.1f %7.2f %6.1f %6.1f\n), (++$i,@$_)); }
		push @evt, [ $wgt, $str ]; }
	my $evt = 0+@evt or die q(no valid event records extracted from file ${fin});
	#print qq(\n   #  typ      eta      phi      pt    jmas   ntrk   btag  had/em   dum1   dum2\n\n);
	print qq(\n# ${evt} EVENT SAMPLE).(($evt==1) ? q() : q(S)).qq( PROCESSED IN TOTAL\n);
	printf qq(\n# %+10.3E PB CROSS SECTION IMPLIES %+10.3E PER PB LUMINOSITY\n), (( &RATIO($xsc,$evt)), ( &RATIO($evt**2,$xsc)));
	for my $i (1..(@evt)) { my ($wgt,$str) = @{$evt[$i-1]};
		printf qq(\n%4i %13d %8d      %10.3E\n), (0,$i,0,( &RATIO($wgt,$evt)));
		print $str; }
	print qq(\n); }

# Returns the validated <MGPythiaCard> tag MG5 library PATH object associated with an open, cued filehandle
sub PYTHIA_PATH { my ($FHI,$pth) = grep { (( ref eq q(GLOB)) or ( &ISA( 1, $_, q(Local::FILE)))) && ((tell $_) >= 0) or (return) } (shift);
	local ($_); while (<$FHI>) {(( m/^\s*#*\s*<\/MGPythiaCard/i ) ? (last) :
		( m/^\s*#*\s*!*\s*DYLD_LIBRARY_PATH=(\/.*?\/)HEPTools\/lib:/i ) && ( $pth = (q(). $1)))}
	( &Local::FILE::PATH($pth)) }

# Returns a Boolean value indicating whether Python and MatPlotLib appear to be suitably configured for piped system calls
sub CAN_MATPLOTLIB { local ($?);
	open my $FHO, q(|-), q(python >/dev/null 2>&1);
	print $FHO <<DONE;
import sys
if (( sys.version_info[0] < 2 ) or (( sys.version_info[0] == 2 ) and ( sys.version_info[1] < 7 ))) : sys.exit( 1 )
import matplotlib as mpl
if (( tuple( map ( int, mpl.__version__.split( '.' ))) + (0,0,0))[0:3] < (1,3,0)) : sys.exit( 1 )
sys.exit( 0 )
DONE
	(( close $FHO ) && (($? >> 8) == 0 )) }

# HERE: Unify with prior
# Returns a Boolean value indicating whether Python and MatPlotLib appear to be suitably configured for piped system calls
sub CAN_XGBOOST { local ($?);
	open my $FHO, q(|-), q(python3 >/dev/null 2>&1);
	print $FHO <<DONE;
import sys
if (( sys.version_info[0] < 2 ) or (( sys.version_info[0] == 2 ) and ( sys.version_info[1] < 7 ))) : sys.exit( 1 )
import matplotlib as mpl
if (( tuple( map ( int, mpl.__version__.split( '.' ))) + (0,0,0))[0:3] < (1,3,0)) : sys.exit( 1 )
import xgboost as xgb
sys.exit( 0 )
DONE
	(( close $FHO ) && (($? >> 8) == 0 )) }

# Returns the sum of a list of values
sub SUM { my ($sum); ($sum += $_) for ( grep { (defined) or (return undef) } (@_)); $sum }

# Returns the product of a list of values
sub PRODUCT { my ($prd) = ((0+@_) ? (1) : (undef)); ($prd *= $_) for ( grep { (defined) or (return undef) } (@_)); $prd }

# Returns the arithmetic mean of a list of values
sub ARITHMETIC { ( &RATIO(( &SUM(@_)),(0+@_))) }

# Returns the geometric mean of a list of values
sub GEOMETRIC { ( exp +( grep { (defined) or (return undef) } ( &ARITHMETIC(
	map {(log)} grep { ($_ > 0) or (return 0) } grep { ((defined) and ($_ >= 0)) or (return undef) } (@_))))[0] ) }

# Returns the harmonic mean of a list of values
sub HARMONIC { ( 1 / +( grep { (defined) or (return undef) } ( &ARITHMETIC(
	map {(1/$_)} grep { ($_ > 0) or (return 0) } grep { ((defined) and ($_ >= 0)) or (return undef) } (@_))))[0] ) }

# Returns the variance of a list of values
sub VARIANCE { my ($avg); (( grep { (wantarray) or (return $_) } map {( &MAX( 0, ( $$_[0] - $$_[1]*$$_[1] )))} [ grep { (defined) or (return undef) }
	(( &ARITHMETIC( map {($_*$_)} grep { (defined) or (return undef) } (@_))), (($avg) = ( &ARITHMETIC ))) ] ), $avg ) }

# Returns the root-mean-square of a list of input values
sub RMS { ( return ((defined) ? (sqrt) : (undef))) for ( &ARITHMETIC( map {($_*$_)} grep { (defined) or (return undef) } (@_))) }

# Returns the cartesian norm of a list of input values
sub NORM { ( return ((defined) ? (sqrt) : (undef))) for ( &SUM( map {($_*$_)} grep { (defined) or (return undef) } (@_))) }

# Returns the minimum of a list of values
sub MIN { ( return ((defined) ? (0+ $_) : (undef))) for ( &CMP( sub ($$) { my ($a,$b) = @_; ($a <=> $b) }, (@_))) }

# Returns the maximum of a list of values
sub MAX { ( return ((defined) ? (0+ $_) : (undef))) for ( &CMP( sub ($$) { my ($a,$b) = @_; ($b <=> $a) }, (@_))) }

# Returns the leading extremal member of a list of values as compared by a user provided subroutine reference
sub CMP { my ($sub,$cmp) = (( grep { ( ref eq q(CODE)) or (return undef) } (shift)),( grep { (defined) or (return undef) } ((shift),(@_))));
	(($cmp) = ( sort $sub ($cmp,$_))) for (@_); $cmp }

# Returns the input or zero if this value is absolutely smaller than a specified threshold
sub FLUSH { ( grep { !(((abs) <= (shift)) && (return 0)) } grep { (defined) or (return undef) } (shift))[0] }

# Returns the first element of a list having a defined value
sub DEFINED { ( undef, ( grep { (defined) && (return $_) } (@_)))[0] }

# Returns the stringwise/defined unique elements of a list of objects
sub UNIQUE { my (%t,$t); map { (wantarray) ? (@$_) : (return $_) } [ grep { !((defined) ? ( $t{$_}++ ) : ( $t++ )) } (@_) ] }

# Returns a Boolean value testing for numerical equality and equality of defined state for a list of objects
sub EQUAL { my ($a) = (shift); for my $b (@_) {((( $a == $b ) && !(( defined $a ) xor ( defined $b ))) or ( return !1 ))}; (1) }

# Returns the safely divided ratio of two values
sub RATIO { my ($n,$d) = grep { (defined) or (return undef) } (shift,shift); ((abs $d) <= (shift)) ? ((abs $n) <= (abs $d)) ? (0) : (shift) : ($n/$d) }

# Rounds numbers to the specified decimal place using "half away from zero"
sub ROUND { my ($val,$dec) = (( grep { (defined) or (return undef) } (shift)),(int shift));
	(( $val <=> 0 )*( int (( abs $val )*( 10**$dec ) + 0.5 ))/( 10**$dec )) }

# Returns the bitwise and, or, and exclusive or of a list of input positive semi-definite integers
sub AND_OR_XOR { my ($a,$o,$x) = (( &MAX(0,( int shift )))x(3)); (($a &= $_),($o |= $_),($x ^= $_)) for ( map {( &MAX(0,(int)))} (@_)); ($a,$o,$x) }

# Returns a list of Boolean flags with optional length specification extracted LSB-first from an input integer
sub INT_2_FLAGS { my ($int,$len,$bas,$flg) = (( map {( &MAX( 0, (int)))} (shift,shift)), 1, [] );
	while (($bas) ? ($len) ? ( @$flg < $len ) : ($bas <= $int ) : ( return undef )) { push @$flg, !( !( $bas & $int )) } 
	continue {( $bas <<= 1 )} ((wantarray) ? (@$flg) : ($flg)) }

# Returns an integer converted from an input LSB-first list of Boolean flags
sub FLAGS_2_INT { my ($flg,$bas,$sum) = (( grep {(( ref eq q(ARRAY)) or ( return undef ))} (shift)), 1, 0 ); 
	for (@$flg) { (($bas) or ( return undef )); (($_) && ($sum += $bas))} continue {( $bas <<= 1 )} ($sum) }

# Returns a verified ordered pair, or undef if unordered, based upon an input [min,max] range
sub ORDERED { my ($min,$max,$eps,$wid) = (( map { ( ref eq q(ARRAY)) ? (@$_)[0,1] : (return) } (shift)),(shift,shift)); (defined $max) && (defined $min) &&
	do { (($max - $min) >= 0) or do { ($eps < 0) or (($min - $max) <= $eps) or (return); ($max) = ($min) = ( &ARITHMETIC($min,$max)) }};
	do { if (defined $max) { $max += $wid } if (defined $min) { $min -= $wid }} if (defined $wid); (wantarray) ? ($min,$max) : [$min,$max] }

# Returns the transformed copy of an input list of values bounded by a user provided [min,max] ordered pair
sub BOUNDED { my ($min,$max) = @{(( &ORDERED(shift)) or (return))}; grep { (wantarray) or (return $_) }
	map { (defined) ? ($min,$_,$max)[ !((defined $min) && ($_ <= $min)) + ((defined $max) && ($_ > $max)) ] : (undef) } ((wantarray) ? (@_) : (shift)) }

# Returns the max of minima and the min of maxima for a list of input [min,max] ordered pairs; Leading input optionally inverts handling
sub MIN_Y_MAX { my ($mode) = (0+ !(!(shift))); map { (wantarray) ? (@$_) : (return $_) } map {[(
	( &MAX ( grep { ($mode) or (defined) } map {($$_[($mode)])} (@$_))),
	( &MIN ( grep { ($mode) or (defined) } map {($$_[(1-$mode)])} (@$_)))
	)[($mode),(1-$mode)]]} [ grep { ( ref eq q(ARRAY)) or (return) } (@_) ] }

# Returns the [min,max] bounding roots (possibly undefined) for an input (physically positive) quadratic polynomial object
sub MAX_O_MIN { my ($q) = map { ( &Local::POLY::OBJECT($_)) || (return) } (shift); map { (wantarray) ? (@$_) : (return $_) } map {[ ((@$_), undef )[
	((@$_==2) ? ($$q[2]<0) ? (0,+1) : (+1,-1) : (0+@$_) ? ($$q[1]<0) ? (-1,0) : (0,-1) : (return)) ]]} [( &QUAD_REAL_ROOTS($q,(shift)))] }

# Returns the integer, ceiling, or floor truncation of a floating point number
sub INT_CEILING_FLOOR { my ($val,$int,$icf) = (( map { (0+($_),(int)) } grep { (defined) or (return undef) } (shift)),
	((shift) <=> 0 )); $int + $icf*(( $val <=> $int ) == ($icf)) }

# Returns the integer quotient and floating point remainder of a pair of numbers
sub INT_QUOTIENT { map { (wantarray) ? ($_,($_[0] - $_*$_[1])) : (return $_) }
	grep { (defined) or (return (undef,undef)) } ( &INT_CEILING_FLOOR(( &RATIO(@_[0,1])),($_[2]))) }

# Returns two raised to a bounded input integer power
{; my ($exp); sub INT_EXP_TWO { $exp ||= do { my ($i) = 1; [ 1, ( map {($i*=2)} (1..( int BMX ))) ] }; $$exp[( int shift )] }}

# Returns the base-two logarithm floor of an input positive semi-definite integer, as well as the ceiling in list context
sub INT_LOG_TWO { my ($val,$flr,$clg,$try) = (( int shift ), -2, 1+( int BMX ));
	do { ($flr,$clg) = @{ ( [$try,$try], [$flr,$try], [$try,$clg] )[ ( &INT_EXP_TWO(
		$try = ( int (( $flr + $clg )/2)))) <=> $val ] }} while ($clg > (1+$flr));
	grep {((wantarray) or (return $_))} map {((($_ < 0) or ($_ > ( int BMX ))) ? (undef) : ($_))} ($flr,$clg) }

# Returns the invariantly ordered real roots of a quadratic polynomial object; Optional parameter imposes positive semi-definite discriminant
sub QUAD_REAL_ROOTS { my ($q) = map { ( &Local::POLY::OBJECT($_)) or (return) } (shift); (@$q == 3) ? (
	map { my ($t) = $_; map { (-1*$$q[1] + $_*$t)/(2*$$q[2]) } (($$q[2] > 0) ? (-1,+1) : (+1,-1)) }
	map { ($_ >= 0) ? (sqrt) : do { my ($t) = (shift); (($t < 0) || ((-1*$_) <= $t)) ? (0) : () }}
		($$q[1]*$$q[1] - 4*$$q[0]*$$q[2])) : (@$q == 2) ? ((-1)*($$q[0]/$$q[1])) : (@$q == 0) ? (undef) : () }

# Returns the reduced discriminant of a set of input polynomial (up to quartic order) coefficients
sub QUARTIC_DISCRIMINANT { ( pop @_ ) while ((@_) && ($_[-1] == 0)); (@_ > 5) ? (undef) : (@_ == 5) ? (
	+ 256*$_[0]**3*$_[4]**3 - 192*$_[0]**2*$_[1]*$_[3]*$_[4]**2 - 128*$_[0]**2*$_[2]**2*$_[4]**2 + 144*$_[0]**2*$_[2]*$_[3]**2*$_[4] - 27*$_[0]**2*$_[3]**4
	+ 144*$_[0]*$_[1]**2*$_[2]*$_[4]**2 - 6*$_[0]*$_[1]**2*$_[3]**2*$_[4] - 80*$_[0]*$_[1]*$_[2]**2*$_[3]*$_[4] + 18*$_[0]*$_[1]*$_[2]*$_[3]**3 + 16*$_[0]*$_[2]**4*$_[4]
	- 4*$_[0]*$_[2]**3*$_[3]**2 - 27*$_[1]**4*$_[4]**2 + 18*$_[1]**3*$_[2]*$_[3]*$_[4] - 4*$_[1]**3*$_[3]**3 - 4*$_[1]**2*$_[2]**3*$_[4] + 1*$_[1]**2*$_[2]**2*$_[3]**2 ) :
	(@_ == 4) ? ( $_[1]**2*$_[2]**2 - 4*$_[0]*$_[2]**3 - 4*$_[1]**3*$_[3] + 18*$_[0]*$_[1]*$_[2]*$_[3] - 27*$_[0]**2*$_[3]**2 ) :
	(@_ == 3) ? ( $_[1]**2 - 4*$_[0]*$_[2] ) : (@_ == 2) ? (1) : (@_ == 1) ? ( &RATIO(1,$_[0]**2)) : (0) }

# (siz,min,max,asc,adj,unq); Returns (ascending) indexed [lists] of length siz with (adjacent) (unique) values from min to max
sub PERMUTE { my ($siz,$min,$max,$asc,$adj,$unq,$sub,$prm,@prm) = @_[0..7]; $siz = (int $siz); $min = (int $min); $max = ((defined $max) ? (int $max) : ($min+$siz-1));
	do { @prm = (( ref $unq eq q(ARRAY)) ? (@$unq[0..($max-$min)]) : ((int $unq)x(1+$max-$min))); (--$prm[$_-$min]) for ((@{$prm||[]}),($min..$max)); } if ($unq);
	map { (@$_ > $siz) ? () : (@$_ == $siz) ? (( ref $sub ne q(CODE)) or ( $sub->($_))) ? ($_) : () : ( &PERMUTE(@_[0..6],$_)) } map {[ @{$prm||[]}, $_ ]}
	grep {(($prm[$_-$min]) >= 0 )} ((($asc) ? ( ${$prm||[$min]}[-1] ) : ($min))..(($adj) ? ( &MIN( $max, 1+( &MAX(@{$prm||[$min-1]})))) : ($max))) }

# Returns every reordered [list] derivable from a single input [list]
sub ORDERINGS { my (@lst) = map {(( ref eq q(ARRAY)) ? @$_ : (return))} (shift); map {[ @lst[@$_]]} ( &PERMUTE((0+@lst),0,(@lst-1),!1,!1,1)) }

# Returns every N-tuple [list] derivable from a single input [list]
sub TUPLES { my ($tup,@lst) = ((shift), ( map {(( ref eq q(ARRAY)) ? (@$_) : (return))} (shift))); map {[ @lst[@$_]]} ( &PERMUTE($tup,0,(@lst-1),1,!1,1)) }

# Returns every partition into N pieces derivable from a single input [list] with minimal redundancy
sub PARTITIONS { my ($prt,@lst) = ((shift), ( map {(( ref eq q(ARRAY)) ? (@$_) : (return))} (shift))); map { my ($prm) = $_; [ map {
	my ($i) = $_; [ @lst[ grep {($$prm[$_] == $i)} (0..(@$prm-1)) ]]} (0..($prt-1)) ]} ( &PERMUTE((0+@lst),0,($prt-1),!1,1,!1)) }

# Returns every grouping into sets of N items derivable from a single input [list] with minimal redundancy
sub SETS { my ($set,@lst) = (( &MAX(0,( int shift ))), ( map {(( ref eq q(ARRAY)) ? (@$_) : (return))} (shift))); map { my (@lst) = @$_;
	map { my ($prm) = $_; [ map { my ($i) = $_; [ @lst[ grep {($$prm[$_] == $i)} (0..(@$prm-1)) ]]} (0..((@lst/$set)-1)) ]}
		( &PERMUTE((0+@lst),0,((@lst/$set)-1),!1,1,$set)) } ( &TUPLES($set*( int &RATIO((0+@lst),$set)),\@lst)) }

# Returns the ordered (optionally transformed) grouping into sets of N items derivable from a single input [list], optionally with ungrouped entries trailing
sub GROUPS { my ($grp,$alt,$lst,$sub) = (( map {(( &MAX(0,( int shift @$_ ))),( int shift @$_ ))} map {(( ref eq q(ARRAY)) ? ($_) : [$_,0] )} (shift)),
	[ map {(( ref eq q(ARRAY)) ? (@$_) : ())} (shift) ], ( map {(( ref eq q(CODE)) ? ($_) : ( sub {(shift)} ))} (shift)));
	(( grep {((wantarray) or (return $_))} map { my ($lst) = []; ( push @$lst, ( scalar $sub->([ splice @$_, 0, $grp ]))) while (@$_); $lst }
		[ splice (@$lst,0,$grp*( int (((@$lst < $grp) and (($alt < 0) or (@$lst == $alt))) or ( &RATIO((0+@$lst),$grp))))) ] ), @$lst ) }

# Returns the skip-counted (optionally transformed) dealing into N sets of items derivable from a single input [list], optionally with unspanned entries trailing
sub SPANS { my ($spn,$res,$lst,$sub) = (( map {(( &MAX(0,( int shift @$_ ))),(( shift @$_ ) <=> 0 ))} map {(( ref eq q(ARRAY)) ? ($_) : [$_,0] )} (shift)),
	[ map {(( ref eq q(ARRAY)) ? (@$_) : ())} (shift) ], ( map {(( ref eq q(CODE)) ? ($_) : ( sub {(shift)} ))} (shift))); my (@i) =
	( &RANGE(0,(@$lst-$spn),$spn,0+($res >= 0))); (( grep {((wantarray) or (return $_))} [ map { my ($i) = $_; ( scalar $sub->(
		[ @$lst[( grep {(( $res == 1 ) or ( $_ < @$lst ))} map {($_+$i)} (@i))]])) } (0..($spn-1)) ] ), @$lst[($i[-1]+$spn)..(@$lst-1)] ) }

# Returns the collection of zipped N-item ordered [lists] (optionally transformed) derivable from a [list] of N input [lists]
sub ZIPS { my ($zip,$clp,$sub) = (( grep {(( ref eq q(ARRAY)) or (return))} (shift)), ((shift) <=> 0 ),
	( map {(( ref eq q(CODE)) ? ($_) : ( sub {(shift)} ))} (shift))); ( map { my ($i) = $_;
	( scalar $sub->([ map {(( ref eq q(ARRAY)) ? ($$_[$i]) : ($_))} (@$zip) ])) } (0..((((( sub {
	my ($n) = (shift); do {(( $_ == $n) or ( return 0 ))} for (@_); 0+($n) }, \&MAX, \&MIN )[$clp] ) ->
	( map {(( ref eq q(ARRAY)) ? ( 0+(@$_)) : ())} (@$zip))) - 1 )))) }

# Returns the ordered assignments of N items derivable from a [list] of input [list] objects, optionally unique and defined
sub ASSIGNMENTS { my (@asn); my ($lst,$unq,$def) = (( grep { ( ref eq q(ARRAY)) or (return) } (shift)), ( map { !(!$_) } (shift,shift)));
	for ( map { ( ref eq q(ARRAY)) ? ($_) : [$_] } (@$lst)) { ((@asn) = ( map { my ($t) = $_; map { my ($h) = ( &CLONE($$_[0]));
		(( defined $t ) ? (($unq) && ( $$h{$t}++ )) : ($def)) ? () : [$h,@$_[1..(@$_-1)],$t] } ((@asn) ? (@asn) : ([{}])) }
		((@$_) ? ( &UNIQUE(@$_)) : (undef)))) or (return); } ( grep {( shift @$_ )} (@asn)) }

# Returns a list of integers from some input minimum up to an input maximum, optionally with a specified stride and with a selected capping strategy
sub RANGE { my ($min,$max,$srd,$sgn,@lst) = (( map {((defined) ? (int) : (return))} (shift,shift)), ( map {( $_, ( $_ <=> 0 ))}
	(( int ( &DEFINED((shift),1))) or (return)))); $max += ((0),($srd-$sgn),(0-$sgn))[((shift) <=> 0 )];
	while (($sgn*$min) <= ($sgn*$max)) { push @lst, $min; $min += $srd; } (@lst) }

# Generate a Poisson-distibuted random number via the algorithm of Knuth; reverts to Gaussian limit for large numbers
sub POISSON_RANDOM { my ($l,$k,$p) = (( map {(($_ > 0) ? ($_ > LPR) ? ( return ( &ROUND( scalar &GAUSS_RANDOM(1*$_,(sqrt))))) : ( exp(-1*$_)) : (return 0))} (shift)),0,1);
	do { $k++; $p *= (rand); } while ($p > $l); ($k-1) }

# Generate a Gauss-distributed random number via the algorithm of Box and Muller in polar form
sub GAUSS_RANDOM { my ($m,$s) = (((@_)?0+(shift):(0)),((@_)?0+(shift):(1))); {; (redo) unless
	(( my ($r)) = map {( sqrt ( -2*(log)/($_)))} grep {(($_ > 0) && ($_ < 1))} ( &SUM( map {($_*$_)} (( my (@r)) = ( map {( 2*(rand) - 1 )} (0,1))))));
	( return ( map { grep {((wantarray) or (return $_))} ( $m + $s*$r*$_ ) } (@r))); }}

# Returns template for extraction of requested kinematic components from a list of input lhco objects
{; my ($n,$s); sub OUTPUT_OBJECT { my (%n); $s ||= do { $n = 1; +{ map {(($_)=>($n++))} ( qw( eta phi ptm mas ep0 ep1 ep2 ep3 dm1 dm2 )) }};
	map { my ($q,$j,$w,$x) = (%$_)[0,1]; map {[ $_, $n{$q}++, $q, $j, $w, $x ]} ((( $j = ( &MAX(0,( int $j )))) && (($$s{$q}) or (( $q =~ /^x(\d{2})$/ ) &&
		do { ($w,$x) = ( 1, { aux => 0+$1 } ); ( $n + $1 ) } ))) or ()) } grep {( ref eq q(HASH))} (@{(shift)||[]}) }}

# Tests whether output of a value has been requested
sub OUTPUT_VALUE { (1,1,!1)[0+(0..2)[ ${ (shift) || [] }[0]]] }

# Tests whether a provided object value matches the specifed range of acceptable values; Optionally caps range and defines UNDEF
sub MATCH_VALUE { my ($min,$max,$mod,$val) = (( map { ( ref eq q(ARRAY)) ? (@$_[0..2]) : (undef,undef,undef) } (shift)), (shift));
	(( my ($ref,$cap)),$val) = (( map { ($_,(undef,$min,$max)[0+(0,+1,-1)[$mod]],$$_) } grep { local ($@);
		eval { ($$_,$$_) = (undef,$$_) }; !$@ } grep {( ref eq q(SCALAR))} ($val)), (undef,undef,$val));
	((defined $min) or (defined $max)) ? ((defined $val) or ( return ( defined ( $$ref = $cap )))) : ( return ((wantarray) ? () : (1))); (
	map { ($$_[1]) or (!($$_[0][0+(-1,0,+1)[0+(0,+1,-1)[$mod]]]) && ( defined ( $$ref = $cap ))) }
	map {[ $_, (((defined $min) && (defined $max) && ($max < $min)) ? ($$_[0] || $$_[1]) : ($$_[0] && $$_[1])) ]}
		[ (!(defined $min) || ($val >= $min)), (!(defined $max) || ($val <= $max)), 1 ])[0] }

# Sorts a pair of [lists] by deep asciibetical comparison of the stored value sequence
sub SORT_LIST_ALPHA ($$) { my ($a,$b) = (shift,shift); for (0..(( &MIN(0+@$a,0+@$b))-1)) { (($_) && (return $_)) for ($$a[$_] cmp $$b[$_]) }; 0 }

# Sorts a pair of [lists] by deep numerical comparison of the stored value sequence
sub SORT_LIST_NUM ($$) { my ($a,$b) = (shift,shift); for (0..(( &MIN(0+@$a,0+@$b))-1)) { (($_) && (return $_)) for ($$a[$_] <=> $$b[$_]) }; 0 }

# Returns code to sort lhco objects according to descending transverse momentum and then ascending pseudorapidity magnitude; Leading input optionally reverses ordering
sub SORT_OBJECT_LORENTZ_CODE { my ($sort) = (0,+1,-1)[(shift)] || -1; sub ($$) { my ($a,$b) = (shift,shift); ($sort) *
	(($$a{ptm} <=> $$b{ptm}) || ((abs $$b{eta}) <=> (abs $$a{eta})) || ($$a{mas} <=> $$b{mas}) || ($$b{eta} <=> $$a{eta}) || ($$a{phi} <=> $$b{phi})) }}

# Sort inserts an input value into a presorted input list reference, optionally according to an input sort code reference
sub SORT_INSERT { my ($lst,$flr,$clg,$val,$srt,$try,$mod) = (
	( map {( $_, 0, (0+@$_))} grep {(( &ISA( 0, ($_), q(ARRAY))) or (return))} (shift)), (shift),
	( map {(( ref eq q(CODE)) ? ($_) : ( sub ($$) { my ($a,$b) = @_; ( $a <=> $b ) } ))} (shift)));
	while ($clg > $flr) {( ${ (\$clg,\$flr)[$mod] } = (( $try = ( int (($flr+$clg)/2))) +
		( $mod = 0+(1,1,0)[(($srt) -> ( our ($a,$b) = ($val,$$lst[$try])))] )))}
	( splice @$lst, $clg, 0, $val ); ((wantarray) ? (@$lst) : ($lst)) }

# Returns a boolean response membership of an object in a list of classes
sub ISA { use Scalar::Util qw(blessed); my ($m,$o,@c) = (((shift) <=> 0 ), (shift),
	( map {(( length ref ) ? ( &DEFINED(( blessed $_ ), (ref))) : qq($_))} (@_)));
	do {((( defined blessed $o ) ? (($o)->isa($_)) : (( $m <= 0 ) && ( &UNIVERSAL::isa(
	(( length ref $o ) ? ($o) : ( $m == -1 ) ? qq($o) : (undef)), $_ )))) and ( return 1 ))} for (@c); !1 }

# Creates a deeply nested independent clone of a recursively hashed or arrayed data structure
sub CLONE { use Scalar::Util qw(blessed); my ($ref) = (shift);
	( map { my ($t) = $_; do {( bless $t, $_ )} for ( grep {(defined)} ( blessed $ref )); ($t) } (
	( &ISA( 0, $ref, q(ARRAY))) ? [ map {( &CLONE($_))} (@$ref) ] :
	( &ISA( 0, $ref, q(HASH))) ? { map {(($_) => ( &CLONE($ref->{$_})))} (keys %$ref) } :
	( &ISA( 0, $ref, qw( SCALAR REF ))) ? \( &CLONE($$ref)) :
	( &ISA( 0, $ref, qw(CODE))) || ( !( length ref $ref )) ? ($ref) :
	(undef)))[0] }

# Implements fixed point Y-Combinator for the Currying of a recursive anonymous closure while avoiding cyclic self-reference
sub Y_COMBINATOR { my ($cry) = (shift);
	sub { my ($sub) = (shift); $cry->( sub { $sub->( $sub )->(@_) }) }->(
	sub { my ($sub) = (shift); $cry->( sub { $sub->( $sub )->(@_) }) }) }

# Eliminates those items within a list of object references that are also present within a leading [list]; Also eliminates redundant entries
sub EXCLUDE_OBJECTS { my (%excl) = map {($_ => 1)} @{ (shift) || [] }; grep { !( $excl{$_}++ ) } (@_) }

# Gathers the list of objects satisfying certain inclusion and exclusion criteria
sub INCLUDE_OBJECTS { my ($inc,$exc,$obj) = (( map {( [ grep {($_ >= 0)} @$_ ], [ grep {($_ < 0)} @$_ ] )} map { my ($lvl) = (shift);
	[ grep { ((abs) < $lvl) } @$_ ] } map { (@$_) ? [ map {(int)} @$_ ] : [0] } [ grep {(defined)} @{(shift)||[]} ]), (shift));
	&EXCLUDE_OBJECTS( [ map { @{ $$obj[(abs)] || [] }} @$exc ], ( map { @{ $$obj[$_] || [] }} @$inc )) }

# Extracts objects from a [[list]] by the specified index; Defaults to empty object list
sub INDEXED_OBJECTS { map {(@$_)} grep {( ref eq q(ARRAY))} map { ${ $$_[1]}[ $$_[0]] }
	grep { (defined $$_[0]) && (ref $$_[1] eq q(ARRAY)) && (( abs ( int ((int $$_[0]) + 0.5))) <= (@{ $$_[1]}-1)) }
	map {[ ${ $_[(2*$_)] || []}[0], $_[(1+2*$_)]]} (0..((@_/2)-1)) }

# Extracts values from a [list] by the specified index; Defaults to undefined value
sub INDEXED_VALUES { grep { (wantarray) || (return $_); 1 }
	map { ((defined $$_[0]) && ( ref $$_[1] eq q(ARRAY)) && (( abs ( int ((int $$_[0]) + 0.5))) <= (@{ $$_[1]}-1))) ?
		${ $$_[1]}[ $$_[0]] : (undef) } map {[ ${ $_[(2*$_)] || []}[0], $_[(2*$_+1)]]} (0..((@_/2)-1)) }

# Helper routines for common invocation of photon, lepton, and jet indexed objects
sub IPHO {( &INDEXED_OBJECTS( $_[0]{pho}, $_[1]{pho} ))}; sub ILEP {( &INDEXED_OBJECTS( $_[0]{lep}, $_[1]{lep} ))}; sub IJET {( &INDEXED_OBJECTS( $_[0]{jet}, $_[1]{jet} ))};

# Helper routine for common invocation of particle indexed objects
sub IOBJ {(( &IPHO ),( &ILEP ),( &IJET ))}

# Helper routine for common invocation of met indexed value
sub IMET {( scalar ( &INDEXED_VALUES( $_[0]{met}, $_[1]{met} )))};

# Clips or pads a list of LHCO objects to a specified count
sub CLIP_PAD_OBJECTS { my ($min,$max,$clp) = map { ( ref eq q(ARRAY)) ? @$_[0..2] : (undef,undef,undef) } (shift);
	( return @_ ) unless ( defined ( my $cnt = (undef,$max,$min)[( $clp = 0+(0,+1,-1)[($clp)])] ));
	($clp < 0) ? ($cnt > @_) ? (@_) : (@_)[(0..($cnt-1))] : ((@_), ( map { +{}} (0..($cnt-(1+@_))))) }

# Returns the principal value (in the range 0 to 2Pi) corresponding to an input angle in radians
sub PRINCIPAL_RAD { ( &INT_QUOTIENT((shift),((shift) ? PIE : 2*PIE),-1))[1] }

# Returns the absolute azimuthal angular separation (in the range 0 to Pi) between two input angles in the range -Pi to 2*Pi
sub DELTA_RAD_ABS { ( map { ($_ <= PIE) ? $_ : ( abs (2*PIE-$_)) } map {( abs ($$_[1]-$$_[0]))} [ grep { (defined) or (return undef) } (shift,shift) ])[0] }

# Returns a {hash} lhco object with unified 4-vector and collider kinematics extracted from an input 4-vector [list]
sub LORENTZ_HASH { my ($tvrs,$msls,$ivrt,$flsh) = map { ( ref eq q(ARRAY)) ? (@$_) : () } (shift);
	map { my ($vctr) = ( &LORENTZ($_,$tvrs,$msls,$ivrt,$flsh)) || (return); grep { (wantarray) or (return $_) }
	grep { @$_{ qw( ep0 ep1 ep2 ep3 eta phi ptm mas )} = ( @$vctr[0..3], @{ ( &ETA_PHI_PTM_MAS($vctr)) || [] }[0..3] ); 1 } +{}} (@_) }

# HERE ... setting UNDEF to 0 will be inconsistent for purely longitudinal states ...
# Returns a [eta,phi,ptm,mas] list reference corresponding to a [4-vector] or lhco object
sub ETA_PHI_PTM_MAS { [ map {(0+$_)} map { ( &ISA( 0, $_, q(HASH))) ? @$_{ qw(eta phi ptm mas)} : ( &ISA( 0, $_, q(ARRAY))) ? do {
	my ($ptm) = ( sqrt (($$_[1]*$$_[1])+($$_[2]*$$_[2]))); my ($ptr) = ( &RATIO((( sqrt (($ptm*$ptm)+($$_[3]*$$_[3]))) - $$_[3] ), $ptm ));
	((( $ptm && $ptr ) ? ( -1*( log ( $ptr )), ( &PRINCIPAL_RAD( atan2 ($$_[2],$$_[1])))) : (undef,undef)),
		($ptm), ( sqrt ( &MAX( 0, (($$_[0]*$$_[0])-($ptm*$ptm)-($$_[3]*$$_[3])))))) } : (return undef) } (shift) ] }

# Returns the pseudorapidity - azimuth separation in radians between two [4-vector] references or lhco objects; longitudinal only is trailing parameter
sub DELTA_RPA { my ($obja,$objb) = grep { ($$_[2] > 0) or (return) } map { ( &ETA_PHI_PTM_MAS($_)) || (return undef) } (shift,shift);
	( sqrt (($$objb[0]-$$obja[0])**2 + ((shift) ? (0) : ( &DELTA_RAD_ABS($$objb[1],$$obja[1]))**2 ))) }

# Returns the spherical angular separation (in the range 0 to Pi) between two [4-vector] references or lhco objects; Transverse only is trailing parameter
sub DELTA_RSA { ( return ( atan2 ( sqrt ( &MAX(0,(1 - $_**2))), $_ ))) for map { ( 1 - ( &RATIO(( &LORENTZ_PRODUCT(@$_[0,1],-1,1)), $$_[0][0]*$$_[1][0] ))) }
	[ grep { ($$_[0] > NIL) or (return undef) } map { ( &LORENTZ($_,!(!$_[0]),1,!1)) or (return undef) } (shift,shift) ] }

# Returns the azimuthal angular separation (in the range 0 to Pi) between two [4-vector] references or lhco objects
sub DELTA_PHI { &DELTA_RSA(shift,shift,1) }

# Returns the magnitude of the longitudinal radian separation in pseudorapidity between two [4-vector] references or lhco objects
sub DELTA_ETA { &DELTA_RPA(shift,shift,1) }

# Cuts a list of lhco objects to enforce a specified separation in delta-R from a leading [list]
sub INTER_OBJECT_RPA { my ($idr,$cmp) = (shift,shift); ( &MATCH_VALUE( $idr, undef )) && ( return @_ );
	grep { my ($t) = $_; !(0+( grep { !( &MATCH_VALUE( $idr, $_ )) } map { &DELTA_RPA($t,$_) } @{$cmp||[]} )) } (@_) }

# Cuts a list of lhco objects to enforce a specified mutual separation in delta-R; Targets greater conflicts or lower sort; Optionally cuts all conflicts
sub INTRA_OBJECT_RPA { my ($iso) = map { 0+(( ref eq q(ARRAY)) && (0,1,0)[$$_[2]]) } ( my ($idr) = (shift));
	( &MATCH_VALUE( $idr, undef )) && ( return @_ ); my (@inc,@idx) = ((1)x(@_)); my (@iso) = ( map {[ ((!1)x(@inc)) ]} (@inc));
	for my $i (0..(@inc-1)) { for my $j (0..($i-1)) { $iso[$i][$j] = $iso[$j][$i] = !( &MATCH_VALUE( $idr, &DELTA_RPA($_[$i],$_[$j]))) }}
	while ((@idx) = ( grep {($inc[$_])} (0..(@inc-1)))) { my (@t) = ( map {[ $_, 0+( grep {($_)} ( @{ $iso[$_] }[ @idx ] )) ]} (@idx));
		if ($iso) { (@idx) = ( map {($$_[0])} grep {( $$_[1] == 0 )} (@t)); (last) } else { $inc[ ( map {($$_[0])} grep {(( $$_[1] > 0 ) or (last))}
			( &CMP( sub ($$) { my ($a,$b) = @_; ($$b[1] <=> $$a[1]) }, ( reverse (@t)))))[0]] = !1 }} ( @_[ @idx ] ) }

# Cuts a list of LHCO objects on specified pseudorapidity and transverse momentum limits; Sorts by transverse momentum
sub SELECT_PTM_PRM { my ($pts,$prs) = map { 0+(( ref eq q(ARRAY)) && (0,+1)[$$_[2]]) } ( my ($ptm,$prm,$prf) = (shift,shift));
	grep { !($prf) && do { ( &MATCH_VALUE( $prm, (abs $$_{eta}))) || do { $prf = ($prs == 1); !1 }}}
	grep {( &MATCH_VALUE( $ptm, $$_{ptm} ))} ( sort ${ \( &SORT_OBJECT_LORENTZ_CODE($pts)) } (@_)) }

# Cuts a list of LHCO objects according to specified subtype flags
sub SELECT_SUB { my (@sub) = ( map {( $$_[0] <=> 0 )} @{(shift)||[]} ); grep { my ($sub,$pas) = ($$_{sub},1);
	for my $i (0..(@sub-1)) {(( $pas = (( $sub[$i] == 0 ) or (( $sub[$i] < 0 ) xor ($$sub[$i])))) or (last))}; ($pas) } (@_) }

# Cuts a list of LHCO objects according to specified lepton flavor mixing
sub SELECT_EMT { my ($emt,$inv) = map { ([1,1,1,1],[!1,1,!1,!1],[!1,!1,1,!1],[!1,!1,!1,1])[0+(0..3)[(abs)]], ($_ < 0) }
	map { 0+(( ref eq q(ARRAY)) && $$_[0]) } (shift); grep { ($$emt[0+(0..3)[$$_{typ}]]) xor ($inv) } (@_) }

# Cuts a list of LHCO objects according to specified electric charge sign
sub SELECT_SGN { my ($sgn) = map {(0+ (( ref eq q(ARRAY)) && (0,+1,-1)[(0+$$_[0])]))} (shift);
	grep { ($sgn == 0) or (($sgn*$$_{sgn}) > 0) } (@_) }

# Cuts a list of LHCO objects according to specified adjacent track transverse momentum
sub SELECT_PTC { my ($ptc) = (shift); grep {(($$_{typ} != 2) || ( &MATCH_VALUE( $ptc, $$_{ptc} )))} (@_) }

# Cuts a list of LHCO objects according to specified transverse calorimeter energy to track momentum ratio
sub SELECT_ETR { my ($etr) = (shift); grep {(($$_{typ} != 2) || ( &MATCH_VALUE( $etr, $$_{etr} )))} (@_) }

# Cuts a list of LHCO objects according to specified heavy flavor tagging
sub SELECT_HFT { my ($hft) = (shift); grep {( &MATCH_VALUE( $hft, $$_{hft} ))} (@_) }

# Cuts a list of LHCO objects according to specified electromagnetic fraction
sub SELECT_FEM { my ($fem) = (shift); grep {( &MATCH_VALUE( $fem, $$_{fem} ))} (@_) }

# Cuts a list of LHCO objects according to specified composite track count
sub SELECT_TRK { my ($trk) = (shift); grep {( &MATCH_VALUE( $trk, $$_{trk} ))} (@_) }

# Cuts a list of LHCO objects according to specified muon integration count
sub SELECT_MUO { my ($muo) = (shift); grep {( &MATCH_VALUE( $muo, $$_{muo} ))} (@_) }

# Cuts a list of LHCO objects according to specified values of the two dummy data columns 
sub SELECT_DUM { my (@dum) = @{(shift)||[]}; grep {(( &MATCH_VALUE( $dum[0], $$_{dm1} )) && ( &MATCH_VALUE( $dum[1], $$_{dm2} )))} (@_) }

# Cuts a list of LHCO objects according to specified values of the auxiliary data columns 
sub SELECT_AUX { my (@aux) = @{(shift)||[]}; grep { my ($o,$p) = ($_,1); for (@aux) {(( $p = ( &MATCH_VALUE( $$_[1], $$o{aux}[$$_[0]] ))) or (last))}; ($p) } (@_) }

# Cuts a list of LHCO objects on specified kinematic and tagging selections; Sorts by transverse momentum
sub SELECT_OBJECTS { my ($obj,$slf,$cmp,$i,$j,$k) = ((shift,shift,shift),( &MAX( 0, ( int shift ))),0+(0,1,-1)[(shift)],( int shift )); (
	map {[ &CLIP_PAD_OBJECTS( $$obj{cut}, (@$_)) ]}
	grep { ( $_ = [ &INTRA_OBJECT_RPA( $$obj{sdr}, (@$_)) ] ) if ($i > 0); 1 }
	grep { ( $_ = [ &HEMISPHERES( 0, $$obj{eff}, (@$_)) ] ) if ($i > 0); 1 }
	grep { ( $_ = [ &HEMISPHERES( 1, $$obj{set}, (@$_)) ] ) if ($i > 0); 1 }
	grep { ( $_ = [ &INTER_OBJECT_RPA( $$obj{cdr}, (($cmp) || [] ), (@$_)) ] ) if ($i > 0); 1 }
	grep { ( $_ = [ &SELECT_AUX( [ map {(( m/^x(\d{2})$/ ) ? [ (0+ $1), $$obj{$_} ] : ())} ( keys %{$obj||+{}} ) ], (@$_)) ] ); 1 }
	grep { ( $_ = [ &SELECT_DUM( [ @$obj{( qw( dm1 dm2 ))} ], (@$_)) ] ); 1 }
	map {[ ($j == 0) ? ( &SELECT_ETR( $$obj{etr}, &SELECT_PTC( $$obj{ptc}, &SELECT_SGN( $$obj{sgn}, &SELECT_SUB( [ @$obj{( qw())} ], (($k == 0) ? ( &SELECT_EMT( $$obj{emt}, (@$_))) : (@$_))))))) :
		($j > 0) ? ( &SELECT_MUO( $$obj{muo}, &SELECT_TRK( $$obj{trk}, &SELECT_FEM( $$obj{fem}, &SELECT_HFT( $$obj{hft}, &SELECT_SUB( [ @$obj{( qw( gen fat ))} ], (@$_))))))) : (@$_) ]}
	map {[ &SELECT_PTM_PRM( @$obj{( qw(ptm prm))}, (@$_)) ]} (($slf) || [] ))[0] }

# Returns a code reference capable of classifying input lepton object pairs according to relative sign and flavor
sub SELECT_DIL_CODE { my ($dls,$dlf) = map {( 0+(0,+1,-1)[(0+$_)] )} (shift,shift); sub { my (@lep) = (shift,shift);
	(($dls == 0) or ( $dls*( $lep[0]{sgn} <=> 0 )*( $lep[1]{sgn} <=> 0 ) > 0 )) &&
	(($dlf == 0) or ( map { ($$_[0] > 0) && ($$_[1] > 0) && (($$_[0] == $$_[1]) xor ($dlf < 0)) } [ map {( 0+(0..3)[ $$_{typ} ] )} (@lep) ] )[0] ) }}

# Returns an array reference format copy of a [4-vector] or lhco object modified for transverse, massless, p-inverted, flush kinematics; With mass in list context
sub LORENTZ { my ($vctr,$tvrs,$msls,$ivrt,$flsh) =
	(( map {[ map {(0+$_)} (( &ISA( 0, $_, q(ARRAY))) ? @$_[0..3] : ( &ISA( 0, $_, q(HASH))) ? @$_{ qw(ep0 ep1 ep2 ep3)} : (return)) ]} (shift)), (@_));
	do { do { $$vctr[0] = ( sqrt ( &MAX( 0, (($$vctr[0]*$$vctr[0])-($$vctr[3]*$$vctr[3]))))); } unless ($msls); $$vctr[3] = 0; } if ($tvrs);
	do { $$vctr[0] = ( sqrt (($$vctr[1]*$$vctr[1])+($$vctr[2]*$$vctr[2])+($$vctr[3]*$$vctr[3]))); } if ($msls);
	do { do { $$vctr[$_] *= -1; } for (1..3) } if ($ivrt); map { ((wantarray) ? ($vctr,$_) : (return $vctr)) } ($msls) ? (0) :
	do { my ($esqr,$psqr) = (($$vctr[0]*$$vctr[0]), (($$vctr[1]*$$vctr[1])+($$vctr[2]*$$vctr[2])+($$vctr[3]*$$vctr[3])));
		map { (((defined) and ($_ < 1) || (($flsh >= 0) && (( &FLUSH(($_-1),$flsh*$flsh)) == 0 ))) ? ($flsh < 0) ? (undef) :
		do { do { my ($s) = (sqrt); do { $$vctr[$_] *= $s; } for (1..3) } unless ($_ == 1); (0) } : ( sqrt ( $esqr - $psqr ))) } ( &RATIO($esqr,$psqr)) }}

# Returns the vector sum over a list of [4-vector] array references or lhco objects
sub LORENTZ_SUM { my ($tvrs,$msls,$ivrt,$flsh) = map { ( ref eq q(ARRAY)) ? (@$_) : () } (shift);
	( my (@vcts) = map { ( &LORENTZ($_,$tvrs,$msls,$ivrt,$flsh)) || (return) } (@_)) or (return);
	( &LORENTZ([ map { my ($i) = $_; &SUM( map { $$_[$i] } (@vcts)) } (0..3) ])) }

# Returns the massless vector difference between a pair of [4-vector] array references or lhco objects
sub LORENTZ_DIFFERENCE { ( &LORENTZ(( scalar &LORENTZ_SUM((undef), (shift), ( scalar &LORENTZ(shift,!1,!1,1)))), (!1,1,!1))) }

# Returns the inner product of a pair of [4-vector] array references or lhco objects; Optional parameters allow Euclidean metric & pos-semi-def filter
sub LORENTZ_PRODUCT { my ($vcta,$vctb,$mtrc,$psdf) = (( map { ( &LORENTZ($_,!1,!1,!1,-1)) or (return undef) } (shift,shift)),
	(((shift) <=> (0)) || (-1)), !(!(shift))); (return $_) for map { ($psdf) ? ( &MAX(0,$_)) : ($_) }
		( $$vcta[0]*$$vctb[0] + ($mtrc)*( $$vcta[1]*$$vctb[1] + $$vcta[2]*$$vctb[2] + $$vcta[3]*$$vctb[3] )) }

# Returns the [3-vector] angle Lorentz-rotated transformation of a [4-vector] array reference or lhco object
sub LORENTZ_ROTATE { my ($vctr,$lrot,$cost,$sint) = (( scalar &LORENTZ(shift)), ( map { my ($t) = $$_[0]; ( [ map {( &RATIO($_,$t))} (@$_) ], (cos $t), (sin $t)) }
	( scalar &LORENTZ([ 0, ( map { @{ ( ref eq q(ARRAY)) ? $_ : [] }[0..2] } (shift)) ], !1, 1, !1 )))); do { ( return ((defined $vctr) ? ( &LORENTZ(
		${ ( &Local::MATRIX::INNER_PRODUCT([$vctr],( scalar &Local::MATRIX::TRANSPOSE($_)))) || [] }[0] )) : ( &Local::MATRIX::OBJECT($_)))) } for
	[ map { my ($i) = $_; [ map { my ($j) = $_; (($i==0)||($j==0)) ? 0+($i==$j) : (($$lrot[$i]*$$lrot[$j]*(1-$cost)) + (($i==$j) ? ($cost) :
		((0,+1,-1)[($i-$j)%(3)]*$sint*$$lrot[(2-$i-$j)%(4)]))) } (0..3) ] } (0..3) ] }

# Returns the [3-vector] velocity Lorentz-boosted transformation of a [4-vector] array reference or lhco object
sub LORENTZ_BOOST { my ($vctr,$lbst,$gmma) = (( scalar &LORENTZ(shift)), ( map {( $_, (( &RATIO( 1, ( &INVARIANT_MASS($_)))) || (return)))}
	[ 1, ( map {(-1*$_)} map { @{ ( ref eq q(ARRAY)) ? $_ : [] }[0..2] } (shift)) ] )); do { ( return ((defined $vctr) ? ( &LORENTZ(
		${ ( &Local::MATRIX::INNER_PRODUCT([$vctr],$_)) || [] }[0] )) : ( &Local::MATRIX::OBJECT($_)))) } for
	[ map { my ($i) = $_; [ map { my ($j) = $_; ($i!=0)*($i==$j) + $gmma*$$lbst[$i]*$$lbst[$j]*(1-($i!=0)*($j!=0)/(1+$gmma)) } (0..3) ] } (0..3) ] }

# Returns a condensed list merged by angular proximity from a list of [4-vector] array references or lhco objects; [size,mode,sort] are leading parameters
sub LORENTZ_MERGE { my ($size,$mode,$sort,@vcts) = (( map { @{ ( ref eq q(ARRAY)) ? ($_) : ( return ( &LORENTZ_CLIP($_,@_))) }[0..2] } (shift)),
	( grep {($$_[0] > 0)} map { ( &LORENTZ($_)) or (return) } (@_))); $sort = 0+(0,0,-1)[$sort]; while (
	((@vcts) = ( &LORENTZ_OBJECT_SORT(-1,@vcts))) && ($size >= 0) && (@vcts > $size)) { my ($vctr) = ( pop @vcts ); do { $vcts[$$_[0]] =
		( scalar &LORENTZ_SUM(undef,$vctr,$vcts[$$_[0]])) } for grep {(defined)} ( &CMP( sub ($$) { my ($a,$b) = @_; ($$a[1] <=> $$b[1])*($sort) },
		( grep {(defined $$_[1])} map {[ $_, ( scalar &{ ( undef, \&DELTA_RPA, \&DELTA_RSA, \&DELTA_PHI, \&DELTA_ETA, sub { undef } )[$mode] ||
			sub { $_ }} ($vctr,$vcts[$_])) ]} (0..(@vcts-1))))); } (@vcts) }

# Returns an identity preserving truncated list of [4-vector] array references or lhco objects; size is leading parameter
sub LORENTZ_CLIP { my ($size,@vcts) = (( map { (defined) ? 0+(int) : (return @_) } (shift)),
	( grep {( ${ ( &LORENTZ($_)) or (return) }[0] > 0 )} (@_)));
	( &LORENTZ_OBJECT_SORT(-1,@vcts))[0..((($size < 0) ? 0+@vcts : ( &MIN($size,0+@vcts))) - 1 )] }

# Physically sorts a list of [4-vector] references or lhco objects; Leading input optionally reverses ordering
sub LORENTZ_OBJECT_SORT { my ($sort) = \( &SORT_OBJECT_LORENTZ_CODE(shift)); map {( $$_{obj} )} sort $$sort
	map { my ($obj) = $_; grep { @$_{ qw(obj eta phi ptm mas)} = ( $obj, @{ ( &ETA_PHI_PTM_MAS($obj)) || (return) } ); 1 } +{}} (@_) }

# Returns the missing transverse energy (met_tot,met_x,met_y,0) components for a list of [4-vector] or lhco objects
sub MET { ( scalar &LORENTZ(( scalar &LORENTZ_SUM(undef,@_)), (1,1,1))) }

# HERE ... deal globally with massless conversion question ... TRANSVERSE MASS, ENERGY, ETC.
# Returns the scalar sum of transverse momenta (Default) or energies for a list of transverse [4-vector] or lhco objects
sub MHT { ${ ( &LORENTZ_SUM([1,!( map { ( ref eq q(ARRAY)) ? ($$_[0]) : ($_) } (shift))[0],!1],@_)) || (return undef) }[0] }

# Returns the transverse energy computed for a list of massless [4-vector] momenta components or lhco objects
sub TRANSVERSE_ENERGY { ${ ( &LORENTZ(( scalar &LORENTZ_SUM(undef,@_)), (1,!1,!1))) || (return undef) }[0] }

# Returns the invariant mass computed from the vector sum over a list of [4-vector] momenta components or lhco objects
sub INVARIANT_MASS { ((undef), ( &LORENTZ_SUM(undef,@_)))[-1] }

# Returns the transverse mass computed for [list]s of massless transverse [4-vector] momenta components or lhco objects for the visible/invisible systems
sub TRANSVERSE_MASS { ((undef), ( &LORENTZ_SUM([1,!1,!1],(shift,shift))))[-1] }

# Returns the s-transverse mass (MT2) computed for a pair of massless transverse [4-vector] momenta components or lhco objects
sub S_TRANSVERSE_MASS { ( sqrt ((2)*( &LORENTZ_PRODUCT(( map { ( &LORENTZ($_,1,1,!1)) or (return undef) } (shift,shift)), +1, 1 )))) }

# Returns the asymmetric s-transverse mass (AMT2) computed for a pair of [4-vector] momenta components and masses; Independent met or undef is leading parameter
sub A_TRANSVERSE_MASS { my ($met,@obj) = (( map { ((defined) ? ( scalar &LORENTZ($_,1,1,!1)) : ( &MET( grep {(defined)} @_[0,1,4,5] ))) or (return undef) } (shift)),(@_));

	do { my ($min,$try,$max,$stp,$rts,$scl,$dsc,$qrt) = my ($min_0,undef,$max_0,$stp_0,$rts_0,undef,undef,undef) = (@$_); LOOP: {; do {
		${ (\$min,\$max)[ ((defined $dsc) && !(defined $stp)) ? do { ($rts) = ( $dsc->($try)); (defined $rts_0) ? 0+($rts < $rts_0) :
			do { ($rts_0) = ($rts); (0) }} : ( $qrt->($try)) ? do { ($dsc,$stp) = (undef,undef); (1) } : (0) ]} = ($try);
		($try) = (defined $stp) ? (($max)*( &RATIO($min_0,$max))**( &RATIO(( &MAX(0,--$stp)),$stp_0))) : ($min + $max)/2; }
	until (((defined $rts) and ($rts == ($rts_0 - 1)) || ($rts_0 == 0)) or (($max - $min)/2 <= EPS*( sqrt ($try/$scl))));
	return (defined $dsc) ? do { ((undef,$min,$try,$max,$stp,$rts) = (undef,$min_0,undef,undef,$stp_0,$rts_0) = (defined $stp) ? ($max == $max_0) ? (!1) :
		(1,$max,$max,$max_0) : (1,$min_0,$min_0,$max,( &ROUND((1 + (1-( &RATIO($min_0,$max)))*TRY)*TRY))))[0] ?
			(redo LOOP) : (undef) } : ($try <= BIG) ? ( sqrt ($try*$scl)) : (undef) }} for # Bisection to AMT2

	map { my (@leg) = @$_; my ($min,$scl) = ($leg[0][2][0],2*$leg[0][0][10]*$leg[1][0][10]); (($min), my ($max,$dsc)) =
		do { my ($max) = ( &MIN( grep {(defined)} ($leg[0][2][1],$leg[1][2][1]))); map { my ($max,$dsc) = (@$_);
			(( map {( &FLUSH(( &MAX(0,$_)), NIL ))} ( &ORDERED([($min),(($dsc) ? (defined $max) ? ( &MIN($max,BIG)) : (BIG) : ($max))],-1))),($dsc)) }
			grep { (( &ORDERED([$min,$max],NIL)) && ($min <= BIG)) or (return undef); 1 }
		0+( map { my ($trg,$prj) = @leg[$_,(1-$_)]; ($$trg[0][0][5]) ? ( grep {((( &MAX(0,$_)) - ($$trg[1][1])) <= NIL )} # Type III unbalanced
			map {(( sqrt ($$trg[0][2]**2+$$_[1]**2+$$_[2]**2)) - ($$trg[0][3][1]*$$_[1]+$$trg[0][3][2]*$$_[2]))}
			( scalar $$prj[5]->($$trg[2][0],1))) : () } (0,1)) ? [$min,!1] :
		do { my ($max) = ( &MIN( map { my ($trg,$prj) = @leg[$_,(1-$_)]; ($$trg[0][0][0]) ? () : ( # Type I unbalanced
			map {( $$trg[0][7]*(( &MAX(0,$_)) - ( $$trg[1][0]->EVALUATE(0))))}
			map {(( sqrt ($$trg[0][2]**2+$$_[1]**2+$$_[2]**2)) - ($$trg[0][1][1]*$$_[1]+$$trg[0][1][2]*$$_[2]))}
				( scalar $$prj[5]->($$prj[2][0],1))) } (0,1))); (defined $max) and ((($max - $min) <= NIL ) && [$min,!1]) ||
				(($max <= BIG) && !($leg[0][0][0][0]) && !($leg[1][0][0][0]) && [$max,!1]) } ||
		( undef, ( map { my ($trg,$prj) = @leg[$_,(1-$_)]; ($$trg[4]) ? do { # Type II unbalanced
			my (@t) = (( $$trg[4]->(undef,undef)), ( $$trg[4]->(undef,( $$prj[3]->(undef,1)))));
			my ($dsc) = -1*(($t[0][2]*$t[1][0]-$t[1][2]*$t[0][0])**2 + ($t[0][2]*$t[1][1]-$t[1][2]*$t[0][1])*($t[0][0]*$t[1][1]-$t[1][0]*$t[0][1]));
			map { ($min,$max) = ( &MIN_Y_MAX(!1,[$min,$max],$_)); (($$prj[4]) || ((defined $$_[0]) &&
				((( 2*($t[0][2]*$t[1][0]+$t[0][0]*$t[1][2]) - $t[0][1]*$t[1][1] )->EVALUATE($$_[0])) <= NIL ))) ? [$min,!1] : [$max,$dsc] }
			map { ( &MAX_O_MIN($_,NIL)) or (return undef) } grep { (pop) if (( abs $$_[2]) <= NIL ); 1 } map {[ @$_[0..2]]}
				(($$prj[4]) ? ($dsc) : ($t[1][1]**2-4*$t[1][0]*$t[1][2])) } : () } (0,1)))[-1] || [$max,1] };
		[ ($min,$max)[0,0,1], (undef,undef,$scl), ( map { my (@u) = map {( &Local::POLY::OBJECT($_,NIL))} (@$_);
			(($dsc) ? ($leg[0][0][0][5] || $leg[1][0][0][5]) ? sub {(0)} : do { my ($sub) = # Closure "dsc": degenerate root count
				map {( &Local::POLY::REAL_ROOTS($_,undef))} ((ref $dsc) ? ($dsc) : ( &Local::POLY::OBJECT(( &QUARTIC_DISCRIMINANT(@u)),NIL)));
				sub { 0+( $sub->([ @{ ( &ORDERED([(shift),($max)],NIL)) || (return 0) }[0,1], 1 ])) }} : (undef)),
			sub { my ($s) = (shift); do { do { (return 1) if (($_ > 0) or !(defined)) } # Closure "qrt": Boolean intersection status
				for ( &Local::POLY::REAL_ROOTS(@$_)) } for map { my ($o) = ( &MIN((NIL)*( &FLUSH( &Local::VECTOR::NORM($$_[0]), ONE )), (EPS)));
				($o > 0) ? ([($$_[0]+$o),$$_[1]],[($$_[0]-$o),$$_[1]]) : ($_) } [ ( &Local::POLY::OBJECT([ map {( $_->EVALUATE($s))} (@u) ],NIL)),
				(( &ORDERED(( scalar &MIN_Y_MAX( !1, ( map {( $leg[$_][6])->($s,$_)} (0,1)))),NIL,EPS)) or (return !1)) ] ; !1 }}
			map { my ($v) = $_; [ ($$v[2]*$$v[10]-$$v[4]**2), ($$v[0]*$$v[10]+$$v[2]*($$v[7]+$$v[9])-2*$$v[3]*$$v[4]), # u_i
				($$v[0]*($$v[7]+$$v[9])+$$v[2]*($$v[6]-$$v[8])-$$v[3]**2-2*$$v[1]*$$v[4]),
				($$v[0]*($$v[6]-$$v[8])+$$v[2]*$$v[5]-2*$$v[1]*$$v[3]), ($$v[0]*$$v[5]-$$v[1]**2) ]}
			map { my ($h,$p) = @$_; [ 2*($$h[0]*$$p[1]-$$p[0]*$$h[1]), 1*($$h[0]*$$p[2]-$$p[0]*$$h[2]), 2*($$h[0]*$$p[3]-$$p[0]*$$h[3]), # v_i
				2*($$h[0]*$$p[4]-$$p[0]*$$h[4]), 1*($$h[0]*$$p[5]-$$p[0]*$$h[5]), 2*($$h[1]*$$p[2]-$$p[1]*$$h[2]), 4*($$h[1]*$$p[4]-$$p[1]*$$h[4]),
				2*($$h[1]*$$p[5]-$$p[1]*$$h[5]), 2*($$h[2]*$$p[3]-$$p[2]*$$h[3]), 4*($$h[3]*$$p[4]-$$p[3]*$$h[4]), 2*($$h[3]*$$p[5]-$$p[3]*$$h[5]) ]}
		[ ( $leg[0][3]->(undef,0)), ( $leg[1][3]->(undef,1)) ] ) ] } map {[ sort { our ($a,$b); ($$b[2][0] <=> $$a[2][0]) } (@$_) ]} # Event leg unification

	map {[ grep { my ($obj,$del) = (@$_); # Event leg preparation 0:obj, 1:del, 2:dot, 3:ell, 4:bar, 5:org, 6:dsh
		my ($shp) = (($$obj[0][3]) || ($$obj[0][5])) ? [ (1-$$obj[3][1]**2), -1*($$obj[3][1]*$$obj[3][2]), (1-$$obj[3][2]**2), # Conic shepherd
			-1*($$del[1]*$$obj[3][1]), -1*($$del[1]*$$obj[3][2]), -1*($$del[1]**2-$$obj[2]**2) ] : (undef);
		my ($dsc) = (($$obj[0][5]) ? (($$obj[0][2]) ? (undef) : ($$obj[6]**2)) : # Conic discriminant
			($$obj[0][0]) ? ($$del[3]*$$del[5]) : (($$obj[0][1]) ? (undef) : ($$obj[5]**2)));
		my ($dot) = grep { (defined $$_[0]) or (return undef) } [ # Scale bounds "dot": [(M_Y^{dot+/-})^2/(2*E_V*E_V')]
			map { (defined) ? $$obj[7]*(( &MAX(0,$_)) - ( $$del[0]->EVALUATE(0))) : (undef) } (
			($$obj[0][5]) ? ($$del[1])[0,0] : ($$obj[0][2] && $$obj[0][4]) ? (0,undef) : !($$obj[0][0]) ? ($$obj[2]*$$obj[5],undef) :
			(( &QUAD_REAL_ROOTS([ -1*($$del[5]*$$obj[2]**2+($$del[1]*$$obj[5])**2), 2*$$del[1]*$$del[6], -1*$$obj[6]**2 ], -1 )), undef, undef )[0,1] ) ];
		my ($ell) = do { my (@ell) = ($$obj[0][5]) ? (@$shp) : do { my ($t) = ($$del[4]*$$del[0]-$$del[1]); # Closure "ell" 0:a, 1:b, 2:c, 3:d, 4:f, 5:g
			($shp) ? ( 0, 0, 0, ($$del[2][1])/2, ($$del[2][2])/2, $t ) : (
			($$del[2][1]**2+$$del[3]*(1-$$obj[1][1]**2)), ($$del[2][1]*$$del[2][2]-$$del[3]*$$obj[1][1]*$$obj[1][2]), ($$del[2][2]**2+$$del[3]*(1-$$obj[1][2]**2)),
			($$del[2][1]*$t-$$del[3]*$$del[0]*$$obj[1][1]), ($$del[2][2]*$t-$$del[3]*$$del[0]*$$obj[1][2]), ($t**2-$$del[3]*($$del[0]**2-$$obj[2]**2))) };
			sub { ( map { (shift) ? [ (@$_[0..2]), -1*$$obj[7]*($$_[0]*$$obj[9][1]+$$_[1]*$$obj[9][2]+$$_[3]), -1*$$obj[7]*($$_[2]*$$obj[9][2]+$$_[1]*$$obj[9][1]+$$_[4]),
				$$obj[7]**2*($$_[0]*$$obj[9][1]**2+2*$$_[1]*$$obj[9][1]*$$obj[9][2]+$$_[2]*$$obj[9][2]**2+2*$$_[3]*$$obj[9][1]+2*$$_[4]*$$obj[9][2]+$$_[5]) ] : ($_) }
				map { my ($s) = (shift); (defined $s) ? [( map {( $_->EVALUATE($s))} (@$_))] : ($_) } [( map {( &Local::POLY::OBJECT($_))} (@ell))] )[0] }};
		my ($bar) = ((!$shp) || ($$obj[0][5])) ? (undef) : do { my ($t) = ($$del[0]-$$del[1]); sub { # Closure "bar": linear intersection finder
			my ($t,$s) = map { (defined) ? ( $t->EVALUATE($_)) : ($t,( grep {( ref eq q(ARRAY))} (shift))) } (shift); my ($p) = (shift); (
			map {[ undef, ( map {( &ARITHMETIC(@$_))} (@$_)), undef ]} [ reverse map { my ($m) = ( pop @$_ ); (ref $t) && (($_) = [( map {( &Local::POLY::OBJECT($_))} (@$_))]);
				grep { (wantarray) ? ( return @$_ ) : 1 } grep { ($p) && (($_) = [( map {($$obj[7]*($m-$_))} ( reverse @$_ ))]); 1 }
				map { (ref $t) ? (return $_) : [(( &QUAD_REAL_ROOTS($_,-1)), undef )[0,1]] } [
					($$_[0]*$t**2-2*$$_[3]*$$_[6]*$t+$$_[5]*$$_[6]**2), 2*(($$_[0]*$$_[7]-$$_[1]*$$_[6])*$t-$$_[3]*$$_[7]*$$_[6]+$$_[4]*$$_[6]**2),
					((ref $s) ? ($$_[0]*$$_[7]**2-2*$$_[1]*$$_[6]*$$_[7]+$$_[2]*$$_[6]**2) : ($$del[5])) ] }
				( [@{$s||$shp}[0..5],@{$$del[2]}[1,2],$$obj[9][2]], [@$shp[2,1,0,4,3,5],@{$$del[2]}[2,1],$$obj[9][1]] ) ] )[0] }};
		my ($org) = sub { my ($s,$p) = (( grep { (defined) or (return undef) } (shift)),(shift)); # Closure "org": conic coordinate center
			( return ((wantarray) ? ($$_[2]) : ($_))) for (($bar) ? ( scalar $bar->($s,$p)) : [ (undef), ((defined $dsc) ?
			( map {( &RATIO($_,$dsc))} map { ($$_[1]*$$_[4]-$$_[2]*$$_[3]), ($$_[1]*$$_[3]-$$_[0]*$$_[4]) } ( $ell->($s,$p))) : do { my ($t) =
				( &RATIO(( map {(($$obj[2]**2-$_**2), 2*($_))} (($$obj[0][0]) ? ($$del[1]) : ( $$del[0]->EVALUATE($s)))),NIL,HUG)); map { my ($i) = $_;
				grep { ($p) && (($_) = ($$obj[7]*($$obj[9][$i]-$_))) ; 1 } ($t)*(($$obj[0][0]) ? ($$obj[3][$i]) : ($$obj[1][$i])) } (1,2) } ), (undef) ] ) };
		my ($dsh) = sub { my ($s,$p) = (( grep { (defined) or (return undef) } (shift)),(shift)); # Closure "dsh": conic coordinate bounds
			($bar) ? [( $bar->($s,$p))[0,1]] : ( map { ((defined $dsc) or 0+( &FLUSH($$_[1],NIL))) ? ( scalar &MAX_O_MIN($_,-1)) : do { my ($i) =
				0+(($p) xor ((($$obj[0][0]) ? ($$obj[3][2]) : ($$obj[1][2])) < 0 )); [ (( $org->($s,$p)), (undef))[$i,(1-$i)]] }}
				map {[ ($$_[3]**2-$$_[0]*$$_[5]), 2*($$_[1]*$$_[3]-$$_[0]*$$_[4]), -1*($dsc) ]} ( $ell->($s,$p)))[0] };
		push @$_, ($dot,$ell,$bar,$org,$dsh); 1 } (@$_) ]}

	grep { my ($t) = map { (defined) ? ( scalar &Local::MATRIX::TRANSPOSE( scalar &LORENTZ_ROTATE(undef,[0,0,((PIE/2)-( &PRINCIPAL_RAD($_,1)))]))) : () } # Optimized rotation
		map { ${ ( &ETA_PHI_PTM_MAS((( &Local::VECTOR::INNER_PRODUCT(@$_)) >= 0 ) ? ( scalar &Local::TENSOR::SUM(@$_)) : ( scalar &Local::TENSOR::DIFFERENCE(@$_)))) || [] }[1] } [
		map {[ 0, ( sqrt ( &MAX(0,(1-2*$$_[0]*$$_[3])))), ((0 <=> $$_[1])||(1))*( sqrt ( &MAX(0,(1-2*$$_[2]*$$_[3])))), 0 ]}
		grep { push @$_, ( &RATIO(1,(( sqrt(($$_[0]-$$_[2])**2+4*$$_[1]**2)) + $$_[0] + $$_[2]))); 1 } map { my ($obj,$del) = (@$_); ($$obj[0][5]) ?
			[ (1-$$obj[3][1]**2), -1*($$obj[3][1]*$$obj[3][2]), (1-$$obj[3][2]**2) ] : [ ($$del[2][1]**2+$$del[3]*(1-$$obj[1][1]**2)),
			($$del[2][1]*$$del[2][2]-$$del[3]*$$obj[1][1]*$$obj[1][2]), ($$del[2][2]**2+$$del[3]*(1-$$obj[1][2]**2)) ] } (@$_) ];
		do { for (@$_) { my ($obj,$del) = @$_; for (@$obj[1,9],(($$obj[0][0])?($$obj[3],$$del[2]):())) {
			($_) = ( scalar &LORENTZ( ${ ( &Local::MATRIX::INNER_PRODUCT([$_],$t)) || (return undef) }[0] )) }}} if (defined $t); 1 }

	map {[ grep { my ($obj) = (@$_); push @$_, (($$obj[0][0]) ? do { my ($t) = grep { push @{$$obj[0]}, ( map {( $$_[0] <= EPS**2 )} (@$_)); 1 } [ # Factors "del"
			( map {[ $$_[3]**2, $_ ]} ( scalar &LORENTZ_DIFFERENCE($$obj[1],$$obj[3]))),
			( map {[ ($_ >= -1*NIL) ? ( &MAX(0,$_)) : (return undef) ]} ($$obj[4]**2-($$obj[2]+$$obj[8]*$$obj[6])**2)/(2*$$obj[8])),
			( map {[ ( &MAX(0,($_-$$obj[5]*$$obj[6]))), ($_+$$obj[5]*$$obj[6]), $_ ]} ( &LORENTZ_PRODUCT($$obj[1],$$obj[3],-1,1))) ];
		[ ( &Local::POLY::OBJECT([ -1*(($$obj[4]**2+$$obj[5]**2)/2+$$obj[8]*$$t[2][2]), +1/$$obj[7]])),
			($$t[1][0]+$$obj[2]*$$obj[6]), $$t[0][1], $$t[0][0], +1, $$t[2][0]*$$t[2][1], $$t[2][2]] } : do {
		[ ( &Local::POLY::OBJECT([ -1*($$obj[2]**2+$$obj[5]**2)/2, +1/$$obj[7]])), 0, ( scalar &LORENTZ([])), +1, 0 ] } ); 1 } (@$_) ]}
		# 0:[Gamma], 1:Delta, 2:[Lambda], 3:Omega, 4:Pi, 5:((P_V.P_S)^2-(M_V*M_S)^2), 6:(P_V.P_S)

	grep { my (@t) = (@$_); do { $t[$_][0][7] /= $t[(1-$_)][0][10] } for (0,1); 1 } map {[ map {[ ( # Objects "obj"
		grep { my ($t) = $$_[10]; ( $_ = ( &MAX(0,( &RATIO($_,$t))))) for (@$_[2,4]); ($$_[9]) = ( scalar &LORENTZ([ map {($_/$t)} (@$met) ])); 1 }
		grep { !($$_[0][0]) or do { ($$_[8]) = ( my $t = $$_[3][0] )/$$_[10]; (($$_[3],$$_[6]) =
			( &LORENTZ([ 1, ( map {($_/$t)} (@{$$_[3]}[1..3])) ],!1,!1,!1,EPS))) && do { ( $$_[0][2] = ($$_[6] == 0)); 1 }}}
		grep { ( $$_[1] = ( scalar &LORENTZ($$_[1],1,!1,!1))) && ($$_[1][0] > 0) && do { ($$_[10]) = ($$_[7]) = ( my $t = $$_[1][0] );
			(($$_[1],$$_[5]) = ( &LORENTZ([ 1, ( map {($_/$t)} (@{$$_[1]}[1,2])), 0 ],!1,!1,!1,EPS))) && do { ( $$_[0][1] = ($$_[5] == 0)); 1 }}}
		grep { !( $$_[0][0] = ((defined $$_[3]) || (undef $$_[4]))) or ( $$_[3] = ( scalar &LORENTZ_BOOST($$_[3],[0,0,($$_[1][3]/$$_[1][0])]))) && ($$_[3][0] > 0) }
		grep { ( $$_[1] = ( scalar &LORENTZ($$_[1]))) && ($$_[1][0] > 0) } ($_))[0] or (return undef) ]} (@$_) ]}
		# 0:[2-Step,M_V~0,M_S~0,S_z~0,H~S,V~S], 1:[P_V]/E_V, 2:M_H/E_V, 3:[P_S]/E_S, 4:M_X/E_V, 5:M_V/E_V, 6:M_S/E_S, 7:E_V/E_V', 8:E_S/E_V, 9:[MET]/E_V, 10:E_V

	[[[], @obj[0,2,4,6]], [[], @obj[1,3,5,7]]]; } # Input kinematics

# Returns the tri-jet invariant mass computed for a list of massless [4-vector] momenta components or lhco objects
sub TRI_JET_MASS { my ($mlim,@vcts) = ((shift),( grep {($$_[0] > 0)} map { ( &LORENTZ($_,!1,1,!1)) or (return undef) } (@_)));
	(( map { &INVARIANT_MASS(@$_[0,1]) } map { my ($dvct) = ( scalar &LORENTZ_SUM(undef,@{ $$_[0]})); ( sort { our ($a,$b); ($$a[2] <=> $$b[2]) }
		grep { (defined $$_[2]) } map {[ $dvct, $_, &DELTA_RPA($dvct,$_) ]} ( &EXCLUDE_OBJECTS($$_[0],@vcts)))[0] }
		( sort { our ($a,$b); ($$a[2] <=> $$b[2]) } grep { ( &MATCH_VALUE($mlim,$$_[1])) && (defined $$_[2]) }
			map {[ $_, &INVARIANT_MASS(@$_), &DELTA_RPA(@$_) ]} ( &TUPLES(2,\@vcts)))[0] ), (undef))[0] }

# Returns the tau-tau invariant mass computed for jet plus lepton pair [4-vector] momenta components or lhco objects
sub TAU_TAU_MASS { my ($llm); my ($a,$b,$j) = map {( scalar &Local::VECTOR::OBJECT($_))} grep { (( shift @{$_||[]} ) > 0 ) or (return undef) } map {( scalar &LORENTZ($_,1,1,!1))}
	(( map { ($llm) = grep { (defined) or (return undef) } ( &INVARIANT_MASS(@$_)); (@$_) } [ @{(shift)||[]}[0,1]] ), ( scalar &LORENTZ_SUM(undef,@{(shift)||[]})));
	($llm)*( map { ( $_ <=> 0 )*( sqrt abs ) } grep { (defined) or (return undef) } ( &RATIO( -1*(($a^$j).($b^$j)), ( map {($_.$_)} ( scalar $a^$b ))[0] )))[0] }

# Returns the jet Z-balance statistic of an ordered pair of [lists] of [4-vector] momenta components or lhco objects
sub JET_Z_BALANCE { ( map { ($$_[1] - $$_[0]) } [ map { ${ ( &MET(@{$_||[]})) or (return undef) }[0] } (shift,shift) ])[0] }

# Returns the alpha_R razor statistic(s) for a pair of massless [4-vector] momenta components or lhco objects; Independent met or undef is leading parameter
sub ALPHA_R { my ($metv,$obja,$objb) = map {(@$_)} grep { (defined $$_[0]) or ( $$_[0] = &MET(@$_[1,2])); 1 } [(shift,shift,shift)]; map {(@$_)}
	grep { (wantarray) or ( return ( &RATIO( map { (defined) ? ($_)**2 : (undef) } (@$_)))) } [
		( map { ((defined $$_[0]) && (defined $$_[1])) ? (( sqrt (($$_[0])**2 + ($$_[1])**2 ))/(2)) : (undef) }
			[ ( &TRANSVERSE_MASS($metv,$obja)), ( &TRANSVERSE_MASS($metv,$objb))] ), ( &TRANSVERSE_ENERGY($obja,$objb)) ] }

# Returns the alpha_T statistic for a pair of transverse [4-vector] momenta components or lhco objects; Independent (met,mht) or undef are leading parameters
sub ALPHA_T { my ($metv,$mhtv) = map {( ${ ((defined $$_[0]) ? ( scalar &LORENTZ($$_[0],1,1,!1)) : ( &MET(@_[0,1]))) || (return undef) }[0],
	((defined $$_[1]) ? ( &MAX(0,$$_[1])) : ( grep { (defined) or (return undef) } ( &MHT( 1, @_[0,1])))))} [(shift,shift)]; ( map { &RATIO(@$_) }
		map {[ ( &MAX( 0, (($mhtv) - ( abs ( $$_[1] - $$_[0] ))))), 2*( sqrt ( &MAX( 0, (($mhtv)**2 - ($metv)**2)))) ]}
			[ grep { (defined) or (return undef) } ( &MHT( 1, (shift)), &MHT( 1, (shift))) ] )[0] }

# Returns the delta phi statistic for a list of [4-vector] momenta components or lhco objects; Independent met [4-vector] or undef is leading parameter
sub MET_DELTA_PHI { my ($metv) = map { ((defined) ? ( scalar &LORENTZ($_,1,1,!1)) : ( &MET(@_))) or (return undef) } (shift);
	( &MIN( map { &DELTA_PHI($metv,$_) } map { ( &LORENTZ($_,1,1,!1)) or (return undef) } (@_))) }

# Returns the biased delta phi statistic for a list of [4-vector] momenta components or lhco objects; Independent met [4-vector] or undef is leading parameter
sub BIASED_DELTA_PHI { my ($metv) = map { ((defined) ? ( scalar &LORENTZ($_,1,1,!1)) : ( &MET(@_))) or (return undef) } (shift);
	( &MIN( map { ( &DELTA_PHI(( scalar &LORENTZ_SUM([1,1,!1],$metv,$_)), $_ )) } map { ( &LORENTZ($_,1,1,!1)) or (return undef) } (@_))) }

# Returns the cosine of the theta-star angle for a pair of [4-vector] momentum components or lhco objects
sub COSINE_THETA_STAR { (( map {( &RATIO(($$_[0]-$$_[1]),($$_[0]+$$_[1])))} map {[ ( exp(+$_)), ( exp(-$_)) ]} grep {(defined)}
	( &DELTA_ETA( map {(( scalar &LORENTZ($_,!1,1,!1)) or (return undef))} ((shift),(shift))) / 2 )), (undef))[0] }

# Returns the lepton W-projection statistic for a massless [4-vector] momentum component set or lhco object and independent met
sub LEP_W_PROJECTION { ( return &RATIO((($$_[0][0]*$$_[1][0]) - ( &LORENTZ_PRODUCT(@$_[0,1],-1,1))), $$_[0][0]**2, NIL )) for
	[ map { ( &LORENTZ($_,1,1,!1)) or (return undef) } (( scalar &LORENTZ_SUM([!1,!1,!1,-1],@_[0,1])), $_[1] ) ]; }

# Returns the transverse thrust shape statistics for a list of [4-vector] momenta components or lhco objects; Independent mht or undef is leading parameter
sub THRUST_SHAPE { my ($mhtv,@vcts) = (( map { (defined) ? ( &MAX(0,$_)) : ( grep { (defined) or (return (undef,undef)) } ( &MHT( !1, @_))) } (shift)),
	( grep {($$_[2] > 0)} map { ( &ETA_PHI_PTM_MAS($_)) || (return (undef,undef)) } (@_))); (@vcts) or (return (undef,undef)); my ($sub) = sub {
		my ($p) = ((shift) + ((shift) && PIE/2)); [ $p, ( &SUM( map {( $$_[2]*( abs ( cos ($$_[1]-$p))))} (@vcts))) ] };
	map { (wantarray) ? (@$_) : (return $$_[0]) } grep { ($$_[0] = (1 - $$_[0])) if (defined $$_[0]); 1 } map {[ map {( &RATIO($$_[1],$mhtv))}
		( $_, ((wantarray) ? ( $sub->($$_[0],1)) : ())) ]} ( sort { our ($a,$b); ($$b[1] <=> $$a[1]) } map { my ($t) = $_; map { my ($min,$try,$max) =
			( undef, ( map {( $sub->($_))} @$t[$_,($_+1)] )); do { ${ (\$min,\$max)[0+((defined $min) && ($$max[1] < $$min[1]))] } = ($try);
			($try) = map {( $sub->($_))} ($$min[0] + $$max[0])/2; } while (($$max[0] - $$min[0])/2 > EPS); $try } (0..(@$t-2)) }
			grep { push @$_, $$_[0]+PIE; 1 } [ ( sort { our ($a,$b); ($a <=> $b) } map {( &PRINCIPAL_RAD(($$_[1]+PIE/2),1))} (@vcts)) ])[0] }

# Returns the transverse spherocity shape statistic for a list of [4-vector] momenta components or lhco objects; Independent mht or undef is leading parameter
sub SPHEROCITY_SHAPE { my ($mhtv,@vcts) = (( map { (defined) ? ( &MAX(0,$_)) : ( grep { (defined) or (return undef) }
	( &MHT( !1, @_))) } (shift)), ( grep {($$_[2] > 0)} map { ( &ETA_PHI_PTM_MAS($_)) || (return undef) } (@_)));
	( return ((defined) ? ((PIE/2)*$_)**2 : (undef))) for map {( &RATIO($_,$mhtv))} (( sort { our ($a,$b); ($a <=> $b) } map {
		my ($p) = $$_[1]; ( &SUM( map {( $$_[2]*( abs ( sin ($$_[1]-$p))))} (@vcts))) } (@vcts)), undef )[0] }

# Returns the transverse sphericity shape statistic for a list of [4-vector] momenta components or lhco objects
sub SPHERICITY_SHAPE { ( return ((defined) ? ( 1 - ( sqrt ( &MAX(0,(1-$_))))) : (undef))) for
	map { ( &RATIO( 4*( &Local::MATRIX::DETERMINANT($_)), ( &Local::MATRIX::TRACE($_))**2 )) }
	map { ( scalar ( &Local::MATRIX::INNER_PRODUCT(( scalar &Local::MATRIX::TRANSPOSE($_)), $_ ))) or (return undef) } [ map {
		($$_[0] > 0) ? [ @$_[1,2]] : () } map { ( &LORENTZ($_,1,1,!1)) or (return undef) } (@_) ] }
# TEST this calculation ...

# Returns the transverse event shape F-statistic for a list of [4-vector] momenta components or lhco objects
sub F_MATRIX_SHAPE { ( return &RATIO( sort { our ($a,$b); ($a <=> $b) } ( &QUAD_REAL_ROOTS([ ( &Local::MATRIX::DETERMINANT($_)), -1*( &Local::MATRIX::TRACE($_)), +1 ])))) for
	map { ( scalar &Local::MATRIX::INNER_PRODUCT(( scalar &Local::MATRIX::TRANSPOSE($_)), $_ )) or (return undef) } [ map { my ($s) = ( sqrt ( &MAX(0,$$_[0])));
		($s > 0) ? [ map {( $_ / $s )} (@$_[1,2]) ] : () } map { ( &LORENTZ($_,1,1,!1)) or (return undef) } (@_) ] }

# Master diobject reconstruction engine for optimization against invariant mass window specification serving lepton and jet subroutines
sub DIOBJECTS { my ($win,$ord,$prs,$obj) = do { my ($t) = ( int shift ); my (@t) = map { (@$_) ? (@$_) : [0,undef] } [
	map {[ 0+(shift @$_), ( map { (defined) ? 0+($_) : (undef) } (shift @$_)), ( grep {( ref eq q(CODE))} (shift @$_)) ]}
	grep { ( ref eq q(ARRAY)) or (return) } map { ( ref eq q(ARRAY)) ? (@$_) : () } (shift) ]; ((\@t),
		(($t > 0) ? (undef,$t) : ($t) ? (undef,undef) : ([ ( &ORDERINGS([0..(@t-1)])) ],0+@t))) };
	my (@vcts) = grep { ( &LORENTZ($_)) || (return) } (@_); ( map { (shift @$_); (@$_) } grep {(defined)}
		( &CMP( sub ($$) { my ($a,$b) = @_; (( @$b <=> @$a ) or ( $$a[0] <=> $$b[0] )) }, ( grep {(defined $$_[0])}
		map { my ($t) = $_; [ ( &MIN( grep {(defined)} map {( &NORM(@$_))} map { my (@t) = @$_;
			[ map {( $$t[$_][$t[$_]] )} (0..(@$t-1)) ]} (@{($ord)||[[(0)x(@$t)]]}))), ( map {[ @$_[-3..-1]]} (@$t)) ]}
		map {[ grep { ($prs) or (defined $$_[0]) } map { $$obj[$$_[0]][$$_[1]] ||= do {
			my ($obj,$mas) = ( &LORENTZ_SUM(undef,( my (@obj) = @vcts[@$_] )));
			my (@t) = map { my ($t) = ( abs ($mas-$$_[0])); &{ ($$_[2]) || sub { 1 }}(@obj) ? (defined $$_[1]) ?
			( map { ((defined) && ($_ <= 1)) ? ($_) : (undef) } ( &RATIO(($t,( abs $$_[1] ))[($$_[1]>=0)?(0,1):(1,0)]))) :
			( sqrt (1+$t**2)) : (undef) } (@$win); [ (($ord) ? (@t) : ( &MIN( grep {(defined)} (@t)))), @obj, $obj ] }} (@$_) ]}
		map {( &SETS(2,$_))} ( &TUPLES(2*(($prs)||( int (@vcts/2))),[0..(@vcts-1)])) )))) }

# Returns a list of pseudo-jets reconstructed from a list of [4-vector] momenta components or lhco objects by specified mode
sub HEMISPHERES { my ($hsph,$mode,$pars) = (!(shift),( map {((shift @$_), $_ )} map {[ ( ref eq q(ARRAY)) ? (@$_) : ($_) ]} (shift)));
	map { (wantarray) ? (@$_) : (return $_) } ((defined $mode) ? do { ($mode) = ( uc ((ref $mode eq q(HASH)) ? (keys %$mode)[0] : ($mode))); ($hsph) } ? do {
		map { (defined) ? [ sort ${ \( &SORT_OBJECT_LORENTZ_CODE(-1)) } grep {( $$_{ep0} > 0 )} ( &LORENTZ_HASH(undef,@$_)) ] : (return) }
		map { my ($l,$m,$c) = @$_; ( scalar (($c) -> ( $pars, [ map { (defined $l) ?
			( grep {( $$_[0] > 0 )} (( &LORENTZ($_,@$l)) or (return))) : ($_) } ( &LORENTZ_MERGE($m,@_)) ] ))) } ${{
	KTJ => [ undef, [TUP,1,+1], sub { # A list of pseudo-jets reconstructed by the KT-Jet family ( +1 => KT , 0 => CA , -1 => ANTI-KT ) of clustering algorithms
		my ($rad,$pow) = map { ((( &MAX(0,0+$$_[0])) || 1 ), 0+(0,+1,-1)[$$_[1]] ) } (shift); my (@vcts,@idx,@del,@jet) =
		( &LORENTZ_HASH(undef,@{(shift)})); while ( ((@idx) = ( grep {( defined $vcts[$_] )} (0..(@vcts-1)))) && ( my ($i,$j,$min) = @{ ( &CMP( sub ($$) {
			my ($a,$b) = @_; ($$a[2] <=> $$b[2]) }, ( map { my ($i,$j) = @$_; (${ $del[$i][$j] ||= [ (( &DELTA_RPA(@vcts[$i,$j]))/($rad))**2 ] }[0] < 1) ?
			( map {[$i,$j,$_]} @{ $del[$i][$j][1] ||= [ map { ($del[$i][$j][0]*($$_{ptm})**(2*$pow)) } (@vcts[$i,$j]) ]} ) : () }
			( &TUPLES(2,\@idx))))) || [] } )) { (($jet[0+@jet],$vcts[$_]) = ($vcts[$_],(undef))) for ( grep {(
				${ $del[$_][$_] ||= \(($vcts[$_]{ptm})**(2*$pow)) } <= $min )} (@idx));
			(@vcts[$i,$j]) = (( scalar &LORENTZ_HASH(undef,( scalar &LORENTZ_SUM(undef,@vcts[$i,$j])))),(undef));
			( undef $del[$$_[0]][$$_[1]] ) for ( map {[ sort { our ($a,$b); ($a <=> $b) } ($i,$_) ]} (@idx)); } [ map {( scalar &LORENTZ($_))} (@jet,@vcts[@idx]) ] } ],
	LND => [[!1,1,!1], [TUP,2,+1], sub { # A pair of massless pseudo-jets reconstructed by the Lund hemisphere algorithm
		(shift); my (@vcts) = map {[0,$_]} @{(shift)}; my ($axis,%axis) = map {(shift @$_)} (( &CMP( sub ($$) { my ($a,$b) = @_;
			($$b[1] <=> $$a[1]) }, ( map {[ $_, &INVARIANT_MASS(@$_) ]} ( &TUPLES( 2, [ map { $$_[1] } (@vcts) ] ))))) or (return));
		while (1) { for my $vctr (@vcts) { $$vctr[0] = ${( &CMP( sub ($$) { my ($a,$b) = @_; ($$a[1] <=> $$b[1]) },
			( map {[$_,(($$axis[$_][0])*(( &INVARIANT_MASS($$vctr[1],$$axis[$_]))/($$vctr[1][0]+$$axis[$_][0]))**2 )]} (0,1))))}[0]; }
			(@$axis) = grep { (defined) or (last) } map { my ($i) = $_; ( scalar &LORENTZ(( scalar &LORENTZ_SUM(undef,( map { $$_[1] }
			grep {($$_[0] == $i)} @vcts ))), (!1,1,!1))) } (0,1); (last) if ( $axis{ join q(), ( map { $$_[0] } (@vcts)) }++ ) }; ($axis) } ],
	MIM => [[!1,1,!1], [PRT,2,+1], sub { # A pair of massless pseudo-jets reconstructed by the minimal invariant mass-square sum
		(shift); my (@vcts) = @{(shift)}; [ map { ( scalar &LORENTZ(( scalar &LORENTZ_SUM(undef,@$_)), (!1,1,!1))) } @{(( &CMP( sub ($$) { my ($a,$b) = @_;
			(($$a[3] <=> $$b[3]) || ($$a[2] <=> $$b[2])) }, ( map {[ @$_, (( &INVARIANT_MASS(@{ $$_[0]}))**2 + ( &INVARIANT_MASS(@{ $$_[1]}))**2 ) ]}
			grep {( $$_[2] < (0+@vcts))} map {[ @$_, ( abs ( @{ $$_[1]} - @{ $$_[0]} )) ]} ( &PARTITIONS(2,\@vcts))))) or (return))}[0,1]] } ],
	MDH => [[1,!1,!1], [PRT,0,-1], sub { # A pair of transverse pseudo-jets reconstructed by the minimal scalar energy difference
		(shift); my (@vcts) = @{(shift)}; [ map {( scalar &LORENTZ_SUM(undef,@$_))} @{(
			( &CMP( sub ($$) { my ($a,$b) = @_; (($$a[3] <=> $$b[3]) || ($$a[2] <=> $$b[2])) },
			( map {[ @$_, ( abs ( &MHT( 1, @{ $$_[1]}) - &MHT( 1, @{ $$_[0]}))) ]} grep {( $$_[2] < (0+@vcts))}
			map {[ @$_, ( abs ( @{ $$_[1]} - @{ $$_[0]} )) ]} ( &PARTITIONS(2,\@vcts))))) or (return))}[0,1]] } ],
	WIN => [ undef, PRT, sub { # A list of pseudo-jets reconstructed by optimization against sets of pair-wise invariant mass window specifications
		[ map {($$_[-1])} ( &DIOBJECTS((( &GROUPS([2,-1],(shift))),(undef))[1,0], @{(shift)} )) ] } ],
	DIL => [ undef, PRT, sub { # A list of dileptons reconstructed by optimization against quartets of sign, flavor, and invariant mass window specifications
		[ map {($$_[-1])} ( &DIOBJECTS((( &GROUPS([4,-1],(shift),( sub { my ($t) = (shift);
			[ (@$t[2,3]), ( &SELECT_DIL_CODE(@$t[0,1])) ] } ))),(undef))[1,0], @{(shift)} )) ] } ],
	SUM => [ undef, undef, sub { # A pseudo-jet reconstructed by an optionally transverse, massless, inverse, or flush Lorentz sum
		my ($tmif) = (shift); my (@vcts) = @{(shift)}; [ ( &LORENTZ_SUM($tmif,@vcts)) or (return) ] } ],
	TMI => [ ($pars), undef, sub { # A list of pseudo-jets modified for transverse, massless, inverse, or flush kinematics
		(shift,shift)[1] } ],
	}}{$mode} || [ undef, undef, sub {} ] } : do { my (%sort) = do { my ($i); map {( $_ => $i++ )} (@_) };
		map { (defined) ? [ sort { our ($a,$b); ($sort{$a} <=> $sort{$b}) } (@$_) ] : (return) } ( scalar (( ${{
	WIN => sub { # A list of objects selected by optimization against sets of pair-wise invariant mass window specifications
		[ map {(@$_[0,1])} ( &DIOBJECTS((( &GROUPS([2,-1],(shift))),(undef))[1,0], @{(shift)} )) ] },
	DIL => sub { # A list of leptons selected by optimization against quartets of sign, flavor, and invariant mass window specifications
		[ map {(@$_[0,1])} ( &DIOBJECTS((( &GROUPS([4,-1],(shift),( sub { my ($t) = (shift);
			[ (@$t[2,3]), ( &SELECT_DIL_CODE(@$t[0,1])) ] } ))),(undef))[1,0], @{(shift)} )) ] },
	VBF => sub { # The inner pair of opposite hemisphere objects, or optionally the highest mass pair; Minimal pseudorapidity gap is leading parameter
		my ($deta,$minv) = map {($_,0+(0,+1)[$$_[2]])} (shift); my (@vcts) = grep { 0+(@$_) or (return) } map {( [ grep {($$_[1] < 0)} (@$_) ],
			[ grep {($$_[1] >= 0)} (@$_) ] )} [ sort { our ($a,$b); ($$a[1] <=> $$b[1]) } map {[ $_, ${ ( &ETA_PHI_PTM_MAS($_)) || (return) }[0]]} @{(shift)} ];
		map {[$$_[0][0],$$_[1][0]]} ( grep {( &MATCH_VALUE($deta,($$_[1][1]-$$_[0][1])))} (($minv) ? ( grep { ( pop @$_ ); 1 } sort { our ($a,$b); ($$b[2] <=> $$a[2]) }
			map { my ($t) = $_; map {[$t,$_,( &INVARIANT_MASS($$t[0],$$_[0]))]} (@{$vcts[1]}) } (@{$vcts[0]})) : [$vcts[0][-1],$vcts[1][0]]))[0] },
	LED => sub { # A list with specified number of objects selected by ranking on a kinematic key, ascending or descending
		my ($len,$end,$sub,@key) = map { (( map { ((abs), ($_ <=> 0)) } (( int ( shift @$_ )) || 1 )),
			( map {(( shift @$_ ), ( map {( lc (( ref eq q(HASH)) ? (keys %$_)[0] : qq($_)))} (@$_)))} map {[
				( ref eq q(ARRAY)) ? (@$_) : (defined) ? ( sub {(shift)} , $_ ) : () ]} ( shift @$_ ))) } (shift);
		[ splice ( @{(($sub) ? [ map {( shift @$_ )} sort { our ($a,$b); $end*( $$b[1] <=> $$a[1] ) }
			map {[ $_, (($sub)->( @{ ( &LORENTZ_HASH(undef,$_)) || +{}}{ @key } )) ]} @{(shift)} ] :
			($end < 0) ? [ reverse @{(shift)} ] : (shift))}, 0, $len ) ] },
	}}{$mode} || sub {} ) -> ( $pars, [ &LORENTZ_CLIP(PRT,@_) ] ))) } : [ @_ ] ) }

# Returns a list of array reference combinatoric roles for AMT2 analysis partioned from a list of [4-vector] momenta components or lhco objects
sub AMT2_ROLES { my ($mode,$pars,$lepx,$leps,$jetx,$jets) =
	(( map {((( defined $$_[0] ) ? ( uc ((ref $$_[0] eq q(HASH)) ? (keys %{$$_[0]})[0] : ($$_[0]))) : (undef)), [ @$_[1..(@$_-1)]] )} ((shift)||[])), (@_));
		( $mode eq q(GEN)) ? do { # General 2-step asymmetric user-defined MT2
	my ($sym,$mas) = map {((( &EQUAL(@$_[0,1])) && ( &EQUAL(@$_[2,3]))), $_ )} [ map {((defined) ? 0+($_) : (undef))} (@$pars[2,3,6,7]) ];
	map {[ (@$_,@$mas)[0,1,4,5,2,3,6,7]]} map { ($$pars[8]) ? ( shift @{( &CMP( sub ($$) { my ($a,$b) = @_; (($$a[1] <=> $$b[1]) or ($$b[2] <=> $$a[2])) },
		( map { my ($t) = $_; [ $t, ( &SUM( map {( &INVARIANT_MASS( grep {(defined)} @$t[@$_] ))**2 } ([0,2],[1,3]))), ( &MHT(@$t[0,1])) ]}
		(@$_))))} ) : (@$_) } do { my (%t); [ grep { !(($sym) && ( $t{ join q(_), @$_[0,1,2,3] }++ + $t{ join q(_), @$_[1,0,3,2] }++ )) }
		( &ASSIGNMENTS( [ map { my ($k,$v) = map { ( ref eq q(HASH)) ? (%$_) : () } ($_); [ &INDEXED_OBJECTS(
			[$v], ${{ LEP => $leps, JET => $jets }}{( uc $k )} ) ] } (@$pars[0,1,4,5]) ],1,!1)) ] }} : do {
	my ($leps,$hfts,$nfts) = map {[ grep {($$_[0] > 0)} map { ( &LORENTZ($_)) or (return) } (@$_) ]} ( [ &INDEXED_OBJECTS( $lepx, $leps )], (
		map { my ($t) = $_; map { [ &EXCLUDE_OBJECTS($_,(@$t)) ], ($_) } [ &SELECT_HFT([0,0],(@$t)) ] } [ &INDEXED_OBJECTS( $jetx, $jets )] ));
		my ($jets,$hdjs) = ( [ ( &LORENTZ_OBJECT_SORT(-1,@$hfts,@$nfts)) ], [ map {( &ORDERINGS($_))}
			(( &TUPLES(2,$hfts)), ( grep { push @$_, (@$hfts); 1 } ( &TUPLES((2-@$hfts),$nfts)))) ] );
		(((defined $mode) && ( $mode ne q(STM))) ? ( ${{
			MBL =>	sub {	(((@$leps) >= 1) && ((@$jets) >= 2)) or (return); # MT2_BL
				map {[ ( scalar &LORENTZ_SUM(undef,$$leps[0],$$_[0])), $$_[1], 0, 80.4 ]} (@$hdjs) },
			MTW =>	sub {	(((@$leps) >= 1) && ((@$jets) >= 2)) or (return); # MT2_W
				map {[ $$_[0], $$_[1], 0, 80.4, $$leps[0], undef, 80.4, undef ]} (@$hdjs) },
			TAU =>	sub {	(((@$leps) >= 1) && ((@$jets) >= 3)) or (return); # MT2_TAU
				map {[ $$_[0], $$_[1], 0, 0, $$leps[0], ( &EXCLUDE_OBJECTS([$$_[0],$$_[1]],@$jets))[0], 80.4, 80.4 ]} (@$hdjs) },
		}}{$mode} or sub {} ) : sub {[ (( &HEMISPHERES( 0, q(LND), @$leps,@$hfts,@$nfts)), undef, undef )[0,1], 0, 0 ]} )->() }}

# Returns a Python-styled raw string, with (LaTeX style \\) newlines and double quotes (") interpolated
sub RAW_STRING { 'r"'.(( map {( join '" "\n" r"', map {( join '" "\"" r"', ( split /"/, $_, -1 ))} ( split /\\\\/, $_, -1 ))} ( qq/${\(shift)}/ ))[0] ).'"' }

# Returns the transformation of a user specified string with escape codes replaced by their corresponding values
{; my ($map); sub UNESCAPE_STRING { my ($str) = map { ((defined) and !(ref)) ? qq($_) : (return undef) } (shift); $map ||= +{
	HSH => q(#),
	DQT => q("),
	LES => q(<),
	GRT => q(>),
	MET => q(${/\!\!\!\!E}_{\rm T}$),
	MHT => q($H_{\rm T}$),
	MEF => q($M_{\rm eff}$),
	RET => q($R\({/\!\!E}_{\rm T}\)$),
	RHR => q(${/\!\!E}_{\rm T}/\sqrt{H_{\rm T}}$),
	REF => q(${/\!\!E}_{\rm T}/M_{\rm eff}$),
	RHH => q($R\(H_{\rm T}\)$),
	DET => q($\Delta\({/\!\!E}_{\rm T}\)$),
	PTM => q(${P}_{\rm T}$),
	SRS => q($S/\sqrt{S+B}$),
	SR1 => q($S/\sqrt{1+B}$),
	SRB => q($S/\sqrt{B}$),
	SO1 => q($S/1+B$),
	SOB => q($S/B$),
	RTS => q($\sqrt{s}$),
	LUM => q($\mathcal{L}$),
	IPB => q(${\rm pb}^{-1}$),
	IFB => q(${\rm fb}^{-1}$),
	IAB => q(${\rm ab}^{-1}$),
	IZB => q(${\rm zb}^{-1}$),
	IYB => q(${\rm yb}^{-1}$),
	DEF => q($d\sigma/d E \div \sigma$),
	DNF => q($\Delta\sigma/\Delta N \div \sigma$),
	}; $str =~ s/<((?i:[A-Z][A-Z\d]{2}))>/$$map{uc $1}/g; ($str) }}

# Returns the initial (or sequential) match of an input string against a list of regular expressions, optionally with interpolation and replacement
sub REX_LIST { my ($mod) = 0+(0..2)[(shift)]; my ($ref,$str) = map {(\($_),( qq($_)))} (shift); do { my ($rex,$val,$rpl) = @$_;
	( return ((wantarray) ? (@$_) : \((ref $val eq q(CODE)) ? ( scalar $val->( map {( qq($_))} (@$_))) : ($val)))) for
	grep { local ($@); eval { $$ref = ( join '', (($$_[-1]),((ref $rpl eq q(CODE)) ? ( scalar $rpl->( map {( qq($_))} (@$_))) : ($rpl)),($$_[-2]))) }
		if ($mod == 2); !($@) or (return) } ((($mod == 1) ? ( $$ref =~ m/\G$rex/gc ) : ( $str =~ m/$rex/ )) ?
		[ ( map {( substr $str, $-[$_], ($+[$_] - $-[$_]))} (0..(@+-1))), ( substr $str, $+[0] ), ( substr $str, 0, $-[0] ) ] :
		()) } for ( map { ( ref eq q(ARRAY)) ? ($_) : [$_] } (@_)); () }

# Returns a closure encapsulating evaluation of an input function at an indexed point in the domain
sub INDEXED_FUNCTIONAL { my ($sub,@vls) = (( map { ( ref eq q(CODE)) ? ($_) : (defined) ? (return undef) : sub {(shift)}} (shift)),(@_));
	sub { ( scalar (($sub) -> ( @{ ( grep { ( ref eq q(ARRAY)) or (return undef) } (shift))[0] }[ @vls ] ))) }}

# Returns a closure encapsulating evaluation of an input function at a hashed point in the domain
sub HASHED_FUNCTIONAL { my ($idx,$sub) = (( grep { ( ref eq q(HASH)) or (return undef) } (shift)), ( grep { ( ref eq q(CODE)) or !(defined) or (return undef) } (shift)));
	( &INDEXED_FUNCTIONAL( $sub, ( grep { (defined) or (return undef) } map { ( ref eq q(HASH)) ? ( $$idx{( sprintf q(%3.3s_%3.3i), %$_ )} ) : (undef) } (@_)))) }

# Returns the compound evaluation or closure encapsulation of an input operation string
{; my ($rex); sub STRING_FUNCTIONAL { my ($str,$vls,@map) = (( grep { s/\s+//g; ( m'(?i:[^-+/*^\$\d.)(,A-Z><=!&|])' ) && (return undef); 1 }
	map {( qq($_))} (shift)), ( map { ( ref eq q(ARRAY)) ? [(undef),(@$_)] : () } (shift)));
		$rex ||= do { my ($val) = qr"(?:(?>${\EXP})|(?>[-+]?[\$@]\d+))"; [
			[ qr"^(${val})$", sub {[1,@_[1..1]]} ],
			[ qr"(?<!(?i:[A-Z]))\((${val})\)", sub {[2,@_[1..1]]} ],
			[ qr"(?<=[-+])(?=[-+])(${val})(?:(?=[-+/*,)><=!&|])|$)", sub {[2,@_[1..1]]} ],
			[ qr"((?i:[A-Z]+))\((${val}(?:,${val})*)?\)", sub {[3,(undef),@_[1..2]]} ],
			[ qr"(?:^|(?<=[-+/*^(,><=&|]))(?![-+])(${val})(\^)(${val})(?:(?=[-+/*,)><=!&|])|$)", sub {[4,@_[1..3]]} ],
			[ qr"(?:^|(?<=[-+(,><=&|]))(?![-+])(${val})([/*])(${val})(?:(?=[-+/*,)><=!&|])|$)", sub {[4,@_[1..3]]} ],
			[ qr"(?:^|(?<=[(,><=&|]))(${val})([-+])(${val})(?:(?=[-+,)><=!&|])|$)", sub {[4,@_[1..3]]} ],
			[ qr"(?:^|(?<=[(,&|]))(${val})(<|<=|<=>|==|!=|>=|>)(${val})(?:(?=[,)&|])|$)", sub {[4,@_[1..3]]} ],
			[ qr"(?:^|(?<=[(,|]))(${val})(&&)(${val})(?:(?=[,)&|])|$)", sub {[4,@_[1..3]]} ],
			[ qr"(?:^|(?<=[(,]))(${val})(\|\|)(${val})(?:(?=[,)|])|$)", sub {[4,@_[1..3]]} ]] }; { my ($mod,$opn,$arg) =
	do { my ($mod,@c) = @{ ${ ( &REX_LIST( 2, $str, ( map {[ @$_[0,1], ( q(@).(0+ @map)) ]} (@$rex)))) || \(!1) } || [] };
		( 0+($mod), $c[1], [ map { my ($s,$o,$n) = /^([-+]?)([\$@]?)(.*)$/; map {(( $s eq q(-)) ? do { my ($sub) = $_;
			sub { -1*( grep { (defined) or (return undef) } ( $sub->()))[0] }} : ($_))} (( length $o ) ? ( $o eq q($)) ?
			sub { ( grep { (defined) or (return undef) } ($$vls[(int $n)]))[0] } : ( splice @map, (int $n), 1, (undef)) :
			sub { (( defined $n ) ? (0+ $n) : (undef)) } ) } map {( split q(,))} (@c[0,2]) ] ) };
	push @map, map { my ($sub) = $_; sub { ( grep { (defined) && !(/inf/) && !(/nan/) or (return undef) } ( $sub->( map {( $_->())} (@$arg))))[0] }}
		($mod == 1) ? do { my ($sub) = (shift @$arg); ( return (($vls) ? ( $sub->()) : ( sub { ($vls) = [(undef),(@_)]; ( $sub->()) } ))) } :
		($mod == 2) ? (shift @$arg) :
		($mod == 3) ? do { my ($sub,$map) = ((((@$arg) ? +{
			q(sin)	=> sub { ( sin (shift)) },
			q(cos)	=> sub { ( cos (shift)) },
			q(tan)	=> sub { ( map { ( scalar ( eval { (( sin ($_)) / ( cos ($_))) } )) } (shift))[0] },
			q(asn)	=> sub { ( map { ( scalar ( eval { atan2 ( $_, sqrt( 1 - $_*$_ )) } )) } (shift))[0] },
			q(acs)	=> sub { ( map { ( scalar ( eval { atan2 ( sqrt( 1 - $_*$_), $_ ) } )) } (shift))[0] },
			q(atn)	=> sub { ( atan2 ((shift), ((@_) ? (shift) : (1)) )) },
			q(exp)	=> sub { ( map { (@_) ? ((shift) ** ($_)) : ( exp ($_)) } (shift))[0] },
			q(log)	=> sub { ( scalar ( eval { ( log (shift)) / ((@_) ? ( log (shift)) : (1)) } )) },
			q(snh)	=> sub { ( map { ( exp (+$_) - exp (-$_)) / 2 } (shift))[0] },
			q(csh)	=> sub { ( map { ( exp (+$_) + exp (-$_)) / 2 } (shift))[0] },
			q(tnh)	=> sub { ( map { (( 1 - exp (-2*$_)) / ( 1 + exp (-2*$_))) } (shift))[0] },
			q(ash)	=> sub { ( map { ( log ( $_ + sqrt( $_*$_ + 1 ))) } (shift))[0] },
			q(ach)	=> sub { ( map { ( scalar ( eval { ( log ( $_ + sqrt( $_*$_ - 1 ))) } )) } (shift))[0] },
			q(ath)	=> sub { ( map { ( scalar ( eval { ( log ( 1 + $_ ) - log ( 1 - $_ )) / 2 } )) } (shift))[0] },
			q(srt)	=> sub { ( scalar ( eval { ( sqrt (shift)) } )) },
			q(abs)	=> sub { ( abs (shift)) },
			q(int)	=> sub { ( &INT_CEILING_FLOOR ((shift), 0 )) },
			q(clg)	=> sub { ( &INT_CEILING_FLOOR ((shift), +1 )) },
			q(flr)	=> sub { ( &INT_CEILING_FLOOR ((shift), -1 )) },
			q(rnd)	=> sub { ( &ROUND (@_)) },
			q(qnt)	=> sub { ( &INT_QUOTIENT ((shift), ((@_) ? (shift) : (1)), -1 ))[0] },
			q(mod)	=> sub { ( &INT_QUOTIENT ((shift), ((@_) ? (shift) : (1)), -1 ))[1] },
			q(ife)	=> sub { ( $_[ 0+( !(shift)) ] ) },
			q(def)	=> sub { 0+( defined (shift)) },
			q(not)	=> sub { 0+( !(shift)) },
			q(les)	=> sub { 0+((shift) < (shift)) },
			q(leq)	=> sub { 0+((shift) <= (shift)) },
			q(cmp)	=> sub { 0+((shift) <=> (shift)) },
			q(eql)	=> sub { 0+((shift) == (shift)) },
			q(neq)	=> sub { 0+((shift) != (shift)) },
			q(geq)	=> sub { 0+((shift) >= (shift)) },
			q(grt)	=> sub { 0+((shift) > (shift)) },
			q(and)	=> sub { 0+((shift) and (shift)) },
			q(orr)	=> sub { 0+((shift) or (shift)) },
			q(xor)	=> sub { 0+((shift) xor (shift)) },
				} : +{
			q(pi)	=> sub {(PIE)},
			q(pie)	=> sub {(PIE)},
			q(udf)	=> sub {(undef)},
			q(undef)=> sub {(undef)},
				} ) -> {( lc $opn )} ) or (( +{
			q(sum)	=> sub { ( &SUM (@_)) },
			q(prd)	=> sub { ( &PRODUCT (@_)) },
			q(nrm)	=> sub { ( &NORM(@_)) },
			q(min)	=> sub { ( &MIN (@_)) },
			q(max)	=> sub { ( &MAX (@_)) },
			q(avg)	=> sub { ( &ARITHMETIC (@_)) },
			q(hrm)	=> sub { ( &HARMONIC (@_)) },
			q(geo)	=> sub { ( &GEOMETRIC (@_)) },
				} -> {( lc $opn )} or (return undef)), 1 )); sub { ( scalar (
			(($map) && (@_ == 1) && ( &ISA( 1, $_[0], q(Local::TENSOR)))) ? (($sub) -> ( map {($$_)} ((shift) -> ELEMENTS()))) :
			((( map { ( &UNIVERSAL::can( $_, q(MAP))) || () } (@_)),(undef))[0] or sub { my ($sub) = (pop); (($sub)->(@_)) } ) -> (@_,$sub))) }} :
		($mod == 4) ? (( +{
			q(^)	=> sub { ( scalar ((shift) ** (shift))) },
			q(/)	=> sub { ( scalar ( eval { ((shift) / (shift)) } )) },
			q(*)	=> sub { ( scalar ((shift) * (shift))) },
			q(-)	=> sub { ( scalar ((shift) - (shift))) },
			q(+)	=> sub { ( scalar ((shift) + (shift))) },
			q(<)	=> sub { 0+((shift) < (shift)) },
			q(<=)	=> sub { 0+((shift) <= (shift)) },
			q(<=>)	=> sub { 0+((shift) <=> (shift)) },
			q(==)	=> sub { 0+((shift) == (shift)) },
			q(!=)	=> sub { 0+((shift) != (shift)) },
			q(>=)	=> sub { 0+((shift) >= (shift)) },
			q(>)	=> sub { 0+((shift) > (shift)) },
			q(&&)	=> sub { 0+((shift) and (shift)) },
			q(||)	=> sub { 0+((shift) or (shift)) },
				} ) -> {( lc $opn )} || (return undef)) :
		(return undef); (redo) }}}

# Returns an anonymous scalar reference, optionally initialized by a user-supplied input
sub SCALAR_REF { my ($s) = (shift); \( $s ) }

# Sleeps for a specified number of seconds (default 1, fractional to milliseconds), optionally up to some maximal value (randomized)
sub SLEEP { my ($min,$max) = ((( &MAX( 0, 0+ (shift))) or (1)), (0+ (shift)));
	select ((undef), (undef), (undef), ( $min + (($max > $min) ? (( $max - $min ) * ( rand )) : (0)))); 1 }

#**********#
# PACKAGES #
#**********#

# Encapsulates path, file, and handle object definition and manipulation
{; package Local::FILE; use overload
	q(*{})		=> sub { ${(shift)}},
	q(fallback)	=> 1;

# Returns a handle corresponding to file opened and locked according to specified protocols with validation of atomicity; cf. http://www.perlmonks.org/?node_id=28996
sub HANDLE { use Fcntl qw(:DEFAULT :flock :seek); my ($wrt,$pth,$fil,$opn,$sek,$trm,$lck,$msk) = ((@_) ? ( map { my ($p,$f,$w,$c,$m) =
		((defined) ? ( ref eq q(ARRAY)) ? (@$_[0,1]) : (/^(.*\/)?([^\/]*)$/) && ((($1) or (undef)),($2)) : (return));
		map { my ($o) = $_; (($w), ( &PATH($p,(($c)?(2,$m):())) or (return)), ( map {
		(($_), ((defined) ? (/\.gz$/) ? ($w) ? (return) : () : (@$o) : (return))) } ( &NAME($f)))) } [
		( map {(((($w,undef) = @{ ([!1,O_RDONLY],[1,O_WRONLY],[1,O_RDWR])[0+(0..2)[ shift @$_ ]] } )[1] |
		(($c,undef) = @{ ([!1,0],[1,O_CREAT],[1,O_CREAT|O_EXCL])[0+(0..2)[(0,( shift @$_ ))[0+$w]]] } )[1] ),
		( @{ ([SEEK_SET,!1],[SEEK_SET,1],[SEEK_END,!1])[0+((0..2)[(0,( shift @$_ ))[0+$w]])] } ),
		((LOCK_SH,LOCK_EX)[0+$w] | (0,LOCK_NB)[0+(0..1)[ shift @$_ ]] ))} map { ( ref eq q(ARRAY)) ? [(@$_)] :
		([0,0,0,0],[1,1,1,0],[1,1,2,0],[2,0,0,0],[2,1,1,0],[2,1,2,0])[0+(0..5)[$_]] } (shift)),
		(( grep {(($_ >= 0) && ($_ <= 04777))} map { (/^0/) ? (oct) : (defined) ? 0+($_) : () }
		(($m,undef) = map { ( ref eq q(ARRAY)) ? (@$_) : (undef,$_) } (shift))[1] ), 00666 )[0]] } (shift)) : (1));
	( map { (wantarray) ? ( $_, ((defined $pth) ? [$pth,$fil] : ())) : (return $_) } map {( bless \($_))} ((( grep {
		( binmode $_, q(:perlio) ) && (!(defined $sek) or ( seek $_, 0, $sek )) && (!($trm) or ( truncate $_, 0 )) && (!($wrt) or ( select $_ )) }
		grep { !(defined $fil) ? ( open $_, q(+<), undef ) : !(defined $opn) ? ( open $_, q(-|), ( q(gzip -cdf ).$pth.$fil )) :
		do { my ($b); ( undef $_ ) while (( $b = (( sysopen $_, $pth.$fil, $opn, $msk ) && ( flock $_, $lck ))) and
		(( stat $_ )[1] ne (( stat $pth.$fil )[1] )) and ( $b = ( close $_ ))); $b }} ( my ($IOH))), (undef))[0] or (return))) }
# path: array ref [path,file] objects OR path string OR undef for temporary file
# mode: int=(0..5): (read,write,append,dual-read,dual-write,dual-append)[mode] OR
#	set explicitly as [priv,make,seek,fast] with:
#	priv : (read,write,read|write)[priv]; default 0
#	make : (no,yes,exclusive)[make]; yes and exclusive require write; default 0
#	seek : (head,trim,tail)[seek]; trim and tail require write; default 0
#	fast : (block,nonblock)[fast]; default 0
# mask: file creation umask in octal OR directory/file mask array ref

# Returns a handle to the next integrally available file name based upon an input directory and file name base
sub NEXT { my ($rpt,$pth,$bas,$idx,$ext) =
	map { ((::RPT), ( &PATH($$_[0],2) or (return)), (@{ &KEYS($$_[1]) or (return) })) }
	map {[ ( ref eq q(ARRAY)) ? (@$_[0,1]) : (/^(.*\/)?([^\/]*)$/) && ((($1) or (undef)),($2)) ]} (shift);
	while (1) { ( return ((wantarray) ? (@$_) : ( shift @$_ ))) for (
		grep {((@$_) or ((($rpt--) <= 0 ) && (return)) or (( &::SLEEP( 0.25, 0.50 )), !1 ))}
		map {[ &HANDLE([$pth,$_],[1,2,0,1]) ]} grep { !( -e $pth.( &NAME($_))) } [$bas,++$idx,$ext] ) }}
# An input array or string is split into the path, which must be writable, and the file object keys
# The file object index is incremented until an open file name is located
# The file is created under an exclusive request and a locked writeable handle is returned
# The file creation is retried on a subsequent index after failure a limited number of times

# Returns a list of files located at some level below the specified path that match a specified case-(in)senstive pattern
sub LIST { my ($pth,$fil,$lvl) = (( map { (( &PATH($$_[0]) or (return)), ( &NAME($$_[1],1) or (return))) }
	map {[ ( ref eq q(ARRAY)) ? (@$_[0,1]) : (/^(.*\/)?([^\/]*)$/) && ((($1) or (undef)),($2)) ]} (shift)), ( int (shift)));
	map { (( map {[$pth,$_]} grep {(( -f ($pth.$_)) && ( $_ =~ $fil ) && ( $_ eq ( &NAME($_))))} (($lvl < 1) ? (@$_) : ())),
		( map {( &LIST([[$pth,$_],$fil],($lvl-1)))} grep { ( -d ($pth.$_)) && !(m/^tmp$/i) } (($lvl) ? (@$_) : ()))) }
	map {[ sort { our ($a,$b); ($a cmp $b) } grep {!(/^\./)} ( readdir $_ ) ]} grep {(( opendir $_, $pth ) or (return))} ( my ($DHI)) }
# An input array or string is split into the path and the file match Regexp
# Hidden files and 'tmp' directories are omitted
# Files must match the pattern and adhere to canonical naming conventions
# Results are alpha sorted as path and file doublets
# A negative level returns files at all directory depths

# Returns the set of file object keys [ key, index, extension ], corresponding to an input file string
sub KEYS { my ($r); (( map { (($r) ? (defined) : (((length) <= 255 ) && !(/\*/))) ? ( grep { ( @$_ == 3 ) && do {
	$$_[1] =~ s/_//g; $$_[1] += 0; ($r) && !(length $$_[2]) && do { ( undef $$_[2] ); $$_[1] ||= (undef); }; 1 }}
		[ /^([*A-Za-z][*\w-]*?)((?:_\d{3})*)((?:\.[A-Za-z0-9]+)*)$/ ] ) : (($r) ? [ ((undef)x(3)) ] : ()) }
	map { $r = !!(shift); (( ref eq q(ARRAY)) ? (( !($r) && ( &NAME($_))) or ()) : ($_)) } (shift)), (undef))[0] }
# Input arrays, if the Regexp flag in the second input is false, are run first through NAME, or failed out
# If the Regexp flag is set but the input is UNDEF then all three keys are set to UNDEF
# Processing continues with the Regexp flag for a defined input, or for strings within the CHAR limit of 255 that have no wildcards
# A pattern match isolates the key (leading alpha and trailing word chars including dashes with wildcards permitted anywhere),
#	index (sets of underscores followed by three digits), and extension (sets of dots followed by alphanumeric chars) or fails out
# The index is converted to a number after joining digits
# For Regexp mode with an empty extension, the extension is undefined and the index is undefined if zero

# Returns the normalized file name string or regular expression corresponding to an input set of file object keys
sub NAME { my ($r); (( map { ($r) ? ( map {( qr($_))} (('(?').(( pop @$_ ) ? ( q(i)) : ( q())).( q(:^)).( join q(), (@$_)).('$)'))) :
	( grep {((length) <= 255 )} ( join q(), (@$_))) } grep {( @$_ == 4 )} map {[
		( map { (($r) ? (defined) : !(/\*/)) ? ( grep { s/\*/.*?/g; 1 } ( /^([*A-Za-z][*\w-]*)$/ )) : (($r) ? ( q(.*?)) : ()) } ($$_[0])),
		( map { (($r) && !(defined)) ? ('(?:_\d{3})*') : (( $_ = ( int )) > 0 ) ? ( join q(), ( map {( sprintf q(_%3.3i), $_ )}
			map {(( grep {(length)} ( substr ( $_, 0, ((length)%(3)), q()))), ( m/(\d{1,3})/g ))} ( qq($_)))) : ( q()) } ($$_[1])),
		( map { (($r) && !(defined)) ? ('(?:\.[A-Za-z0-9]+)*') : ( grep { s/^\.?/\./; ($r) && ( s/\./\\\./g ); 1 }
			( /^((?:\.?[A-Za-z0-9]+)*)$/ )) } ($$_[2])), (($r) && !!($$_[3])) ]}
	map { (( $r = !!(shift)) && ( ref eq q(Regexp))) ? (return $_) : ( ref eq q(ARRAY)) ? [ @$_[0..3]] : (( &KEYS($_,$r)) or ()) } (shift)), (undef))[0] }
# If second input (RegExp flag) is true and leading input is a Regexp object, then that object is returned directly
# If input is not an array it is passed first through KEYS or failed out
# A key matching the wildcard '*' that is not in Regexp mode is failed out
# In Regexp mode an undefined key is converted to the non-greedy arbitrary match '.*?'
# Otherwise the call fails out unless the key is leading alpha and trailing word chars including
#	dashes with wildcards permitted anywhere; Wildcards are converted to arbitrary match '.*?'
# In Regexp mode an undefined index is converted to match sets of unscores followed by three digits
# Non-positive index values are set to the empty string ''
# The index is joined on underscore every three characters with leading zeroes as needed
# In Regexp mode an undefined extension is converted to match sets of dots followed by alphanumeric chars
# Otherwise the call fails out unless extension is sets of dots (leading dot optional) followed by alphanumeric chars
# A leading dot is forced; In Regexp mode dots are escaped
# A fourth entry indicating case insensitivity retains its truth value in Regexp mode, or is set false
# In Regexp mode the string is joined, and padded with suitable endcaps, and quoted into a Regexp object
# Otherwise, the string is joined and failed out if it is longer than the CHAR limit of 255

	# length here, elsewhere combined to total length ?
# Returns the normalized path string corresponding to an input path object, optionally creating path or enforcing writeability
sub PATH { my ($i); my ($pth,$mod,$msk,$str) = ( [ map {(($_).q(/))} map { (/^(?:|\.|(\.\.)|([\w-]+)|(~))$/) or (return undef);
	(($i++) ? (($3) ? (return undef) : (($2) or ($1)) ? ($_) : ()) : ((($2) ? ( q(.)) : ()), ($_))) }
		(( map {( split /\//, $_, -1 )} map { ( ref eq q(ARRAY)) ? (@$_) : ($_) } (shift)), ( q(.))) ],
	0+(0,1,2)[(shift)], (( grep {(($_ >= 0) && ($_ <= 04777))} map { (/^0/) ? (oct) : (defined) ? 0+($_) : () } (shift)), 00777 )[0] );
	while (@$pth) { $str .= ( shift @$pth ); (return undef) unless (( -d $str ) ? ((@$pth) || (($mod) ? ( -w $str ) : ( -r $str ))) :
		(($mod > 1) && !( -e $str ) && ( mkdir $str, $msk ) && ( -w $str ))) } ($str) }
# Acceptable path elements are '', '.', '..', '~', corresponding to root, current, parent, and home, and also word characters with dashes
# '~' may be leading only; '' and '.' are omitted if not leading; '.' is prefixed if word characters are leading
# An undefined or empty path canonicalizes to './'; All paths have a trailing slash

# Deselects, closes, and frees extraneous references to associated filehandle object upon exit of scope
sub DESTROY { do { (($_) eq (select)) && ( select STDOUT ); ( return ( close $_ )); } for ( grep {(defined)} ${((shift) or (return))} ); () }

}	# End of package Local::FILE

# Encapsulates HISTOGRAM object definition and manipulation
{; package Local::HISTOGRAM; BEGIN { our (@ISA) = ( q(Local::TENSOR), q(ARRAY)) }; use overload
	q(@{})		=> sub { ${(shift)}{TNSR}},
	q(bool)		=> sub { 1 },
	q(fallback)	=> 1;

# Returns a new multi-dimensional HISTOGRAM object initialized for left offset, span size, and bin count
sub NEW { my ($obj) = ( bless +{} ); @$obj{( qw( HDIM BINC LEFT SPEC TNSR ))} = ( map {((@$_), ( scalar &Local::TENSOR::OBJECT((undef), 0, @{$$_[1]} )))}
	(( &::ISA( 1, $_[0], __PACKAGE__ )) ? [ map {( &::CLONE($_))} @{$_[0]}{( qw( HDIM BINC LEFT SPEC ))} ] : [
	map { my ($t) = $_; ( 0+(@$t), ( map { my ($i) = $_; [ map {($$t[$_][$i])} (0..(@$t-1)) ] } (0..2)))} [ map { my ($i);
	my ($bnc,$lft,$rgt,$spn,$bns,@spc) = ((2), ( grep { (@$_) ? (( @$_ == 1 ) or (( $i ||= @$_ ) == @$_ ) or (return undef)) :
	( push @$_, (undef)) } map { ( ref eq q(ARRAY)) ? ($_) : [$_] } (@$_[0..3]))); (((( my $l ) = (@$lft)) == 1 ) or (return undef));
	for my $j (0..(($i)&&($i-1))) { my ($r,$s,$b) = map {( $$_[(@$_-1)&&($j)] )} ($rgt,$spn,$bns);
	($s,$b) = map {( &::MAX(0,$_))} ((0+ $s ),( int $b )); (defined $r) && ($r > $l) &&
	(($s) ? ( $b = ( int &::RATIO(($r-$l),$s))) : ( $s = (0+ &::RATIO(($r-$l),$b)))); ($s) && ($b) or (next);
	push @spc, [ $s, $b, 0+($l), ($l += $s*$b), 0+($bnc), ($bnc += $b) ] } $bnc++; [ 0+$bnc, 0+( shift @$lft ), \@spc ] }
	( grep {(( ref eq q(ARRAY)) or (return undef))} (( ref $_[0] ) ? (@_) : ( map {($_[1])} (1..( int $_[0] ))))) ]] )); ($obj) }

# Method for the population of object bins by an input list of values
sub BIN { my ($obj) = (shift); my ($dim,$bnc,$lft,$spc,$tns) = @$obj{( qw( HDIM BINC LEFT SPEC TNSR ))};
	my ($val) = ( grep {(( @$_ == $dim ) or (return))} map {(( ref eq q(ARRAY)) ? ($_) : [$_] )} (shift));
 	${ (((($tns) -> ELEMENTS( map { my ($idx,$val) = ($_,$$val[$_]); !(defined $val) ? (0) : do {
		(( $val += ((::NIL)*($$spc[$idx][0][0]))) < $$lft[$idx] ) ? (1) : do {
		my ($s,$b,$l,$r,$i,$j); for (@{$$spc[$idx]}) { ($s,$b,$l,$r,$i,$j) = @$_; (last) if ($val < $r); }
		( &::MIN( $j, ( $i + ( &::INT_CEILING_FLOOR( &::RATIO(($val-$l),$s), -1 ))))) }}}
		(0..($dim-1)))),(undef))[0] or (next)) } += ((@_) ? (0+ (shift)) : (1)); ($obj) }

# Method for the generation of an array holding the bin widths
sub WIDTHS { my ($obj) = (shift); my ($dim,$bnc,$lft,$spc,$tns) = @$obj{( qw( HDIM BINC LEFT SPEC TNSR ))};
	grep {((wantarray) or (return $_))} map {[ map { my ($s,$b) = @$_; (($s)x($b)) } (@{$$spc[$_]}) ]} (0..($dim-1)) }

# Method for the generation of an array holding the bin edge locations
sub EDGES { my ($obj) = (shift); my ($dim,$bnc,$lft,$spc,$tns) = @$obj{( qw( HDIM BINC LEFT SPEC TNSR ))};
	grep {((wantarray) or (return $_))} map {[ 0+($$lft[$_]), ( map { my ($s,$b,$l) = @$_; map {($l + $s*$_)} (1..$b) } (@{$$spc[$_]})) ]} (0..($dim-1)) }

# Method for the generation of an array holding the bin center locations
sub CENTERS { my ($obj) = (shift); grep {((wantarray) or (return $_))} map { my ($e) = $_; [ map {( &::ARITHMETIC(@$e[$_-1,$_]))} (1..(@$e-1)) ] } ( $obj->EDGES()) }

}	# End of package Local::HISTOGRAM

# Encapsulates VECTOR object definition and manipulation
{; package Local::VECTOR; BEGIN { our (@ISA) = q(Local::TENSOR) }; use overload
	q(.)		=> sub { ( &::ISA( 1, $_[1], q(Local::MATRIX))) ? ( &MATRIX_PRODUCT(shift,shift)) : ( &INNER_PRODUCT ) },
	q(^)		=> \&CROSS_PRODUCT,
	q(x)		=> \&OUTER_PRODUCT,
	q(bool)		=> sub { 1 },
	q(fallback)	=> 1;

# Validates and returns a one-dimensional VECTOR object generated from input nested ARRAY reference; includes dimension in list context
sub OBJECT { ( &Local::TENSOR::OBJECT( __PACKAGE__, (shift), ( map {((defined) ? ((int) or (return)) : (undef))} (shift)), 0 )) }

# Returns the cartesian scalar product of two input VECTOR objects
sub INNER_PRODUCT { ( my ($a,$e) = ( &OBJECT(shift))) or (return undef); my ($b) = ( &OBJECT((shift),$e) || (return undef));
	( &::SUM( map { $$a[$_]*$$b[$_] } (0..($e-1)))) }

# Returns the row vector times matrix inner product of input vector and MATRIX objects; includes dimension in list context
sub MATRIX_PRODUCT { my ($a,$b) = (shift,shift); ( return ( &Local::MATRIX::VECTOR_PRODUCT($b,$a))) if (shift);
	( &OBJECT( ${( &Local::MATRIX::INNER_PRODUCT( [ ( &OBJECT($a) or (return)) ], ($b)) or (return))}[0] )) }

# Returns the cartesian cross product of two input three-dimensional VECTOR objects
sub CROSS_PRODUCT { my ($a,$b,$i) = (( map { ( &OBJECT($_,3)) || (return) } (shift,shift)), ((shift) ? -1 : +1 ));
	( &OBJECT([ map { $i*( $$a[$_]*$$b[$_+1] - $$a[$_+1]*$$b[$_] ) } (-2..0) ])) }

# Returns the outer matrix product of two input VECTOR objects
sub OUTER_PRODUCT { my ($a,$b) = map {( &OBJECT($_) or (return))} (shift,shift);
	(( \&Local::MATRIX::OBJECT, \&Local::MATRIX::TRANSPOSE )[ 0+(0,1)[(shift)]] ) ->
	([ map { my ($t) = $_; [ map {($t*$_)} (@$b) ] } (@$a) ]) }

# Returns the cartesian norm of a list of input VECTOR objects
sub NORM { ( return ((defined) ? (sqrt) : (undef))) for ( &INNER_PRODUCT(( &Local::TENSOR::SUM(@_))[0,0])) }

# Returns an array reference with specified count of ascending input coordinate powers; includes dimension in list context
sub POWERS { my ($x,$e,$a) = (( map { (defined) or (return); (0+$_) } (shift)), (( &::MAX(0,( int shift ))) || (return)), 0+(shift));
	( &OBJECT([ 0+( $a ||= 1 ), ( map { $a *= $x } (1..($e-1))) ])) }

}	# End of package Local::VECTOR

# Encapsulates MATRIX object definition and manipulation
{; package Local::MATRIX; BEGIN { our (@ISA) = q(Local::TENSOR) }; use overload
	q(.)		=> sub { ( &::ISA( 1, $_[1], q(Local::VECTOR))) ? ( &VECTOR_PRODUCT(shift,shift)) : ( &INNER_PRODUCT ) },
	q("")		=> \&STRING,
	q(bool)		=> sub { 1 },
	q(fallback)	=> 1;

# Validates and returns a two-dimensional MATRIX object generated from input nested ARRAY reference; includes dimension in list context
sub OBJECT { ( &Local::TENSOR::OBJECT( __PACKAGE__, (shift), ( map {((defined) ? ((int) or (return)) : (undef))} (shift,shift)), 0 )) }

# Returns a unit MATRIX object of a specified size, optionally with a specified scale factor
sub UNIT { my ($r,$s) = ((( &::MAX(0,( int shift ))) or (return)), ( map {((defined) ? 0+($_) : (1))} (shift)));
	( &OBJECT( [ map { my ($t) = $_; ( grep { $$_[$t] = $s; 1 } [((0)x($r))] ) } (0..($r-1)) ] )) }

# Returns the matrix inner product of two input MATRIX objects; includes dimension in list context
sub INNER_PRODUCT { my (@m) = @_[($_[2])?(1,0):(0,1)]; (( my ($a,$r,$x) = ( &OBJECT($m[0]))) or (return));
	(( my ($b,$c,$l,(undef)) = ( &TRANSPOSE($m[1],$x,(undef)))) or (return));
	( &OBJECT( [ map { my ($i) = $_; [ map { my ($j) = $_; ( &Local::VECTOR::INNER_PRODUCT($$a[$i],$$b[$j])) } (0..($c-1)) ]} (0..($r-1)) ] )) }

# Returns the matrix times column vector inner product of input MATRIX and VECTOR objects; includes dimension in list context
sub VECTOR_PRODUCT { my ($a,$b) = (shift,shift); ( return ( &Local::VECTOR::MATRIX_PRODUCT($b,$a))) if (shift);
	( &Local::VECTOR::OBJECT( [ map {( shift @$_ )} @{( &INNER_PRODUCT( $a, [ map {[ $_ ]}
		@{( &Local::VECTOR::OBJECT($b) or (return))} ] ) or (return))} ] )) }

# Returns a transposed nested array reference copy of input MATRIX object; includes dimension in list context
sub TRANSPOSE { ( my ($m,$r,$c) = ( &OBJECT(shift))) or (return); ( &OBJECT([ map { my ($i) = $_; [ map { $$m[$_][$i] } (0..($r-1)) ] } (0..($c-1)) ])) }

# Returns the inverse of an input square MATRIX object; includes dimension in list context
sub INVERSE { ( my ($m,$e) = ( &OBJECT((shift),(undef),-1))) or (return); ( &SOLVE( $m, [ map {[ ((0)x($_)), (1), ((0)x($e-1-$_)) ]} (0..($e-1)) ] )) }

# Returns the determinant of an input square MATRIX object
sub DETERMINANT { ( scalar &SOLVE(shift)) }

# Returns the trace of an input square MATRIX object
sub TRACE { ( my ($m,$e) = ( &OBJECT((shift),(undef),-1))) or (return undef); ( &::SUM( map { $$m[$_][$_] } (0..($e-1)))) }

# Returns coordinate inversions of input square coefficient matrix and augmenting solution MATRIX objects; includes dimension in list context
sub SOLVE { (( my ($m,$r,$c) = ( &TRANSPOSE([ @{(( &TRANSPOSE(( &::CLONE(shift)),(undef),-1)) or (return))}, ( map {
	(defined) ? @{(( &TRANSPOSE($_)) or (return))} : () } ( &::CLONE(shift))) ]))) || (return)); my ($p,@c) = (1,(0..($c-1)));
	for my $i (0..($r-1)) { do { my ($t) = $_; ($$t[2] > 0) or ( return (($c == $r) ? 0 : ())); do { ($$t[$_] == $i) or do {
		(@{($m,\@c)[$_]}[$$t[$_],$i]) = (@{($m,\@c)[$_]}[$i,$$t[$_]]); $p *= -1 }} for (0,1) } for ( sort { our ($a,$b); ($$b[2] <=> $$a[2]) }
			map { my ($k) = $_; map {[ $k, $_, ( abs ($$m[$k][$c[$_]])) ]} (($i)..($r-1)) } (($i)..($r-1)))[0];
		for my $j (($i+1)..($r-1)) { my ($s) = $$m[$j][$c[$i]]/$$m[$i][$c[$i]]; do { $$m[$j][$c[$_]] -= $s*$$m[$i][$c[$_]] }
			for (($i+1)..($c-1)); }} ( return ( &::PRODUCT( $p, ( map { $$m[$_][$c[$_]] } (0..($r-1)))))) if ($c == $r);
	( &TRANSPOSE([ map { my ($j,$x) = ($_,[]); for my $i ( reverse (0..($r-1))) { $$x[$c[$i]] = ( $$m[$i][$j] -
		( &::SUM( map { $$x[$c[$_]]*$$m[$i][$c[$_]] } (($i+1)..($r-1)))))/($$m[$i][$c[$i]]) }; ($x) } (($r)..($c-1)) ])) }

# Returns a stringified representation of an input MATRIX object
sub STRING { my ($s) = ( scalar &Local::TENSOR::STRING(( scalar &OBJECT(shift)) or ( return q()))); ( $s =~ s/],\[/],\n [/g ); ($s) }

}	# End of package Local::MATRIX

# Encapsulates tensor object definition and manipulation
{; package Local::TENSOR; use Scalar::Util qw(blessed); use overload
	q(+)		=> sub { ( &SUM(shift,shift)) },
	q(*)		=> sub { ( &PRODUCT(shift,shift)) },
	q(-)		=> sub { ( &DIFFERENCE((shift,shift)[ (pop) ? (1,0) : (0,1) ])) },
	q(/)		=> sub { ( &QUOTIENT((shift,shift)[ (pop) ? (1,0) : (0,1) ])) },
	q(**)		=> sub { ( &EXPONENTIATION((shift,shift)[ (pop) ? (1,0) : (0,1) ])) },
	q(x)		=> sub { ( &OUTER_PRODUCT((shift,shift)[ (pop) ? (1,0) : (0,1) ])) },
	q("")		=> \&STRING,
	q(bool)		=> sub { 1 },
	q(fallback)	=> 1;

# Validates and returns a multi-dimensional TENSOR object generated from input nested ARRAY reference; includes dimension in list context
sub OBJECT { my ($b,$o,$s,$f,@d) = (( map {(( &::DEFINED(( blessed $_ ), (( length ref ) ? () : (($_),( blessed $_[0] ))), __PACKAGE__ )),
	(( length ref ) ? \($_) : \(shift)))} (shift)), (( ref $_[0] eq q(CODE)) ? ((shift),1) : (( sub { ${(shift)}} ),!1)), (@_)); do { $b = ( &::DEFINED(
	( blessed ( ${ $o = \( $$$o{TNSR} ) } )), __PACKAGE__ )) } if ( &::ISA( 0, $$o, qw(HASH))); if (( &::ISA( 1, $$o, __PACKAGE__ )) && !($f)) { my ($o,$i) = ($o,0);
	{; (((@$$o) == (( $d[$i++] = ( &::DEFINED($d[$i],0+@$$o))) or (return))) or (return)); do { ((( int $d[$i] ) < 0 ) and ( $d[$i] = $d[$i-1] )); redo }
	if ( &::ISA( 1, ${ $o = \( $$$o[0] ) }, __PACKAGE__ )); }; ( 0+(( splice @d, $i, (@d-$i), (0))[0] ) && (return)); } else { my ($d) = ( shift @d ); ($o) = \( map {
		((( int $d[0] ) < 0 ) and ( $d[0] = $d )); [ map { grep {((defined) or (return))} ((((undef),@d) = ( &OBJECT((undef), $_, (($f)?($s):()), (@d)))),(undef))[0] }
		(@$_) ] } map { (( &::ISA( 0, $_, q(ARRAY))) ? (((@$_) == (( $d = ( &::DEFINED($d,0+@$_))) or (return))) ? ($_) : (return)) :
		((((defined) && !(/inf/) && !(/nan/)) or (return)) && ((( int $d ) > 0 ) ? ((length ref) ? (return) : [ ($_) x ( $d = ( int $d )) ] ) :
		( 0+$d ) ? (return) : do { do {((wantarray) ? ( return ($_,0)) : ( return $_ ))} for ( map {((((defined) && !(/inf/) && !(/nan/)) or (return)) &&
		(( length ref ) ? ( ref =~ /^(?:SCALAR|REF)$/ ) ? ( &::ISA( 0, $_, qw( ARRAY HASH ))) ? (return) : ($_) : (return) : 0+($_)))} ( scalar $s->($o))) } ))) }
	($$o))[0]; ( unshift @d, $d ); } (( &::ISA( 1, ( bless $$o, $b ), __PACKAGE__ )) or (return)); ((wantarray) ? ($$o,@d) : ($$o)) }

# Returns references to the requested elements of an input multi-dimensional TENSOR object; array references span and undef slurps
{; my ($ele); sub ELEMENTS { $ele ||= ( &::Y_COMBINATOR( sub { my ($sub) = (shift); sub {( map {(( &::ISA( 1, $_, __PACKAGE__ )) ?
	( map {(($sub)->($_))} (@$_)) : ($_))} (shift))}} )); ( grep {((wantarray) or (return $_))} (($ele)->(
	( scalar &SLICE((shift), (( ref $_[0] eq q(CODE)) ? (shift) : ( sub {(shift)} )), (@_))) or (return)))) }}

# Returns a TENSOR object slice of an input multi-dimensional TENSOR object; array references span and undef slurps
sub SLICE { (( my ($o,$d) = ( map { my ($s) = (( ref $_[0] eq q(CODE)) && (shift));
	((( &::ISA( 1, $_, __PACKAGE__ )) && !($s)) ? ( $_, 0+(@$_)) : ( &OBJECT((undef), $_, (($s)||())))) } (shift))) or (return));
	( &OBJECT(( blessed $o ), (( $d == 0 ) ? ((@_) ? (return) : ($o,0)) : [
		grep {((defined) or (return))} map {( scalar &SLICE( $$o[$_], (@_)))} map {((defined) ? ( map {(( ref eq q(ARRAY)) ?
		do { ((( my ($n,$x) = grep {((defined) or (return))} map { (0..($d-1))[$_] } (@$_)) == 2 ) or (return)); (($n)..($x)) } :
		( length ref ) ? (return) : ( grep {((exists $$o[$_]) or (return))} (int)))} (( ref eq q(ARRAY)) ?
		( grep {((defined) or (return))} (@$_)) : ($_))) : (0..($d-1)))} (shift) ] ))) }

# Returns a TENSOR object after functional transformation and substitution of a user specified slice
sub SPLICE { (( my ($c) = ( &OBJECT((undef), ( &::CLONE(shift))))) or (return)); my ($s) = (shift); (( my ($o,@d) = (($c)->SLICE(( sub {(shift)} ), (@_)))) or (return));
	do { $${( shift @$_ )} = ${( shift @$_ )}} for ( &::ZIPS( [[ (( scalar &OBJECT((undef), (( ref $s eq q(CODE)) ? ( scalar (($s) ->
	( scalar (($o) -> OBJECT( sub { $${(shift)}} ))))) : ($s)), (@d))) or (return)) -> ELEMENTS() ], [ ($o) -> ELEMENTS() ]], 0 )); (($c) -> OBJECT()) }

# Returns a TENSOR object after functional transformation of specified multi-dimensional projections
sub PROJECT { (( my ($o,@d) = ( &OBJECT((undef), ( &::CLONE(shift))))) or (return)); my ($s) = (( ref $_[0] eq q(CODE)) ? (shift) : ( \&::SUM )); ( pop @d );
	my (@p) = ( map {(!(!$_))} map {(( ref eq q(ARRAY)) ? ((@$_) > (@d)) ? (return) : ( @$_[(0..(@d-1))] ) : (($_)x(@d)))} (shift));
	do { my (@e) = (($o) -> ELEMENTS(@$_)); my ($e) = ( scalar (($s) -> ( map {($$_)} (@e)))); do {( $$_ = $e )} for (@e); } for
		( &::ASSIGNMENTS([ map {(($p[$_]) ? (undef) : [0..($d[$_]-1)])} (0..(@d-1)) ])); ( &OBJECT(( blessed $o ), $o, (@d,0))) }

# Returns the multi-dimensional compounding in user specified directions for values of an input TENSOR object
sub COMPOUND { (( my ($o,@d) = ( &OBJECT((undef), ( &::CLONE(shift))))) or (return)); my ($s) = (( ref $_[0] eq q(CODE)) ? (shift) : ( \&::SUM )); ( pop @d );
	my (@m) = ( map {( $_ <=> 0 )} map {(( ref eq q(ARRAY)) ? ((@$_) > (@d)) ? (return) : ( @$_[(0..(@d-1))] ) : (($_)x(@d)))} (shift));
	do { my ($i,$m,$d,@e) = ($_,$m[$_],$d[$_]); do { my ($j) = $_; do { ${( shift @$_ )} = ( scalar (($s) -> ( map {($$_)} (@$_)))); } for
	( &::ZIPS( [ map {( [ ($o) -> ELEMENTS(((undef)x($i)), ($_), ((undef)x(@d-1-$i))) ] )} (($m*($j+(1-$m)/2)), ($m*($j-(1+$m)/2))) ], 0 )) } for
	( @{ ([],[1..($d-1)])[$m] } ); } for (0..(@d-1)); ( &OBJECT(( blessed $o ), $o, (@d,0))) }

# Returns the rolling average over a specified multi-dimensional correlation length for values of an input TENSOR object
sub SMOOTH { (( my ($o,@d) = ( &OBJECT((undef), (shift)))) or (return)); my ($s) = (( ref $_[0] eq q(CODE)) ? (shift) : ( \&::ARITHMETIC )); ( pop @d );
	my (@w) = do { my ($i); map {( &::BOUNDED([0,($d[$i++]-1)],(int)))} map {(( ref eq q(ARRAY)) ? ((@$_) > (@d)) ? (return) :
	( @$_[(0..(@d-1))] ) : (($_)x(@d)))} (shift) }; my (@i,$r); do { my ($i,$r) = ($_,\($r)); do { $r = \( ${ $$r ||= [] }[$_] ) }
	for (@$i); $$r = ( scalar (($s) -> ( map {($$_)} (($o) -> ELEMENTS( map { my ($j,$i) = ($_,$$i[$_]); $i[$j][$i] ||= do { my ($w,$d) =
	($w[$j],$d[$j]); [ (($i)x($w-$i)), (( &::MAX(0,($i-$w)))..( &::MIN(($i+$w),($d-1)))), (($i)x($i+$w-($d-1))) ] }} (0..(@d-1))))))); }
	for ( &::ASSIGNMENTS([ map {[ 0..($_-1) ]} (@d) ])); ( &OBJECT(( blessed $o ), (($r) or ($o)), (@d,0))) }

# Returns the tensor sum over a list of input TENSOR objects; includes dimension in list context
sub SUM { ( &MAP((@_), \&::SUM )) }

# Returns the element by element (Hadamard) product over a list of input TENSOR objects; includes dimension in list context
sub PRODUCT { ( &MAP((@_), \&::PRODUCT )) }

# Returns the tensor difference between a pair of input TENSOR objects; includes dimension in list context
sub DIFFERENCE { ( &MAP((shift,shift), sub {((shift)-(shift))} )) }

# Returns the element by element quotient of two input TENSOR objects; includes dimension in list context
sub QUOTIENT { ( &MAP((shift,shift), \&::RATIO )) }

# Returns the element by element exponentiation of two input TENSOR objects; includes dimension in list context
sub EXPONENTIATION { ( &MAP((shift,shift), sub {((shift)**(shift))} )) }

# Returns the outer product of a list of input TENSOR objects
sub OUTER_PRODUCT { my ($a,@b) = map {(( &OBJECT((undef), ( &::CLONE($_)))) or (return))} ((shift),(@_));
	for my $b (@b) { for my $e (($a) -> ELEMENTS()) {( $$e *= $b )}} (($a) -> OBJECT()) }

# Returns a stringified representation of an input TENSOR object
sub STRING { ( my ($o,$d) = ( map {( &OBJECT(( blessed $_ ), ($_)))} (shift))) or (return q());
	(($d) ? ( q([).( join qq(,), map {( &STRING($_))} (@$o)).q(])) : ( ref =~ /^(SCALAR|REF)$/ ) ? qq(*${1}*) : ( sprintf q(%+.3E), ( 0+$o ))) }

# Returns the TENSOR object element-by-element map of a user specified function over input TENSOR objects
{; my ($map); sub MAP { $map ||= ( &::Y_COMBINATOR( sub { my ($sub) = (shift); sub { my ($d,$s,@o) = ((pop,pop),(@_)); (($$d[0]) ? [ map { my ($i) = $_;
	(($sub)->(( map {(( &::ISA( 1, $_, __PACKAGE__ )) ? ($$_[$i]) : ($_))} (@o)), $s, [ @$d[1..(@$d-1)]] )) } (0..($$d[0]-1)) ] : ( scalar (($s)->(@o)))) }} ));
	my ($b,$d); my ($s,@o) = (( grep {(( ref eq q(CODE)) or (return))} (pop)), ( map {((( length ref ) && ( ref !~ /^(?:SCALAR|REF)$/ )) ? do {
	grep { ((defined) or (return)); ( $b = ( &::DEFINED( $b, ( blessed $_ )))); 1 } ((((undef),@$d) = ( &OBJECT( $b, $_, @{$d||=[]} ))),(undef))[0] } :
	((((defined) && !(/inf/) && !(/nan/)) or (return)) && (( length ref ) ? ( &::ISA( 0, $_, qw( ARRAY HASH ))) ? (return) : ($_) : 0+($_))))} (@_)));
	( &OBJECT( $b, (($map)->(@o,$s,$d)), @{$d||[0]} )) }}

}	# End of package Local::TENSOR

# Encapsulates polynomial object definition and manipulation
{; package Local::POLY; use overload
	q(+)		=> sub { ( &SUM(shift,shift)) },
	q(*)		=> sub { ( &PRODUCT(shift,shift)) },
	q(-)		=> sub { ( &DIFFERENCE((shift,shift)[ (pop) ? (1,0) : (0,1) ])) },
	q(/)		=> sub { ( &QUOTIENT((shift,shift)[ (pop) ? (1,0) : (0,1) ]))[0] },
	q(%)		=> sub { ( &QUOTIENT((shift,shift)[ (pop) ? (1,0) : (0,1) ]))[1] },
	q(**)		=> sub { my ($b,$n) = (shift,shift)[ (pop) ? (1,0) : (0,1) ];
		(ref $n) ? (undef) : (($n = (int $n)) < 0) ? 1/($b**(-$n)) : ($n == 0) ? ( &OBJECT(1)) : ((1-($n%2))||$b)*($b**($n/2))*($b**($n/2)) },
	q(==)		=> sub { ( return ( $$_[1] eq $$_[0] )) for [ map {( join '_', (@$_))} map { ( &OBJECT($_)) or (return !1) } (shift,shift) ] },
	q(fallback)	=> 1;

# Returns an array reference format copy of polynomial object normalized for numerical entries and nonzero leading coefficent
sub OBJECT { ( return ( bless $_ )) for grep { my ($t) = 0+(shift); ( pop @$_ ) while ((@$_) && (( abs $$_[-1] ) <= ($t))); 1 }
	map {[ ( &::ISA( 0, $_, q(ARRAY))) ? ( map {(0+$_)} (@$_)) : (defined) ? (0+$_) : (return undef) ]} (shift); }

# Returns the derivative of a polynomial object
sub DERIVATIVE { (return $_) for ( grep { my ($i); (shift @$_); ($_ *= ++$i) for (@$_); 1 } (( &OBJECT(shift,shift)) or (return undef))) }

# Returns the sum of a list of polynomial objects
sub SUM { my (@p) = map { ( &OBJECT($_)) or (return undef) } (@_); &OBJECT([ map { my ($i) = $_;
	&::SUM( map { 0+$$_[$i] } (@p)) } (0..( &::MAX( map {(@$_-1)} (@p)))) ]) }

# Returns the difference of a pair of polynomial objects
sub DIFFERENCE { &SUM((shift),( &PRODUCT((shift),[-1]))) }

# Returns the product of a pair of polynomial objects
sub PRODUCT { my ($a,$b) = map { ( &OBJECT($_)) or (return undef) } (shift,shift);
	&SUM( map { my ($t) = $_; ([ map { $_*$t } (@$a) ],( unshift @$a, 0 ))[0] } (@$b)) }

# Returns the quotient and remainder of a pair of polynomial objects
sub QUOTIENT { my ($n,$d,@q) = map { ( &OBJECT($_)) or ( return (undef,undef)) } (shift,shift);
	(0+@$d) or (return (( &OBJECT([])),( &OBJECT([])))); while (@$n >= @$d) { $n = &DIFFERENCE( grep { $$_[@$n-1] = 0; 1 } (
		$n, ( &PRODUCT($d,(($q[@q]) = grep { $$_[@$n-@$d] = ($$n[-1]/$$d[-1]); 1 } [] ))))); }
	map { (wantarray) ? ($_,$n) : (return $_) } ( &SUM(@q)) }

# Returns the numerical evaluation of a polynomial object at a specified coordinate
sub EVALUATE { my ($y,$x) = ((( &OBJECT(shift)) or (return undef)),( map { (defined) or (return undef); (0+$_) } (shift)));
	0+( &Local::VECTOR::INNER_PRODUCT(( scalar &Local::VECTOR::POWERS($x,0+@$y)), $y )) }

# Returns the count of distinct (optionally all) real polynomial object roots within specified (default infinite) inclusive bounds by the Sturm sequence method
sub REAL_ROOTS { my ($p,$b,@s) = (( map { ( &OBJECT($_)) or (return undef) } (shift)),(shift)); do { push @s,
	grep { while ((0+@{$$_[-1]||[]}) > 1) { push @$_, ( grep { (0+@$_) or (last) } ( &PRODUCT(( &QUOTIENT( @$_[-2,-1]))[1],[-1]))) }; 1 }
	map {[ grep {(0+@$_)} ( $_, ( &DERIVATIVE($_))) ]} ($p) } while ((0+@{$p=$s[-1][-1]||[]}) > 1); do { ( return ((defined $b) ? ( $_->($b)) : ($_))) } for
	( &::Y_COMBINATOR( sub { my ($sub) = (shift); sub { my ($b) = ( grep { ( ref eq q(ARRAY)) or (return undef) } (shift));
		(defined $$b[0]) && (defined $$b[1]) && ($$b[0] > $$b[1]) && ( return &::SUM( map {( $sub->([(@$b[0..2],undef)[@$_]]))} ([-1,1,2],[0,-1,2])));
		( &::MAX( 0, ( &::SUM( map { my (@r); map { my ($p) = $_; map { my ($j) = $_; ((+1,-1)[$j]) * ((0+$r[$j]) != 0) * ((0+$r[$j]) != ( $r[$j] =
		(( 0+((defined $$b[$j]) ? do { my ($p,$i,$e) = ($p,0); while ((( $e = ( &EVALUATE($p,$$b[$j]))) == 0) && ((0+@$p) > 1)) { $p =
		( &DERIVATIVE($p)); $i++; }; ($e)*( $j || (1 - 2*($i%2))) } : ($$p[-1])*( $j || (2*((0+@$p)%2) - 1))) <=> 0 ) || (0+$r[$j])))) }
			(0,1) } (@$_) } ( @s[(0..0+(($$b[2]) && (@s-1)))] ))))) }} )) }

# Returns a least squares polynomial with specified limit of terms for input range, domain, and weight ~ 1/Sqrt(Var) factors
sub LEAST_SQUARES { my ($x,$y) = map { my ($t) = $_; my ($e) = ( map {( &::MIN((( &::MAX( 0, ( int (shift)))) || $_ ), $_ ))}
	1+( grep {($$t[$_][0] > $$t[($_-1)][0])} (1..(@$t-1)))); ([ map { ( scalar &Local::VECTOR::POWERS($$_[0],$e,$$_[2])) } (@$t) ],
		[ map {[ $$_[1] ]} (@$t) ]) } map {[ sort { our ($a,$b); ($$a[0] <=> $$b[0]) } (@$_) ]} map { ($$_[2] == 3) ? [ grep { ($$_[2] > 0)
			&& do { $$_[1] *= $$_[2]; 1 }} (@{ $$_[0] }) ] : ($$_[2] == 2) ? ($$_[0]) : () } [ &Local::MATRIX::OBJECT( &::CLONE(shift)) ];
	( &OBJECT( shift @{( &Local::MATRIX::TRANSPOSE( scalar &Local::MATRIX::SOLVE( do { my ($t) = ( scalar &Local::MATRIX::TRANSPOSE($x));
		( map { ( scalar &Local::MATRIX::INNER_PRODUCT($t,$_)) } ($x,$y)) }))) || [] })) }

}	# End of package Local::POLY

# Encapsulates Balanced KD-type TREE object definition and manipulation.
{; my ($srt,$gft,$nbr,$nhd,$rsq); package Local::TREE; use Scalar::Util qw(weaken);

# Returns a new TREE object data structure initialized for dimension. After Procopiuc, Agarwal, Arge, and Vitter (2003)
sub NEW { ( grep { @{ $_ }{( qw( GFT PRN NBR NDX DIM QAD BLK LID FAR PRL NOD LEF CRD MRG MTR RNK BUF BKD CSH NHD ))} = ( do {
	$srt ||= ( sub ($$) { my ($a,$b) = @_; (( $$a[1] <=> $$b[1] ) or ( $$b[2] <=> $$a[2] )) } );
	my (@dim,@far,@prd); do { my ($min,$max,$mod,$def,$prd) = (( ref eq q(ARRAY)) ? (@$_) : ($_))[0..2];
		push @far, (!1,!1,1)[ $mod = 0+(0,+1,-1)[$mod]]; push @dim, (( $def = ( grep {((defined) && (($_+=0),1))} ($min,$max))) ? ( do {
			(((0,1,0)[$mod] ) && ((( $def == 2 ) && (( $prd = ( $max - $min )) > 0 )) or (( return undef )))); (
			sub { do {( return (( &::MATCH_VALUE([$min,$max],$_)) ? ($_) : (undef)))} for (shift) },
			sub { do {( return [ ( $min + $_ ), ( $min + $_ + ((( 2 * $_ ) <= ($prd)) ? ($prd) : ( -1 * $prd ))) ] )} for
				( &::INT_QUOTIENT(((shift) - $min ),($prd),(-1)))[1] } )[ (0,1,0)[$mod]] } ) : ( sub {(shift)} )); push @prd, $prd; }
		for ( map {(( ref eq q(ARRAY)) ? ((@$_) ? (@$_) : (( return undef ))) : ((undef)x( &::MAX( 1, (int)))))} (shift));
	my ($dim,$qad,$blk,@nod,@lef) = (( 0+@dim ), ( map {(( $_ < 0 ) or ( !1, $_ ))} map {( int ((defined) ? ($_) : ( ::BLK )))} (shift)));
	my ($lid,@prl) = ( map {(( ref eq q(CODE)) ? ( do { push @lef, $_; !1 } ) : ((defined) && ( do { my ($i) = 0+((0..($dim-1)),-1)[$_];
		(( $i < 0 ) or ( !1, ( do { (( $prd[$i] ) && ( return undef )); push @nod, ( sub {( ${(shift)}[$i] >= ${(shift)}[$i] )} );
		(($far[$i]) ? ( map { my ($prl) = $_; [ map {(( $_ == $i ) ? ($prl) : (0))} (0..($dim-1)) ] } (-1,+1)) : ()) } ))) } )))} ((@_) ? (shift) : (-1)));
	push @nod, ( map { my ($i,$prd) = ($_,$prd[$_]); (( defined $prd ) ? ( sub {(( 2 * ( abs ( ${(shift)}[$i] - ${(shift)}[$i] ))) <= $prd )} ) : ()) } (0..($dim-1)));
	push @lef, ((($lid) ? ( sub {( ${(shift)}{LID} >= ${(shift)}{LID} )} ) : ()), ( sub {((shift) != (shift))} ));
	my ($crd) = ( map {(( ref eq q(CODE)) ? ($_) : ( sub {(shift)} ))} (shift)); (((0)x(4)), $dim, $qad, $blk, $lid, \@far, \@prl, \@nod, \@lef,
	( sub { ((( my (@crd) = ( map {(0+$_)} map {(( &::ISA( 0, $_, q(ARRAY))) ? (@$_) : ($_))} (($crd) -> (shift)))) == ($dim)) or
		(( return undef ) )); ( &::ASSIGNMENTS([ map {( scalar (($dim[$_]) -> ($crd[$_])))} (0..(@crd-1)) ], !1, 1 )) } ),
	( map {(( ref eq q(CODE)) ? ($_) : ( sub {[ map { my ($i) = $_; ( &::SUM( map {($$_[$i])} (@_))) } (0..($dim-1)) ]} ))} (shift)),
	( map {(( ref eq q(CODE)) ? ($_) : ( sub { my ($a,$b) = (shift,shift); ( &::SUM( map {(($$b[$_]-$$a[$_])**2)} (0..($dim-1)))) } ))} (shift)),
	( map {(( ref eq q(CODE)) ? ($_) : ( sub { ${(shift)||{}}{RSQ}} ))} (shift))) } ); 1 } ( bless +{} ))[0] }

# Populates a TREE object with trunks, branches, and leaves of data
sub GRAFT { my ($obj) = (shift); do { my ($prm); for my $slf ( $_, @{$$_{IMG}||[]} ) {
	( @{ $slf }{(( qw( KDX IDX )), (($prm) ? () : ( q(LID))))} = ((-1), (( push @{$$obj{BUF}||=[]}, ($slf)) - 1 ), (($prm) or ( ++$$obj{GFT} )))) }
		continue {( $prm ||= $slf )}} for ( my (@lvs) = ( map {((($obj) -> LEAF($_)) or (return))} (@_)));
	my ($flr,$clg) = (($$obj{QAD}) ? () : ( &::INT_LOG_TWO(0+@{$$obj{BUF}||[]}))); if (( defined $flr ) && ($flr >= $$obj{BLK})) { my ($buf,$min,$max) =
		( [ grep {(defined)} ( splice @{ $$obj{BUF}} ) ], ($$obj{BLK}), ( grep {((defined) or ( return ((wantarray) ? () : (-1))))} ($clg)));
	do { while ( defined $$obj{BKD}[($max-$$obj{BLK})] ) {((( ++$max ) <= ( int ::BMX )) or ( return ((wantarray) ? () : (-1))))} do { do {
		push @$buf, ( grep {(defined)} @{$$_{BUF}||[]} ); ( undef $_ ); } for ($$obj{BKD}[($_-$$obj{BLK})]) } for (($min)..($max-1)); }
	while ((@$buf) && (($min,$max) = ( map {(($_ > $max) ? ($max+1,$_) : ($_ < $$obj{BLK}) ?
		do { my ($i) = 0; for ( @{$$obj{BUF}} = @$buf ) {( @{ $_ }{( qw( PAR HND KDX IDX ))} = (undef,undef,-1,$i++))} () } :
		do { package Local::TREE::TRUNK; use Scalar::Util qw(weaken); my ($dim,$lvl,$kdx) = ($$obj{DIM},$_,($_-$$obj{BLK}));
		($$obj{BKD}[$kdx]) = ( grep { my ($slf) = $_; ( weaken +( @{ $slf }{( qw( OBJ KDX BUF NOD ))} = ( $obj, $kdx, $buf,
	(( $gft ||= ( &::Y_COMBINATOR( sub { package Local::TREE::BRANCH; use Scalar::Util qw(weaken); my ($sub) = (shift); sub {
		my ($dim,$lvl,$par,$hnd,$kdx,$buf,$idx,$slf,$key,$val,@lft) = ((@_[0..6]), ( bless +{} )); do { do {
		( weaken +( @{ $$buf[$_] }{( qw( PAR HND KDX IDX ))} = ( $par, $hnd, $kdx, $_ ))[0] );
			( return $_ ) } for ( $$idx[0][0] ) } if ( @{$$idx[0]} == 1 ); ( $lft[ $val = $_ ] = 1 ) for
			( @{ $$idx[ $key = ( $lvl % $dim ) ] }[0..(( &::INT_EXP_TWO((( &::INT_LOG_TWO(0+@{$$idx[0]}))[1] ) - 1 )) - 1 )] );
		( weaken +( @{ $slf }{( qw( PAR HND KEY VAL BND RSQ LFT RGT ))} = ( $par, $hnd, $key, $$buf[$val][$key], [ map { my ($i) = $_;
			[ map {( $$buf[$_][$i] )} @{ $$idx[$i] }[(0,-1)] ]} (0..($dim-1)) ], ( map {((($$obj{LID}) ? (undef) : ( &::MAX( grep {(defined)}
			map {( ${(( length ref ) ? ($_) : ($$buf[$_]))}{RSQ} )} (@$_)))), (@$_))} [ map { my ($rgt) = $_; (($sub) ->
		( $dim, ($lvl-1), $slf, ( 0+($rgt)), $kdx, $buf, [ map {[ grep {( $rgt xor $lft[$_] )} @{ $$idx[$_] } ]} (0..($dim-1)) ] )) }
			((undef),1) ] )))[0] ); ($slf) }} ))) ->
		( $dim, $lvl, $slf, (-1), $kdx, $buf, [ map { my ($i) = $_; [ map {($$_[0])} sort { our ($a,$b); ( $$a[1] <=> $$b[1] ) }
			map {[ $_, $$buf[$_][$i]]} (0..(@$buf-1)) ] } (0..($dim-1)) ] ))))[0] ); 1 } ( bless +{} )); () } )}
		grep {((defined) or ( return ((wantarray) ? () : (-1))))} ( &::INT_LOG_TWO(0+@$buf))[1] ))) }
	if ($$obj{NBR}) { for my $slf ( map {( $_, @{$$_{IMG}||[]} )} (($$obj{LID}) ? () : (@lvs))) {
		for my $kdx ((-1)..(@{$$obj{BKD}||[]}-1)) { for my $idx (($kdx < 0) ? (0..(@{$$obj{BUF}||[]}-1)) : (-1)) { (($nbr) ->
			( $obj, $slf, (($idx < 0) ? ( @{(($$obj{BKD}[$kdx]) or (next))}{( qw( BUF NOD ))} ) : ($$obj{BUF},$idx)), (-1), (-1), $$obj{PRL}[1] )) }}}
		if ($$obj{NHD}) { for my $slf (@lvs) {( &::SORT_INSERT( $$obj{NHD}, ((($nhd) -> ($slf)) or (next)), $srt ))}}}
	((wantarray) ? (@lvs) : (0+@lvs)) }

# Vacates a TREE object with trunks, branches, and leaves of data
sub CLEAR { my ($obj) = (shift); for (($obj) -> LEAVES()) {( undef $$_{OBJ} )}; @{ $obj }{( qw( GFT PRN NBR NDX BUF BKD CSH NHD ))} = ((0)x(4)); }

# Validates and returns a TREE::LEAF object generated from input dimension and data point ARRAY reference
sub LEAF { package Local::TREE::LEAF; BEGIN { our (@ISA) = ( q(ARRAY)) }; use Scalar::Util qw(weaken); use overload
	q(@{}) => sub { ${(shift)}{CRD}}, q(fallback) => 1; {; do { ( return ((wantarray) ? (@$_) : ( shift @$_ ))) }
		for do { my ($obj,$raw,$prm) = (shift,shift); [ map { my ($crd) = $_; ((defined) ? ( grep { ( weaken +( @{
	(($prm) ? ( do { package Local::TREE::LEAF::IMAGE; BEGIN { our (@ISA) = ( q(Local::TREE::LEAF)) }; ( bless $_ ) } ) : ( bless $_ )) }{ (
	(($prm) ? ( q(PRM)) : ( qw( OBJ RAW ))), ( qw( CRD PAR HND KDX IDX )), (($prm) ? () : ( qw( LID IMG NBR RSQ AUX NDX )))) } = (( $prm or ($obj,$raw)), $crd ))[0] );
	(($prm) ? ( push @{$$prm{IMG}||=[]}, $_ ) : ( $prm = $_ )) } ( +{} )) : (undef)) } (($$obj{CRD}) -> ($raw)) ] }}

# Returns the TREE:LEAF object(s) nearest to an input TREE::LEAF by local metric consistent with input count and radius
sub NEIGHBOR { my ($buf,$bkd,$cnt,$dom,@nbr) = (( map { @{ $_ }{( qw( BUF BKD ))}}
	(( my $obj = ((( my ($slf) = ( map {(($$_{PRM}) or ($_))} (shift)))[0] ) -> {OBJ} )) or (return))),
	((@_) ? ( map {((defined) ? ($_ >= 0) ? (0,0+($_)) : (( int abs ),-1) : (0,-1))} (shift)) : (1,-1))); (( $cnt == 1 ) &&
	( map { my ($nbr) = $_; ( return ((wantarray) ? ($nbr,$$slf{RSQ}) : ($nbr))) } grep {($$_{OBJ})} (($$slf{NBR}) or ()))) or do {
	for my $kdx ((-1)..(@{$bkd||[]}-1)) { for my $idx (($kdx < 0) ? (0..(@{$buf||[]}-1)) : (-1)) {
		(( $nbr ||= ( &::Y_COMBINATOR( sub { my ($sub) = (shift); sub { my ($obj,$slf,$buf,$nod,$cnt,$dom,$prl,$nbr) = @_;
			my ($sek) = ($cnt >= 0); my ($prm,$lef) = ( map {((($$slf{PRM}) or ($slf)),$_)}
				map {((( grep {(($_) && (($sek) or ( return undef )))} ($$_{PRM})), ($_))[0] )}
				(( length ( ref $nod )) ? () : (( $nod = $$buf[$nod] ) or ( return undef ))));
			my ($ret) = (($sek) ? (undef) : ( grep {((defined) or ( return undef ))} ($$nod{RSQ})));
			for my $sub (($lef) ? (@{$$obj{LEF}}) : ()) {((($sub) -> (($sek)?($prm,$lef):($lef,$prm))) or ( return $ret ))}
			my (@prl,$VTX); for my $i (0..(((@prl) = (($lef) ? () : ( @{$prl||[]}[0..($$obj{DIM}-1)] ))) - 1 )) { my ($far) = $$obj{FAR}[$i]; $prl[$i] ||=
				(( $$slf[$i] <= $$nod{BND}[$i][0] ) ? (($far) ? (+1) : (-1)) : ( $$slf[$i] <= $$nod{BND}[$i][1] ) ? (0) : (($far) ? (-1) : (+1))); }
			my ($vtx) = sub {( $VTX ||= (($lef) ? ($nod) : [ map { my ($prl,$p,$r,$l) = ( $prl[$_], $$slf[$_], @{$$nod{BND}[$_]}[(1,0)] );
				(($prl) ? (((undef), $r, $l )[$prl] ) : ($$obj{FAR}[$_]) ? ((( 2 * $p ) <= ( $l + $r )) ? ($r) : ($l)) : ($p)) } (0..(@prl-1)) ] ))};
			for my $sub (@{$$obj{NOD}}) {((($sub) -> (($sek)?($slf,(($vtx)->())):((($vtx)->()),$slf))) or ( return $ret ))} my ($min,@aux);
			for my $cmp (($sek) ? ((($cnt) && ( @$nbr == $cnt )) ? ($$nbr[-1][1]) : (($dom < 0) ? (($lef) ? (undef) : ()) : ($dom))) : ($ret)) {
				($min,@aux) = @{ ${(($lef) ? ( do { my ($i,$j) = sort { our ($a,$b); ( $a <=> $b ) } map {($$_{LID})} ($prm,$lef);
					\( ${$$obj{CSH}||=[]}[$i][$j] ) } ) : ( do { \( my $tmp ) } ))} ||= [ ($$obj{MTR}) -> ($slf,(($vtx)->())) ] };
				$min += 0; ((defined $cmp) && (( $min <= $cmp ) or ( return $ret ))); }
			if ($lef) { if ($sek) { my (%nbr); @$nbr = ( grep {( !( $nbr{( int $$_[0] )}++ ))}
				( sort { our ($a,$b); ( $$a[1] <=> $$b[1] ) } ( [($lef,$min,\@aux)], @$nbr )));
				(($cnt) && ( splice @$nbr, $cnt )); ( return undef ) } else { ( weaken +( @{ $lef }{( qw( NBR RSQ AUX ))} = ($prm,$min,\@aux))[0] );
				do {( &::SORT_INSERT( $_, ((($nhd) -> ($lef)) or (next)), $srt ))} for (($$obj{NHD}) or ()); ( return $min ) }}
			else { my ($key) = $$nod{KEY}; do {( return (($sek) ? (undef) : (($$nod{RSQ}) = ( &::MAX( grep {(defined)} (@$_))))))} for [
				map { my ($hnd,$prl) = (@$_); my ($nod) = $$nod{(( qw( LFT RGT ))[$hnd] )}; (($prl) && ( $prl[$key] = $prl ));
				(($sub) -> ( $obj, $slf, $buf, $nod, $cnt, $dom, [ @prl ], $nbr )) } ( do { my ($prl) = $prl[$key];
				(($prl) ? ([ (undef,+1,0)[$prl]],[ (undef,0,+1)[$prl]]) : ( do { my ($far) = $$obj{FAR}[$key];
				((( $$slf[$key] <= $$nod{VAL} ) xor ($far)) ? (($far) ? ([0,-1],[+1,0]) : ([0,0],[+1,-1])) :
				(($far) ? ([+1,+1],[0,0]) : ([+1,0],[0,+1]))) } )) } ) ]; }
		}} ))) -> ( $obj, $slf, (($idx < 0) ? ( @{(($$bkd[$kdx]) or (next))}{( qw( BUF NOD ))} ) : ($buf,$idx)), $cnt, $dom, $$obj{PRL}[0], \@nbr ))
	}} (( $cnt == 1 ) ? ( map {(@$_)} grep { (( defined $$slf{RSQ} ) or ( $$obj{NBR}++ ));
		( weaken +(( @{ $slf }{( qw( NBR RSQ AUX ))} ) = ( my ($nbr,(undef,undef)) = (@$_)))[0] );
		(( $rsq ||= sub { my ($obj,$slf,$buf,$hnd,$rsq) = @_[0..2]; (($$obj{LID}) && (return));
			while ( &::ISA( 1, (($slf,$hnd,$rsq) = @{ $slf }{( qw( PAR HND RSQ ))} )[0], q(Local::TREE::BRANCH))) { (($$slf{RSQ}) = ( &::MAX(
			grep {(defined)} (($rsq), ( map {( ${(( length ref ) ? ($_) : ($$buf[$_]))}{RSQ} )} ( $$slf{(( qw( RGT LFT ))[$hnd] )} )))))) }} ) ->
			($obj,$slf,( map {(($_ < 0) ? ($buf) : ($$bkd[$_]{BUF}))} ($$slf{KDX}))));
		((wantarray) or ( return $nbr )) } (( shift @nbr ) or [((undef),( ::INF ))] )) : ((wantarray) ? (@nbr) : (\@nbr)))
}}

# Deletes a leaf of data from a TREE object
sub PRUNE { (((( my $obj = ((( my ($slf) = ( map {(($$_{PRM}) or ($_))} (shift)))[0] ) -> {OBJ} )) or (( return undef ))) -> {PRN} )++ );
	( undef $$slf{OBJ} ); (( defined $$slf{RSQ} ) && ( do {(( --$$obj{NBR} ) or ( undef $$obj{NHD} ))} ));
	my ($prm); for my $slf ( $slf, @{$$slf{IMG}||[]} ) { my ($buf) = ${(( @{$$obj{BKD}||[]}, $obj )[$$slf{KDX}] )}{BUF};
		my ($par); if ( $par = $$slf{PAR} ) { my ($sib); (( &::ISA( 1, $par, q(Local::TREE::TRUNK))) && ( do { ( undef ( $$obj{BKD}[$$slf{KDX}] )); (next) } ));
		( weaken +((( @{(( length ( ref $sib )) ? ($sib) : ( $sib = $$buf[$sib] ))}{( qw( PAR HND ))} ), ( $$par{PAR}{(( qw( LFT RGT NOD ))[$$par{HND}] )} )) =
			(( @{ $par }{( qw( PAR HND ))} ), ((undef) = ( $sib = $$par{(( qw( RGT LFT ))[$$slf{HND}] )} ))))[0] );
		(( !($prm)) && ( defined $$slf{RSQ} ) && (($rsq) -> ($obj,$sib,$buf))); do { my ($slf,@bdx) =
			( $sib, ( map {([$_,0],[$_,-1])} (0..($$obj{DIM}-1)))); while ( &::ISA( 1, ( $slf = $$slf{PAR} ), q(Local::TREE::BRANCH))) {
			my ($bnd,@bnd) = ( $$slf{BND}, ( map {(( length ref ) ? ($$_{BND}) : ($$buf[$_]{CRD}))} ( @{ $slf }{( qw( LFT RGT ))} )));
			(((@bdx) = ( grep { my ($i,$j) = @$_; ( map { (( $$bnd[$i][$j] != $_ ) && ( 1, ( $$bnd[$i][$j] = $_ )))}
			( sort { our ($a,$b); ( $a <=> $b ) } ( map {(( length ref ) ? ($$_[$j]) : ($_))} map {($$_[$i])} (@bnd)))[$j] )[0] }
			(@bdx))) or (last)) }}} ( undef ( $$buf[$$slf{IDX}] )) } continue {( $prm ||= $slf )} 1 }

# Merges multiple input TREE::LEAF objects into a composite data structure
sub MERGE { (( my $obj = ((( my ($slf,@slf) = ( map {( {((int) => (1))}, $_ )}
		map {(($$_{PRM}) or ($_))} (shift)))[-1] ) -> {OBJ} )) or (( return undef )));
	do {( push @slf, $_ )} for ( grep {(((($$_{OBJ}) == ($obj)) && ( !( $$slf{(int)}++ ))) or ( return undef ))}
		map {(($$_{PRM}) or ($_))} grep {(( &::ISA( 1, $_, q(Local::TREE::LEAF))) or (( return undef )))} (@_));
	do {(($_) -> PRUNE())} for (@slf); ((($obj) -> GRAFT( scalar (($$obj{MRG}) -> (@slf)))), (undef))[0] }

# Returns the raw input from which the data array of a TREE::LEAF object was derived
sub RAW {((( map {(($$_{PRM}) or ($_))} (shift))[0] ) -> {RAW} )}}

# Returns all TREE::LEAF objects attached to the TREE object or the count of such in scalar context
sub LEAVES {( map {((wantarray) ? ( grep {((defined) && !($$_{PRM}))} map { @{${$_||{}}{BUF}||[]}}
	(($_), @{$$_{BKD}||[]} )) : ( return ( $$_{GFT} - $$_{PRN} )))} (shift))}

# Returns the most proximal TREE::LEAF object pair grafted into an input TREE object
sub NEIGHBORHOOD { my ($obj,$slf,$rnk,$ndx) = (shift); $$obj{NHD} ||= [ ( sort $srt ( map {(( $nhd ||=
	( sub {( grep {((wantarray) or ( return $_ ))} map {[ $_, 0+(($$_{OBJ}{RNK}) -> ($_)), ( $$_{NDX} = ++$$_{OBJ}{NDX} ) ]}
		grep {((($_) -> NEIGHBOR()) or (( !(wantarray)) && ( return undef )))} (shift))} )) -> ($_))} (($obj) -> LEAVES()))) ];
	while (( @{ $$obj{NHD}} ) ? (($slf,$rnk,$ndx) = @{ $$obj{NHD}[0] } ) : (( undef $$obj{NHD} ) or (return))) {
		(((($$slf{OBJ}) && ( $$slf{NDX} == $ndx )) or ( undef $slf )) && (${$$slf{NBR}||{}}{OBJ}) && (last));
		( shift @{( $$obj{NHD} )} ); (($slf) && ( &::SORT_INSERT( $$obj{NHD}, ((($nhd) -> ($slf)) or (next)), $srt ))); }
	((wantarray) ? ( $slf, @{ $slf }{( qw( NBR RSQ ))}, $rnk, @{ $$slf{AUX}} ) : ($slf)) }

# Returns the dimension of an input TREE object
sub DIMENSION {( ${(shift)}{DIM} )}

}	# End of package Local::TREE

1

