#!/usr/bin/perl

#*******************************#
# rhadamanthus.pl Version 1.8	#
# September '14 - July '20	#
# Joel W. Walker		#
# Sam Houston State University	#
# jwalker@shsu.edu		#
# arXiv:XXXX.XXXX		#
# Copy: GNU Public License V3	#
# Student Assistants:		#
# Kebur Fantahun		#
# B. Ashen Fernando		#
# Nicolle Schachtner		#
#*******************************#

# Apply a strict coding pragma and locate filesystem
use strict; use sort q(stable); use FindBin qw($Bin); use lib qq($Bin);

# Import AEACuS subroutine library and perform version compatibility check
require q(aeacus.pl); ( &UNIVERSAL::VERSION(q(Local::AEACuS),3.30)); our ($OPT);

# Read event plotting specifications from cardfile
my ($PLT) = map { (/^(.*\/)?([^\/]*?)(?:\.dat)?$/); my ($crd,$err,$fil) =
	( &LOAD_CARD([ (($1) or ( q(./Cards/))), [ ((defined) ? ($2) : ( q(plt_card))), q(), q(.dat) ]]));
	($crd) or ((defined) ? ( die 'Cannot open card file for read' ) : (exit));
	( die ( join "\n", ( 'Malformed instruction(s) in card '.($$fil[0].$$fil[1]),
		( map {( "\t".'* Line '.$$_[0].':'."\t".'>> '.$$_[1].' <<' )} (@$err)), q()))) if (@$err);
	@$crd{( qw( plt ))}} ( &$OPT( q(crd)));

# Establish whether Python 2.7/3.X and MatPlotLib 1.3.0+ are suitably configured for piped system calls
my ($cpl) = ( &CAN_MATPLOTLIB());

# Generate histograms
for my $dim (1..3) { my ($hky) = ((undef), qw( hst h2d h3d ))[$dim]; my ($def) = ( ${$$PLT{$hky}||[]}[0] || {} ); HIST: for my $i (1..(@{$$PLT{$hky}||[]}-1)) {

	my ($hst) = (( $$PLT{$hky}[$i] ) || (next HIST)); do {(( exists $$hst{$_} ) or ( $$hst{$_} = $$def{$_} ))} for ( keys %{$def} );
	my ($ipb,$fix) = do { my ($lum,$ipb,$fix) = [ @$hst{( qw( ipb ifb iab izb iyb ))} ]; for my $i (0..4) {
		($ipb,$fix) = @{(($$lum[$i]) or (next))}; $ipb *= (10)**(3*$i); (last) } ( map {( $_, ((defined) && ($fix < 0)))} ($ipb)) };
	my ($obj,$eql,$wdt,$edg,$cnt) = map {( $_, ( map {( [ map {( &EQUAL(@$_))} (@$_) ], $_ )} [ ($_) -> WIDTHS() ] ),[ ($_) -> EDGES() ], [ ($_) -> CENTERS() ] )}
		grep { (( grep {($_ > 3)} (( &Local::TENSOR::OBJECT((undef), $_ )), (undef,undef))[1..$dim] ) == $dim ) or
			do { print STDERR 'INVALID BINNING SPECIFICATION IN HISTOGRAM '.$i."\n"; (next HIST) }}
		do { my (@t) = map {((@{$_||[]} == 1) ? do { my ($t) = @$_; [ map {[$t]} (1..$dim) ] } : ( &SPANS($dim,$_)))} @$hst{( qw( lft rgt spn bns ))};
			( &Local::HISTOGRAM::NEW( map {[ map {[ grep {(defined)} @{( shift @$_ )} ]} (@t) ]} (1..$dim))) };
	my ($sum,$nrm,$per,$avg) = map {[ ((@{$_||[]}) ? ( map {((defined) ? (0+ $_) : (undef))} (@$_)) : (undef)) ]} ( @$hst{( qw( sum nrm per avg ))} );
	my ($out,$nam,$fmt) = ( @$hst{( qw( out nam fmt ))} );
	my ($ttl,$lbl,$lgd) = map {[ map { ((defined) && !(ref)) ? qq($_) : (undef) } (@{$_||[]}) ]} ( @$hst{( qw( ttl lbl lgd ))} );
	my (@vls) = do { my ($chn,$set) = []; map { my ($sub,@cid) = ((ref eq 'ARRAY') ? (@$_) : ((ref) or !(defined)) ? (undef) : ( sub {(shift)} , $_ ));

		# Merge, transform, normalize, integrate, and average a dataset and its threaded channels
		map { my ($s,$n,$p,$a) = map {[ (($dim == 1) ? ($$_[(@$_-1)&&($set)]) : (@$_-1) ? (@$_[0..($dim-1)]) : (($$_[0])x($dim))) ]} ($sum,$nrm,$per,$avg);

			# Apply indicated noise averaging along each dimensional axis
			map { $a = [ map {( &MAX(0,(int)))} (@$a) ]; (( grep {($_)} (@$a)) ? ( scalar (($_) -> SMOOTH($a))) : ($_)) }

			# Correct normalization for differential binning with variable width along each dimensional axis
			map { my ($msr) = 1; my (@msr) = ( map { my ($i) = $_; (($$s[$i]) or ($$p[$i] < 0) or ( do { $msr *=
				(((0+ $$p[$i]) or ($$n[$i] != 0) or ($$wdt[$i][0])) / (($$eql[$i]) ? ($$wdt[$i][0]) : (1))); ($$eql[$i]) } )) } (0..($dim-1)));
				((($_)*($msr)) / (( grep {(!($_))} (@msr)) ? ( &Local::TENSOR::OUTER_PRODUCT(
					map {(($$msr[$_]) ? [(1)x(@{$$wdt[$_]})] : ($$wdt[$_]))} (0..($dim-1)))) : (1))) }

			# Extract portion of data intended for visualization, discarding undefined (0), underflow (1), and overflow (-1) bins
			map {( scalar (($_) -> SLICE( map {[[2,-2]]} (0..($dim-1)))))}

			# Apply indicated post-normalization along each dimensional axis
			map { my ($p,$z) = ([ map {(($$n[$_] > 0) && (! ($$s[$_])))} (0..($dim-1)) ]);
				(( grep {($_)} (@$p)) ? (($_) / ( map {(($z) ? ( do { print STDERR 'INDETERMINATE NORMALIZATION IN HISTOGRAM '.$i."\n"; 1 } ) : ($_))}
					( scalar (($_) -> PROJECT(( sub { ( grep {((0+ $_) or (++$z))} ( &::SUM ))[0] } ), $p ))))[0] ) : ($_)) }

			# Apply and validate binwise functional transformation and recombination to datasets threaded across channels
			grep { (( &ISA( 1, $_, q(Local::TENSOR))) && (++$set)) or
				do { print STDERR 'INVALID CHANNEL OBJECT GENERATED IN HISTOGRAM '.$i."\n"; !1 }} ( scalar (($sub) -> (

			# Apply indicated left/right binwise integration to each threaded channel along each dimensional axis
			map { (( grep {($_)} (@$s)) ? ( scalar (($_) -> SPLICE(( sub {((shift) -> COMPOUND($s))} ), ( map {[[1,-1]]} (0..($dim-1)))))) : ($_)) }

			# Apply indicated pre-normalization to each threaded channel along each dimensional axis
			map { ($s) = [ map {($_ <=> 0)} (@$s) ]; my ($p,$z) = ([ map {(($$n[$_] < 0) or (($$s[$_]) && ($$n[$_] > 0)))} (0..($dim-1)) ]);
				(( grep {($_)} (@$p)) ? (($_) / ( map {(($z) ? ( do { print STDERR 'INDETERMINATE NORMALIZATION IN HISTOGRAM '.$i."\n"; 1 } ) : ($_))}
					( scalar (($_) -> PROJECT(( sub { ( grep {((0+ $_) or (++$z))} ( &::SUM ))[0] } ), $p ))))[0] ) : ($_)) }

			# Make an independent copy of each threaded channel in a dataset
			map {( scalar &Local::TENSOR::OBJECT((undef), ( &CLONE($_))))} (@$_)))) }

		# Read and bin data files into channels, combining like samples by luminosity and discrete samples by cross section
		do { my (@set); CHN: {; my ($j); my (@chn) = grep { ((@$_ == 1) or ((( $j ||= @$_ ) == @$_ ) && ($j))) or
				do { print STDERR 'DATA SET MULTIPLICITY MISMATCH IN HISTOGRAM '.$i."\n"; (last CHN) }}

			map { $$chn[$_] ||= do { my ($dat,$key,$esc,$wgt) = @{ $$PLT{chn}[$_] or 
				do { print STDERR 'CHANNEL '.$_.' IS NOT DEFINED'."\n"; (last CHN) }}{( qw( dat key esc wgt ))}; [

				map { my ($dir,$fil,%bin,%ipb) = @{ $$PLT{dat}[$_] or
						do { print STDERR 'DATA SET '.$_.' IS NOT DEFINED'."\n"; (last CHN) }}{( qw( dir fil ))};
					my ($dir) = ( map { ((defined) ? qq($_) : q(./Cuts/)) } ($$dir[0]));

					FIL: for my $fil ( sort SORT_LIST_ALPHA ( values %{{
							map {(($$_[0].$$_[1]) => ($_))} grep {( $$_[1] =~ /^[\w-]+\.cut$/ )}
							map {( &Local::FILE::LIST([$dir,$_]))} grep {(defined)} (@{$fil||[]}) }} )) {

						( my $tag = $$fil[1] ) =~ s/(?:_\d+)*\.cut$//; ( my ($FHI) = ( &Local::FILE::HANDLE($fil))) or
							do { print STDERR 'CANNOT READ FROM FILE '.$$fil[0].$$fil[1]."\n"; (last CHN) };
						my ($xsc,$idx) = ( &AEACUS_XSEC( $FHI )); my ($e,$x,$l,$s,$r) = (@$xsc); my ($l,$w) = (( &RATIO($e,$x)), ( &RATIO($x,$e)));
						(defined $e) or do { print STDERR 'CANNOT ESTABLISH EVENT COUNT FOR FILE '.$$fil[0].$$fil[1]."\n"; (last CHN) };
						(defined $x) or do { print STDERR 'CANNOT ESTABLISH EVENT CROSS SECTION FOR FILE '.$$fil[0].$$fil[1]."\n"; (last CHN) };
						($s > 0) or do { print STDERR 'NO SURVIVING EVENTS IN FILE '.$$fil[0].$$fil[1]."\n"; (next FIL) };
						(defined $idx) or do { print STDERR 'CANNOT ESTABLISH STATISTICS INDEX FOR FILE '.$$fil[0].$$fil[1]."\n"; (last CHN) };

						my ($key) = [ map {(( &HASHED_FUNCTIONAL( $idx, ((ref eq 'ARRAY') ? (@$_) : ((undef),$_)))) or (
							do { print STDERR 'INVALID CHANNEL KEY SPECIFICATION IN HISTOGRAM '.$i."\n"; (last CHN) } ))} @{$key||[]}[0..($dim-1)]];
						my (@cut) = grep {(( $$_[2] = ( &HASHED_FUNCTIONAL( $idx, ( map {((ref eq 'ARRAY') ? (@$_) : ((undef),$_))} ($$_[2][0]))))) or (
								do { print STDERR 'INVALID KEY IN SELECTION '.$$_[1].' FOR HISTOGRAM '.$i.' ON FILE '.$$fil[0].$$fil[1]."\n"; !1 } ))}
							grep {(! ( &MATCH_VALUE( $$_[3], undef )))} map {[ ($_ < 0), 0+( $_ = ( int abs )), ( @{ (($_) && ( $$PLT{esc}[$_] )) or
								do { print STDERR 'INVALID EVENT SELECTION CUT SPECIFICATION IN HISTOGRAM '.$i."\n"; +{}}}{( qw( key cut ))} ) ]}
							grep {(defined)} (@{$esc||[]});
						my ($wgt) = map {((defined) ? (( &HASHED_FUNCTIONAL( $idx, ((ref eq 'ARRAY') ? (@$_) :
								((undef),((ref eq 'HASH') ? ($_) : +{ wgt => (0+ $_) } ))))) or
							do { print STDERR 'INVALID CHANNEL WEIGHT SPECIFICATION IN HISTOGRAM '.$i."\n"; (last CHN) } ) :
							(( &HASHED_FUNCTIONAL( $idx, (undef), +{ wgt => 0 } )) or ( sub {($w)} )))} (${$wgt||[]}[0]);
						my ($nmx) = (( do { my ($f); ((($fix) && ( &MATCH_VALUE( [0,1], (($f) = (($w)*($ipb-$ipb{$tag})))))) ?
							(( &ROUND(($e)*($f))), ( $ipb{$tag} = $ipb )) : (($e), ( $ipb{$tag} += $l )))[0] } ) or (next));

						push @{ $bin{$tag} ||= [] }, [ $l, ( my ($bin) = (($obj) -> NEW())) ]; local ($_); while (<$FHI>) {
							((/^\s*$/) and (next)); ((/^\s*(\d+)/) && ($1 <= $nmx)) or (last);
							do { my ($val) = $_; (($bin) -> BIN( [ map {(($_) -> ($val))} (@$key) ], (($wgt) -> ($val)))) } for
							grep { my ($val,$mch) = ($_,1); for (@cut) { my ($inv,$eid,$key,$cut) = (@$_);
								(( $mch = (($inv) xor ( &MATCH_VALUE( $cut, (($key) -> ($val)))))) or (last)) }; ($mch) }
							[ map { (/^UNDEF$/) ? (undef) : (0+ $_) } ( split ) ]; }}

					( scalar (($obj) -> SUM( map { my ($tag) = $_; my ($scl) = ( &RATIO(( &DEFINED($ipb,1)), $ipb{$tag} ));
						if ((defined $ipb) && (($scl < 0) or ($scl > 1))) { print STDERR 'RESCALING BY '.( sprintf '%+10.3E', ($scl)) .
							' TO TARGET LUMINOSITY OF '.( sprintf '%+10.3E', ($ipb)).' PER PB IN CHANNEL '.($tag)."\n"; }
 						map {(($scl)*($$_[0])*($$_[1]))} (@{ $bin{$tag}}) } (keys %bin)))) }

				map { ( int ( &MAX(0,$_))) or do { print STDERR 'INVALID DATA SET SPECIFICATION'.$_."\n"; (last CHN) }} grep {(defined)} (@{$dat||[]}) ] }}

			map { ( int ( &MAX(0,$_))) or do { print STDERR 'INVALID CHANNEL SPECIFICATION'.$_."\n"; (last CHN) }} grep {(defined)} (@cid);

			# Pass through datasets one at a time, threading over channels indicated for merging and functional transformation
			(@set) = map { my ($i) = $_; [ map {( $$_[(@$_-1)&&($i)] )} (@chn) ] } (0..(($j)&&($j-1))); } (@set) }} (@{$$hst{chn}||[]}) };

	do { print STDERR 'NO BINNED EVENTS ESTABLISHED FOR HISTOGRAM '.$i."\n"; (next HIST) } unless (@vls);

	my ($log,$min,$fpo); my (%dat) = (
		DIM	=> $dim,
		VAL	=> '['."\n".( join ','."\n", map {( qq($_))} map {(($dim == 1) ? ($_) : ( scalar &Local::MATRIX::TRANSPOSE($_)))} (@vls)).']',
		BIN	=> (($dim == 1) ? ( scalar &Local::TENSOR::STRING($$edg[0])) : ( '['.( join ',', map {( &Local::TENSOR::STRING($_))} (@$edg[0..($dim-1)])).']' )),
		LOG	=> (( $log = ($$hst{log}[0] > 0)) ? q(True) : q(False)),
		STK	=> (($$hst{stk}[0] > 0) ? q(True) : q(False)),
		MIN	=> ( map { ((defined) && (($_ > 0) or (!($log)))) ? ( $min = $_ ) : q(None) } ($$hst{min}[0]))[0],
		MAX	=> ( map { ((defined) && ((defined $min) ? ($_ > $min) : (($_ > 0) or (!($log))))) ? ($_) : q(None) } ($$hst{max}[0]))[0],
		LOC	=> ( map { ${{ TR => 1, TL => 2, BL => 3, BR => 4, RR => 5, ML => 6, MR => 7, BC => 8, TC => 9, MC => 10 }}{( uc $_ )} or
				0+(0..10)[$_] } ($$hst{loc}[0]))[0],
		CLR	=> (( @{ $$hst{clr}||[]} < 2 ) ? ((0..8)[ $$hst{clr}[0]] || (1)) :
				('('.( join ', ', map { m/^(?:([a-z]{2,})|([0-9a-f]{6})|(\d+\.\d*|\.\d+)|(r|g|b|w|c|m|y|k))$/i ?
				(($1)?q({):q()).q(").(($2)?q(#):q()).(($3) ? ( sprintf "%.2f", 0+( &BOUNDED([0,1],$_))) : ( lc $_ )).q(").(($1)?q(:1}):q()) :
				q("") } (@{ $$hst{clr}||[]})).')')),
		HCH	=> (( @{ $$hst{hch}||[]} < 2 ) ? ((0..3)[ $$hst{hch}[0]] || (1)) : ('('.( join ', ', map { m/^[\\\/|xoO.*+-]{1,3}$/ ? q(").(
				do { ( my $t = $_ ) =~ s/\\/\\\\/g; $t } ).q(") : q("") } (@{ $$hst{hch}||[]})).')')),
		FLL	=> '('.(( join ' ', map { ( 0+(1,1,0)[0+(0..2)[$_]]).q(,) } (@{ $$hst{fll}||[]})) or q(1,)).')',
		BLD	=> '('.(( join ' ', map { ( 0+(0,1,0)[0+(0..2)[$_]]).q(,) } (@{ $$hst{bld}||[]})) or q(0,)).')',
		TTL	=> '['.( join ',', ( map {( &RAW_STRING($_))} @$ttl[0,1] )).']',
		XLB	=> ( q()).( &RAW_STRING($$lbl[0])),
		YLB	=> ( q()).( &RAW_STRING($$lbl[1])),
		LGD	=> (( grep {(defined)} (@$lgd)) ? '['.( join ',', map {( q()).( &RAW_STRING($_))} (@$lgd[0..(@vls-1)])).']' : q(False)),
		OUT	=> do { my ($d,$f,$e) = ((( &Local::FILE::PATH(( $out = q().( &DEFINED($$out[0],q(./Plots/)))), 2 )) or ( die 'Cannot write to directory '.$out )),
				( map { (/^[\w:~-]+$/) ? ($_) : ( sprintf "HST_%3.3i", $i ) } qq($$nam[0])),
				( map { s/^\.//; ( ${{ map {( $_ => 1 )} ( qw( pdf eps svg png jpg )) }}{$_} ) ? ( q(.).($_)) : ( q(.pdf)) } ( lc $$fmt[0] )));
				(((undef),$fpo) =	map {( join q(), (@$_))} ((($$fmt[1] > 0) or ( !($cpl) &&
					do { print STDERR 'CANNOT VERIFY PYTHON 2.7/3.X WITH MATPLOTLIB 1.3.0+ (DELIVERING SCRIPT LITERAL)'."\n"; 1 } )) ?
					([q(./),($f),($e)],[($d),($f),q(.py)]) : ([($d),($f),($e)])))[0] },
	);

	do { use Fcntl qw(:seek); local ($.,$?); my ($t,$FHO) = ( tell DATA ); ( defined $fpo ) ?
		do { ( $FHO = ( &Local::FILE::HANDLE($fpo,1))) or ( die 'Cannot write to file '.($fpo)) } :
		do { ( open $FHO, q(|-), q(python 2>&1)) or ( die 'Cannot open pipe to Python' ) };
		local ($_); while (<DATA>) { s/<\[(\w+)]>/$dat{$1}/g; ( print $FHO $_ ) }
		( close $FHO ) && (($? >> 8) == 0 ) && ( defined $fpo ) && ( chmod 0755, $fpo ); ( seek DATA, $t, SEEK_SET ) }}}

1

__DATA__
#!/usr/bin/env python

import sys
if ((sys.version_info[0] < 2) or ((sys.version_info[0] == 2) and (sys.version_info[1] < 7))) :
	sys.exit( 'RHADAManTHUS requires Python versions 2.7 or 3.X' )

import matplotlib as mpl
if (( tuple( map ( int, mpl.__version__.split("."))) + (0,0,0))[0:3] < (1,3,0)) :
	sys.exit( 'RHADAManTHUS requires MatPlotLib version 1.3.0 or Greater' )

import warnings as wrn; wrn.filterwarnings("ignore")

import matplotlib.font_manager as mfm
font = mfm.FontProperties( fname=( list( filter( lambda x: "/times new roman.ttf" in x.lower(), mfm.findSystemFonts(fontext="ttf")))
	or [ mfm.findfont( mfm.FontProperties(family="serif")) ] )[0], size=12 )
mpl.rcParams['mathtext.fontset'] = 'stix'; mpl.rcParams['font.family'] = 'STIXGeneral'

import matplotlib.pyplot as plt
dim = <[DIM]>
fig = plt.figure( figsize=(( 7.5 if dim == 1 else 6.5 ),5), tight_layout=True )
ax = fig.add_subplot(1,1,1)

import math

import numpy as np

color = <[CLR]>
clist = set( k.lower() for k in  mpl.colors.cnames )
clrs = ((),
	( "#CE0000", "#1500AA", "#107C00", "#EF6C00", "#9500AF", "#00B2DB", "#3DD100", "#FFCE00" ),
	( "#E41A1C", "#377EB8", "#4DAF4A", "#FF7F00", "#984EA3", "#FFFF33", "#A65628", "#F781BF" ),
	( "#1F78B4", "#33A02C", "#E31A1C", "#FF7F00", "#A6CEE3", "#B2DF8A", "#FB9A99", "#FDBF6F" ),
	( "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#E6F598", "#ABDDA4", "#66C2A5", "#3288BD" ),
	( "#B10026", "#E31A1C", "#FC4E2A", "#FD8D3C", "#FEB24C", "#FED976", "#FFEDA0", "#FFFFCC" ),
	( "#004529", "#006837", "#238443", "#41AB5D", "#78C679", "#ADDD8E", "#D9F0A3", "#F7FCB9" ),
	( "#084081", "#0868AC", "#2B8CBE", "#4EB3D3", "#7BCCC4", "#A8DDB5", "#CCEBC5", "#E0F3DB" ),
	( "#4D004B", "#810F7C", "#88419D", "#8C6BB1", "#8C96C6", "#9EBCDA", "#BFD3E6", "#E0ECF4" ),
	) # http://colorbrewer2.org/
color = clrs[color] if isinstance(color, int) else tuple( map( lambda x:
	(( lambda x: x if x in clist else "blue" )( list( x.keys())[0].lower())) if isinstance(x,dict)
	else ( x if len(x)>0 else "blue" ), color ))
c = len(color)

hatch = <[HCH]>
hchs = ((),
	( "/"*3, "\\"*3, "|"*2, "-"*2, "x"*3, "+"*2, "."*2, "O"*2 ),
	( "/", "\\", "|", "-", "x", "+", ".", "O" ),
	(),
	)
if isinstance(hatch, int): hatch = hchs[hatch]
h = len(hatch)

fill = <[FLL]>
f = len(fill)

bold = <[BLD]>
b = len(bold)

stack = <[STK]>

log = <[LOG]>
if log: ax.set_yscale("log")

wght = <[VAL]>
n = len(wght)
if n == 0: sys.exit(0)

bins = <[BIN]>

(ymin, ymax) = (<[MIN]>, <[MAX]>)

ttl = <[TTL]>

if dim == 1:

	if stack:
		for i in range(1,n):
			wght[i] = [ sum(l) for l in zip(wght[i],wght[i-1]) ]

	if ymin is not None:
		ymin = float(ymin)
	if log and ( ymin is None or ymin <= 0.0 ):
		ymin = math.pow( 10, math.floor( math.log10( min( [ k for j in wght for k in j if k > 0.0 ] or [1.0] ))))
	if ymax is not None:
		ymax = float(ymax)
		if ymax <= ymin: ymax = None

	edges = [ k for j in zip(bins[:-1],bins[1:]) for k in j ]
	widths = [ j[1] - j[0] for j in zip(bins[:-1],bins[1:]) ]
	lower = [ ymin/2.0 if log else 0.0 ]*(len(bins)-1)*2

	for i in range(n):
		upper = [ k if ( ymin is None or k >= ymin ) else ymin/2.0 if log else ymin-1.0 for j in zip(wght[i],wght[i]) for k in j ]
		if fill[i%f]: ax.fill_between(edges, upper, lower, linewidth=0.0, alpha=0.6, edgecolor="none", facecolor=color[i%c], hatch="", zorder=i-3*n )
		if h: ax.fill_between(edges, upper, lower, linewidth=0.0, alpha=1.0, edgecolor="0.25", facecolor="none", hatch=hatch[i%h], zorder=i-2*n )
		if bold[i%b]: ax.plot( edges, upper, color="#DDDDDD", linewidth=2.0, linestyle="solid", zorder=i-1*n )
		ax.plot( edges, upper, color=color[i%c], linewidth=(2.0,2.0)[bold[i%b]], linestyle=("solid","dashed")[bold[i%b]], zorder=i-1*n )
		if stack: lower = upper

	if ( ymin is None or ymin <= 0.0 ):
		ax.plot( edges, [0.0]*(len(bins)-1)*2, color="black", linewidth=0.8, linestyle="solid", zorder=0 )

	lgd = None
	lgnd = <[LGD]>
	if lgnd:
		patches = [(
			( mpl.patches.Rectangle((0,0), 1, 1, fill=True, facecolor=color[i%c], alpha=0.6, linewidth=0.0 ) if fill[i%f] else ()),
			( mpl.patches.Rectangle((0,0), 1, 1, fill=None, color="0.25", hatch=hatch[i%h], linewidth=0.0 ) if h else ()),
			( mpl.patches.Rectangle((0,0), 1, 1, fill=None, linestyle="solid", edgecolor="#DDDDDD", linewidth=1.4 ) if bold[i%b] else ()),
			mpl.patches.Rectangle((0,0), 1, 1, fill=None, linestyle=("solid","dashed")[bold[i%b]], edgecolor=color[i%c], linewidth=1.4 ))
			for i in range(n) ]
		lgd = ax.legend( patches, lgnd, loc=<[LOC]>, prop=font )
		# lgd = ax.legend( patches, lgnd, loc="center right", bbox_to_anchor=(1.22, 0.5), prop=font )
		lgd.get_frame().set( facecolor="white", linewidth=0.8, edgecolor="black", alpha=1.0 ); lgd.set(zorder=0)

	for (_,i) in ax.spines.items(): i.set( linewidth=0.8, color="black", alpha=1.0, zorder=0 )
	for i in (ax.get_xticklabels() + ax.get_yticklabels()): i.set_fontproperties(font)
	ax.tick_params( axis="both", which="major", zorder=0 ); ax.minorticks_on()

	if ( len(widths) <= 12 and sum( x!=1.0 for x in widths ) == 0 ):
		ax.tick_params(axis="both", which="minor", left="on", right="on", top="off", bottom="off");
		ax.set_xticks(bins, minor=False); ax.set_xticklabels(map((lambda x:""), bins), minor=False)
		xt = range(int(math.ceil(bins[0])),1+int(math.floor(bins[-1])))
		ax.set_xticks(xt, minor=True); ax.set_xticklabels(map(str,xt), minor=True)

	ax.set_xlim([bins[0],bins[-1]]); ax.set_ylim([ymin,ymax])

	ax.set_xlabel( <[XLB]>, fontproperties=font, size=14 )
	ax.set_ylabel( <[YLB]>, fontproperties=font, size=14 )
	ax.set_title( ttl[0], fontproperties=font, size=17, verticalalignment="bottom" )

	fig.savefig( "<[OUT]>", facecolor="white" )
	# fig.savefig( "<[OUT]>", facecolor="white", bbox_extra_artists=(lgd,), bbox_inches="tight")

elif dim == 2:

	x = bins[0]
	y = bins[1]
	intensity = wght[0]
	cmap = plt.get_cmap('RdYlGn')
	(x,y) = np.meshgrid(x, y)
	c= np.array(intensity)
	im = ax.pcolormesh( x, y, c, vmin=ymin, vmax=ymax, cmap=cmap, facecolor="black", edgecolor="black" )

	cbar = plt.colorbar( im )
	cbar.ax.get_yaxis().labelpad = 35
	cbar.ax.set_ylabel( ttl[1], rotation=270, size=16 )

	for (_,i) in ax.spines.items(): i.set( linewidth=0.8, color="black", alpha=1.0, zorder=0 )
	for i in (ax.get_xticklabels() + ax.get_yticklabels()): i.set_fontproperties(font)
	ax.tick_params( axis="both", which="major", zorder=0 ); ax.minorticks_on()

	for i in (cbar.ax.get_yticklabels()):
		i.set_fontproperties(font)
		i.set_horizontalalignment("right")
		i.set_x(2.5)

	x0,x1 = ax.get_xlim()
	y0,y1 = ax.get_ylim()
	ax.set_aspect((x1-x0)/(y1-y0))

	ax.set_xlabel( <[XLB]>, fontproperties=font, size=16 )
	ax.set_ylabel( <[YLB]>, fontproperties=font, size=16 )
	ax.set_title( ttl[0], fontproperties=font, size=19, verticalalignment="bottom" )

	fig.savefig( "<[OUT]>", facecolor="white" )

sys.exit(0)

