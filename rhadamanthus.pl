#!/usr/bin/perl

#*******************************#
# rhadamanthus.pl Version 1.7	#
# September '14 - December '17	#
# Joel W. Walker		#
# Sam Houston State University	#
# jwalker@shsu.edu		#
# www.joelwalker.net		#
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
require q(aeacus.pl); ( &UNIVERSAL::VERSION(q(Local::AEACuS),3.26)); our ($OPT);

# Read event plotting specifications from cardfile
my ($PLT) = map { (/^(.*\/)?([^\/]*?)(?:\.dat)?$/); my ($crd,$err,$fil) =
	( &LOAD_CARD([ (($1) or ( q(./Cards/))), [ ((defined) ? ($2) : ( q(plt_card))), q(), q(.dat) ]]));
	($crd) or ((defined) ? ( die 'Cannot open card file for read' ) : (exit));
	( die ( join "\n", ( 'Malformed instruction(s) in card '.($$fil[0].$$fil[1]),
		( map {( "\t".'* Line '.$$_[0].':'."\t".'>> '.$$_[1].' <<' )} (@$err)), q()))) if (@$err);
	@$crd{( qw( plt ))}} ( &$OPT( q(crd)));

# Establish whether Python 2.6/2.7 and MatPlotLib 1.3.0+ are suitably configured for piped system calls
my ($cpl) = ( &CAN_MATPLOTLIB());

# Generate histograms
for my $dim (1..3) { my ($hky) = ((undef), qw( hst h2d h3d ))[$dim]; HIST: for my $i (1..(@{$$PLT{$hky}||[]}-1)) { my ($hst) = (( $$PLT{$hky}[$i] ) || (next HIST));
	my ($ipb,$ovr) = do { my ($lum,$ipb,$ovr) = [ @$hst{( qw( ipb ifb iab izb iyb ))} ]; for my $i (0..4) { ($ipb,$ovr) = @{(($$lum[$i]) or (next))};
		$ipb *= (10)**(3*$i); (last) } ( map {( $_, !((defined) && ($ovr < 0)))} (( &MAX(0,$ipb)) or (undef))) };
	my ($obj,$wdt,$edg,$cnt) = map {( $_, [ $_->WIDTHS() ], [ $_->EDGES() ], [ $_->CENTERS() ] )} grep {
		(( grep {($_ > 3)} (( &Local::TENSOR::OBJECT((undef), $_ )), (undef,undef))[1..$dim] ) == $dim ) or
		do { print STDERR 'INVALID BINNING SPECIFICATION IN HISTOGRAM '.$i."\n"; (next HIST) }}
do { my (@t) = map {((@{$_||[]} == 1) ? do { my ($t) = @$_; [ map {[$t]} (1..$dim) ] } : ( &SPANS($dim,$_)))} @$hst{( qw( lft rgt spn bns ))};
( &Local::HISTOGRAM::NEW( map {[ map {[ grep {(defined)} @{( shift @$_ )} ]} (@t) ]} (1..$dim))) };
# think carefully about how the old iteration over channels translates to the new iteration over dimensions
	my ($sum,$nrm,$per,$avg) = map { (@{$_||[]}) ? ($_) : [(undef)] } ( @$hst{( qw( sum nrm per avg ))} ); my ($out,$nam,$fmt) = ( @$hst{( qw( out nam fmt ))} );
	my ($ttl,$lbl,$lgd) = map {[ map { ((defined) && !(ref)) ? qq($_) : (undef) } (@{$_||[]}) ]} ( @$hst{( qw( ttl lbl lgd ))} );
	my (@vls) = do { my ($chn,$set) = []; map { my ($sub,@cid) = ((ref eq 'ARRAY') ? (@$_) : ((ref) or !(defined)) ? (undef) : ( sub {(shift)} , $_ ));
		map { my ($s,$n,$p,$a) = map {( 0+ $$_[(@$_-1)&&($set)] )} ($sum,$nrm,$per,$avg);

			map {( scalar (($_)->SMOOTH($a)))} # set up for SPANS ? do option to skip/return if all are null? 
				# this/others get -1 early options? full order specification via "layers"? 2nd arg for (geom,arith,harm)[+1,0,-1]?
				# speed up this / COMPOUND, etc., if no action ... correct the multi-parameters to span correctly ...

			map {(($dim == 1) ? do { ($s) ? ($_) : (( &MAX(0,$p)) or ($n != 0) or ($$wdt[0]))*($_/$$wdt[0]) } : ($_))}
				# figure out how to do multi-dim norms and scaling

			map {( scalar (($_)->SLICE( map {[[1,-2]]} (0..($dim-1)))))}

			map {(($dim == 1) ? do { my ($u) = ( shift @$_); ($n > 0) ? ($_)*(( &RATIO( $n, (($u) +
				(($s) ? ($$_[((0<=>$s)-1)/2]) : ( &SUM( map {($$_)} (($_)->ELEMENTS()))))))) or
				do { print STDERR 'INDETERMINATE NORMALIZATION FACTOR IN HISTOGRAM '.$i."\n"; 1 } ) : ($_) } :
				( scalar (($_)->SLICE( map {[[1,-1]]} (0..($dim-1))))))}   # combine strategy with $n < 0 ?

			map {(($dim == 1) ? ( scalar (($_)->MAP( sub {( &MAX(0,(shift)))} ))) : ($_))}

			grep { (( &ISA( 1, $_, q(Local::TENSOR))) && (++$set)) or
				do { print STDERR 'INVALID CHANNEL OBJECT GENERATED IN HISTOGRAM '.$i."\n"; !1 }} ( scalar (($sub) -> (

			map {( scalar &Local::TENSOR::SPLICE( sub {( &Local::TENSOR::COMPOUND((shift), $s ))}, $_, ( map {[[1,-1]]} (0..($dim-1)))))}

			map {(($n < 0) ? ($_)*(( &RATIO(( abs $n ), ( &SUM( map {($$_)} (($_)->ELEMENTS()))))) or # GENERALIZE across (separately?) TENSOR & single pass average ...

				do { print STDERR 'INDETERMINATE NORMALIZATION FACTOR IN HISTOGRAM '.$i."\n"; 1 } ) : ($_))}
			map {( scalar &Local::TENSOR::OBJECT((undef), ( &CLONE($_))))} (@$_)))) } # is clone necessary? or rebless? do elsewhere also?
		do { my (@set); CHN: {; my ($j); my (@chn) = grep { ((@$_ == 1) or ((( $j ||= @$_ ) == @$_ ) && ($j))) or
			do { print STDERR 'DATA SET MULTIPLICITY MISMATCH IN HISTOGRAM '.$i."\n"; (last CHN) }} map { $$chn[$_] ||= do { my ($dat,$key,$esc) =
				@{ $$PLT{chn}[$_] || do { print STDERR 'CHANNEL '.$_.' IS NOT DEFINED'."\n"; (last CHN) }}{( qw( dat key esc ))}; [ map {
				my ($dir,$fil,%bin,%ipb) = @{ $$PLT{dat}[$_] or do { print STDERR 'DATA SET '.$_.' IS NOT DEFINED'."\n"; (last CHN) }}{( qw( dir fil ))};
				my ($dir) = ( map { ((defined) ? qq($_) : q(./Cuts/)) } ($$dir[0])); FIL: for my $fil ( sort SORT_LIST_ALPHA ( values %{{
				map {(($$_[0].$$_[1]) => ($_))} grep {( $$_[1] =~ /^[\w-]+\.cut$/ )} map {( &Local::FILE::LIST([$dir,$_]))}
				grep {(defined)} (@{$fil||[]}) }} )) { ( my $crp = $$fil[1] ) =~ s/(?:_\d+)*\.cut$//; ( my ($FHI) = ( &Local::FILE::HANDLE($fil))) or
				do { print STDERR 'CANNOT READ FROM FILE '.$$fil[0].$$fil[1]."\n"; (last CHN) }; my ($xsc,$idx) = ( &AEACUS_XSEC($FHI));
					(defined $$xsc[0]) or do { print STDERR 'CANNOT ESTABLISH EVENT COUNT FOR FILE '.$$fil[0].$$fil[1]."\n"; (last CHN) };
					(defined $$xsc[1]) or do { print STDERR 'CANNOT ESTABLISH EVENT CROSS SECTION FOR FILE '.$$fil[0].$$fil[1]."\n"; (last CHN) };
					($$xsc[3] > 0) or do { print STDERR 'NO SURVIVING EVENTS IN FILE '.$$fil[0].$$fil[1]."\n"; (next FIL) };
					(defined $idx) or do { print STDERR 'CANNOT ESTABLISH STATISTICS INDEX FOR FILE '.$$fil[0].$$fil[1]."\n"; (last CHN) };
				my ($nmx) = ( map { (($ovr) or (($ipb-$ipb{$crp}) > $_ )) ? do { $ipb{$crp} += $_; ($$xsc[0]) } :
					(( &ROUND($$xsc[0]*($ipb-$ipb{$crp})/$_)), ($ipb{$crp} = $ipb))[0] } ( &RATIO($$xsc[0],$$xsc[1])))[0];
				my ($key) = [ map { ( &HASHED_FUNCTIONAL( $idx, ((ref eq 'ARRAY') ? (@$_) : (undef,$_)))) or
					do { print STDERR 'INVALID CHANNEL KEY SPECIFICATION IN HISTOGRAM '.$i."\n"; (last CHN) }} @$key[0..($dim-1)]];
				my (@cut) = grep { ( $$_[2] = ( &HASHED_FUNCTIONAL( $idx, ( map { (ref eq 'ARRAY') ? (@$_) : (undef,$_) } ($$_[2][0]))))) or
					do { print STDERR 'INVALID KEY IN SELECTION '.$$_[1].' FOR HISTOGRAM '.$i.' ON FILE '.$$fil[0].$$fil[1]."\n"; !1 }}
					grep { !( &MATCH_VALUE( $$_[3], undef )) } map {[ ($_ < 0), 0+( $_ = (int abs)), ( @{ (($_) && ( $$PLT{esc}[$_] )) or
					do { print STDERR 'INVALID EVENT SELECTION CUT SPECIFICATION IN HISTOGRAM '.$i."\n"; +{}}}{( qw( key cut ))} ) ]}
					grep {(defined)} (@{$esc||[]});
				local ($_); while (<$FHI>) { ((/^\s*(\d+)/) && ($1 <= $nmx)) or (last); do { my ($val) = $_; ( $bin{$crp} ||=
					( $obj->NEW()))->BIN([ map {(($_)->($val))} (@$key) ]) } for grep { my ($val,$mch) = ($_,1); do { my ($inv,$eid,$key,$cut) = (@$_);
						( $mch = (($inv) xor ( &MATCH_VALUE( $cut, ( $key->($val)))))) or (last); } for (@cut); ($mch) }
					[ map { (/^UNDEF$/) ? (undef) : (0+ $_) } ( split ) ]; }}
				( scalar &Local::TENSOR::SUM( $obj, ( map { my ($crp) = $_; $bin{$crp}*( grep {
					( print STDERR 'RESCALING BY '.( sprintf '%9.3e', ($_)).' TO TARGET LUMINOSITY OF ' .
					( sprintf '%9.3e', ($ipb)).' PER PB IN CHANNEL '.($crp)."\n" ) if ((defined $ipb) && ($_ > 1)); 1 }
					( &RATIO((($ipb)||(1000)),$ipb{$_})))[0] } (keys %bin)))) }
				map { ( int ( &MAX(0,$_))) or do { print STDERR 'INVALID DATA SET SPECIFICATION'.$_."\n"; (last CHN) }} grep {(defined)} (@{$dat||[]}) ] }}
			map { ( int ( &MAX(0,$_))) or do { print STDERR 'INVALID CHANNEL SPECIFICATION'.$_."\n"; (last CHN) }} grep {(defined)} (@cid);
			(@set) = map { my ($i) = $_; [ map {( $$_[(@$_-1)&&($i)] )} (@chn) ] } (0..(($j)&&($j-1))); } (@set) }} (@{$$hst{chn}||[]}) };
	my ($log,$min,$fpo); my (%dat) = (
		DIM	=> $dim,
		VAL	=> '['."\n".( join ','."\n", map {( qq($_))} map {(($dim == 1) ? ($_) : ( scalar &Local::MATRIX::TRANSPOSE($_)))} (@vls)).']', # 2-d is x-y swapped ... generalize for n-d
		BIN	=> (($dim == 1) ? ( scalar &Local::TENSOR::STRING($$edg[0])) : ( '['.( join ',', map {( &Local::TENSOR::STRING($_))} (@$edg[0..($dim-1)])).']' )),
		LOG	=> (( $log = ($$hst{log}[0] > 0)) ? q(True) : q(False)),
		STK	=> (($$hst{stk}[0] > 0) ? q(True) : q(False)),
		MIN	=> ( map { ((defined) && (($log) ? ($_ > 0) : (($dim > 1) or ($_ >= 0)))) ? ( $min = $_ ) : q(None) } ($$hst{min}[0]))[0],
		MAX	=> ( map { ((defined) && ((defined $min) ? ($_ > $min) : ($log) ? ($_ > 0) : (($dim > 1) or ($_ >= 0)))) ? ($_) : q(None) } ($$hst{max}[0]))[0],
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
				((undef,$fpo) =	map {( join q(), (@$_))} ((($$fmt[1] > 0) or ( !($cpl) &&
					do { print STDERR 'CANNOT VERIFY PYTHON 2.6/2.7/3.X WITH MATPLOTLIB 1.3.0+ (DELIVERING SCRIPT LITERAL)'."\n"; 1 } )) ?
					([q(./),($f),($e)],[($d),($f),q(.py)]) : ([($d),($f),($e)])))[0] },
	); do { use Fcntl qw(:seek); local ($.,$?); my ($t,$FHO) = ( tell DATA ); ( defined $fpo ) ?
		do { ( $FHO = ( &Local::FILE::HANDLE($fpo,1))) or ( die 'Cannot write to file '.($fpo)) } :
		do { ( open $FHO, q(|-), q(python 2>&1)) or ( die 'Cannot open pipe to Python' ) };
		local ($_); while (<DATA>) { s/<\[(\w+)]>/$dat{$1}/g; ( print $FHO $_ ) }
		( close $FHO ) && (($? >> 8) == 0 ) && ( defined $fpo ) && ( chmod 0755, $fpo ); ( seek DATA, $t, SEEK_SET ) }}}

1

# allow some command line, e.g. $ipb overrides?
# unify setup and configuration with aeacus.pl ... can call rhad from aeacus
# extend to function plotting, best fit spline, overlay with PLT_SHW_001
# add functions to CHN = DAT:{} ? NO; multiple level stacking -- two types pre/post of channel merging? -- use for multiplots!?
# consider zeroth channel, where undef is implied => wildcard ... invoked by the absence of a channel designation; numbered channels require "*" to wildcard ?
# intercept zero content errors propagated to MPL & log => 1 (empty)
# automatic axis labeling, bin configuration & scaling
# log binning
# default L is the "natural" one ... or the one from aeacus .cut ? ... second entry of $$xsc
# if channel is skipped bc of failure readig event count, etx. don't generate plot ... currently generates empty plot
# heat maps, exclusion plots , triangle heat maps ... convolve two 1-d histograms VS each cell counted.  S/B w/ wout compounding ... optimization, descent, correlation.
# can do heat map of correlation ... how to scale ...
# import the smuon upgrades & generalize ...
# are the errors, eg. indet norm firing in right sequence ... allow no event count &/or no xsec? .. for external applications
# use HST_000 to set global parameters
# assume chn matches hst # if empty? autoplot hst from chn if 000 defaults?
# think carefully through logic of pre-norm (-1) and confirm good idea & implementation
# for multi-dim with additional plots .. output them too, in a grid or in sequence, numbered
# settle new dashed / ext legend into card format
# allow for external data with no xsec, etc ... just warn
# font size and family; lines dashing / multiple types / thickness / bg line
# legend float right
# finish 2D norm, etc. & color schemes

__DATA__
#!/usr/bin/env python

import sys
if ((sys.version_info[0] < 2) or ((sys.version_info[0] == 2) and (sys.version_info[1] < 6))) :
	sys.exit( 'RHADAManTHUS requires Python versions 2.6, 2.7, or 3.X' )

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

import numpy as np # wanted to avoid this

color = <[CLR]>
clist = set( k.lower() for k in  mpl.colors.cnames )
clrs = ((),
	( "#CE0000", "#1500AA", "#107C00", "#EF6C00", "#9500AF", "#00B2DB", "#3DD100", "#FFCE00" ),
	( "#e41a1c", "#377eb8", "#4daf4a", "#ff7f00", "#984ea3", "#ffff33", "#a65628", "#f781bf" ),
	( "#1f78b4", "#33a02c", "#e31a1c", "#ff7f00", "#a6cee3", "#b2df8a", "#fb9a99", "#fdbf6f" ),
	( "#d53e4f", "#f46d43", "#fdae61", "#fee08b", "#e6f598", "#abdda4", "#66c2a5", "#3288bd" ),
	( "#b10026", "#e31a1c", "#fc4e2a", "#fd8d3c", "#feb24c", "#fed976", "#ffeda0", "#ffffcc" ),
	( "#004529", "#006837", "#238443", "#41ab5d", "#78c679", "#addd8e", "#d9f0a3", "#f7fcb9" ),
	( "#084081", "#0868ac", "#2b8cbe", "#4eb3d3", "#7bccc4", "#a8ddb5", "#ccebc5", "#e0f3db" ),
	( "#4d004b", "#810f7c", "#88419d", "#8c6bb1", "#8c96c6", "#9ebcda", "#bfd3e6", "#e0ecf4" ),
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

bins = <[BIN]>

(ymin, ymax) = (<[MIN]>, <[MAX]>)

ttl = <[TTL]>

if dim == 1:

	ymin = max( 0.0 if ymin is None else float(ymin), 0.0 )
	if log and ymin == 0.0:
		ymin = math.pow( 10, math.floor( math.log10( min( [ k for j in wght for k in j if k > 0.0 ] or [1.0] ))))
	if ymax is not None:
		ymax = float(ymax)
		if ymax <= ymin: ymax = None

	edges = [ k for j in zip(bins[:-1],bins[1:]) for k in j ]
	widths = [ (j[1]-j[0]) for j in zip(bins[:-1],bins[1:]) ]
	value = empty = [0]*2*(len(bins)-1); lower = [ymin/2.0]*2*(len(bins)-1)

	for i in range(n):
		value = [ sum(l) for l in zip ([ k for j in zip(wght[i],wght[i]) for k in j ],( value if stack else empty )) ]
		upper = [ j if j>=ymin else ymin/2.0 for j in value ]
		if fill[i%f]: ax.fill_between(edges, upper, lower, linewidth=0.0, alpha=0.6, edgecolor="none", facecolor=color[i%c], hatch="", zorder=i-3*n )
		if h: ax.fill_between(edges, upper, lower, linewidth=0.0, alpha=1.0, edgecolor="0.25", facecolor="none", hatch=hatch[i%h], zorder=i-2*n )
		if bold[i%b]: ax.plot( edges, upper, color="#dddddd", linewidth=2.0, linestyle="solid", zorder=i-1*n )
		ax.plot( edges, upper, color=color[i%c], linewidth=(2.0,2.0)[bold[i%b]], linestyle=("solid","--")[bold[i%b]], zorder=i-1*n )
		if stack: lower = upper

	lgd = None
	lgnd = <[LGD]>
	if lgnd:
		patches = [(
			( mpl.patches.Rectangle((0,0), 1, 1, fill=True, facecolor=color[i%c], alpha=0.6, linewidth=0.0 ) if fill[i%f] else ()),
			( mpl.patches.Rectangle((0,0), 1, 1, fill=None, color="0.25", hatch=hatch[i%h], linewidth=0.0 ) if h else ()),
			( mpl.patches.Rectangle((0,0), 1, 1, fill=None, linestyle="solid", edgecolor="#dddddd", linewidth=1.4 ) if bold[i%b] else ()),
			mpl.patches.Rectangle((0,0), 1, 1, fill=None, linestyle=("solid","--")[bold[i%b]], edgecolor=color[i%c], linewidth=1.4 ))
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
# change if ymin < 0 ...

	x0,x1 = ax.get_xlim()
	y0,y1 = ax.get_ylim()
	ax.set_aspect((x1-x0)/(y1-y0))

	ax.set_xlabel( <[XLB]>, fontproperties=font, size=16 )
	ax.set_ylabel( <[YLB]>, fontproperties=font, size=16 )
	ax.set_title( ttl[0], fontproperties=font, size=19, verticalalignment="bottom" )

	fig.savefig( "<[OUT]>", facecolor="white" )

sys.exit(0)

