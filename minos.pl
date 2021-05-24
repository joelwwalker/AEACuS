#!/usr/bin/perl

#*******************************#
# minos.pl Version 0.1 (BETA)	#
# August '20 - May '21		#
# Joel W. Walker		#
# Sam Houston State University	#
# jwalker@shsu.edu		#
# arXiv:XXXX.XXXX		#
# Copy: GNU Public License V3	#
#*******************************#

# Apply a strict coding pragma and locate filesystem
use strict; use sort q(stable); use FindBin qw($Bin); use lib qq($Bin);

# Import AEACuS subroutine library and perform version compatibility check
BEGIN { require q(aeacus.pl); ( &UNIVERSAL::VERSION( q(Local::AEACuS), 3.032 )); }

# Read event plotting specifications from cardfile
our ($OPT); my ($MIN) = map { (/^(.*\/)?([^\/]*?)(?:\.dat)?$/); my ($crd,$err,$fil) =
	( &LOAD_CARD([ (($1) or ( q(./Cards/))), [ ((defined) ? ($2) : ( q(min_card))), q(), q(.dat) ]]));
	($crd) or ((defined) ? ( die 'Cannot open card file for read' ) : (exit));
	( die ( join "\n", ( 'Malformed instruction(s) in card '.($$fil[0].$$fil[1]),
		( map {( "\t".'* Line '.$$_[0].':'."\t".'>> '.$$_[1].' <<' )} (@$err)), q()))) if (@$err);
	@$crd{( qw( min ))}} ( &$OPT( q(crd)));

# Establish whether Python 2.7/3.X and MatPlotLib 1.3.0+ are suitably configured for piped system calls
my ($cxg) = ( &CAN_XGBOOST());

# Perform machine learning training
{; my ($def) = ( ${$$MIN{trn}||[]}[0] || {} ); my ($cdf) = ( ${$$MIN{chn}||[]}[0] || {} );

	TRN: for my $i (1..(@{$$MIN{trn}||[]}-1)) {

	my ($trn) = (( $$MIN{trn}[$i] ) || (next TRN)); do {(( exists $$trn{$_} ) or ( $$trn{$_} = $$def{$_} ))} for ( keys %{$def} );

	my ($ipb,$fix) = do { my ($lum,$ipb,$fix) = [ @$trn{( qw( ipb ifb iab izb iyb ))} ]; for my $i (0..4) {
		($ipb,$fix) = @{(($$lum[$i]) or (next))}; $ipb *= (10)**(3*$i); (last) } ( map {( $_, ((defined) && ($fix < 0)))} ($ipb)) };
#	my ($sum,$nrm,$per,$avg) = map {[ ((@{$_||[]}) ? ( map {((defined) ? (0+ $_) : (undef))} (@$_)) : (undef)) ]} ( @$trn{( qw( sum nrm per avg ))} );
	my ($out,$nam,$fmt,$key) = ( @$trn{( qw( out nam fmt ))} );

	$out = (( &Local::FILE::PATH( [ ( $out = q().( &DEFINED( $$out[0], q(./Models/)))), ( sprintf "TRN_%3.3i", $i ) ], 2 )) or ( die 'Cannot write to directory '.$out ));
	my (@vls) = do { my ($chn) = []; map { my ($cid,@set) = $_; CHN: {; (@set) =

		# Read and bin data files into channels, combining like samples by luminosity and discrete samples by cross section
		map { $$chn[$_] ||= do {

			my ($chn) = (( $$MIN{chn}[$_] ) or do { print STDERR 'CHANNEL '.$_.' IS NOT DEFINED'."\n"; (last CHN) } );
			do {(( exists $$chn{$_} ) or ( $$chn{$_} = $$cdf{$_} ))} for ( keys %{$cdf} );

			my ($dat,$inc,$exc,$ftr,$esc,$wgt,$lbl) = @{ $chn }{( qw( dat inc exc ftr esc wgt lbl ))};
			push @{$exc||=[]}, { wgt => (undef) };

			[

			map { my ($dir,$fil,%ipb,%FHT) = @{ $$MIN{dat}[$_] or
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

					my (@key) = do { my ($k,@key,%idx,%inc);
						for ( keys %$idx ) { ((/^([a-z][a-z\d]{2})_(\d{3})$/) or (next)); push @{$idx{( qq($1))}||=[]}, (0+ $2); }
						for (@{$inc||=[]}) { my ($k,$v) = (%$_); push @{$inc{$k}||=[]}, ((defined $v) ? ($v) : (@{$idx{$k}||[]})); } ((@$inc) or ((%inc) = (%idx)));
						for (@$exc) { my ($k,$v) = (%$_); if ( defined $v ) { $inc{$k} = [ grep {( $_ != $v )} @{$inc{$k}||[]} ] } else { delete $inc{$k}}}
						do {(( not defined $key ) ? ( $key = $_ ) : (( $key eq $_ ) or (
							do { print STDERR 'INCONSISTENT STATISTIC INDEXES IN TRAINING CYCLE '.$i."\n"; (last CHN) } )))} for ( join q(,), (
							map {( sprintf q('%3.3s_%3.3i'), ((ref eq 'ARRAY') ? ( q(FTR), (++$k)) : ( map {(uc)} (%$_))))}
							grep { ( push @key, (( &HASHED_FUNCTIONAL( $idx, ((ref eq 'ARRAY') ? (@$_) : ((undef),$_)))) or (
								do { print STDERR 'INVALID CHANNEL KEY SPECIFICATION IN TRAINING CYCLE '.$i."\n"; (last CHN) } ))) } ((
							map { my ($k) = $_; map { +{ $k => (0+ $_) }} sort { our ($a,$b); ( $a <=> $b ) } ( &UNIQUE( @{$inc{$k}||[]} )) }
							sort { our ($a,$b); ( $a cmp $b ) } ( keys %inc )), (@{$ftr||[]}))));
						((@key) or do { print STDERR 'EMPTY FEATURE KEY LIST IN TRAINING CYCLE '.$i."\n"; (last CHN) } ); (@key) };
					my (@cut) = grep {(( $$_[2] = ( &HASHED_FUNCTIONAL( $idx, ( map {((ref eq 'ARRAY') ? (@$_) : ((undef),$_))} ($$_[2][0]))))) or (
							do { print STDERR 'INVALID KEY IN SELECTION '.$$_[1].' FOR TRAINING CYCLE '.$i.' ON FILE '.$$fil[0].$$fil[1]."\n"; !1 } ))}
						grep {(! ( &MATCH_VALUE( $$_[3], undef )))} map {[ ($_ < 0), ( abs ), ( @{ (($_) && ( $$MIN{esc}[( abs )] )) or
							do { print STDERR 'INVALID EVENT SELECTION CUT SPECIFICATION IN TRAINING CYCLE '.$i."\n"; +{}}}{( qw( key cut ))} ) ]}
						map {((defined) ? ( int ) : ())} (@{$esc||[]});
					my ($wgt) = map {((defined) ? (( &HASHED_FUNCTIONAL( $idx, ((ref eq 'ARRAY') ? (@$_) :
							((undef),((ref eq 'HASH') ? ($_) : +{ wgt => (0+ $_) } ))))) or
						do { print STDERR 'INVALID CHANNEL WEIGHT SPECIFICATION IN TRAINING CYCLE '.$i."\n"; (last CHN) } ) :
						(( &HASHED_FUNCTIONAL( $idx, (undef), +{ wgt => 0 } )) or ( sub {($w)} )))} (${$wgt||[]}[0]);
					my ($nmx) = (( do { my ($f); ((($fix) && ( &MATCH_VALUE( [0,1], (($f) = (($w)*($ipb-$ipb{$tag})))))) ?
						(( &ROUND(($e)*($f))), ( $ipb{$tag} = $ipb )) : (($e), ( $ipb{$tag} += $l )))[0] } ) or (next));

					push @{ $FHT{$tag} ||= [] }, [ $l, ( grep {((defined) or ( die 'Cannot open temporary file for read/write' ))} ( &Local::FILE::HANDLE())) ];
						local ($_); while (<$FHI>) { ((/^\s*$/) and (next)); ((/^\s*(\d+)/) && ($1 <= $nmx)) or (last);
							do { my ($val) = $_; print (( join q(,), ( map {((defined) ? ($_) : ( q(NAN)))}
								((($wgt) -> ($val)), ( map {(($_) -> ($val))} (@key))))),( qq(\n))); } for
							grep { my ($val,$mch) = ($_,1); for (@cut) { my ($inv,$eid,$key,$cut) = (@$_);
								(( $mch = (($inv) xor ( &MATCH_VALUE( $cut, (($key) -> ($val)))))) or (last)) }; ($mch) }
							[ map { (/^UNDEF$/) ? (undef) : (0+ $_) } ( split ) ]; }}

				my ($FHO); do { my ($tag) = $_; my ($scl) = ( &RATIO(( &DEFINED($ipb,1)), $ipb{$tag} ));
					if ((defined $ipb) && (($scl < 0) or ($scl > 1))) { print STDERR 'RESCALING BY '.( sprintf '%+10.3E', ($scl)) .
						' TO TARGET LUMINOSITY OF '.( sprintf '%+10.3E', ($ipb)).' PER PB IN CHANNEL '.($tag)."\n"; }
					do { my ($l,$FHT) = @$_; (( seek $FHT, 0, SEEK_SET ) or ( die 'Cannot rewind temporary file' )); local ($_); while (<$FHT>) {
						my ($wgt,$lin) = split (( q(,)), $_, 2 );  (( $wgt *= ( $l * $scl )) or (next));
						$FHO ||= ( grep {((@$_) or ( die 'Cannot open file in directory '.$out.' for write' ))}
							[ ( &Local::FILE::NEXT([[ $out, q(CSV) ], [ (($$lbl[0]) ? q(POS) : q(NEG)), 0, q(csv) ]])) ] )[0];
						print +( $wgt, ( q(,)), $lin ); }} for (@{ $FHT{$tag}}) } for (keys %FHT);

				() }

			map { ( int ( &MAX(0,$_))) or do { print STDERR 'INVALID DATA SET SPECIFICATION'.$_."\n"; (last CHN) }} grep {(defined)} (@{$dat||[]}) ] }}

		map { ( int ( &MAX(0,$_))) or do { print STDERR 'INVALID CHANNEL SPECIFICATION'.$_."\n"; (last CHN) }} grep {((defined) && !(ref))} ($cid) } (@set) } (@{$$trn{chn}||[]}) };

	if (1) {
# ??? how to include/exclude? +/- or sep inputs? do target luminosity &  allow rescale/ absolute? how to file/merge ... memory? group by DAT or by CHN?
	do { print STDERR 'NO BINNED EVENTS ESTABLISHED FOR TRAINING CYCLE '.$i."\n"; (next TRN) } unless (@vls);
	my ($fpo) = $out . q(minos.py); my (%dat) = ();
	do { use Fcntl qw(:seek); local ($.,$?); my ($t,$FHO) = ( tell DATA );
		( $FHO = ( &Local::FILE::HANDLE($fpo,1))) or ( die 'Cannot write to file '.($fpo));
		local ($_); while (<DATA>) { s/<\[(\w+)]>/$dat{$1}/g; ( print $FHO $_ ) }
		( close $FHO ) && ( chmod 0755, $fpo ); ( seek DATA, $t, SEEK_SET );
		if ($cxg) { system( qq( cd ${out} && ./minos.py )) } else { print STDERR 'CANNOT VERIFY PYTHON 2.7/3.X WITH MATPLOTLIB 1.3.0+ AND XGBOOST (DELIVERING SCRIPT LITERAL)'."\n"; }}}}};

1

__DATA__
#!/usr/bin/env python3

import sys
if ((sys.version_info[0] < 2) or ((sys.version_info[0] == 2) and (sys.version_info[1] < 7))) :
        sys.exit( 'MInOS requires Python versions 2.7 or 3.X' )

import matplotlib as mpl
if (( tuple( map ( int, mpl.__version__.split("."))) + (0,0,0))[0:3] < (1,3,0)) :
        sys.exit( 'MInOS requires MatPlotLib version 1.3.0 or Greater' )

import warnings as wrn; wrn.filterwarnings("ignore")

import matplotlib.pyplot as plt

import math, os, glob

import numpy as np

import xgboost as xgb

mpl.rcParams['mathtext.fontset'] = 'cm'; mpl.rcParams['font.family'] = 'STIXGeneral'

def main () :
    lum = 300000 # events per pb
    xx = np.linspace( 0., 1., num=( 1 + 2500 ), endpoint=True )
    bgs = tuple( x for x in ( loadCSV( x ) for x in glob.glob( "./CSV/NEG_*.csv" )) if x is not None )
    sig = tuple( x for x in ( loadCSV( x ) for x in glob.glob( "./CSV/POS_*.csv" )) if x is not None )
    dmx_all = tuple( getDMatrix(( mergeCSV( bgs ), mergeCSV( sig )), 3, False ))[0]
    bdt_all = trainBDT( dmx_all )
    scr_all = scoreBDT( bdt_all, dmx_all )
    f = tuple( {
        "density":interpolateTrapezoid(
            *( pointsFromDensityEdge(
            *( binDensity( c["score"], c["weight"], bins=25, b2=500, smooth=4. ))))),
        "xsec":c["xsec"] } for c in scr_all )
    ff = tuple( tuple((( 1. - _["density"]( xx, -1 )), ( lum * _["xsec"] ))) for _ in f )
    if not os.path.exists( "./Plots/" ) : os.mkdir( "./Plots/" )
    plotDistribution( xx, f[0]["density"]( xx ), f[1]["density"]( xx ))
    plotROC ( ff[0][0], ff[1][0] )
    plotSignificance( xx, *( tuple( _[0] * _[1] for _ in ff )))
    plotImportance( bdt_all )
    return

# load CSV data record into numpy array
def loadCSV ( f ) :
    text = np.loadtxt( f, delimiter=",", dtype=np.float64 )
    if len( text.shape ) == 1 :
        if text.shape[0] == 0 : return None
        text = np.array([ text ])
    return text

# merge numpy arrays generated from CSV data records
def mergeCSV ( csv ) : return np.concatenate( tuple( csv ), axis=0 )

# convert CSV data records into DMatrix objects for training and testing
def getDMatrix ( csv, prt=2, fld=False ) :
    # shuffle csv data and split into prt partitions for each supervisory class
    prt = max( 2, int( prt )); csv = tuple( np.array_split(
        ( np.random.default_rng( seed=0 ).shuffle( x ) or x ), prt, axis=0 )
        for x in ( np.copy( x ) for x in csv ))
    # internal method for construction of dmx objects
    def dmx ( csv, prt, idx ) :
        # merge selected csv partitions and separate weights from features
        csv = tuple(( lambda x : ( x[1], x[0].flatten()))( np.split(
            np.concatenate( tuple( csv[i][j] for j in idx )), [1], axis=1 ))
            for i in range( len( csv )))
        # count residual events in each class
        events = tuple( len( x[0] ) for x in csv )
        # sum cross section of residual events in each class
        xsec = tuple( np.sum( x[1] ) for x in csv )
        # return data structure with merged dmatrix, event counts, and cross sections
        return {
            "dmatrix":xgb.DMatrix(
                data=np.concatenate( tuple( x[0] for x in csv ), axis=0 ),
                label=np.concatenate( tuple( np.full( x, i ) for i,x in enumerate( events ))),
                # normalize total cross section weight to unity for each class
                weight=np.concatenate( tuple(( x[1] / xsec[i] ) for i,x in enumerate( csv ))),
                feature_names=feature_names ),
            "events":np.array( events ),
            # scale up to represent the inclusive physical sample for each class
            "xsec":np.array( tuple(( x * prt / len( idx )) for x in xsec )) }
    # generator function yields tuple of (train,test) dmx objects for each fold
    for i in range( prt if fld else 1 ) : yield tuple( dmx( csv, prt, idx )
        for idx in ( tuple( j for j in range( prt ) if j != i ), ( i, )))

# Train a dmx object and generate a bdt object
def trainBDT ( dmx ) :
    ( train, test ) = ( x["dmatrix"] for x in dmx ); progress = dict()
    return xgb.train( param, train,
        num_boost_round=40, early_stopping_rounds=5,
	evals=[ ( train, "train" ), ( test, "test" ) ],
        verbose_eval=False, evals_result=progress )
# HERE print ( bst.best_score , bst.best_iteration )

# Score a dmx object against a trained bdt
def scoreBDT ( bdt, dmx, test=1 ) :
    ( dmatrix, events, xsec ) = ( dmx[test][x] for x in ( "dmatrix", "events", "xsec" ))
    ( score, weight ) = ( bdt.predict( dmatrix ), dmatrix.get_weight())
    return tuple((( lambda i,j,k : (( lambda i,x : {
        # generate sorted data structure with scores, weights, and cross section
        "score":score[x], "weight":weight[x], "xsec":xsec[i] } )(
        # sort class indices according to classification score
        i, ( j + np.argsort( score[j:k], kind="stable" )))))(
        # establish boundary indices of events in ith class
        i, events[0:i].sum(), events[0:i+1].sum()))
        # loop over each supervisory class
        for i in range( len( events )))

# Integrates histogram densities and edges
def binIntegrate ( d, e, leftward=False ) :
    d = np.asarray( d ); e = np.asarray( e ); a = d * ( e[1:] - e[:-1] )
    return np.flip( np.cumsum( np.flip( a ))) if leftward else np.cumsum( a )

# Bins samples and weights into densities and edges by edge specification
def binEdge ( s, w, bins=None ) :
    return np.histogram( s, weights=w, range=( 0., 1. ), density=True,
        bins=( 2 * binsSturges( len( s )) if bins is None else bins ))

# Bins samples and weights into densities and edges by density specification
def binDensity ( x, y, bins=None, b2=500, smooth=1. ) :
    ( p, s, a, e, w, n, t ) = ( 0., 0., 0., [ 0. ], [ 0. ], len( x ),
        ( 1. / ( 2 * binsSturges( len( x )) if bins is None else bins )))
    for i in range( n ) :
        p += y[i]
        if i < n - 1 and x[i+1] == x[i] : continue
        ( q, r ) = np.divmod(( s + p ), t )
        if i == n - 1 and r > t / 2. : q += 1.; r = 0.
        if q :
            if s :
                e.append(( a + ( t - s ) * x[i] ) / t )
                w.append( t ); q -= 1.
            if q :
                e.append( x[i] )
                w.append( q * t )
            s = r; a = r * x[i]
        else : s += p; a += p * x[i]
        p = 0.
    e.append( 1. ); e = np.array( e )
    w.append( 0. ); w = np.array( w )
    edg = (( w[1:] * e[:-1] + w[:-1] * e[1:] ) / ( w[:-1] + w[1:] ))
    wdt = ( edg[1:] - edg[:-1] ); e = e[1:-1]; w = w[1:-1]; w /= w.sum(); w /= wdt
    if b2 is None : return tuple(( w, edg ))
    smooth /= 2.; wdt *= smooth; w /= smooth
    edg = np.linspace( 0., 1., num=( 1 + b2 ), endpoint=True )
    cen = ( edg[1:] + edg[:-1] ) / 2.
    val = np.add.reduce( tuple ( w[i] * np.exp( np.square(( cen - e[i] ) / wdt[i] ) / ( -2. )) for i in range( len( e ))), axis=0 )
    val /= ( val.sum() / b2 )
    return tuple(( val, edg ))
# HERE : consider non-flipped association? or even sharing / other? see minos_break_ ...

def evaluatePolynomial ( p, x ) :
    v = 0.
    for c in p[::-1] :
        v = c + ( v * x )
    return v

def polynomialIntegral ( p, a=0. ) :
    p = np.asarray( p )
    return np.concatenate(( [ a ], ( p / ( 1. + np.arange( len( p ))))))

def polynomialDerivative ( p ) :
    p = np.asarray( p )
    return (( p * np.arange( len( p )))[1:] )

def pointsFromDensityEdge ( d, e ) :
    d = np.asarray( d ); e = np.asarray( e ); h = ( e[1:] - e[:-1] )
    return tuple(( e, np.concatenate(
        ( [ 0. ], (( d[:-1] * h[:-1] + d[1:] * h[1:] ) / ( h[:-1] + h[1:] )), [ 0. ] ))))

# Cubic Spline interpolation, after Burden and Faires
def interpolateSpline ( x, y, slope_left=None, slope_right=None ) :
    x = np.asarray( x ); y = np.asarray( y ); h = ( x[1:] - x[:-1] ); n = len( h )
    a = np.concatenate((
        [[ 0, 1, 0 ] if slope_left is None else [ 0, 2*h[0], h[0]]],
        [[ h[i], 2*( h[i] + h[i+1] ), h[i+1]] for i in range( n - 1 ) ],
        [[ 0, 1, 0 ] if slope_right is None else [ h[-1], 2*h[-1], 0 ]] ), axis=0 )
    b = np.concatenate((
        [ (( y[1] - y[0] ) / h[0] ) if slope_left is None else slope_left ],
        [ (( y[i+1] - y[i] ) / h[i] ) for i in range( n ) ],
        [ (( y[-1] - y[-2] ) / h[-1] ) if slope_right is None else slope_right ] ), axis=0 )
    ( l, u, z, c ) = ( np.zeros( n + 1 ) for _ in range( 4 ))
    for i in range( n + 1 ) : l[i] = a[i,1] - a[i,0] * u[i-1]; u[i] = a[i,2] / l[i]
    for i in range( n + 1 ) : z[i] = ( 3*( b[i+1] - b[i] ) - a[i,0] * z[i-1] ) / l[i]
    for i in range( n + 1 ) : c[-(i+1)] = z[-(i+1)] - u[-(i+1)] * c[-i]
    return polynomialAccessor( x, np.array(
        [[ y[i], b[i+1] - h[i] * ( c[i+1] + 2*c[i] ) / 3, c[i],
        ( c[i+1] - c[i] ) / ( 3*h[i] ) ] for i in range( n ) ] ))

# Linear Trapezoidal interpolation
def interpolateTrapezoid ( x, y ) :
    x = np.asarray( x ); y = np.asarray( y ); h = ( x[1:] - x[:-1] ); n = len( h )
    return polynomialAccessor( x, np.array(
        [[ y[i], (( y[i+1] - y[i] ) / h[i] ) ] for i in range( n ) ] ))

def polynomialAccessor ( e, p ) :
    e = np.asarray( e ); p = np.asarray( p ); ( n, o ) = np.shape( p )
    ( t, d ) = ( 0., ( p, np.empty(( n, ( o - 1 ))), np.empty(( n, ( o + 1 )))))
    for i in range( n ) :
        d[1][i] = polynomialDerivative( d[0][i] )
        d[-1][i] = polynomialIntegral( d[0][i], t )
        t = evaluatePolynomial( d[-1][i], ( e[i+1] - e[i] ))
    for _ in d : _ /= t
    return np.vectorize( lambda x, i=0 : (
        lambda j : evaluatePolynomial( d[i][j], ( x - e[j] )))(
        np.clip(( np.searchsorted( e, x ) - 1 ), 0, ( n - 1 ))))

def binsSturges ( n ) : return ( 1 + math.ceil( math.log( n , 2 )))

# Generate donut plot of feature importance to total gain
def plotImportance ( bdt, pth="./Plots/", idx=None ) :
    scores = bdt.get_score( importance_type="total_gain" )
    keys = sorted( scores, key=scores.get )
    values = np.array([ scores[k] for k in keys ]); values /= values.sum()
    ( t, l, d ) = ( 0.0, 0.0, [] )
    for k, v in zip( keys, values ) :
        p = d[-1][0] if d else 0.03
        if p >= l and t <= 1 - l : l = p; d.append( [ 0.0, []] )
        t += v; d[-1][0] += v; d[-1][1].append( k )
    fig = plt.figure( figsize=( 7.5, 5 ), tight_layout=False )
    ax = fig.add_axes([ 0., .075, 1., .75 ])
    labels = tuple((( lambda x : ", ".join( tuple( feature_dict[k] for k in x[:3] ) +
        (( r"$\ldots$", ) if ( len( x ) > 3 ) else ())))( tuple( reversed( r[1] )))) for r in reversed( d ))
    wedges, texts, autotexts = ax.pie( [ _[0] for _ in reversed( d ) ],
        normalize=True, radius=1.0, startangle=90.0, autopct='%1.1f%%', pctdistance=0.75,
        wedgeprops={ "width":0.5, "edgecolor":"black", "linewidth":0.8 }, colors=( plt.get_cmap("tab10")( range( 10 ))))
    plt.setp( autotexts, size=12, color="white" )
    kw = { "zorder": 0, "verticalalignment":"center", "arrowprops":{ "arrowstyle":"-", "linewidth":0.8 },
        "bbox":{ "boxstyle":"round,pad=0.3,rounding_size=0.3", "facecolor":"white", "edgecolor":"black", "linewidth":0.8 }}
    for i, p in enumerate(wedges) :
        ang = (( p.theta2 + p.theta1 ) / 2.0 ); y = np.sin( np.deg2rad( ang )); x = np.cos( np.deg2rad( ang ))
        horizontalalignment = { -1: "right", 1: "left"}[ int( np.sign(x)) ]
        connectionstyle = "angle,angleA=0,angleB={}".format(ang)
        kw["arrowprops"].update( { "connectionstyle":connectionstyle } )
        ax.annotate( labels[i], xy=( x, y ), size=12, xytext=( 1.25*np.sign(x), 1.25*y ),
            horizontalalignment=horizontalalignment, **kw )
    ax.set_title( "Feature Importance to Total Gain", size=17, verticalalignment="bottom", pad=20 )
    fig.savefig( pth + "donut_plot" + ( "" if idx is None else ( "_" + str(idx))) + ".pdf", facecolor="white" )

# Generate plot of Receiver Operating Characteristic evolution with bdt score
def plotROC ( xx, yy, pth="./Plots/" ) :
    auc = ((( yy[:-1] + yy[1:] ) * ( xx[1:] - xx[:-1] )).sum() / ( -2. ))
    fig = plt.figure( figsize=( 5., 5. ), tight_layout=False )
    ax = fig.add_axes([ 0.05, .12, 0.95, 0.75 ], aspect="equal" )
    color=( plt.get_cmap("tab10")( range( 1 )))
    ax.set_xlim([ 0.0, 1.0 ])
    ax.set_ylim([ 0.0, 1.0 ])
    ax.set_xlabel( r"False Positive Rate", size=14, color="black" )
    ax.set_ylabel( "True Positive Rate", size=14, color="black" )
    ax.xaxis.grid( True, linewidth=0.8, color="black" )
    ax.yaxis.grid( True, linewidth=0.8, color="black" )
    ax.text( 0.53, 0.05, "Area Under Curve: {:4.2f}".format( auc ), fontsize=12,
        bbox={ "boxstyle":"round,pad=0.3,rounding_size=0.3", "facecolor":"white", "edgecolor":"black", "linewidth":0.8 } )
    tkw = dict( size=4, width=0.8 )
    ax.tick_params( axis='x', **tkw )
    ax.plot( xx, yy, color=color[0], linewidth=1.4, linestyle="solid", zorder=1 )
    ax.plot( xx, xx, color="red", linewidth=1.4, linestyle="dashed", zorder=1 )
    ax.fill_between( xx, yy, 0.0, linewidth=0.0, alpha=0.1, edgecolor="none", facecolor=color[0], hatch="", zorder=-1 )
    ax.tick_params( axis='y', colors="black", **tkw )
    for x in ( ax.get_xticklabels() + ax.get_yticklabels()): x.set_fontsize( 12 )
    ax.tick_params( axis="y", which="both", zorder=0, color="black" ); ax.minorticks_on()
    for (_,i) in ax.spines.items(): i.set( linewidth=0.8, color="black", alpha=1.0, zorder=0 )
    fig.suptitle( r"Receiver Operating Characteristic", x=0.5, y=0.9, size=17, verticalalignment="bottom" )
    fig.savefig( pth + "roc_plot.pdf", facecolor="white" )

# Generate plot of bdt score distributions for signal and background
def plotDistribution ( xx, bb, ss, pth="./Plots/" ) :
    fig = plt.figure( figsize=( 7.5, 5. ), tight_layout=False )
    ax = fig.add_axes([ 0.1, .14, 0.8, 0.75 ] )
    color=( plt.get_cmap("tab10")( range( 2 )))
    ax.set_xlim([ 0.0, 1.0 ])
    label = ( r"Background", r"Signal" )
    ax.set_xlabel( r"Signal Classification Score", size=14, color="black" )
    ax.set_ylabel( "Probability Density", size=14, color="black" )
    ax.xaxis.grid( True, linewidth=0.8, color="black" )
    tkw = dict( size=4, width=0.8 )
    ax.tick_params( axis='x', **tkw )
    for i, dd in enumerate(( bb, ss )) : 
        ax.plot( xx, dd, color=color[i], linewidth=1.4, linestyle="solid", zorder=1 )
        ax.fill_between( xx, dd, 0.0, linewidth=0.0, alpha=0.1, edgecolor="none", facecolor=color[i], hatch="", zorder=-1 )
    patches = [(
        mpl.patches.Rectangle((0,0), 1, 1, fill=True, facecolor=color[i], alpha=0.1, linewidth=0.0 ),
        mpl.patches.Rectangle((0,0), 1, 1, fill=None, linestyle="solid", edgecolor=color[i], linewidth=1.4 ))
       for i in range(2) ]
    lgd = ax.legend( patches, label, loc="best", fontsize=12 )
    lgd.get_frame().set( facecolor="white", linewidth=0.8, edgecolor="black", alpha=1.0 ); lgd.set(zorder=2)
    ax.set_ylim([ 0.0, None ])
    ax.tick_params( axis='y', colors="black", **tkw )
    for x in ( ax.get_xticklabels() + ax.get_yticklabels()): x.set_fontsize( 12 )
    ax.tick_params( axis="y", which="both", zorder=0, color="black" ); ax.minorticks_on()
    for (_,i) in ax.spines.items(): i.set( linewidth=0.8, color="black", alpha=1.0, zorder=0 )
    fig.suptitle( r"Signal and Background Score Distribution", x=0.5, y=0.9, size=17, verticalalignment="bottom" )
    fig.savefig( pth + "distribution_plot.pdf", facecolor="white" )

# Generate plot of significance as a function of bdt score
def plotSignificance ( xx, bb, ss, pth="./Plots/" ) :
    fig = plt.figure( figsize=( 7.5, 5 ), tight_layout=False )
    ax = fig.add_axes([ 0.11, .12, 0.65, 0.75 ])
    twin1 = ax.twinx()
    twin2 = ax.twinx()
    twin2.spines["right"].set_position(( "axes", 1.2 ))
    yy = [ ss, ( ss / ( 1. + bb )), ( ss / np.sqrt( 1. + bb )) ]
    label = ( r"$S$", r"$S\div(1+B)$", r"$S\div\sqrt{(1+B)}$" ) 
    color=( plt.get_cmap("tab10")( range( 3 )))
    axes = ( ax, twin1, twin2 )
    side = ( "left", "right", "right" )
    ax.set_xlim([ 0.0, 1.0 ])
    ax.set_xlabel( r"Signal Classification Threshold", size=14 )
    ax.xaxis.grid( True, linewidth=0.8, color="black" )
    tkw = dict( size=4, width=0.8 )
    ax.tick_params( axis='x', **tkw )
    for i in range(3) :
        axes[i].plot( xx, yy[i], color=color[i], linewidth=1.4, linestyle="solid", label=label[i], zorder=i+1 )
        axes[i].fill_between( xx, yy[i], 0.0, linewidth=0.0, alpha=0.1, edgecolor="none", facecolor=color[i], hatch="", zorder=i-3 )
        axes[i].set_ylabel( label[i], size=14, color=color[i] )
        axes[i].set_ylim([ 0.0, None ])
        axes[i].tick_params( axis='y', colors=color[i], **tkw )
        for x in ( axes[i].get_xticklabels() + axes[i].get_yticklabels()): x.set_fontsize( 12 );
        axes[i].tick_params( axis="y", which="both", zorder=0, color=color[i] ); axes[i].minorticks_on()
        axes[i].spines[ side[i]].set( linewidth=1.4, color=color[i], alpha=1.0, zorder=9 )
    fig.suptitle( r"Signal vs. Background Significance", x=0.5, y=0.9, size=17, verticalalignment="bottom" )
    fig.savefig( pth + "sig_plot.pdf", facecolor="white" )

# NOTE: Various lists / parameters below are temporarily hard coded for BETA distribution

# Ordered list of training feature keys
feature_names = [
    "ATM_001","ATM_002","CTS_001","ETA_001","ETA_002","ETA_003","MAS_001","MAS_003",
    "MDP_001","MDP_002","MDP_003","MEF_000","MET_000","MHT_000","ODP_001","ODP_002","ODP_004",
    "PTM_001","PTM_002","PTM_003","TTM_001","VAR_001","VAR_002","VAR_004","VAR_011","VAR_012","VAR_013","VAR_021" ]

# Ordered list of training feature TeX symbols
feature_tex = [
    r"$M_{\rm T2}^{100}$", r"$M_{\rm T2}^0$", r"$\cos \theta^*$", r"$\eta_{\ell_1}$", r"$\eta_{\ell_2}$",
    r"$\eta_{j_1}$", r"$M_{\ell\ell}$", r"$M_{j}$", r"$\Delta \phi_{\ell_1 {/\!\!\!\!E}_{\rm T}}$", r"$\Delta \phi_{\ell_2 {/\!\!\!\!E}_{\rm T}}$",
    r"$\Delta \phi_{j_1 {/\!\!\!\!E}_{\rm T}}$", r"$M_{\rm eff}$", r"${/\!\!\!\!E}_{\rm T}$", r"$H_{\rm T}$", r"$\Delta \phi_{\ell_1 \ell_2}$",
    r"$\Delta \phi_{\ell_1 j_1}$", r"$\Delta \phi_{\ell_2 j_1}$", r"${P}_{\rm T}^{\ell_1}$", r"${P}_{\rm T}^{\ell_2}$", r"${P}_{\rm T}^{j_1}$",
    r"$M_{\tau\tau}$", r"$\tanh \vert \Delta \eta_{\ell_1 \ell_2} \vert$",
    r"$\tanh \vert \Delta \eta_{\ell_1 j_1} \vert$", r"$\tanh \vert \Delta \eta_{\ell_2 j_1} \vert$",
    r"${P}_{\rm T}^{\ell_1}\div {/\!\!\!\!E}_{\rm T}$", r"${P}_{\rm T}^{\ell_2}\div {/\!\!\!\!E}_{\rm T}$",
    r"${P}_{\rm T}^{j_1}\div {/\!\!\!\!E}_{\rm T}$", r"$(M_{\rm T2}^{100}-100) \div M_{\rm T2}^0$" ]

# Construct feature dictionary
feature_dict = { feature_names[i]:feature_tex[i] for i in range( len( feature_names )) }

# Configure hyperparameters
param = {
    "booster":"gbtree",
    "max_depth":8,
    "objective":"binary:logistic",
    "eval_metric":["auc","logloss"],
    "disable_default_eval_metric":1,
    "lambda":0.05,
    "alpha":0.,
    "gamma":0.,
    "eta":0.2,
    "base_score":0.5,
    "min_child_weight":0.05,
    "max_delta_step":0.,
    "tree_method":"auto",
    "scale_pos_weight":1.,
    "subsample":1.,
    "colsample_bytree":1.,
    "colsample_bylevel":1.,
    "colsample_bynode":1.,
    "verbosity":1,
    }

if __name__ == "__main__" : main()

sys.exit(0)
