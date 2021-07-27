#!/usr/bin/perl

#*******************************#
# minos.pl Version 0.3 (BETA)	#
# August '20 - July '21		#
# Joel W. Walker		#
# Sam Houston State University	#
# jwalker@shsu.edu		#
# arXiv:XXXX.XXXX		#
# Copy: GNU Public License V3	#
#*******************************#

# Apply a strict coding pragma and locate filesystem
use strict; use sort q(stable); use FindBin qw($Bin); use lib qq($Bin);

# Import AEACuS subroutine library and perform version compatibility check
BEGIN { require q(aeacus.pl); ( &UNIVERSAL::VERSION( q(Local::AEACuS), 3.033 )); }

# Read event plotting specifications from cardfile
our ($OPT); my ($crd) = map { (/^(.*\/)?([^\/]*?)(?:\.dat)?$/); my ($crd,$err,$fil) =
	( &LOAD_CARD([ (($1) or ( q(./Cards/))), [ ((defined) ? ($2) : ( q(min_card))), q(), q(.dat) ]]));
	($crd) or ((defined) ? ( die 'Cannot open card file for read' ) : (exit));
	( die ( join "\n", ( 'Malformed instruction(s) in card '.($$fil[0].$$fil[1]),
		( map {( "\t".'* Line '.$$_[0].':'."\t".'>> '.$$_[1].' <<' )} (@$err)), q()))) if (@$err);
	($crd) } ( &$OPT( q(crd)));

# Establish whether Python 2.7/3.X and MatPlotLib 1.3.0+ are suitably configured for piped system calls
my ($cxg) = ( &CAN_XGBOOST());

# Perform machine learning training
{; my ($def) = ( ${$$crd{trn}||[]}[0] || {} ); my ($cdf) = ( ${$$crd{chn}||[]}[0] || {} );

	TRN: for my $i (1..(@{$$crd{trn}||[]}-1)) {

	my ($trn) = (( $$crd{trn}[$i] ) || (next TRN)); do {(( exists $$trn{$_} ) or ( $$trn{$_} = $$def{$_} ))} for ( keys %{$def} );

	my ($ipb,$fix) = do { my ($lum,$ipb,$fix) = [ @$trn{( qw( ipb ifb iab izb iyb ))} ]; for my $i (0..4) {
		($ipb,$fix) = @{(($$lum[$i]) or (next))}; $ipb *= (10)**(3*$i); (last) } ( map {( $_, ((defined) && ($fix < 0)))} ($ipb)) };
#	my ($sum,$nrm,$per,$avg) = map {[ ((@{$_||[]}) ? ( map {((defined) ? (0+ $_) : (undef))} (@$_)) : (undef)) ]} ( @$trn{( qw( sum nrm per avg ))} );
	my ($out,$nam,$fmt,$key) = ( @$trn{( qw( out nam fmt ))} );

	$out = (( &Local::FILE::PATH( [ ( $out = q().( &DEFINED( $$out[0], q(./Models/)))), ( sprintf "TRN_%3.3i", $i ) ], 2 )) or ( die 'Cannot write to directory '.$out ));
	my (@vls) = do { my ($chn) = []; map { my ($cid,@set) = $_; CHN: {; (@set) =

		# Read and bin data files into channels, combining like samples by luminosity and discrete samples by cross section
		map { $$chn[$_] ||= do {

			my ($chn) = (( $$crd{chn}[$_] ) or do { print STDERR 'CHANNEL '.$_.' IS NOT DEFINED'."\n"; (last CHN) } );
			do {(( exists $$chn{$_} ) or ( $$chn{$_} = $$cdf{$_} ))} for ( keys %{$cdf} );

			my ($dat,$inc,$exc,$ftr,$esc,$wgt,$lbl) = @{ $chn }{( qw( dat inc exc ftr esc wgt lbl ))};
			push @{$exc||=[]}, { wgt => (undef) };

			[

			map { my ($dir,$fil,%ipb,%FHT) = @{ $$crd{dat}[$_] or
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
						grep {(! ( &MATCH_VALUE( $$_[3], undef )))} map {[ ($_ < 0), ( abs ), ( @{ (($_) && ( $$crd{esc}[( abs )] )) or
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
    sys.exit( "MInOS requires Python versions 2.7 or 3.X" )

import matplotlib as mpl
if (( tuple( map ( int, mpl.__version__.split("."))) + (0,0,0))[0:3] < (1,3,0)) :
    sys.exit( "MInOS requires MatPlotLib version 1.3.0 or Greater" )

import warnings; warnings.filterwarnings("ignore")

import math, os, glob

import numpy as np

import matplotlib.pyplot as plt

import xgboost as xgb

mpl.rcParams["mathtext.fontset"] = "cm"; mpl.rcParams["font.family"] = "STIXGeneral"

"""
cross-validate ... distinct backgrounds & diff signals
plot names and labels / diff versions distinct
truth histogram validation under interpolation
composite trainings ...
    need the mdl for each & their "own-dmx" xsec / score interpolations
    generate a .predict object composite ...
    apply it to the MERGED dmx & generate an interpolant
    do the merger auto/internal or manual/external?
    general tools can be bundled into safe/compact macros
card structure / logic / calling methodology
card options & parameters
one script vs many scripts
"""

def main () :
    # do merged training
    csv_all = tuple( tuple( CSV( csv ) for csv in glob.glob( cls )) for cls in ( "./CSV/NEG_*.csv", "./CSV/POS_*.csv" ))
    dmx_all = next( DMX( tuple( mergeCSV( *x ) for x in csv_all )))
    bdt_all = BDT( dmx_all )
    doPlots( bdt_all, lum=300000., idx=1 )
    # do split training with recombination
    bdt = tuple( BDT( next( DMX(( x, csv_all[1][0] )))) for x in csv_all[0] )
    bdt_merge = BDT( dmx_all, mergeMDL( *bdt ))
    doPlots( bdt_merge, lum=300000., idx=2 )
    return

# instantiate object with method for the weighted recombination of a training ensemble
class mergeMDL :
    def __init__ ( self, *bdt ) :
        self.predict = lambda dmx : (( lambda s :
            s * np.array( tuple(( lambda a,b : ( a / b ) if b else np.full( len( a ), ( 1 / ( len( a ) or 1 )), dtype=np.float64 ))( x, x.sum())
            for x in np.array( tuple( np.sum(( x["xsec"] * x["density"](s) for x in bdt[i]["density"] ), axis=0 )
            for i,s in enumerate( s )), dtype=np.float64, ndmin=2 ).transpose())).transpose())
            ( np.array( tuple( x["model"].predict( dmx ) for x in bdt ), dtype=np.float64, ndmin=2 ))).sum( axis=0 )

# generate BDT dictionary object with test dmatrix, model, and density
def BDT ( dmx, mdl=None ) :
    if mdl is None : mdl = MDL( dmx )
    return { "dmatrix":dmx, "model":mdl, "density":DNS( SCR( mdl, dmx )) }

# generate collection of plots corresponding to an input dmx object
def doPlots ( bdt, lum=1., idx=None ) :
    if not os.path.exists( "./Plots/" ) : os.mkdir( "./Plots/" )
    plotImportance( bdt["model"], idx=idx )
    plotDistribution( bdt["density"], idx=idx )
    plotROC ( bdt["density"], idx=idx )
    plotSignificance( bdt["density"], lum=lum, idx=idx )
    return

# load CSV data into tuple of numpy arrays (events) of arrays (weight & features)
def CSV ( *fil, prt=3 ) :
    # minimal partition is 2 for separation of training and test samples
    prt = max( 2, int( prt ))
    # numpy import utility generates uniform 2x2 matrices of floating point values to be joined across files
    csv = (( lambda x : np.concatenate( x, axis=0 ) if len(x) > 1 else x[0] if x else np.empty((0,1), dtype=np.float64 ))(
        tuple( filter(( lambda x : x.shape[0] ), ( np.loadtxt( x, delimiter=",", dtype=np.float64, ndmin=2 ) for x in fil )))))
    # fail out records with no features or fewer events than the partition
    if csv.shape[1] < 2 or csv.shape[0] < prt : return tuple()
    # shuffle event records in place using fixed seed
    np.random.default_rng( seed=0 ).shuffle( csv, axis=0 )
    # return tuple of event partitions
    return tuple( np.array_split( csv, prt, axis=0 ))

# merge sets of CSV data partitions into a single tuple of numpy arrays of arrays
def mergeCSV ( *csv ) : return tuple( np.concatenate( x, axis=0 ) if len(x) > 1 else x[0] for x in zip( *filter( None, csv )))

# yield dictionaries of dmx objects for each fold from input CSV partitions for each class
def DMX ( cls, fld=False ) :
    # count available partitions for each class
    prt = min( tuple( len( x ) for x in cls ) or (0,))
    # internal method for construction of dmx objects
    def dmx ( cls, prt, idx ) :
        # merge selected csv partitions and separate weights from features
        cls = tuple(( lambda x : ( x[1], x[0].flatten()))( np.split(
            np.concatenate( tuple( cls[i][j] for j in idx )), [1], axis=1 ))
            for i in range( len( cls )))
        # count residual events in each class
        events = tuple( len( x[0] ) for x in cls )
        # sum cross section of residual events in each class
        xsec = tuple( np.sum( x[1] ) for x in cls )
        # return data structure with merged dmatrix, event counts, and cross sections
        return {
            # instance of core XGBoost data matrix object
            "dmatrix":xgb.DMatrix(
                # merge event samples from each class into a list of feature lists
                data=np.concatenate( tuple( x[0] for x in cls ), axis=0 ),
                # generate a list with the class label for each merged event
                label=np.concatenate( tuple( np.full( x, i ) for i,x in enumerate( events ))),
                # normalize total cross section weight to unity for each class
                weight=np.concatenate( tuple(( x[1] / xsec[i] ) for i,x in enumerate( cls ))),
                # assign plain text names for each ordered training feature
                feature_names=feature_names ),
            # array with number of sequential event samples belonging to each class
            "events":np.array( events ),
            # array with scaled inclusive physical cross sections for each class
            "xsec":np.array( tuple(( x * prt / len( idx )) for x in xsec )) }
    # generator function yields dictionary of (train,test) dmx objects for each fold
    for i in range( prt if fld else prt and 1 ) : yield { key : dmx( cls, prt, idx )
        for key, idx in (( "train", tuple( j for j in range( prt ) if j != i )), ( "test", ( i, ))) }

# generate a mdl object via training and validation of dmx object pair
def MDL ( dmx, trs=40, stp=5 ) :
    prg = dict(); return xgb.train( param, dmx["train"]["dmatrix"],
        num_boost_round=trs, early_stopping_rounds=stp,
        evals=[ ( dmx["train"]["dmatrix"], "train" ), ( dmx["test"]["dmatrix"], "test" ) ],
        verbose_eval=False, evals_result=prg )
# HERE print ( bst.best_score , bst.best_iteration )

# generate score object by projecting member (default:test) of dmx object pair onto mdl object
def SCR ( mdl, dmx ) :
    ( dmatrix, events, xsec ) = ( dmx["test"][x] for x in ( "dmatrix", "events", "xsec" ))
    ( score, weight ) = ( mdl.predict( dmatrix ), dmatrix.get_weight())
    # return tuple with dictionary of sorted data for each supervisory class
    return tuple((( lambda i,j,k : (( lambda i,x : {
        # generate sorted data structure with scores, weights, and cross section
        "score":score[x], "weight":weight[x], "xsec":xsec[i] } )(
        # sort class indices according to classification score
        i, ( j + np.argsort( score[j:k], kind="stable" )))))(
        # establish boundary indices of events in ith class
        i, events[0:i].sum(), events[0:i+1].sum()))
        # perform outer loop over each supervisory class
        for i in range( len( events )))

# generate an interpolated score density object dns from a score object
def DNS ( scr, bins=None, b2=1000, smooth=4. ) : return tuple( {
    "density":interpolateTrapezoid(
        *( pointsFromDensityEdge(
        *( binDensity( x["score"], x["weight"], bins=bins, b2=b2, smooth=smooth ))))),
    "xsec":x["xsec"] } for x in scr )

# integrate histogram densities and edges
def binIntegrate ( d, e, left=False ) :
    d = np.asarray( d ); e = np.asarray( e ); a = d * ( e[1:] - e[:-1] )
    return np.flip( np.cumsum( np.flip( a ))) if left else np.cumsum( a )

# bin samples and weights into densities and edges by edge specification
def binEdge ( s, w, bins=None ) : return np.histogram(
    s, weights=w, range=( 0., 1. ), density=True,
    bins=( 2 * binsSturges( len( s )) if bins is None else bins ))

# bin samples and weights into densities and edges by density specification
def binDensity ( x, y, bins=None, b2=None, smooth=1. ) :
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
    for c in p[::-1] : v = c + ( v * x )
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

# perform Cubic Spline interpolation, after Burden and Faires
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

# perform Linear Trapezoidal interpolation
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
    return np.vectorize(( lambda x, i=0 : (
        lambda j : evaluatePolynomial( d[i][j], ( x - e[j] )))(
        np.clip(( np.searchsorted( e, x ) - 1 ), 0, ( n - 1 )))),
        otypes=[np.float64], excluded=[1] )

# implement the Sturges bin count model
def binsSturges ( n ) : return ( 1 + math.ceil( math.log( n , 2 )))

# generate donut plot of feature importance to total gain
def plotImportance ( mdl, pth="./Plots/", idx=None ) :
    if not isinstance( mdl, xgb.core.Booster ) : return
    # HERE ... setup mdl for multiple and composite print ...
    scores = mdl.get_score( importance_type="total_gain" )
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
        normalize=True, radius=1.0, startangle=90.0, autopct="%1.1f%%", pctdistance=0.75,
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
    fig.savefig( pth + "donut_plot" + indexString(idx) + ".pdf", facecolor="white" )

# generate plot of mdl score distributions for signal and background
def plotDistribution ( dns, pth="./Plots/", idx=None ) :
    xx = np.linspace( 0., 1., num=( 1 + 2500 ), endpoint=True )
    dd = tuple( _["density"]( xx ) for _ in dns[0:2] )
    fig = plt.figure( figsize=( 7.5, 5. ), tight_layout=False )
    ax = fig.add_axes([ 0.1, .14, 0.8, 0.75 ] )
    color=( plt.get_cmap("tab10")( range( 2 )))
    ax.set_xlim([ 0.0, 1.0 ])
    label = ( "Background", "Signal" )
    ax.set_xlabel( "Signal Classification Score", size=14, color="black" )
    ax.set_ylabel( "Probability Density", size=14, color="black" )
    ax.xaxis.grid( True, linewidth=0.8, color="black" )
    tkw = dict( size=4, width=0.8 )
    ax.tick_params( axis="x", **tkw )
    for i in range( 2 ) :
        ax.plot( xx, dd[i], color=color[i], linewidth=1.4, linestyle="solid", zorder=3+i )
        ax.fill_between( xx, dd[i], 0.0, linewidth=0.0, alpha=0.1, edgecolor="none", facecolor=color[i], hatch="", zorder=i-2 )
    patches = [(
        mpl.patches.Rectangle((0,0), 1, 1, fill=True, facecolor=color[i], alpha=0.1, linewidth=0.0 ),
        mpl.patches.Rectangle((0,0), 1, 1, fill=None, linestyle="solid", edgecolor=color[i], linewidth=1.4 )) for i in range( 2 ) ]
    lgd = ax.legend( patches, label, loc="best", fontsize=12 )
    lgd.get_frame().set( facecolor="white", linewidth=0.8, edgecolor="black", alpha=1.0 ); lgd.set(zorder=5)
    ax.set_ylim([ 0.0, None ])
    ax.tick_params( axis="y", colors="black", **tkw )
    for x in ( ax.get_xticklabels() + ax.get_yticklabels()) : x.set_fontsize( 12 )
    ax.tick_params( axis="y", which="both", zorder=0, color="black" ); ax.minorticks_on()
    for (_,i) in ax.spines.items(): i.set( linewidth=0.8, color="black", alpha=1.0, zorder=0 )
    fig.suptitle( "Signal and Background Score Distribution", x=0.5, y=0.9, size=17, verticalalignment="bottom" )
    fig.savefig( pth + "density_plot" + indexString(idx) + ".pdf", facecolor="white" )

# generate plot of Receiver Operating Characteristic evolution with mdl score
def plotROC ( dns, pth="./Plots/", idx=None ) :
    xx = np.linspace( 0., 1., num=( 1 + 2500 ), endpoint=True )
    ( fp, tp ) = tuple(( 1. - _["density"]( xx, -1 )) for _ in dns[0:2] )
    auc = ((( tp[:-1] + tp[1:] ) * ( fp[1:] - fp[:-1] )).sum() / ( -2. ))
    fig = plt.figure( figsize=( 5., 5. ), tight_layout=False )
    ax = fig.add_axes([ 0.05, .12, 0.95, 0.75 ], aspect="equal" )
    color=( plt.get_cmap("tab10")( range( 1 )))
    ax.set_xlim([ 0.0, 1.0 ]); ax.set_ylim([ 0.0, 1.0 ])
    ax.set_xlabel( "False Positive Rate", size=14, color="black" )
    ax.set_ylabel( "True Positive Rate", size=14, color="black" )
    ax.xaxis.grid( True, linewidth=0.8, color="black" )
    ax.yaxis.grid( True, linewidth=0.8, color="black" )
    ax.text( 0.53, 0.05, "Area Under Curve: {:4.2f}".format( auc ), fontsize=12, zorder=5,
        bbox={ "boxstyle":"round,pad=0.3,rounding_size=0.3", "facecolor":"white", "edgecolor":"black", "linewidth":0.8 } )
    tkw = dict( size=4, width=0.8 ); ax.tick_params( axis="x", **tkw )
    ax.plot( xx, xx, color="red", linewidth=1.4, linestyle="dashed", zorder=3 )
    ax.plot( fp, tp, color=color[0], linewidth=1.4, linestyle="solid", zorder=4 )
    ax.fill_between( fp, tp, 0.0, linewidth=0.0, alpha=0.1, edgecolor="none", facecolor=color[0], hatch="", zorder=-1 )
    ax.tick_params( axis="y", colors="black", **tkw )
    for x in ( ax.get_xticklabels() + ax.get_yticklabels()): x.set_fontsize( 12 )
    ax.tick_params( axis="y", which="both", zorder=0, color="black" ); ax.minorticks_on()
    for (_,i) in ax.spines.items(): i.set( linewidth=0.8, color="black", alpha=1.0, zorder=0 )
    fig.suptitle( "Receiver Operating Characteristic", x=0.5, y=0.9, size=17, verticalalignment="bottom" )
    fig.savefig( pth + "roc_plot" + indexString(idx) + ".pdf", facecolor="white" )

# generate plot of significance as a function of mdl score
def plotSignificance ( dns, lum=1.0, pth="./Plots/", idx=None ) :
    xx = np.linspace( 0., 1., num=( 1 + 2500 ), endpoint=True )
    ( bb, ss ) = tuple((( lum * _["xsec"] ) * ( 1. - _["density"]( xx, -1 ))) for _ in dns[0:2] )
    yy = ( ss, ( ss / ( 1. + bb )), ( ss / np.sqrt( 1. + bb )))
    fig = plt.figure( figsize=( 7.5, 5 ), tight_layout=False )
    ax = fig.add_axes([ 0.11, .12, 0.65, 0.75 ])
    ax.set_zorder( 0 ); ax.set_xlim([ 0.0, 1.0 ]); ax.set_ylim([ 0.0, 1.0 ])
    ax.set_xlabel( "Signal Classification Threshold", size=14 )
    ax.xaxis.grid( True, linewidth=0.8, color="black" )
    for x in ax.get_xticklabels() : x.set_fontsize( 12 )
    tkw = dict( size=4, width=0.8 ); ax.tick_params( axis="x", **tkw )
    ax.spines["left"].set_visible(False); ax.spines["right"].set_visible(False)
    ax.yaxis.set_visible(False); ax.set_facecolor("none")
    axes = tuple( ax.twinx() for _ in range( 9 ))
    axes[5].spines["right"].set_position(( "axes", 1.2 ))
    label = ( r"$S$", r"$S\div(1+B)$", r"$S\div\sqrt{(1+B)}$" )
    side = ( "left", "right", "right" )
    color=( plt.get_cmap("tab10")( range( 3 )))
    for i in range( 3 ) :
        axes[i].set_frame_on(False)
        axes[i].yaxis.set_visible(False)
        axes[i].set_zorder( i-3 )
        axes[i].fill_between( xx, yy[i], 0.0, linewidth=0.0, alpha=0.1, edgecolor="none", facecolor=color[i], hatch="" )
        axes[i].set_ylim([ 0.0, None ])
    for i in range( 3 ) :
        axes[3+i].set_zorder( 1+i )
        axes[3+i].set_ylim( axes[i].get_ylim())
        for sp in axes[3+i].spines.values() : sp.set_visible(False)
        axes[3+i].spines[ side[i]].set_visible(True)
        axes[3+i].spines[ side[i]].set( linewidth=0.8, color=color[i], alpha=1.0 )
        axes[3+i].yaxis.set_label_position( side[i] )
        axes[3+i].yaxis.set_ticks_position( side[i] )
        axes[3+i].set_ylabel( label[i], size=14, color=color[i] )
        axes[3+i].tick_params( axis="y", colors=color[i], **tkw )
        axes[3+i].tick_params( axis="y", which="both", color=color[i] ); axes[3+i].minorticks_on()
        for x in axes[3+i].get_yticklabels() : x.set_fontsize( 12 )
    for i in range( 3 ) :
        axes[6+i].set_frame_on(False)
        axes[6+i].yaxis.set_visible(False)
        axes[6+i].set_zorder( 4+i )
        axes[6+i].set_ylim( axes[i].get_ylim())
        axes[6+i].plot( xx, yy[i], color=color[i], linewidth=1.4, linestyle="solid" )
    fig.suptitle( "Signal vs. Background Significance", x=0.5, y=0.9, size=17, verticalalignment="bottom" )
    fig.savefig( pth + "significance_plot" + indexString(idx) + ".pdf", facecolor="white" )

# canonical three-digit numerical identifier string in range "_000" to "_999"
def indexString( idx=None ) : return "" if idx is None else "_{0:03d}".format( min( max( 0, int(idx)), 999 )) 

# NOTE: Various lists / parameters below are temporarily hard coded for BETA distribution

# ordered list of training feature keys
feature_names = [
    "ATM_001","ATM_002","CTS_001","ETA_001","ETA_002","ETA_003","MAS_001","MAS_003",
    "MDP_001","MDP_002","MDP_003","MEF_000","MET_000","MHT_000","ODP_001","ODP_002","ODP_004",
    "PTM_001","PTM_002","PTM_003","TTM_001","VAR_001","VAR_002","VAR_004","VAR_011","VAR_012","VAR_013","VAR_021" ]

# ordered list of training feature TeX symbols
feature_tex = [
    r"$M_{\rm T2}^{100}$", r"$M_{\rm T2}^0$", r"$\cos \theta^*$", r"$\eta_{\ell_1}$", r"$\eta_{\ell_2}$",
    r"$\eta_{j_1}$", r"$M_{\ell\ell}$", r"$M_{j}$", r"$\Delta \phi_{\ell_1 {/\!\!\!\!E}_{\rm T}}$", r"$\Delta \phi_{\ell_2 {/\!\!\!\!E}_{\rm T}}$",
    r"$\Delta \phi_{j_1 {/\!\!\!\!E}_{\rm T}}$", r"$M_{\rm eff}$", r"${/\!\!\!\!E}_{\rm T}$", r"$H_{\rm T}$", r"$\Delta \phi_{\ell_1 \ell_2}$",
    r"$\Delta \phi_{\ell_1 j_1}$", r"$\Delta \phi_{\ell_2 j_1}$", r"${P}_{\rm T}^{\ell_1}$", r"${P}_{\rm T}^{\ell_2}$", r"${P}_{\rm T}^{j_1}$",
    r"$M_{\tau\tau}$", r"$\tanh \vert \Delta \eta_{\ell_1 \ell_2} \vert$",
    r"$\tanh \vert \Delta \eta_{\ell_1 j_1} \vert$", r"$\tanh \vert \Delta \eta_{\ell_2 j_1} \vert$",
    r"${P}_{\rm T}^{\ell_1}\div {/\!\!\!\!E}_{\rm T}$", r"${P}_{\rm T}^{\ell_2}\div {/\!\!\!\!E}_{\rm T}$",
    r"${P}_{\rm T}^{j_1}\div {/\!\!\!\!E}_{\rm T}$", r"$(M_{\rm T2}^{100}-100) \div M_{\rm T2}^0$" ]

# construct feature dictionary
feature_dict = { feature_names[i]:feature_tex[i] for i in range( len( feature_names )) }

# configure hyperparameters
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

