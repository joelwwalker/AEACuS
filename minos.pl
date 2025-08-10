#!/usr/bin/perl

#*******************************#
# minos.pl Version 0.7 ALPHA	#
# August 2020 - January 2023	#
# Joel W. Walker		#
# Sam Houston State University	#
# jwalker@shsu.edu		#
# arXiv:XXXX.XXXX		#
# Copy: GNU Public License V3	#
#*******************************#

# Apply a strict coding pragma and locate filesystem
use strict; use sort q(stable); use FindBin qw($Bin); use lib qq($Bin);

# Import AEACuS subroutine library and perform version compatibility check
BEGIN { require q(aeacus.pl); ( &UNIVERSAL::VERSION( q(Local::AEACuS), 4.000 )); }

# Read event machine learning specifications from cardfile
our ($OPT); my ($crd) = ( map { my ($crd,$err,$fil) = ( &LOAD_CARD(
	map {[ $$_[0], [ (( $$_[1] =~ /^(.*?)(?:\.dat)?$/ ) && qq($1)), q(), q(.dat) ]]}
		( scalar &Local::FILE::SPLIT( $_, q(./Cards/), q(min_card)))));
	($crd) or ((length) ? ( die 'Cannot open card file for read' ) : (exit));
	( die ( join "\n", ( 'Malformed instruction(s) in card '.($$fil[0].$$fil[1]), ( map {((ref) ?
		( "\t".'* Line '.$$_[0].':'."\t".'>> '.( grep { s/^\s+//; s/\s+/ /g; 1 } ($$_[1]))[0].' <<' ) :
		( "\t".'* Duplicative use of shelf '.$_ ))} (@$err)), q()))) if (@$err); ($crd) } ( &$OPT( q(crd))));

# Perform machine learning training
{; my ($def) = ( ${$$crd{bdt}||$$crd{trn}||[]}[0] || {} ); my ($cdf) = ( ${$$crd{chn}||[]}[0] || {} );

	BDT: for my $i (1..(@{$$crd{bdt}||$$crd{trn}||[]}-1)) {

	my ($bdt) = (( ${$$crd{bdt}||$$crd{trn}}[$i] ) || (next BDT)); do {(( exists $$bdt{$_} ) or ( $$bdt{$_} = $$def{$_} ))} for ( keys %{$def} );

	my ($ipb,$fix) = do { my ($lum,$ipb,$fix) = [ @$bdt{( qw( ipb ifb iab izb iyb ))} ]; for my $i (0..4) {
		($ipb,$fix) = @{(($$lum[$i]) or (next))}; $ipb *= (10)**(3*$i); (last) } ( map {( $_, ((defined) && ($fix < 0)))} ($ipb)) };
	# my ($sum,$nrm,$per,$avg) = map {[ ((@{$_||[]}) ? ( map {((defined) ? (0+ $_) : (undef))} (@$_)) : (undef)) ]} ( @$bdt{( qw( sum nrm per avg ))} );
	my ($key,$inc,$exc,$tex,$out) = ( @$bdt{( qw( key inc exc tex out ))} ); my ($pyt,$key_str,$tex_str,@lbl_str,%tex);
	for (0..(( int ( @{$tex||[]} / 2 )) - 1 )) { my ($k,$t) = ( $$tex[2*$_], $$tex[1+2*$_] );
		(( ref $k eq q(HASH)) or (next)); my ($k,$v) = ( @{(( &PAIR_KEY_IDX( $k ))||[])} ); $tex{( uc $k )}[$v] = ( &RAW_STRING( $t )) }
	$out = (( &Local::FILE::PATH( [ ( $out = q().( &DEFINED(( map {((length) ? qq($_) : ())} ($$out[0])),
		q(./Models/)))) ], 2 )) or ( die 'Cannot write to directory '.$out ));
	(( &Local::FILE::PATH(( $out = [ $out, ( uc sprintf "BDT_%3.3u", $i ) ] ), 0 )) and do { print STDERR 'TRAINING PATH ' .
		( join q(), ( @$out )).' ALREADY EXISTS ... DELETE PRIOR TO REGENERATION'."\n"; (next BDT) });
	$out = (( &Local::FILE::PATH( $out, 2 )) or ( die 'Cannot write to directory '.( join q(), ( @$out ))));
        print STDOUT 'Generating Analysis for '.( $out )."\n";
	# Read and bin data files into channels, combining like samples by luminosity and discrete samples by cross section
	CHN: for my $cdx (0..(@{$$bdt{chn}||[]}-1)) { my ($chn) = (
		map {(( $$crd{chn}[$_] ) or do { print STDERR 'CHANNEL '.$_.' IS NOT DEFINED'."\n"; (last CHN) } )}
		map { ( int ( &MAX(0,$_))) or do { print STDERR 'INVALID CHANNEL SPECIFICATION '.$_."\n"; (last CHN) }}
		grep {(((defined) && !(ref)) or (next CHN))} ( ${$$bdt{chn}||[]}[$cdx] ));
		do {(( exists $$chn{$_} ) or ( $$chn{$_} = $$cdf{$_} ))} for ( keys %{$cdf} );
		my ($dat,$esc,$wgt,$lbl) = @{ $chn }{( qw( dat esc wgt lbl ))}; push @lbl_str, [];
		$lbl = (( q(X)) . ( sprintf  q(%02.2u), ( &BOUNDED( [0,99], ( int ( &DEFINED( ${$lbl||[]}[0], $cdx )))))));

		DAT: for (@{$dat||[]}) { my ($dir,$fil,$tex,%ipb,%FHT) = (
			map {( @{ $$crd{dat}[$_] or do { print STDERR 'DATA SET '.$_.' IS NOT DEFINED'."\n"; (last CHN) }}{( qw( dir fil tex ))} )}
			(((defined) && !(ref)) ? (( int ( &MAX(0,$_))) or do { print STDERR 'INVALID DATA SET SPECIFICATION '.$_."\n"; (last CHN) } ) : (next DAT)));
			($dir) = (( map {((length) ? qq($_) : ())} (${$dir||[]}[0])), q(./Cuts/));
			($tex) = (( map {((length) ? ( &RAW_STRING( $_ )) : ())} (${$tex||[]}[0])), q(None));

			FIL: for my $fil ( sort SORT_LIST_ALPHA ( values %{{ map {((( &Local::FILE::DEVICE_INODE( $_ )) or
					( die 'Invalid Device/Inode for file '.$$_[0].$$_[1] )) => ($_))}
				grep {( $$_[1] =~ /^[\w-]+\.cut$/ )} map {( &Local::FILE::LIST( @$_[0,1], 0 ))}
				grep {(( length $$_[1] ) or ( die 'Invalid file name specification in card file' ))}
				map {( scalar &Local::FILE::SPLIT( $_, $dir ))} grep {(length)} (@{$fil||[]}) }} )) {

				( my $tag = $$fil[1] ) =~ s/(?:_\d+)*\.cut$//; ( my ($FHI) = ( &Local::FILE::HANDLE($fil))) or
					( do { print STDERR 'CANNOT READ FROM FILE '.$$fil[0].$$fil[1]."\n"; (last CHN) } );
				my (undef,$nnn,undef,undef,$idx) = ( &IMPORT_HEADER($FHI)); my ($pre,$pst) = ( map {[ grep { (defined) or
					do { print STDERR 'CANNOT ESTABLISH EVENT COUNT AND CROSS SECTION FOR FILE '.$$fil[0].$$fil[1]."\n"; (last CHN) }}
					( @{${$nnn||[]}[$_]||{}}{( qw( epw enw ezw abs ))} ) ]} (0,-1));
				my ($l,$w) = (( &RATIO(( &SUM( @$pre[0,1] )), $$pre[3] )), ( &RATIO( $$pre[3], ( &SUM( @$pre[0,1] )))));
				my ($nmx) = (( do { my ($f); ((( $fix ) && ( &MATCH_VALUE( [0,1], ( $f = (( $w ) * ( $ipb - $ipb{$tag} )))))) ?
					(( &ROUND(( &SUM( @$pre[0..2] )) * ( $f ))), ( $ipb{$tag} = $ipb )) :
					(( &SUM( @$pre[0..2] )), ( $ipb{$tag} += $l )))[0] } ) or (next FIL)); (( &SUM( @$pst[0..2] )) ? (( %{$idx||{}} ) or
					( do { print STDERR 'CANNOT ESTABLISH STATISTICS INDEX FOR FILE '.$$fil[0].$$fil[1]."\n"; (last CHN) } )) : (next FIL));

				my (@key) = ( do { my ($j,@key,%idx,%inc); # move some of this up ...? dont have to protect if read once ... #HERE
					for ( keys %$idx ) { (( m/^([a-z][a-z\d]{2})_(\d{3})$/ ) or (next)); push @{$idx{( qq($1))}||=[]}, (0+ $2); }
					for ( grep {(( not $key ) or ( not $$_[2] ))} map {( &PAIR_KEY_IDX( $_, !!( $key ), 1 ))} ( @{$key||$inc||[]} )) { my ($k,$v) = ( @$_ );
						push @{$inc{$k}||=[]}, (( defined $v ) ? ($v) : ( @{$idx{$k}||[]} )) } ((%inc) or ((%inc) = (%idx)));
					for ( grep {(( not $key ) or ( $$_[2] ))} map {( &PAIR_KEY_IDX( $_, !!( $key ), 1 ))} ( @{$key||$exc||[]} )) { my ($k,$v) = ( @$_ );
						if ( defined $v ) { $inc{$k} = [ grep {( $_ != $v )} ( @{$inc{$k}||[]} ) ] } else { delete $inc{$k}}} (( delete $inc{eid} ), ( delete $inc{wgt} ));
					do { my ($k) = ( q([ ) . ( join q(, ), ( map {($$_[0])} (@$_))) . q( ])); if ( not defined $key_str ) {
						($key_str,$tex_str) = ( $k, ( q([ ) . ( join q(, ), ( map {($$_[1])} (@$_))) . q( ]))); } else {
						( $key_str eq $k ) or do { print STDERR 'INCONSISTENT STATISTIC INDEXES IN TRAINING CYCLE '.$i."\n"; (last CHN) }}} for [
						map { my ($k,$v) = ((ref eq q(ARRAY)) ? (( q(KEY)), (++$j)) : (( uc ((%$_)[0] )), (0+ ((%$_)[1] ))));
							[ ( uc sprintf q("%3.3s_%3.3u"), ($k,$v)), (( defined $key_str ) ? () : ( &DEFINED(( ${$tex{$k}||[]}[$v] ), ( &LATEX_KEY_IDX($k,$v))))) ] }
						grep { ( push @key, (( &HASHED_FUNCTIONAL( $idx, (( ref eq q(ARRAY)) ? ( @$_ ) : ((undef), $_ )))) or (
							do { print STDERR 'INVALID CHANNEL KEY SPECIFICATION IN TRAINING CYCLE '.$i."\n"; (last CHN) } ))) } ((
						map { my ($k) = $_; map { +{ $k => (0+ $_) }} sort { our ($a,$b); ( $a <=> $b ) } ( &UNIQUE( @{$inc{$k}||[]} )) }
						sort { our ($a,$b); ( $a cmp $b ) } ( keys %inc )), ( grep {( ref eq q(ARRAY))} ( @{$key||[]} ))) ];
					((@key) or do { print STDERR 'EMPTY FEATURE KEY LIST IN TRAINING CYCLE '.$i."\n"; (last CHN) } ); (@key) } );
				my (@cut) = ( grep {(( $$_[2] = ( &HASHED_FUNCTIONAL( $idx, (
					map {((ref eq 'ARRAY') ? ( @$_ ) : ((undef), $_ ))} ($$_[2][0]))))) or (
						do { print STDERR 'INVALID KEY IN SELECTION '.$$_[1].' FOR TRAINING CYCLE '.$i.' ON FILE '.$$fil[0].$$fil[1]."\n"; !1 } ))}
					grep {(! ( &MATCH_VALUE( $$_[3], undef )))} map {[ ($_ < 0), ( abs ), ( @{ (($_) && ( $$crd{esc}[( abs )] )) or
						do { print STDERR 'INVALID EVENT SELECTION CUT SPECIFICATION IN TRAINING CYCLE '.$i."\n"; +{}}}{( qw( key cut ))} ) ]}
					map {((defined) ? ( int ) : ())} ( @{$esc||[]} ));
				my ($wgt) = ( map {((defined) ? (( &HASHED_FUNCTIONAL( $idx,
					((ref eq 'ARRAY') ? ( @$_ ) : ((undef),((ref eq 'HASH') ? ( $_ ) : +{ wgt => (0+ $_) } ))))) or
					do { print STDERR 'INVALID CHANNEL WEIGHT SPECIFICATION IN TRAINING CYCLE '.$i."\n"; (last CHN) } ) :
					(( &HASHED_FUNCTIONAL( $idx, (undef), +{ wgt => 0 } )) or ( sub {( $w )} )))} (${$wgt||[]}[0]));

				push @{ $FHT{$tag} ||= [] }, [ $l, (( &Local::FILE::HANDLE()) or ( die 'Cannot open temporary file for read/write' )) ];
				local ($_); while (<$FHI>) { ((/^\s*$/) and (next)); ((/^\s*(\d+)/) && ($1 <= $nmx)) or (last);
					do { my ($val) = $_; print (( join q(,), ( map {((defined) ? ($_) : ( q(NAN)))}
						((($wgt) -> ($val)), ( map {(($_) -> ($val))} (@key))))),( qq(\n))); } for
					grep { my ($val,$mch) = ($_,1); for (@cut) { my ($inv,$eid,$key,$cut) = (@$_);
						(( $mch = (($inv) xor ( &MATCH_VALUE( $cut, (($key) -> ($val)))))) or (last)) }; ($mch) }
					[ map { (/^UNDEF$/) ? (undef) : (0+ $_) } ( split ) ]; }}

			my ($FHO); do { my ($tag) = $_; if ((defined $ipb) && ((( $ipb * $ipb{$tag} ) <= 0 ) or ( $ipb > $ipb{$tag} ))) {
				print STDERR 'RESCALING BY '.( uc sprintf '%+10.3e', ( &RATIO( $ipb, $ipb{$tag} ))) .
					' TO TARGET LUMINOSITY OF '.( uc sprintf '%+10.3e', ($ipb)).' PER PB IN CHANNEL '.($tag)."\n"; }
				do { my ($l,$FHT) = @$_; (( seek $FHT, 0, SEEK_SET ) or ( die 'Cannot rewind temporary file' )); local ($_); while (<$FHT>) {
					my ($wgt,$lin) = split (( q(,)), $_, 2 ); (( $wgt *= ( &RATIO( $l, $ipb{$tag} ))) or (next));
					$FHO ||= ( grep {((( @$_ ) && ( push @{ $lbl_str[-1] }, $tex )) or ( die 'Cannot open file in directory '.$out.' for write' ))}
						[ ( &Local::FILE::NEXT([[ $out, q(CSV) ], [ $lbl, 0, q(csv) ]])) ] )[0];
					print +( $wgt, ( q(,)), $lin ); }} for (@{ $FHT{$tag}}) } for ( sort { our ($a,$b); ( $a cmp $b ) } keys %FHT ) }}

	if (1) {
	# ??? how to include/exclude? +/- or sep inputs? do target luminosity &  allow rescale/ absolute? how to file/merge ... memory? group by DAT or by CHN?
	# do { print STDERR 'NO BINNED EVENTS ESTABLISHED FOR TRAINING CYCLE '.$i."\n"; (next BDT) } unless (1);
	my ($fpo) = $out . q(minos.py); my (%dat) = (
		PYT	=> ( $pyt = ( scalar &PYTHON_THREE($$bdt{py3}))),
		VRB	=> ( q(False)),
		KEY	=> ( $key_str ),
		TEX	=> ( $tex_str ),
		LBL	=> ( q([) . ( join q(, ), ( map {( q([ ) . ( join q(, ), ( @$_ )) . q( ]))} ( @lbl_str ))) . q(]) ),
		PRT	=> ( 3 ),
		FLD	=> ( q(True)),
		MIX	=> ( uc sprintf q(%+5.2f), 0. ),
		TPR	=> ( uc sprintf q(%+5.2f), 0. ),
		BAL	=> ( uc sprintf q(%+5.2f), 0. ),
		LUM	=> ( &DEFINED( $ipb, 1 )),
		MOD	=> ( 0 ),
		BNS	=> ( 100 ),
		LVL	=> ( 5 ),
		ALP	=> ( 0. ),
		GAM	=> ( 0. ),
		ETA	=> ( 0.5 ),
		LAM	=> ( 0.01 ),
		MCW	=> ( 0. ),
		MDS	=> ( 0. ),
		BAS	=> ( 0.5 ),
		SPW	=> ( 1. ),
		SMP	=> ( 1. ),
		CST	=> ( 1. ),
		CSL	=> ( 1. ),
		CSN	=> ( 1. ),
		);
	do { use Fcntl qw(:seek); local ($.,$?); my ($t) = ( tell DATA );
		( my ($FHO) = ( &Local::FILE::HANDLE($fpo,1))) or ( die 'Cannot write to file '.($fpo));
		local ($_); while (<DATA>) { s/<\[(\w+)]>/$dat{$1}/g; ( print $FHO $_ ) }
		( close $FHO ) && ( chmod 0755, $fpo ); ( seek DATA, $t, SEEK_SET );
		if ( &CAN_MATPLOTLIB( $pyt, 1 )) { system( qq( cd ${out} && ./minos.py )) }
		else { print STDERR 'CANNOT VERIFY PYTHON 2.7/3.X WITH MATPLOTLIB 1.3.0+ AND XGBOOST (DELIVERING SCRIPT LITERAL)'."\n"; }}}}};

1

__DATA__
#!/usr/bin/env <[PYT]>

#####################################
# import modules and validate runtime
#####################################

import sys
if ((sys.version_info[0] < 2) or ((sys.version_info[0] == 2) and (sys.version_info[1] < 7))) :
    sys.exit( "MInOS requires Python versions 2.7 or 3.X" )

import matplotlib as mpl
if (( tuple( map ( int, mpl.__version__.split("."))) + (0,0,0))[0:3] < (1,3,0)) :
    sys.exit( "MInOS requires MatPlotLib version 1.3.0 or Greater" )

import warnings; warnings.filterwarnings("ignore")

import math, os, glob, re

import numpy as np

import matplotlib.pyplot as plt

import xgboost as xgb

#####################
# define main routine
#####################

def main () :
    # analysis parameters
    vrb = <[VRB]>; prt = <[PRT]>; fld = <[FLD]>; mix = <[MIX]>; tpr = <[TPR]>; bal = <[BAL]>; lum = <[LUM]>; mod = <[MOD]>; bns = <[BNS]>
    # tuple per class of tuple per valid channel of dictionary with channel index and tuple of csv partitions
    csv = tuple(( lambda x : x if len( x ) else sys.exit( "Insufficient Events Available for Training" ))(
        tuple( { "chn":chn, "csv":csv } for chn, csv in ((( lambda y : int( y.group( 1 )) if y else 0 )(
            re.search( "\/(?:X01|X00)_(\d{3})\.csv$", fil )), CSV( fil, prt=prt ))
            for fil in sorted( glob.glob( cls ))) if ( chn > 0 ) and len( csv )))
        for cls in ( "./CSV/X00_*.csv", "./CSV/X01_*.csv" ))
    for i, dmx in enumerate( DMX( csv, fld=fld, mix=mix, tpr=tpr, bal=bal )) :
        if vrb and i == 0 : print( "Initial Cross Sections [fb]:\n" + str( fbXSC( dmx, prt=prt )))
        doPlots( BDT( dmx, typ=1, fld=1+i ), vrb=vrb, pth="./Plots", lum=lum, mod=mod, bns=bns, typ=1, fld=1+i )

    #dmx_all = tuple( DMX( tuple( { "chn":None, "csv":mergeCSV( *( x["csv"]
    #    for x in cls )) } for cls in csv ), fld=fld, mix=mix, tpr=tpr, bal=bal ))
    #bdt_all = tuple( BDT( f ) for f in dmx_all )
    #for i, f in enumerate( bdt_all ) :
    #    doPlots( f, lum=lum, typ=1, fld=1+i,
    #        dat=tuple( "Joint" if len( c ) > 1 else label_tex[j][c[0]["chn"]-1] for j, c in enumerate( csv )))
    # do split training with recombination
    #if len( csv[0] ) > 1 :
        #dmx = tuple( tuple( DMX(( c, csv[1][0] ), fld=fld, mix=mix, tpr=tpr, bal=bal )) for c in csv[0] )
        #bdt = tuple( tuple( BDT( f ) for f in c ) for c in dmx )
        # test each background separately against the joint training
        #for c in dmx :
        #    for i, f in enumerate( c ) :
        #        cid = f["test"]["chn"][0]
        #        doPlots( BDT( f, mdl=bdt_all[i]["mdl"] ),
        #            lum=lum, typ=100+cid, fld=1+i, dat=( label_tex[0][cid-1], label_tex[1][0] ))
        # train and test against each background separately
        #for c in bdt :
        #    for i, f in enumerate( c ) :
        #        cid = f["dmx"]["test"]["chn"][0]
        #        doPlots( f, lum=lum, typ=200+cid, fld=1+i, dat=( label_tex[0][cid-1], label_tex[1][0] ))
        # merge separate trainings and test against the joint background
        #for i, f in enumerate( dmx_all ) :
        #    doPlots( BDT( f, mdl=mergeMDL( *( c[i] for c in bdt ))),
        #        lum=lum, typ=301, fld=1+i, dat=( "Joint", label_tex[1][0] ))
    return

#################################################
# specify global settings, inputs, and parameters
#################################################

# configure LaTeX and fonts settings
mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'
mpl.rcParams['mathtext.fontset'] = 'cm';
mpl.rcParams['font.family'] = 'STIXGeneral'

# select color palette
colors = plt.get_cmap( "tab10", lut=10 )

# ordered list of training feature keys
feature_names = <[KEY]>

# ordered list of training feature TeX symbols
feature_tex = <[TEX]>

# ordered list of data set TeX labels
label_tex = <[LBL]>

# configure hyperparameters
param = {
    "booster":"gbtree",
    "tree_method":"auto",
    "eval_metric":["auc","logloss"],
    "disable_default_eval_metric":1,
    "max_depth":<[LVL]>,
    "alpha":<[ALP]>,
    "gamma":<[GAM]>,
    "eta":<[ETA]>,
    "lambda":<[LAM]>,
    "min_child_weight":<[MCW]>,
    "max_delta_step":<[MDS]>,
    "objective":"binary:logistic",
    "base_score":<[BAS]>,
    "scale_pos_weight":<[SPW]>,
    "subsample":<[SMP]>,
    "colsample_bytree":<[CST]>,
    "colsample_bylevel":<[CSL]>,
    "colsample_bynode":<[CSN]>,
    "verbosity":1
    }

####################
# define subroutines
####################

# generate collection of plots corresponding to an input dmx object
def doPlots ( bdt, vrb=False, pth="./Plots", lum=1., mod=0, bns=0, typ=None, fld=None ) :
    try : os.makedirs( pth )
    except FileExistsError : pass
    except : sys.exit( "Unable to Create File Path {0:s}".format( pth ))
    plotImportance( bdt, pth=pth, lum=lum, typ=typ, fld=fld )
    plotROC ( bdt, pth=pth, lum=lum, typ=typ, fld=fld )
    plotScatter( bdt, pth=pth, lum=lum, mod=mod, typ=typ, fld=fld )
    plotDistribution( bdt, pth=pth, lum=lum, mod=mod, bns=bns, typ=typ, fld=fld )
    plotSigma( bdt, pth=pth, lum=lum, mod=mod, bns=bns, typ=typ, fld=fld )
    plotSignificance( bdt, vrb=vrb, pth=pth, lum=lum, bns=bns, typ=typ, fld=fld )
    return

# rounded fb cross sections averaged over training and validation samples
def fbXSC ( dmx, prt=3 ) :
    ( xsc, train, test ) = ( [], dmx["train"]["xsc"], dmx["test"]["xsc"] )
    for i, cls in enumerate( train ) :
        xsc.append( [] )
        for j, chn in enumerate( cls ) :
            xsc[-1].append( float( '%.3g' %
                ( 1000. * (( prt - 1) * train[i][j] + test[i][j] ) / prt )))
    return xsc

# partition shuffled CSV data into tuple of numpy arrays (events) of arrays (weight & features)
def CSV ( *fil, prt=3, rns=0 ) :
    # minimal partition is 2 for separation of training and test samples
    prt = max( 2, int( prt ))
    # numpy import utility generates uniform 2x2 matrices of floating point values to be joined across files
    csv = (( lambda x : np.concatenate( x, axis=0 ) if len(x) > 1 else x[0] if x else np.empty((0,1), dtype=np.float64 ))(
        tuple( filter(( lambda x : x.shape[0] ), ( np.loadtxt( x, delimiter=",", dtype=np.float64, ndmin=2 ) for x in fil )))))
    # fail out records with no features or fewer events than the partition
    if csv.shape[1] < 2 or csv.shape[0] < prt : return tuple()
    # shuffle event records in place using fixed seed
    np.random.default_rng( seed=rns ).shuffle( csv, axis=0 )
    # return tuple of event partitions
    return tuple( np.array_split( csv, prt, axis=0 ))

# merge sets of CSV data partitions into a single tuple of numpy arrays of arrays
def mergeCSV ( *csv ) :
    return tuple( np.concatenate( x, axis=0 ) if len(x) > 1 else x[0]
        for x in zip( *filter( None, csv )))

# yield dictionaries of dmx objects for each fold from input CSV partitions for each class
def DMX ( csv, fld=False, mix=0., tpr=0., bal=0. ) :
    bal = max( -0.99, min( +0.99, float( bal ))) if len( csv ) == 2 else 0.
    # count available partitions for each class
    prt = min( tuple( len( chn["csv"] ) for cls in csv for chn in cls ) or (0,))
    # internal method for construction of dmx objects
    def dmx ( csv, prt, inc ) :
        # fold selected csv partitions and separate weights from features
        fld = tuple( tuple(( lambda chn, csv : { "chn":chn, "wgt":csv[0].flatten(), "ftr":csv[1] } )(
            chn["chn"], np.split( np.concatenate( tuple( chn["csv"][j] for j in inc )), [1], axis=1 ))
            for chn in cls ) for cls in csv )
        # copy per-event weights for each class
        wgt = tuple( tuple( np.absolute( chn["wgt"] ) for chn in cls ) for cls in fld )
        # count residual events in each class
        evt = tuple( tuple( len( chn ) for chn in cls ) for cls in wgt )
        # sum cross section of residual events in each class
        xsc = tuple( tuple( chn.sum() for chn in cls ) for cls in wgt )
        # return data structure with merged dmatrix, event counts, and cross sections
        return {
            # instance of core XGBoost data matrix object
            "dmx":xgb.DMatrix(
                # merge event samples from each class into a list of feature lists
                data=np.concatenate( tuple( chn["ftr"] for cls in fld for chn in cls ), axis=0 ),
                # generate a list with the class label for each merged event
                label=np.concatenate( tuple( np.full( sum( chn for chn in cls ), i, dtype=np.int32 )
                    for i, cls in enumerate( evt )), axis=0 ),
                # normalize total cross section weight to unity for each class
                weight=np.concatenate( tuple(( lambda i, x : x * (( 1. - bal if i == 0 else 1. + bal ) / x.sum()))( i,
                    taper( np.concatenate( tuple( c * ( m / x if x > 0. else 0. ) for c, x, m in
                    zip( cls, xsc[i], taper( xsc[i], mix ))), axis=0 ), tpr )) for i, cls in enumerate( wgt )), axis=0 ),
                # assign plain text names for each ordered training feature
                feature_names=feature_names ),
            # array with normalized event weights belonging to each class
            "wgt":np.concatenate( tuple( chn / ( xsc[i][j] or 1. )
                for i, cls in enumerate( wgt ) for j, chn in enumerate( cls )), axis=0 ),
            # array with number of sequential event samples belonging to each class
            "evt":evt,
            # array with scaled inclusive physical cross sections for each class
            "xsc":tuple( tuple( chn * prt / len( inc ) for chn in cls ) for cls in xsc ),
            # tuple with channel identification indices for each class
            "chn":tuple( tuple( chn["chn"] for chn in cls ) for cls in fld ) }
    # generator function yields dictionary of (train,test) dmx objects for each fold
    for i in range( prt if fld else prt and 1 ) : yield { key : dmx( csv, prt, inc )
        for key, inc in (( "train", tuple( j for j in range( prt ) if j != i )), ( "test", ( i, ))) }
#HERE ... generalize for negative weights ...

# generate BDT dictionary object with test dmatrix, model, score, and density
def BDT ( dmx, mdl=None, typ=None, fld=None ) :
    if mdl is None :
        if not os.path.exists( "./JSON/" ) :
            os.mkdir( "./JSON/" )
        fil = os.path.join( "./JSON/", "model" +
            indexString( typ ) + indexString( fld ) + ".json" )
        if os.path.isfile( fil ) :
            mdl = xgb.Booster(); mdl.load_model( fil )
        else :
            mdl = MDL( dmx ); mdl.save_model( fil )
    if not isinstance( mdl, xgb.core.Booster ) : return
    return { "dmx":dmx, "mdl":mdl, "scr":SCR( dmx, mdl ) }
#HERE generalize save and load to work for mergeMDL class via overload

# generate a mdl object via training and validation of dmx object pair
def MDL ( dmx, trs=50, stp=5 ) :
    prg = dict(); return xgb.train( param, dmx["train"]["dmx"],
        num_boost_round=trs, early_stopping_rounds=stp,
        evals=[ ( dmx["train"]["dmx"], "train" ), ( dmx["test"]["dmx"], "test" ) ],
        verbose_eval=False, evals_result=prg )
#HERE print ( bst.best_score , bst.best_iteration )

# generate score object by projecting member (default:test) of dmx object pair onto mdl object
def SCR ( dmx, mdl ) :
    ( dmx, wgt, evt, xsc, chn ) = ( dmx["test"][x]
        for x in ( "dmx", "wgt", "evt", "xsc", "chn" ))
    scr = mdl.predict( dmx ); cnt = counter( 0 )
    # return tuple with dictionary of sorted data for each supervisory class
    return tuple( tuple(
        # generate sorted data structure with scores, weights, and cross section
        ( lambda i, j, idx : { "scr":scr[idx], "wgt":wgt[idx], "dns":DNS( scr[idx], wgt[idx] ), "xsc":xsc[i][j], "chn":chn[i][j] } )(
            # sort class indices according to classification score
            i, j, cnt( 0 ) + np.argsort( scr[ cnt( 0 ) : cnt( evt[i][j] ) ], kind="stable" ))
            # establish boundary indices of events in ith class
            # perform outer loop over each supervisory class
        for j in range( len( evt[i] ))) for i in range( len( evt )))

# generate an interpolated score density object dns from input scores and weights
def DNS ( scr, wgt, bins=None, b2=1000, smooth=1. ) : return interpolateTrapezoid(
    *( pointsFromDensityEdge( *( binDensity( scr, wgt, bins=bins, b2=b2, smooth=smooth )))))

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
    wdt *= smooth; w /= smooth
    edg = np.linspace( 0., 1., num=( 1 + b2 ), endpoint=True )
    cen = ( edg[1:] + edg[:-1] ) / 2.
    val = np.add.reduce( tuple ( w[i] * np.exp( np.square(( cen - e[i] ) / wdt[i] ) / ( -2. )) for i in range( len( e ))), axis=0 )
    val /= ( val.sum() / b2 )
    return tuple(( val, edg ))
#HERE consider non-flipped association? or even sharing / other? see minos_break_ ...

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

# vectorized accessor method for polynomial objects
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

# normalize ndarray from 0. to 1. and apply taper
def taper ( x, tpr=0. ) :
    x = np.absolute( x ); m = np.max( x ); return (
        np.zeros( x.shape, dtype=x.dtype ) if m == 0. else
        np.vectorize( lambda x : 0. if x == 0. else 1. )( x ) if tpr <= -0.99 else
        x / m if tpr == 0. else
        np.power( x / m, math.tan(( 1. + tpr ) * math.pi / 4. )) if tpr < +0.99 else
        np.vectorize( lambda x : 1. if x == m else 0. )( x ))

# construct closure for accumulating counts
def counter ( init=0 ) :
    value = [ init ]
    def closure ( step=1 ) :
        value[0] += step
        return value[0]
    return closure

# generate donut plot of feature importance to total gain
def plotImportance ( bdt, lum=1., pth="./Plots", typ=None, fld=None ) :
    if not isinstance( bdt["mdl"], xgb.core.Booster ) : return
    scores = bdt["mdl"].get_score( importance_type="total_gain" )
    keys = sorted( scores, key=scores.get, reverse=True )
    values = np.array( [ scores[k] for k in keys ] )
    importance = []; values /= values.sum()
    for k, v in zip( keys, values ) :
        if not importance or importance[-1][0] + v > 0.075 : importance.append( [ 0., []] )
        importance[-1][0] += v; importance[-1][1].append( k )
    fig = plt.figure( figsize=( 7.5, 5. ), tight_layout=False )
    ax = fig.add_axes( [ 0., .05, 1., .75 ] )
    labels = tuple((( lambda x : ", ".join( tuple( feature_dict[k] for k in x[:3] ) +
        (( r"$\ldots$", ) if ( len( x ) > 3 ) else ())))( tuple( _[1] ))) for _ in importance )
    wedges, texts, autotexts = ax.pie( [ _[0] for _ in importance ],
        normalize=True, radius=1., startangle=90., autopct="%1.1F%%", pctdistance=0.75,
        wedgeprops={ "width":0.5, "edgecolor":"black", "linewidth":0.8 }, colors=colors.colors )
    plt.setp( autotexts, size=12, color="white" )
    kw = { "zorder": 0, "verticalalignment":"center", "arrowprops":{ "arrowstyle":"-", "linewidth":0.8 },
        "bbox":{ "boxstyle":"round,pad=0.3,rounding_size=0.3", "facecolor":"white", "edgecolor":"black", "linewidth":0.8 }}
    for i, p in reversed( list( enumerate( wedges ))) :
        ang = (( p.theta2 + p.theta1 ) / 2. ); y = np.sin( np.deg2rad( ang )); x = np.cos( np.deg2rad( ang ))
        horizontalalignment = { -1: "right", 1: "left"}[ int( np.sign( x )) ]
        connectionstyle = "angle,angleA=0,angleB={0:G}".format( ang )
        kw["arrowprops"].update( { "connectionstyle":connectionstyle } )
        ax.annotate( labels[i], xy=( x, y ), size=12, xytext=( 1.25*np.sign( x ), 1.25*y ),
            horizontalalignment=horizontalalignment, **kw )
    fig.suptitle( "Feature Importance to Total Gain in Training Fold {0:d}\n".format( fld ) + " ".join( tuple(
        (( lambda x : "" if x is None or len( x ) == 0 else x + " " )( label_tex[i][ cls[0]["chn"] - 1 ] )
        if len( cls ) == 1 else "Joint " ) + ( "Background vs.", "Signal" )[i] for i, cls in enumerate( bdt["scr"][0:2] ))),
        x=0.5, y=0.9, size=17, verticalalignment="center" )
    fig.savefig( os.path.join( pth, "donut_plot" + indexString( typ ) + indexString( fld ) + ".pdf" ), facecolor="white" )

# generate plot of Receiver Operating Characteristic evolution with mdl score
def plotROC ( bdt, lum=1., pth="./Plots", typ=None, fld=None ) :
    fig = plt.figure( figsize=( 7.5, 5. ), tight_layout=False )
    ax = fig.add_axes([ 0.16, .14, 0.68, 0.68 ], aspect="equal" )
    ax.set_xlim([ 0., 1. ]); ax.set_ylim([ 0., 1. ])
    ax.set_xlabel( "False Signal Rate", size=14, color="black" )
    ax.set_ylabel( "True Signal Rate", size=14, color="black" )
    ax.xaxis.grid( True, linewidth=0.8, color="black" )
    ax.yaxis.grid( True, linewidth=0.8, color="black" )
    xxx = np.linspace( 0., 1., num=( 1 + 2500 ), endpoint=True ); yyy = []
    for i, cls in enumerate( bdt["scr"][0:2] ) :
        if len( cls ) > 1 : xsc = np.array( tuple(( chn["xsc"] ,) for chn in cls ))
        yyy.append( 1. - (( xsc * np.array( tuple( chn["dns"]( xxx, -1 )
            for chn in cls ))).sum( axis=0 ) / xsc.sum( axis=None )
            if len( cls ) > 1 else cls[0]["dns"]( xxx, -1 )))
    ax.plot( yyy[0], yyy[1], color=colors(0), linewidth=1.4, linestyle="solid", zorder=3 )
    ax.fill_between( yyy[0], yyy[1], 0., linewidth=0., alpha=0.1, edgecolor="none", facecolor=colors(0), hatch="", zorder=-1 )
    ax.plot( xxx, xxx, color="red", linewidth=1.4, linestyle="dashed", zorder=4 )
    ax.text( 0.45, 0.055, "Area Under Curve: {0:4.2F}".format(
        (( yyy[1][:-1] + yyy[1][+1:] ) * ( yyy[0][+1:] - yyy[0][:-1] )).sum() / ( -2. )), fontsize=12, zorder=5,
        bbox={ "boxstyle":"round,pad=0.3,rounding_size=0.3", "facecolor":"white", "edgecolor":"black", "linewidth":0.8 } )
    tkw = dict( size=4, width=0.8 ); ax.tick_params( axis="x", **tkw )
    ax.tick_params( axis="y", colors="black", **tkw )
    for x in ( ax.get_xticklabels() + ax.get_yticklabels()): x.set_fontsize( 12 )
    ax.tick_params( axis="y", which="both", zorder=0, color="black" ); ax.minorticks_on()
    for (_,i) in ax.spines.items(): i.set( linewidth=0.8, color="black", alpha=1., zorder=0 )
    fig.suptitle( "Receiver Operating Characteristic for Validation Fold {0:d}\n".format( fld ) + " ".join( tuple(
        (( lambda x : "" if x is None or len( x ) == 0 else x + " " )( label_tex[i][ cls[0]["chn"] - 1 ] )
        if len( cls ) == 1 else "Joint " ) + ( "Background vs.", "Signal" )[i] for i, cls in enumerate( bdt["scr"][0:2] ))),
        x=0.5, y=0.9, size=17, verticalalignment="center" )
    fig.savefig( os.path.join( pth, "roc_plot" + indexString(typ) + indexString(fld) + ".pdf" ), facecolor="white" )

# generate scatter plot of signal/bg event weights as a function of mdl score
def plotScatter ( bdt, mod=0, lum=1., pth="./Plots", typ=None, fld=None ) :
    fig = plt.figure( figsize=( 7.5, 5. ), tight_layout=False )
    ax = fig.add_axes([ 0.125, .14, 0.8, 0.725 ] )
    ax.set_xlim([ 0., 1. ])
    ax.set_xlabel( "Signal Classification Score", size=14, color="black" )
    ax.set_ylabel( "Fractional Cross Section", size=14, color="black" )
    ax.xaxis.grid( True, linewidth=0.8, color="black" )
    tkw = dict( size=4, width=0.8 )
    ax.tick_params( axis="x", **tkw )
    labels = []; patches = []
    for i, cls in enumerate( bdt["scr"][0:2] ) :
        for j, chn in ((( 0, cls[0] ),) if mod >= 0 and len( cls ) == 1 else ((( 0, None ),) if mod >= 0 else tuple()) +
            ( tuple(( j+1, chn ) for j, chn in enumerate( cls )) if mod <= 0 else tuple())) :
            if chn is None : xsc = np.array( tuple( chn["xsc"] for chn in cls ))
            ax.scatter( np.concatenate( tuple( chn["scr"] for chn in cls ), axis=0 ) if chn is None else chn["scr"],
                np.concatenate( tuple( xsc[i] * chn["wgt"] for i, chn in enumerate( cls ))) / xsc.sum( axis=0 ) if chn is None else chn["wgt"],
                facecolor=( colors(i) if j == 0 else str( 0.85 - 0.1*j )), edgecolor=shades[i](j), alpha=0.1, zorder=( 2 if j == 0 else 3 ))
            patches.append(( mpl.patches.Rectangle(( 0, 0 ), 1, 1, linewidth=0., facecolor=( colors(i) if j == 0 else str( 0.85 - 0.1*j )), alpha=0.5, fill=True ),
                mpl.patches.Rectangle(( 0, 0 ), 1, 1, linewidth=1.4, linestyle="solid", edgecolor=shades[i](j), fill=None )))
            labels.append( " ".join( filter( None, (( label_tex[i][ chn["chn"] - 1 ] or None ) if chn is not None else "Joint" ,) +
                (( "Background", "Signal" )[i] if j == 0 else None ,))))
    lgd = ax.legend( patches, labels, loc="best", fontsize=12 )
    lgd.get_frame().set( facecolor="white", linewidth=0.8, edgecolor="black", alpha=1. ); lgd.set(zorder=4)
    ax.set_ylim([ 1.0E-12, 1.0E+00 ])
    ax.set_yscale("log")
    ax.tick_params( axis="y", colors="black", **tkw )
    for x in ( ax.get_xticklabels() + ax.get_yticklabels()) : x.set_fontsize( 12 )
    ax.tick_params( axis="y", which="both", zorder=0, color="black" ); ax.minorticks_on()
    for (_,i) in ax.spines.items(): i.set( linewidth=0.8, color="black", alpha=1., zorder=0 )
    fig.suptitle( "Normalized Sample Scatter for Validation Fold {0:d}".format( fld ), x=0.5, y=0.9, size=17, verticalalignment="bottom" )
    fig.savefig( os.path.join( pth, "scatter_plot" + indexString(typ) + indexString(fld) + ".pdf" ), facecolor="white" )

# generate plot of mdl score distributions for signal and background
def plotDistribution ( bdt, mod=0, bns=0, lum=1., pth="./Plots", typ=None, fld=None ) :
    fig = plt.figure( figsize=( 7.5, 5. ), tight_layout=False )
    ax = fig.add_axes([ 0.1, .14, 0.8, 0.725 ] ); ax.set_xlim([ 0., 1. ])
    ax.set_xlabel( "Signal Classification Score", size=14, color="black" )
    ax.set_ylabel( "Probability Density", size=14, color="black" )
    ax.xaxis.grid( True, linewidth=0.8, color="black" )
    tkw = dict( size=4, width=0.8 )
    ax.tick_params( axis="x", **tkw )
    labels = []; patches = []
    xxx = np.linspace( 0., 1., num=( 1 + 2500 ), endpoint=True )
    for i, cls in enumerate( bdt["scr"][0:2] ) :
        for j, chn in ((( 0, cls[0] ),) if mod >= 0 and len( cls ) == 1 else ((( 0, None ),) if mod >= 0 else tuple()) +
            ( tuple(( j+1, chn ) for j, chn in enumerate( cls )) if mod <= 0 else tuple())) :
            if chn is None : xsc = np.array( tuple(( chn["xsc"] ,) for chn in cls ))
            yyy = ( xsc * np.array( tuple( chn["dns"]( xxx ) for chn in cls ))).sum( axis=0 ) / xsc.sum( axis=None ) if chn is None else chn["dns"]( xxx )
            ax.plot( xxx, yyy, linewidth=1.4, linestyle=( 0, (( 8 - j )/1.5, j/1.5 )),
                color=shades[i](j), gapcolor=str( 0.85 - 0.1*j ), zorder=( 2 if j == 0 else 3 ))
            if j == 0 : ax.fill_between( xxx, yyy, 0., linewidth=0., edgecolor="none", facecolor=colors(i), alpha=0.1, hatch="", zorder=-1 )
            patches.append((( mpl.patches.Rectangle(( 0, 0 ), 1, 1, linewidth=0., facecolor=colors(i), alpha=0.1, fill=True ), ) if j == 0 else
                ( mpl.patches.Rectangle(( 0, 0 ), 1, 1, linewidth=1.4, linestyle="solid", edgecolor=str( 0.85 - 0.1*j ), fill=None ), )) +
                ( mpl.patches.Rectangle(( 0, 0 ), 1, 1, linewidth=1.4, linestyle=( 0, (( 8 - j )/1.5, j/1.5 )), edgecolor=shades[i](j), fill=None ), ))
            labels.append( " ".join( filter( None, (( label_tex[i][ chn["chn"] - 1 ] or None ) if chn is not None else "Joint" ,) +
                (( "Background", "Signal" )[i] if j == 0 else None ,))))
            if j == 0 and int( bns ) > 0 : ax.hist(
                np.concatenate( tuple( chn["scr"] for chn in cls ), axis=0 ) if chn is None else chn["scr"],
                weights=(( int( bns ) / xsc.sum( axis=None )) * np.concatenate( tuple( xsc[i][0] * chn["wgt"]
                    for i, chn in enumerate( cls )), axis=0 )) if chn is None else int( bns ) * chn["wgt"],
                bins=int( bns ), range=( 0., 1. ), color=colors(i), alpha=0.3, histtype="stepfilled", zorder=1 )
    lgd = ax.legend( patches, labels, loc="best", fontsize=12 )
    lgd.get_frame().set( facecolor="white", linewidth=0.8, edgecolor="black", alpha=1. ); lgd.set(zorder=4)
    ax.set_ylim([ 0., 20. ]) #HERE temporarily insert hard upper bound - make user adjustable # ax.set_ylim([ 0., None ])
    ax.tick_params( axis="y", colors="black", **tkw )
    for x in ( ax.get_xticklabels() + ax.get_yticklabels()) : x.set_fontsize( 12 )
    ax.tick_params( axis="y", which="both", zorder=0, color="black" ); ax.minorticks_on()
    for (_,i) in ax.spines.items(): i.set( linewidth=0.8, color="black", alpha=1., zorder=0 )
    fig.suptitle( "Normalized Event Distribution for Validation Fold {0:d}".format( fld ), x=0.5, y=0.9, size=17, verticalalignment="bottom" )
    fig.savefig( os.path.join( pth, "density_plot" + indexString(typ) + indexString(fld) + ".pdf" ), facecolor="white" )

# generate plot of signal/bg events as a function of mdl score
def plotSigma ( bdt, mod=0, bns=0, lum=1., pth="./Plots", typ=None, fld=None ) :
    fig = plt.figure( figsize=( 7.5, 5. ), tight_layout=False )
    ax = fig.add_axes([ 0.1, .14, 0.8, 0.725 ] ); ax.set_xlim([ 0., 1. ])
    ax.set_xlabel( "Signal Classification Threshold", size=14, color="black" )
    ax.set_ylabel( r"Residual Cross Section (fb)", size=14, color="black" )
    ax.xaxis.grid( True, linewidth=0.8, color="black" )
    tkw = dict( size=4, width=0.8 )
    ax.tick_params( axis="x", **tkw )
    labels = []; patches = []
    xxx = np.linspace( 0., 1., num=( 1 + 2500 ), endpoint=True )
    for i, cls in enumerate( bdt["scr"][0:2] ) :
        for j, chn in ((( 0, cls[0] ),) if mod >= 0 and len( cls ) == 1 else ((( 0, None ),) if mod >= 0 else tuple()) +
            ( tuple(( j+1, chn ) for j, chn in enumerate( cls )) if mod <= 0 else tuple())) :
            yyy = 1000. * sum( chn["xsc"] * ( 1. - chn["dns"]( xxx, -1 )) for chn in ( cls if chn is None else ( chn, )))
            ax.plot( xxx, yyy, linewidth=1.4, linestyle=( 0, (( 8 - j )/1.5, j/1.5 )),
                color=shades[i](j), gapcolor=str( 0.85 - 0.1*j ), zorder=( 2 if j == 0 else 3 ))
            if j == 0 : ax.fill_between( xxx, yyy, 0., linewidth=0., edgecolor="none", facecolor=colors(i), alpha=0.1, hatch="", zorder=-1 )
            patches.append((( mpl.patches.Rectangle(( 0, 0 ), 1, 1, linewidth=0., facecolor=colors(i), alpha=0.1, fill=True ), ) if j == 0 else
                ( mpl.patches.Rectangle(( 0, 0 ), 1, 1, linewidth=1.4, linestyle="solid", edgecolor=str( 0.85 - 0.1*j ), fill=None ), )) +
                ( mpl.patches.Rectangle(( 0, 0 ), 1, 1, linewidth=1.4, linestyle=( 0, (( 8 - j )/1.5, j/1.5 )), edgecolor=shades[i](j), fill=None ), ))
            labels.append( " ".join( filter( None, (( label_tex[i][ chn["chn"] - 1 ] or None ) if chn is not None else "Joint" ,) +
                (( "Background", "Signal" )[i] if j == 0 else None ,))))
            if j == 0 and int( bns ) > 0 :
                yyy, edg = np.histogram( np.concatenate( tuple( chn["scr"] for chn in cls ), axis=0 ) if chn is None else chn["scr"],
                    weights=1000.*(( np.concatenate( tuple( chn["xsc"] * chn["wgt"] for chn in cls ), axis=0 )) if chn is None else chn["xsc"] * chn["wgt"] ),
                    bins=int( bns ), range=( 0., 1. ))
                for k in range( len( yyy ) -1, 0, -1 ) : yyy[ k - 1 ] += yyy[k]
                ax.hist( edg[:-1], bins=edg, weights=yyy, color=colors(i), alpha=0.3, histtype="stepfilled", zorder=1 )
    lgd = ax.legend( patches, labels, loc="best", fontsize=12 )
    lgd.get_frame().set( facecolor="white", linewidth=0.8, edgecolor="black", alpha=1. ); lgd.set(zorder=4)
    ax.set_ylim([ 0.01, None ]); ax.set_yscale("log")
    ax.tick_params( axis="y", colors="black", **tkw )
    for x in ( ax.get_xticklabels() + ax.get_yticklabels()) : x.set_fontsize( 12 )
    ax.tick_params( axis="y", which="both", zorder=0, color="black" ); ax.minorticks_on()
    for (_,i) in ax.spines.items(): i.set( linewidth=0.8, color="black", alpha=1., zorder=0 )
    fig.suptitle( "Integrated Event Distribution for Validation Fold {0:d}".format( fld ), x=0.5, y=0.9, size=17, verticalalignment="bottom" )
    fig.savefig( os.path.join( pth, "sigma_plot" + indexString(typ) + indexString(fld) + ".pdf" ), facecolor="white" )

# generate plot of significance as a function of mdl score
def plotSignificance ( bdt, bns=0, lum=1., vrb=False, pth="./Plots", typ=None, fld=None ) :
    fig = plt.figure( figsize=( 7.5, 5. ), tight_layout=False )
    ax = fig.add_axes([ 0.11, .14, 0.65, 0.68 ])
    ax.set_zorder( 0 ); ax.set_xlim([ 0., 1. ]); ax.set_ylim([ 0., 1. ])
    ax.set_xlabel( "Signal Classification Threshold", size=14 )
    ax.xaxis.grid( True, linewidth=0.8, color="black" )
    for x in ax.get_xticklabels() : x.set_fontsize( 12 )
    tkw = dict( size=4, width=0.8 ); ax.tick_params( axis="x", **tkw )
    ax.spines["left"].set_visible(False); ax.spines["right"].set_visible(False)
    ax.yaxis.set_visible(False); ax.set_facecolor("none")
    axes = tuple( ax.twinx() for _ in range( 9 ))
    axes[5].spines["right"].set_position(( "axes", 1.2 ))
    label = ( r"$S$", r"$S\div(1+B)$", r"$\sqrt{\,2\,[\,(S+B)\,\log\,(1+S/B)-S\,]\,}$" )
    #HERE ... replace with Asimov formula here and below ... allow user selection # label = ( r"$S$", r"$S\div(1+B)$", r"$S\div\sqrt{(1+B)}$" )
    side = ( "left", "right", "right" )
    xxx = np.linspace( 0., 1., num=( 1 + 2500 ), endpoint=True ); sbi = []; sbs = []; edg = None; zzz = []
    for cls in bdt["scr"][0:2] :
        sbi.append( lum * sum( chn["xsc"] * ( 1. - chn["dns"]( xxx, -1 )) for chn in cls ))
        if int( bns ) > 0 :
            one, two = np.histogram( np.concatenate( tuple( chn["scr"] for chn in cls ), axis=0 ),
                weights=( lum * np.concatenate( tuple( chn["xsc"] * chn["wgt"] for chn in cls ), axis=0 )),
                bins=int( bns ), range=( 0., 1. ))
            for k in range( len( one ) -1, 0, -1 ) : one[ k - 1 ] += one[k]
            sbs.append( one ); edg = two
    yyy = ( sbi[1], ( sbi[1] / ( 1 + sbi[0] )), ( np.sqrt( 2. * (( sbi[0] + sbi[1] ) * np.log( 1. + ( sbi[1] / ( 0.001 + sbi[0] ))) - sbi[1] ))))
    if sbs : zzz = ( sbs[1], ( sbs[1] / ( 1 + sbs[0] )), ( np.sqrt( 2. * (( sbs[0] + sbs[1] ) * np.log( 1. + ( sbs[1] / ( 0.001 + sbs[0] ))) - sbs[1] ))))
    # yyy = ( sbi[1], ( sbi[1] / ( 1. + sbi[0] )), ( sbi[1] / np.sqrt( 1. + sbi[0] )))
    # if sbs : zzz = ( sbs[1], ( sbs[1] / ( 1. + sbs[0] )), ( sbs[1] / np.sqrt( 1. + sbs[0] )))
    for i in range( 3 ) :
        axes[i].set_frame_on(False)
        axes[i].yaxis.set_visible(False)
        axes[i].set_zorder( i-3 )
        axes[i].fill_between( xxx, yyy[i], 0., linewidth=0., alpha=0.1, edgecolor="none", facecolor=colors(i), hatch="" )
        if zzz :
            if vrb: print(( "Signal Events", "Systematics Ratio", "Statistical Significance")[i] + " in Fold " + str(fld) + ":\n" + str( zzz[i] ))
            axes[i].hist( edg[:-1], bins=edg, weights=zzz[i], color=colors(i), alpha=0.3, histtype="stepfilled" )
        axes[i].set_ylim([ 0., None ])
    for i in range( 3 ) :
        axes[3+i].set_zorder( 1+i )
        axes[3+i].set_ylim( axes[i].get_ylim())
        for sp in axes[3+i].spines.values() : sp.set_visible(False)
        axes[3+i].spines[ side[i]].set_visible(True)
        axes[3+i].spines[ side[i]].set( linewidth=0.8, color=colors(i), alpha=1. )
        axes[3+i].yaxis.set_label_position( side[i] )
        axes[3+i].yaxis.set_ticks_position( side[i] )
        axes[3+i].set_ylabel( label[i], size=14, color=colors(i) )
        axes[3+i].tick_params( axis="y", colors=colors(i), **tkw )
        axes[3+i].tick_params( axis="y", which="both", color=colors(i) ); axes[3+i].minorticks_on()
        for x in axes[3+i].get_yticklabels() : x.set_fontsize( 12 )
    for i in range( 3 ) :
        axes[6+i].set_frame_on(False)
        axes[6+i].yaxis.set_visible(False)
        axes[6+i].set_zorder( 4+i )
        axes[6+i].set_ylim( axes[i].get_ylim())
        axes[6+i].plot( xxx, yyy[i], color=colors(i), linewidth=1.4, linestyle="solid" )
    fig.suptitle( r"Event Significance at $\mathcal{{L}} = {0:G}~{{\rm {1:s}}}^{{-1}}$ for Validation Fold {2:d}".format(
        *(( lum, "pb" ) if lum < 1000. else ( lum / 1000., "fb" )), fld ) + "\n" + " ".join( tuple(
        (( lambda x : "" if x is None or len( x ) == 0 else x + " " )( label_tex[i][ cls[0]["chn"] - 1 ] )
        if len( cls ) == 1 else "Joint " ) + ( "Background vs.", "Signal" )[i] for i, cls in enumerate( bdt["scr"][0:2] ))),
        x=0.5, y=0.9, size=17, verticalalignment="center" )
    fig.savefig( os.path.join( pth, "significance_plot" + indexString(typ) + indexString(fld) + ".pdf" ), facecolor="white" )

# construct an color map of lightened shades associated with a provided base color
def shadeMap( clr, nsc=7, dup=True ) :
    ( rgb, alp, nsc ) = ( clr[0:3], clr[3] if len( clr ) > 3 else 1., 1 + max( int( nsc ), 0 ))
    hsv = np.resize( mpl.colors.rgb_to_hsv( rgb ), ( nsc, 3 ))
    hsv[:,1] = np.linspace( hsv[0][1], 0.25, nsc ); hsv[:,2] = np.linspace( hsv[0][2], 1., nsc )
    return mpl.colors.ListedColormap( tuple( np.concatenate(( _, [ alp ] ))
        for _ in mpl.colors.hsv_to_rgb( hsv )[ 0 if dup else 1 : ] ))
# https://stackoverflow.com/questions/47222585/matplotlib-generic-colormap-from-tab10

# canonical three-digit numerical identifier string in range "_000" to "_999"
def indexString( idx=None ) : return "" if idx is None else "_{0:03d}".format( min( max( 0, int(idx)), 999 ))

################
# define classes
################

# instantiate object with method for the weighted recombination of a training ensemble
class mergeMDL :
    def __init__ ( self, *bdt ) :
        self.predict = lambda dmx : (( lambda s :
            s * np.array( tuple(( lambda a,b : ( a / b ) if b else np.full( len( a ), ( 1 / ( len( a ) or 1 )), dtype=np.float64 ))( x, x.sum())
            for x in np.array( tuple( np.sum(( x["xsc"] * x["density"](s) for x in bdt[i]["density"] ), axis=0 )
            for i,s in enumerate( s )), dtype=np.float64, ndmin=2 ).transpose())).transpose())
            ( np.array( tuple( x["mdl"].predict( dmx ) for x in bdt ), dtype=np.float64, ndmin=2 ))).sum( axis=0 )

#####################################################
# construct dependent globals and invoke main routine
#####################################################

# construct associated shades of each base color
shades = tuple( shadeMap( c ) for c in colors.colors )
# revert to distinct color bases
k = 0; shades = []
for i in range( len( label_tex )) :
    shades.append( mpl.colors.ListedColormap(( colors(i), ) + tuple(
        colors( j + k + len( label_tex )) for j in range( len( label_tex[i] )))))
    k += len( label_tex[i] )

# construct feature dictionary
feature_dict = { feature_names[i]:feature_tex[i] for i in range( len( feature_names )) }

if __name__ == "__main__" : main()

sys.exit(0)

