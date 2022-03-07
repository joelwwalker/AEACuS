#!/usr/bin/env python

import sys
if ((sys.version_info[0] < 2) or ((sys.version_info[0] == 2) and (sys.version_info[1] < 6))) :
	print( "RHADAManTHUS requires Python versions 2.6, 2.7, or 3.X" )
	sys.exit(1)

import matplotlib as mpl
if (( tuple( map ( int, mpl.__version__.split("."))) + (0,0,0))[0:3] < (1,3,0)) :
	print( "RHADAManTHUS requires MatPlotLib version 1.3.0 or Greater" )
	sys.exit(1)

import warnings as wrn; wrn.filterwarnings("ignore")

import matplotlib.font_manager as mfm
font = mfm.FontProperties( fname=( list( filter( lambda x: "/times new roman.ttf" in x.lower(), mfm.findSystemFonts(fontext="ttf")))
	or [ mfm.findfont( mfm.FontProperties(family="serif")) ] )[0], size=12 )
mpl.rcParams['mathtext.fontset'] = 'stix'; mpl.rcParams['font.family'] = 'STIXGeneral'

import matplotlib.pyplot as plt
dim = 1
fig = plt.figure( figsize=(( 7.5 if dim == 1 else 6.5 ),5), tight_layout=True )
ax = fig.add_subplot(1,1,1)

import math

import numpy as np # wanted to avoid this

color = ("#e41a1c", "#377eb8", "#4daf4a", "#ff7f00", "#984ea3", "#79e6ff", "#ce0000", "#1500aa", "#107c00", "#ef6c00", "#9500af", "#00b2db")
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

hatch = 3
hchs = ((),
	( "/"*3, "\\"*3, "|"*2, "-"*2, "x"*3, "+"*2, "."*2, "O"*2 ),
	( "/", "\\", "|", "-", "x", "+", ".", "O" ),
	(),
	)
if isinstance(hatch, int): hatch = hchs[hatch]
h = len(hatch)

fill = (0,)
f = len(fill)

bold = (1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0,)
b = len(bold)

stack = False

log = True
if log: ax.set_yscale("log")

wght = [
[+0.000e+00,+9.540e-05,+9.063e-04,+2.278e-03,+3.470e-03,+5.295e-03,+5.426e-03,+6.237e-03,+6.845e-03,+7.191e-03,+7.859e-03,+7.692e-03,+7.990e-03,+8.157e-03,+7.513e-03,+6.941e-03,+7.775e-03,+6.750e-03,+6.440e-03,+5.998e-03,+5.867e-03,+5.271e-03,+4.651e-03,+5.319e-03,+4.961e-03,+4.603e-03,+4.055e-03,+3.804e-03,+3.458e-03,+3.148e-03,+3.017e-03,+2.910e-03,+2.194e-03,+2.480e-03,+2.194e-03,+2.003e-03,+2.015e-03,+1.789e-03,+1.717e-03,+1.610e-03],
[+0.000e+00,+0.000e+00,+2.784e-03,+7.532e-03,+1.801e-02,+2.966e-02,+2.779e-02,+3.777e-02,+2.327e-02,+2.029e-02,+1.330e-02,+8.337e-03,+4.651e-03,+1.533e-03,+1.083e-03,+3.386e-04,+1.419e-05,+0.000e+00,+0.000e+00,+0.000e+00,+0.000e+00,+0.000e+00,+0.000e+00,+0.000e+00,+2.027e-06,+8.110e-06,+7.299e-05,+1.419e-04,+2.433e-04,+2.798e-04,+2.413e-04,+1.257e-04,+7.096e-05,+1.419e-05,+8.110e-06,+2.027e-06,+4.055e-05,+5.210e-04,+8.576e-04,+3.325e-04],
[+0.000e+00,+0.000e+00,+0.000e+00,+1.603e-04,+7.212e-04,+1.603e-03,+1.843e-03,+1.522e-03,+1.522e-03,+8.013e-04,+1.042e-03,+7.212e-04,+5.609e-04,+7.212e-04,+4.006e-04,+1.122e-03,+5.449e-03,+5.649e-02,+1.064e-01,+1.771e-02,+1.122e-03,+0.000e+00,+0.000e+00,+8.013e-05,+0.000e+00,+0.000e+00,+0.000e+00,+0.000e+00,+0.000e+00,+0.000e+00,+0.000e+00,+0.000e+00,+0.000e+00,+0.000e+00,+0.000e+00,+0.000e+00,+0.000e+00,+0.000e+00,+0.000e+00,+0.000e+00],
[+0.000e+00,+6.526e-06,+6.526e-05,+1.468e-04,+3.981e-04,+6.069e-04,+8.157e-04,+8.875e-04,+9.169e-04,+8.190e-04,+8.059e-04,+5.743e-04,+7.276e-04,+7.602e-04,+8.908e-04,+1.361e-03,+6.186e-03,+5.833e-02,+1.088e-01,+1.501e-02,+1.132e-03,+1.436e-04,+3.589e-05,+5.221e-05,+2.284e-05,+4.242e-05,+2.610e-05,+3.915e-05,+1.631e-05,+1.958e-05,+2.610e-05,+2.610e-05,+1.631e-05,+1.631e-05,+3.915e-05,+2.937e-05,+2.610e-05,+1.958e-05,+1.631e-05,+6.526e-06],
[+0.000e+00,+5.852e-06,+2.282e-04,+5.618e-04,+9.480e-04,+1.311e-03,+1.715e-03,+2.013e-03,+1.855e-03,+2.089e-03,+1.902e-03,+2.101e-03,+1.978e-03,+1.843e-03,+2.130e-03,+2.657e-03,+6.660e-03,+5.183e-02,+9.706e-02,+1.361e-02,+1.498e-03,+4.038e-04,+4.682e-04,+2.692e-04,+3.628e-04,+2.458e-04,+2.692e-04,+2.692e-04,+3.394e-04,+1.873e-04,+2.048e-04,+1.990e-04,+1.639e-04,+2.165e-04,+1.814e-04,+1.229e-04,+1.112e-04,+1.463e-04,+1.580e-04,+9.363e-05],
[+0.000e+00,+1.157e-04,+8.724e-04,+2.332e-03,+4.095e-03,+5.350e-03,+6.410e-03,+7.060e-03,+8.306e-03,+8.324e-03,+8.671e-03,+8.591e-03,+8.822e-03,+8.288e-03,+8.155e-03,+7.336e-03,+7.175e-03,+6.890e-03,+6.979e-03,+5.973e-03,+5.751e-03,+5.057e-03,+4.861e-03,+4.807e-03,+4.157e-03,+3.979e-03,+3.383e-03,+3.196e-03,+3.098e-03,+2.715e-03,+2.680e-03,+2.573e-03,+2.350e-03,+2.012e-03,+2.154e-03,+1.807e-03,+1.451e-03,+1.513e-03,+1.415e-03,+1.433e-03],
[+0.000e+00,+7.596e-04,+2.324e-03,+5.541e-03,+1.452e-02,+2.033e-02,+2.542e-02,+2.668e-02,+2.046e-02,+1.729e-02,+1.318e-02,+1.077e-02,+9.115e-03,+6.122e-03,+4.960e-03,+3.932e-03,+3.262e-03,+2.234e-03,+2.234e-03,+1.877e-03,+1.475e-03,+1.251e-03,+1.206e-03,+7.596e-04,+8.937e-04,+4.021e-04,+4.915e-04,+4.021e-04,+4.468e-04,+2.681e-04,+1.787e-04,+1.340e-04,+1.340e-04,+2.234e-04,+4.468e-05,+1.340e-04,+8.937e-05,+8.937e-05,+1.340e-04,+0.000e+00],
[+0.000e+00,+7.063e-04,+2.961e-03,+5.216e-03,+6.466e-03,+9.128e-03,+1.144e-02,+1.250e-02,+1.361e-02,+1.299e-02,+1.407e-02,+1.290e-02,+1.130e-02,+1.035e-02,+8.612e-03,+8.340e-03,+6.574e-03,+6.248e-03,+5.162e-03,+5.080e-03,+4.374e-03,+3.423e-03,+3.776e-03,+2.662e-03,+2.363e-03,+1.739e-03,+2.065e-03,+1.711e-03,+1.739e-03,+1.358e-03,+1.005e-03,+1.250e-03,+9.237e-04,+6.520e-04,+7.878e-04,+7.878e-04,+5.433e-04,+4.347e-04,+4.890e-04,+3.260e-04],
[+0.000e+00,+5.115e-04,+2.491e-03,+3.691e-03,+5.048e-03,+6.604e-03,+7.494e-03,+8.450e-03,+8.962e-03,+8.939e-03,+9.984e-03,+9.962e-03,+9.384e-03,+9.273e-03,+8.895e-03,+8.450e-03,+7.650e-03,+7.538e-03,+6.048e-03,+6.226e-03,+6.404e-03,+5.226e-03,+4.759e-03,+4.447e-03,+3.736e-03,+3.447e-03,+3.291e-03,+2.935e-03,+3.180e-03,+2.491e-03,+2.424e-03,+1.957e-03,+1.601e-03,+1.757e-03,+1.223e-03,+1.112e-03,+1.112e-03,+8.450e-04,+1.290e-03,+6.893e-04],
[+0.000e+00,+9.798e-05,+9.406e-04,+2.606e-03,+3.958e-03,+4.507e-03,+5.546e-03,+6.506e-03,+5.703e-03,+6.722e-03,+7.015e-03,+6.898e-03,+7.564e-03,+6.898e-03,+7.603e-03,+7.721e-03,+7.153e-03,+7.015e-03,+6.349e-03,+6.349e-03,+6.290e-03,+6.232e-03,+5.017e-03,+5.193e-03,+4.860e-03,+5.017e-03,+3.665e-03,+3.880e-03,+3.135e-03,+3.508e-03,+3.233e-03,+2.665e-03,+2.097e-03,+2.802e-03,+2.685e-03,+2.234e-03,+2.097e-03,+1.862e-03,+1.881e-03,+1.333e-03],
[+0.000e+00,+3.698e-05,+7.766e-04,+2.052e-03,+2.958e-03,+3.587e-03,+4.955e-03,+5.436e-03,+4.992e-03,+5.473e-03,+5.750e-03,+5.528e-03,+5.972e-03,+6.268e-03,+6.897e-03,+6.675e-03,+7.026e-03,+5.935e-03,+6.046e-03,+5.880e-03,+6.138e-03,+5.362e-03,+5.806e-03,+5.306e-03,+5.029e-03,+4.752e-03,+4.123e-03,+4.419e-03,+3.957e-03,+3.494e-03,+3.273e-03,+3.291e-03,+3.550e-03,+3.310e-03,+2.810e-03,+2.552e-03,+2.441e-03,+2.237e-03,+2.071e-03,+1.886e-03],
[+0.000e+00,+2.919e-05,+2.482e-04,+1.168e-03,+2.233e-03,+2.963e-03,+3.795e-03,+4.087e-03,+4.905e-03,+4.423e-03,+5.197e-03,+5.124e-03,+5.445e-03,+5.781e-03,+5.649e-03,+5.022e-03,+5.985e-03,+6.043e-03,+5.343e-03,+5.839e-03,+5.036e-03,+5.693e-03,+4.861e-03,+5.051e-03,+4.773e-03,+4.832e-03,+4.890e-03,+4.642e-03,+3.883e-03,+4.014e-03,+4.408e-03,+3.445e-03,+3.518e-03,+3.328e-03,+3.036e-03,+2.890e-03,+2.701e-03,+2.978e-03,+2.496e-03,+2.219e-03]]
n = len(wght)

bins = [+0.000e+00,+5.000e+00,+1.000e+01,+1.500e+01,+2.000e+01,+2.500e+01,+3.000e+01,+3.500e+01,+4.000e+01,+4.500e+01,+5.000e+01,+5.500e+01,+6.000e+01,+6.500e+01,+7.000e+01,+7.500e+01,+8.000e+01,+8.500e+01,+9.000e+01,+9.500e+01,+1.000e+02,+1.050e+02,+1.100e+02,+1.150e+02,+1.200e+02,+1.250e+02,+1.300e+02,+1.350e+02,+1.400e+02,+1.450e+02,+1.500e+02,+1.550e+02,+1.600e+02,+1.650e+02,+1.700e+02,+1.750e+02,+1.800e+02,+1.850e+02,+1.900e+02,+1.950e+02,+2.000e+02]

(ymin, ymax) = (0.0001, 0.5)

ttl = [r"$\sqrt{s} = 14$ TeV",r""]

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
		if fill[i%f]: ax.fill_between(edges, upper, lower, linewidth=0.0, alpha=0.6, edgecolor="none", color=color[i%c], hatch="", zorder=i-3*n )
		if h: ax.fill_between(edges, upper, lower, linewidth=0.0, alpha=1.0, edgecolor="0.25", color="none", hatch=hatch[i%h], zorder=i-2*n )
		if bold[i%b]: ax.plot( edges, upper, color="#dddddd", linewidth=2.0, linestyle="solid", zorder=i-1*n )
		ax.plot( edges, upper, color=color[i%c], linewidth=(2.0,2.0)[bold[i%b]], linestyle=("solid","--")[bold[i%b]], zorder=i-1*n )
		if stack: lower = upper

	lgd = None
	lgnd = [r"$t\bar{t}jj$",r"$\tau\tau jj$",r"Zjjjj",r"$ZZjj$",r"$WZjj$",r"$WWjj$",r"$S_{10}^{110}$",r"$S_{20}^{110}$",r"$S_{30}^{110}$",r"$S_{40}^{110}$",r"$S_{50}^{110}$",r"$S_{60}^{110}$"]
	if lgnd:
		patches = [(
			( mpl.patches.Rectangle((0,0), 1, 1, fill=True, facecolor=color[i%c], alpha=0.6, linewidth=0.0 ) if fill[i%f] else ()),
			( mpl.patches.Rectangle((0,0), 1, 1, fill=None, color="0.25", hatch=hatch[i%h], linewidth=0.0 ) if h else ()),
			( mpl.patches.Rectangle((0,0), 1, 1, fill=None, linestyle="solid", edgecolor="#dddddd", linewidth=1.4 ) if bold[i%b] else ()),
			mpl.patches.Rectangle((0,0), 1, 1, fill=None, linestyle=("solid","--")[bold[i%b]], edgecolor=color[i%c], linewidth=1.4 ))
			for i in range(n) ]
		lgd = ax.legend( patches, lgnd, loc="center right", bbox_to_anchor=(1.22, 0.5), prop=font )
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

	ax.set_xlabel( r"$M_{\ell \ell}$ [GeV]", fontproperties=font, size=14 )
	ax.set_ylabel( r"Event Fraction per GeV  ($d\sigma/d E \div \sigma$)", fontproperties=font, size=14 )
	ax.set_title( ttl[0], fontproperties=font, size=17, verticalalignment="bottom" )

	fig.savefig( "./shape_MAS_001_n1.pdf", facecolor="white", bbox_extra_artists=(lgd,), bbox_inches="tight")

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

	ax.set_xlabel( r"$M_{\ell \ell}$ [GeV]", fontproperties=font, size=16 )
	ax.set_ylabel( r"Event Fraction per GeV  ($d\sigma/d E \div \sigma$)", fontproperties=font, size=16 )
	ax.set_title( ttl[0], fontproperties=font, size=19, verticalalignment="bottom" )

	fig.savefig( "./shape_MAS_001_n1.pdf", facecolor="white" )

sys.exit(0)

