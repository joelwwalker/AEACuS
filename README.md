# AEACuS and RHADAManTHUS
## Consumer-Level Tools for General Collider Data Analysis, Display, and Optimization
<p align="center"><img src="LOGO/Judges_420.png"/></p>

AEACuS (Algorithmic Event Arbiter and Cut Selector)
computes event statistics
and performs event selection cuts
on event files in the .lhco format.

RHADAManTHUS (Recursively Heuristic Analysis, Display,
and Manipulation: The Histogram Utility Suite)
does automated plotting
of statistics computed by AEACuS.

MInOS (Machine Intelligent Optimization of Statistics)
performs automated signal and background discrimination
using a Boosted Decision Tree (XGBoost backend).

Enter EXAMPLE directories for
five quick-start tutorials, including
sample card files for all three programs.

Please contact the author directly
( Joel Walker <jwalker@shsu.edu> )
with any questions, comments, requests,
or invitations to present a talk.

Both programs are included in this
directory, and will run without
installation or use of any external
libraries on any Linux based system
(including MacOS) with Perl 5.8+.
Windows users may also have success,
but this is not officially supported.

A manual for an older version
of AEACuS is available here:
http://arxiv.org/abs/1207.3383
Note that several newer features are
not described, and that certain aspects
of the interface have changed.  Most
critically, program invocation is now:
./aeacus.pl card_file lhco_file cross_section_pb

Slides from a talk on MInOS, AEACuS, and RHADAManTHUS at PHENO 2021
workshop are included with the package distribution
Full documentation will be forthcoming on the arXiv at a later date.

