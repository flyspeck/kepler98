# The 1998 Proof of the Kepler Conjecture

The Kepler conjecture asserts that no packing of congruent balls in Euclidean 3-space has
density greater than the familiar pyramid-shaped packing used to stack oranges at the market.

This repository contains the computer code and other documentation for the
1998 proof of the Kepler Conjecture by Sam Ferguson and Tom Hales.  This code
is not regularly maintained, but it has been deposited at github as a historical record.

## The non-computer portions of proof are given in a series of eight papers available from the arXiv.

* [An overview of the Kepler Conjecture](http://arxiv.org/abs/math/9811071)
* [A Formulation of the Kepler Conjecture](http://arxiv.org/abs/math/9811072)
* [Sphere Packings I](http://arxiv.org/abs/math/9811073)
* [Sphere Packings II](http://arxiv.org/abs/math/9811074)
* [Sphere Packings III](http://arxiv.org/abs/math/9811075)
* [Sphere Packings IV](http://arxiv.org/abs/math/9811076)
* [Sphere Packings V](http://arxiv.org/abs/math/9811077)
* [the Kepler Conjecture](http://arxiv.org/abs/math/9811078)

## Other repositories

Portions of this code have been deposited at the [arXiv](http://arxiv.org/abs/math/9811078) and with the 
[Annals of Mathematics](http://annals.math.princeton.edu/2005/162-3/p01).

## Computer portions of the proof

The computer portions of the proof fall into three main bodies of code, as reflected below in the directory
structure.

* The verification of nonlinear inequalities by [`interval`](interval) arithmetic.
* The verification of [`linear`](linear) programs.  
* The generation of all planar [`graphs`](graph) satisfying a technical list of properties.

## Directory structure

[`documents/index.html`](documents/index.html) This is the entry point for all html documentation.

[`documents`](documents) Html code, graphics, and other documentation.

[`dodec`](dodec) Materials related to Sean McLaughlin's proof of the dodecahedral conjecture.	

[`graph`](graph) A java program to generate tame planar graphs.

[`honey`](honey) Materials related to the proof of the honeycomb conjecture.

[`interval`](interval) Interval arithmetic verification of nonlinear inequalities.

[`linear`](linear) Linear programming.  This directory contains the bulk of the project code in the form
of 

[`numerical`](numerical) Code for the numerical verification of nonlinear inequalities.  This part is not
mathematically rigorous.

[`samf`](samf)  Sam Ferguson's code related to the Kepler Conjecture.

## Changes

This repository of code is essentially the same as the collection of files available from Hales's website
at the University of Michigan and the University of Pittsburgh, from around 1997-2013.  Html files have been
moved to a new documents subdirectory, updating some links.

To accommodate github's 100MB size limits on files, we have split the
zip file in [`linear/SHORT/PENT`](linear/SHORT/PENT) in two.

A more recent project, github.com/flyspeck/flyspeck, gives a formal proof of the Kepler conjecture.
