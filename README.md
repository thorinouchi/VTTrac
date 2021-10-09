# VTTrac
VTTrac: Velocimetry by Template Tracking
===============

Copyright 2021 (C) Takeshi Horinouchi. All rights reserved.

LICENSE: see LICENSE.txt

This library provides the core C-language implementation and its Ruby
interface. The core C library is independent from Ruby, and it should
be good for making wrappers for other languages. It does not use
module variables, so it should be good for parallel execution.

The algorithm used in this library is the simple template matching of
PIV (particle image velocimetry) for monochromatic image-like data,
but the matching is conducted multiple times in a Lagrangian manner as
in PTV (particle tracking velocimetry) over a number of times
specified by the parameter named ntrac. The default scoring method for
template matching is the cross correlation coefficient, as in the basic
PIV.

Available scoring methods
* XCOR: cross-correlation, cov(x',y')/sig(x)/sig(y), where
  x represents the template sub-image and y represents
  the target sub-image that are slid.
* NCOV: normalized covariance, cov(x',y')/sig(x)^2:
  covariance normalized the variance of the template sub-image x.

Both forward and backward tracking is available. Use the parameter
itstep; tracking is backward along time sequence, if it is negative.

A check (result screening) based on velocity change along trajectory
is available (by using the threshold parameter named vxch and vych),
so it is recommended to always set ntrac >= 2. Further screening is
available for initial templates (e.g., in terms of the complexity and
contrast) and the quality of the results (e.g., score threshold);
see the C source code.

Spatial coordinates are based on array indices, with the distance
between adjacent grid points always being 1, so they are
non-dimensional. The velocities are based on non-dimensional spacial
displacement over time difference, where, time can either be
dimensional or non-dimensional.
