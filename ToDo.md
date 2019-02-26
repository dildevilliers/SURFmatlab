# ToDo list for new code

## @Farfield
- [ ] Field rotations (DdV)
- [ ] Generate FarField from specified gain pattern (Brandt)
- [ ] Make empty constructor (Brandt)
- [ ] readFITS (DdV)
- [ ] readGRASPcut (DdV)
- [ ] Overlap integral calculator (DdV)
- [ ] Weighted power integral - for antenna noise (DdV)
- [ ] CBFP expansion (Fahmi)
- [ ] SWE of a given field (Fahmi + Brandt)
- [ ] Multiple frequency concat
- [ ] Shift field in 3D space - phase change
- [ ] Typical pattern parameters calculator
- [ ] writeCSTffs
- [ ] writeFEKOfft
- [ ] writeGRASPcut
- [ ] Gaussian pattern fitter
- [ ] Test the sym/pos and 180/360 plotting order rules.  Should be X and then Y shifts always - force this in the code somehow.
- [ ] Fix AzEl and ElAz poles in getELudwig2EA and getELudwig2AE: should not be 0
- [ ] Array pattern adder
- [ ] Fix 3D plot for negative y-(th)axis cases
- [ ] readMeasurements
- [ ] Jones getter
- [ ] Jones plotter (fix up)
- [ ] Stokes getter
- [ ] Stokes plotter
- [ ] Resample of SWE results on different grid
- [ ] Use above resample for interpolation when plotting
- [ ] CBFP/SWE/Zernike interpolation in freq
- [ ] General CBFP/SWE/Zernike interpolation over parameters
- [ ] Zernike expansion
- [ ] Major Rework: Change angle base to degrees and not radians (for angular grids - not DirCos type grids)
- [x] readCSTffs
- [x] readFEKOffe
- [x] Calculate power through field integration
- [x] Subtraction of fields (for error calculations and comparisons)
- [x] Normalise to isotropic
- [x] Norm of fields (for error calculations and comparisons)

## ReflectorGeometry
# General (class not assigned yet)
- [ ] Masking of a FarField by a reflector system given the feed position and orientation
- [ ] basic ray tracing of a reflector system
- [ ] full GO of a reflector system
- [ ] add offsetParaboloid class
- [ ] add offset Gregorian/Cassegrain class
- [ ] add offsetGregorianShaped class
- [ ] add symmetrical Gregorian/Cassegrain class
- [ ] add symmetrical Gregorian/Cassegrian shaped class

# @reflector
- [ ] Calculate projected area
- [ ] Calculate actual area
- [ ] Write GRASP outputs
- [ ] getMaskFunction
- [ ] add more functionality for point cloud grids (polar, thinned polar)

# @symmetricParaboloid
- [ ] calculate a path length structure
- [ ] calculate the rho-th mapping
- [ ] sort out the 2D plotting
- [ ] include a feed picture for plotting

## utils
- [x] rotation of spherical coordinates