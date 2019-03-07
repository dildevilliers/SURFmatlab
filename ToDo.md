# ToDo list for new code

## @Farfield
- [ ] Rename coorSys to coorType to avoid confusion with coordinateSystem class
- [ ] Speed up interpolateGrid - dont add the full grid left and right of the x-axis, only a section is required (DdV)
- [ ] testScript for projections into Coor/Grid/symmetry and plotting (3,2,1)D
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
- [ ] plot on a spherical surface
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
- [x] Field rotations (DdV)
- [x] Generate FarField from specified gain pattern (Brandt)
- [x] Make empty constructor (Brandt)
- [x] Plot Grid (DdV)


## ReflectorGeometry
# General (class not assigned yet)
- [ ] Masking of a FarField by a reflector system given the feed position and orientation
- [ ] basic ray tracing of a reflector system
- [ ] add offsetParaboloid class (William)
- [ ] add offset Gregorian/Cassegrain class
- [ ] add offsetGregorianShaped class
- [ ] add symmetrical Gregorian/Cassegrain class
- [ ] add symmetrical Gregorian/Cassegrain shaped class
- [ ] full GO of a reflector system

# @reflector
- [ ] Calculate actual area (William)
- [ ] Write GRASP outputs
- [ ] getMaskFunction
- [ ] add more functionality for point cloud grids (polar, thinned polar)

# @symmetricParaboloid
- [ ] calculate a path length structure
- [ ] include a feed picture for plotting
- [x] Calculate projected area
- [x] calculate the rho-th mapping
- [x] sort out the 2D plotting

## utils
- [x] rotation of spherical coordinates