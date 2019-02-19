# ToDo list for new code

## @Farfield
- [ ] Generate FarField from specified gain pattern (Brandt)
- [ ] Make empty constructor (Brandt)
- readFITS (DdV)
- readGRASPcut (DdV)
- Overlap integral calculator (DdV)
- Weighted power integral - for antenna noise (DdV)
- CBFP expansion (Fahmi)
- SWE of a given field (Fahmi + Brandt)
- Field rotations
- Multiple frequency concat
- Shift field in 3D space - phase change
- Typical pattern parameters calculator
- writeCSTffs
- writeFEKOfft
- writeGRASPcut
- Gaussian pattern fitter
- Test the sym/pos and 180/360 plotting order rules.  Should be X and then Y shifts always - force this in the code somehow.
- Fix AzEl and ElAz poles in getELudwig2EA and getELudwig2AE: should not be 0
- Array pattern adder
- Fix 3D plot for negative y-(th)axis cases
- readMeasurements
- Jones getter
- Jones plotter (fix up)
- Stokes getter
- Stokes plotter
- Resample of SWE results on different grid
- Use above resample for interpolation when plotting
- CBFP/SWE/Zernike interpolation in freq
- General CBFP/SWE/Zernike interpolation over parameters
- Zernike expansion
- Major Rework: Change angle base to degrees and not radians (for angular grids - not DirCos type grids)
~~- readCSTffs~~
~~- readFEKOffe~~
~~- Calculate power through field integration~~
~~- Subtraction of fields (for error calculations and comparisons)~~
- [x] Normalise to isotropic
- [x] Norm of fields (for error calculations and comparisons)
