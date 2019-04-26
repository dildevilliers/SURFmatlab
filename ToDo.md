# ToDo list for new code

## @Farfield
- [ ] readMeasurements
- [ ] Check readGRASPgrd. Not sure of E1 and E2 order for all cases. Preallocate the matrices for speed like in readGRASPcut
- [ ] Fix FarField.rotate field components poles (DdV)
- [ ] Field symmetries: XY plane (DdV)
- [ ] Gaussian/cosn pattern fitter
- [ ] Rework FarFieldFromPower pattern workflow.  Can use BOR1 functions here to shorten the code. Merge powerPattern in...
- [ ] readFITS (DdV)
- [ ] Overlap integral calculator (DdV)
- [ ] CBFP expansion (Fahmi)
- [ ] SWE of a given field (Fahmi + Brandt)
- [ ] Multiple frequency concat
- [ ] Typical pattern parameters calculator: SLL, XP, Beamwidth, etc.
- [ ] writeCSTffs
- [ ] writeFEKOfft
- [ ] Test the sym/pos and 180/360 plotting order rules.  Should be X and then Y shifts always - force this in the code somehow.
- [ ] Fix AzEl and ElAz poles in getELudwig2EA and getELudwig2AE: should not be 0
- [ ] Array pattern adder
- [ ] plot on a spherical surface
- [ ] Fix 3D plot for negative y-(th)axis cases
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
- [x] Weighted power integral - for antenna noise (DdV)
- [x] Fix FarField th rotation direction (DdV)  
- [x] Field symmetries: XZ and YZ planes.
- [x] Speed up interpolateGrid - dont add the full grid left and right of the x-axis, only a section is required (DdV)
- [x] testScript for projections into Coor/Grid/symmetry and plotting (3,2,1)D (Ridalise and William)
- [x] Get only the BOR1 components of the field and return as FarField object (DdV)
- [x] Conjugate fields overload
- [x] Shift field in 3D space - phase change
- [x] Made a FarField.rms function to check RMS field values over angle/freq
- [x] readGRASPcut (DdV)
- [x] writeGRASPcut
- [x] setFreq should not be private. Use to set both freq and freqHz, with varargin for backwards compatibility.
- [x] Rename coorSys to coorType to avoid confusion with coordinateSystem class
- [x] Move all parameters to setAccess private, and make setters as required


# ReflectorGeometry
## General (class not assigned yet)
- [ ] Change all class names to be Capital letter first 
- [ ] Antenna temperature calculator
- [ ] include a feed picture for plotting
- [ ] include a detail level selector in the plotting of single/dualReflector classes. Level 1 only dish, 2 rays, 3 points, etc.
- [ ] add numeric surface class (DdV)
- [ ] add numeric rim class (DdV)
- [ ] add offsetGregorianShaped class
- [ ] add symmetrical Gregorian/Cassegrain shaped class
- [ ] full GO of a reflector system
- [x] add hyperboloid class (DdV)
- [x] add offsetParaboloid class (William)
- [x] add symmetrical Gregorian/Cassegrain class
- [x] add offset Gregorian/Cassegrain class
- [x] Aperture efficiency calculator

## coordinateSystem
- [x] make x_axis and y_axis private SetAccess properties 

## @reflector
- [ ] Calculate actual area (William)
- [ ] Write GRASP outputs (William)
- [x] getMaskFunction
- [x] add more functionality for point cloud grids (polar, thinned polar)
- [x] basic ray tracing of a reflector system

## @hyperboloid
- [ ] Sort out the negative e concave case like in GRASP

## @ellipsoid
- [x] Make rotation as in GRASP possible

## @dualReflector
- [x] SR extensions
- [x] masking (both SR and PR masking)
- [x] ray tracing
- [x] path length structure
- [x] ray trace plots
- [x] constructor function for different symmetrical design options
- [x] constructor function for different offset design options
- [x] Re-implement the Gregorian case completely. Do as in legacy work with a rotated ellipsoid. Fixed legacy issues with different extended and non-extended OG systems 
- [x] Fixed Bug: testScript_dualReflector, exNumber = 1|3; th_ext = 20; The DR.SR.surface.F1 point is not plotted correctly for Cassegrains
- [x] Fixed Bug: testScript_dualReflector, exNumber = 1; th_ext = 20 or 10 deg; symFact_ext = 1. Inner SR edge (Q1) wrong...

## @pnt3D
- [x] plot line between 2 points
- [x] plot a ray from a point in a given direction

# utils
- [x] rotation of spherical coordinates

# Arrays
## General
- [ ] Make all classes frequency dependent

## PlaneWaveSignal
- [ ] freq, sigPower, sigPhase frequency dependent
- [ ] include polarization information
- [ ] include plot for the direction and polarization vectors

## ArrayElements
- [ ] channelErrors frequency dependent
- [ ] different field patterns for each element
- [ ] include plot for the element positions

## ArrayReceiver
- [ ] Full noise coupling matrix
- [ ] LNA noise parameters specification method
- [ ] Gains frequency dependent
- [ ] More receiver (amp) characteristics - IP1/3, 3dB comp, etc

## ArrayDAC
- [ ] Signed binary
- [ ] 2s complement

## ArraySystem
- [ ] Expand plotting option for line styles etc.
- [ ] Automate the time scale to Engineering units

## ArrayBeamformer
- [ ] Implement for analog arrays
- [ ] Include transmit array pattern calculator

## ArrayDBE
- [ ] Sort out power levels in the PSD plots (Units)
- [ ] plotScanBeam must be extended to handle 1D and 2D scan plots
