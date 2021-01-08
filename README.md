# WhittThomas2013raytracing
ray tracing and numerical solution. top level script is fig10_withnumericalsolution.m <br>
WARNING: This is not a black box. Tuning is required to apply this to other problems.
<br><br>
Please cite:
<br>
Whitt, Daniel B., and Leif N. Thomas. "Near-inertial waves in strongly baroclinic currents." Journal of Physical Oceanography 43.4 (2013): 706-725.<br>
Whitt, Daniel B., et al. "Interaction of Superinertial Waves with Submesoscale Cyclonic Filaments in the North Wall of the Gulf Stream." Journal of Physical Oceanography 48.1 (2018): 81-99.

8 January 2021 - To address some questions that I received via email, here are some additional notes.

Run the example fig10_withnumericalsolution.m first. It illustrates how to use the various underlying tools.
Since the construction of the inputs for raytraceR.m is not fully provided, here are some additional comments for raytraceR.m:
- The first four inputs (Z,dz,Y,dy) are coordinates, vertical and horizontal, e.g. vertical/depth grid Z and dz is grid spacing.
- startz and starty are the starting coordinates of the ray to trace
- Lz and Ly are the length of the domain to ray trace in
- tsteps,dt are timestepping parameters for the ray tracing (number of steps, time step in seconds)
- thresh defines a threshold (small) vertical group velocity at which to internally reflect
- m0 is an initial vertical wavenumber of the wave packet (the vertical wavenumber evolves evolves through the ray tracing)
- omega is the intrinsic frequency of the wave packet (the frequency remains fixed through the ray tracing)
- F2,S2,N2 ("F squared, S squared and N squared" are as defined in Whitt and Thomas 2013)
- v0 is an initial wave velocity amplitude (e.g. 0.1 m/s), perpendicular to the mean current u (this evolves through the ray tracing)
- f is the coriolis frequency
- chstart (1 or 2) indicates which of two characteristics the initial wave packet is on; see Whitt and Thomas 2013
- s_M (=F2/S2) is the slope of absolute geostrophic momentum surfaces Mg=ug-fy
- ug is the background geostrophic flow field (function of y,z)
- various second derivatives of the background geostrophic velocity and buoyancy field: d2udy2,d2udz2,d2bdz2,d2bdy2,d2bdzdy,d2udzdy;
these get set to zero in the test case.

MIT License

Copyright (c) 2018 Daniel Whitt

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
