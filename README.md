# SWANS

<hr/>
Authors: Thomas Beuler (1), Sam Delamere (1), Denys Dutykh (2), Alexei Rybkin (1), Alex Suleimani (1).
<br/>
<br/>
Affiliations: (1) Columbia University, (2) Bates College, (3) University of Savoie, (4) University of Alaska Fairbanks, (5) Arizona State University.
<hr/>

### Introduction
SWANS (Shallow Water Analytical Numerical Solver) provides a fast, direct comparison of a 1-D general finite volume solution to the 1+1 shallow water equations (SWE) with a 1-D analytical solution presented by Nicolsky et al (2018).

### Installation

Users must download Chebfun. To do so, go to the Chebfun github repo ([here](https://github.com/chebfun/chebfun)). Follow the instructions in the Chebfun readme to install.

Note: If running this code through a command line interface, make sure chebfun is in your MATLAB path.

Download `swe_runup` from github [here](https://github.com/twbf/swe_runup) and note which directory the folder is saved to.

For users operating within MATLAB gui:
1. Navigate to `swe_runup`.
2. Open `run_swe_runup` and run script from editor.
3. Alternatively, if you are not running from the `swe_runup` folder, then you must add `swe_runup` to path from the gui. This can be done by right clicking on the folder and selecting `Add to Path`.

For users operating MATLAB from terminal window (macOS/linux):
1. Add MATLAB to your PATH.
2. Change directory path/to/swe_runup
3. matlab -nosplash -nodesktop
4. From the MATLAB CLI, type run_swe_runup  

### Usage

All simulations are produced through the `run_swe_runup` executable. All wave/simulation parameters are also located in this script. These include
- Simulation resolution (`t_res` and `x_res`)
- Plots (toggled 'on' or 'off')
- Order of data projection (`n`)
- Spatial and temporal domains (`x0`,`Xf`,`t0`,`Tf`)
- Beach slope (`td`)
- Initial wave function parameters (`H1`,`H2`,`c2`,`c2`,`x1`,`x2`).
- Initial wave functions (`eta_0` and `u_0`)

Saving figures can be done through MATLAB's gui.

### Third party packages.
Third party files included:
- [genpath2](https://www.mathworks.com/matlabcentral/fileexchange/72791-genpath2): Adds all subfolders to path (excluding hidden files with `.git` extension).
- [Hankel Transform](https://www.mathworks.com/matlabcentral/fileexchange/13371-hankel-transform): Computes fast/inverse fast Hankel transform.
- [Chebfun](https://www.chebfun.org/download/)

### License
This package is free and open source. See [LGPL-3.0](https://opensource.org/licenses/LGPL-3.0) for more licensing information.

### More information
A full account of software functionalities and implemented methods can be found [here](link-to-paper) (link to paper).

Please file any comments or concerns through the github issue tracker [here](/https://github.com/twbf/swe_runup/issues).

### References
Nicolsky, Dmitry, et al. "General initial value problem for the nonlinear shallow water equations: Runup of long waves on sloping beaches and bays." Physics Letters A 382.38 (2018): 2738-2743.
