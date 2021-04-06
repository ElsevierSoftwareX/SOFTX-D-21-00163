# SWE Runup
<hr/>
Authors:

Latest Update:
<hr/>

### Introduction
This software provides a fast, direct comparison of a 1-D general finite volume solution to the 1+1 shallow water equations (SWE) with a 1-D analytical solution presented by Nicolsky et al.

### Installation

Download `swe_runup` [here](https://github.com/twbf/swe_runup) and note which directory the folder is saved to.

For users operating within MATLAB gui:
1. Navigate to `swe_runup` in current folder window.
2. Right click on `swe_runup` and select `Add to Path --> Selected Folers and Subfolders`.

For users operating MATLAB from command window:
1. Locate directory in which `swe_runup` is located.
2. Input
```Matlab
addpath('Path/to/swe_runup/')
savepath
```

Users must also download Chebfun. To do so, input
```bash
unzip('https://github.com/chebfun/chebfun/archive/master.zip')
movefile('chebfun-master', 'chebfun'), addpath(fullfile(cd,'chebfun')), savepath
```
into bash command line or follow this [link](https://www.chebfun.org) for direct download.

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
- [Hankel Transform](https://www.mathworks.com/matlabcentral/fileexchange/13371-hankel-transform): Computes fast/inverse fast hankel transform.
- [Chebfun](https://www.chebfun.org/download/)

### License
This package is free and open source. See [LGPL-3.0](https://opensource.org/licenses/LGPL-3.0) for more licensing information.

### More information
A full account of software functionalities and implemented methods can be found [here](link-to-paper) (link to paper). If you have any problems or questions contact (email).
