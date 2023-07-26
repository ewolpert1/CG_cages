# CG_cages
This repository contains code to simulate the packing of cages based off of the paper [Coarse-grained modelling to predict the packing of porous organic cages](https://pubs.rsc.org/en/content/articlelanding/2022/sc/d2sc04511g)

Before running the code you will need to make sure you are using [HOOMD v2.9.6](https://hoomd-blue.readthedocs.io/en/v2.9.6/) specifically with jit installed.

This code assumes there is a folder in the working directory named `output` where your output files (a log file and a gsd file) are written to.

To run the code you will need to specify the patch width, $\sigma_{ang}$, for both the window-to-window patches and window-to-arene, as well as the interaction strength, $J$, for each interaction type. For example, for simulations where the window-to-window(arene) interactions have a patch width where $\sigma_{ang}=0.3(0.4)$ and $J=30(50)$ the simulation would be run as

```
python patchy_temp_decrease.py 0.3 0.4 30 50
```

An example of output files for some simulations are given in the output directory. The form of the titles is phase_ "$\sigma_{tor}$(window-to-window)"_ "$\sigma_{tor}$(window-to-arene)"_ "$\sigma_{ang}$(window-to-window)"_ "$\sigma_{ang}$(window-to-arene)"_ "$J$(window-to-window)"_ "$J$(window-to-arene)".log/gsd

