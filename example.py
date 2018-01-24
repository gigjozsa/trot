#! /usr/bin/env python
import trot
import numpy as np

indeffile = 'in.def'
outdeffile = 'out.def'
outcube = 'cube.fits'
outdefprefix = 'out'
outcubeprefix = 'cube'
tirexec = 'tirific'

rot = trot.Trot(seed = 1, tirexec = tirexec)
    
# Just one well-defined rotation, take input def, rotate about z-axis by alpha, then about x-axis by beta, then about z-axis by gamma. Here we only choose gamma, to see if only position angle changes in output.indeffile is the input deffile, outdeffile is the output deffile, outcube is the output cube.
alpha = 0.
beta = 0.
gamma = np.pi*90./180.
rot.rotatesingle(indeffile, outdeffile, outcube, alpha, beta, gamma)

# 10 random rotations with actual generation of cubes
rot.rotatemulti(indeffile, outdefprefix, outcubeprefix, 5)
