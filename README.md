# ROINets_Yokogawa

This repository will hold scripts to compute MEG source-space connectivity between ROIs using data from the Yokogawa MEG160 system. The code is heavily based on MEG ROINets from OHBA, Oxford, UK; but uses the Fieldtrip toolbox for source analysis (LCMV beamformer).

![Imgur](https://i.imgur.com/LjxbKFF.png)

ROInets: A Matlab package for performing leakage-robust network inference 
between ROIs in MEG data

The methodology used in this pipeline is set out in

Colclough, G. L., Brookes, M., Smith, S. M. and Woolrich, M. W., 
"A symmetric multivariate leakage correction for MEG connectomes," 
NeuroImage 117, pp. 439-448 (2015)

We are grateful for Colclough and colleagues at OHBA for making their
code openly available for re-use. Please check copyright information in 
OSL for more information.

These scripts use the MEG ROINets code, but performs source analysis in
Fieldtrip.

Users will need to download OSL (https://github.com/OHBA-analysis).
However please don't initiate OSL - this script will add the relevent
bits of code for you.
