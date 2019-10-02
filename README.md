# Resting-State Analysis @ Macquarie University

A collection of scripts to analyse "resting-state" and/or naturalistic viewing MEG data acquired from the KIT-Macquarie Research Laboratory

## Source Analysis:
Scripts will include instructions to perform:
- Data cleaning
- Beamforming
- Parcellation

## Stationary Whole-brain Connectivity

This repository will hold scripts to compute MEG source-space connectivity between ROIs using data from the Yokogawa MEG160 system. The code is heavily based on MEG ROINets from OHBA, Oxford, UK; but uses the Fieldtrip toolbox for source analysis (LCMV beamformer).

![Imgur](https://i.imgur.com/LjxbKFF.png)

ROInets: A Matlab package for performing leakage-robust network inference 
between ROIs in MEG data

The methodology used in this pipeline is set out in

Colclough, G. L., Brookes, M., Smith, S. M. and Woolrich, M. W., 
"A symmetric multivariate leakage correction for MEG connectomes," 
NeuroImage 117, pp. 439-448 (2015)

## Non-stationary Connectivity using Hidden Markov Models (HMM)


