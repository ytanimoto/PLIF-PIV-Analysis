# PLIF-PIV-Analysis

## Motivation
Simultaneous density and velocity measurements are key in understanding any stratified fluid flow. Gettingm full field (x,y) resolved measurements gives much more insight than from a singe ADV, especially if your flow has spatial gradients. As a graduate student I found very little information about how to actually sync up the PIV and PLIF measurements, so I'm hoping that this might help others who are thinking about implementing a similar system. This is a the second part of the measurement system (like [this one](https://github.com/ytanimoto/PLIF-PIV-Trigger) for example), walking through the steps to get the images into real data!

## Requirements
The scripts utilize functions in the Matlab Computer Vision Toolbox, as well as MATLAB's parallel processing tool. The code has been tested only on Matlab 9.8.0.1417392 (R2020a) Update 4. The user will also need to choose their own PIV code, such as MATPIV 1.7 by J.K. Sveen (a version such as https://github.com/alexlib/matpiv will work). 

## How should I use this?
The workflow in `main.m` is meant to be used for each experimental set. Naturally you'll want to write an external loop to process different experiments sequentially, but the processing steps are the same for each experiment! 

## What it does
* specifies what files and inputs are needed 
* prepare the output folders 
* use images from two cameras to find out how to match the two
* build the dark response, flat-field image, and calibrate the PLIF
* run the PIV on the images to get velocity fields
* uses the PIV outputs and collocates them with the corrected PLIF fields to write out collocated, simultaneous density and velocity fields.

