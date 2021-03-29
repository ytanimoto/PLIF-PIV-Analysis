% declare the files/folders where the images are
% this file can be replaced with a longer file/and or loop to automate 
% the analysis to a sequence of experiments

% 0. Give a unique ID to this run 
casenum=1;

% 1. Transform (TF) images, which are images to spatially
% calibrate and join the two cameras. this version uses a 2D grid
TF_PIV='Example_images/TF_PIV';
TF_PLIF='Example_images/TF_PLIF';

% 2. PLIF calibration images (see PLIF calibration notebook for more details!)
% FF_FLG = the PIV/PLIF system gets alternarting PLIF images, so need to
% specify whether odd (img 1,3,5), or even, should be used
% FFDIR = flat field image
% DDIR = dark response image
% B  = calibration image, with the intial stratification shown. This is
% usually the first n images of an experiment, but here I've included it
% for reference
FF_FLG='even'; 
FFDIR='Example_images/PLIF_CAL/FF';
DDIR='Example_images/PLIF_CAL/D';
BDIR='Example_images/PLIF_CAL/B';

% 3. PIV/PLIF experimental images
PIV_FLG='odd';
PLIF_FLG='even';
PIVDIR='Example_images/PIV_EXP';
PLIFDIR='Example_images/PLIF_EXP';

%4. top/bottom concentrations of dye in ppb 
C.ctop=3.33;
C.cbot=18.57;

%5. PIV properties
PIV_DT=9.5; % miliseconds
