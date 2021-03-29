clear all; clc; close all;

addpath('plif_piv_functions');
addpath('plif_piv_scripts');

%% Check output paths and mkdir as needed
clear all;
run('output_paths');
cv=who;
for i=1:numel(cv)
 if(~exist(eval(cv{i}),'dir'))
  mkdir(eval(cv{i}))
 end
end

%% Check all paths and avariables of inputs
run('inputs');

check_list={'TF_PIV','TF_PLIF','FFDIR','DDIR','PIVDIR','PLIFDIR'};
for i=1:numel(check_list)
 if(~exist(eval(check_list{i}),'dir'))
  error(['directory not found for ',check_list{i}])
 end
end

check_var_list={'FF_FLG','PIV_FLG','PLIF_FLG','C'};
for i=1:numel(check_var_list)
 if (~exist(check_var_list{i},'var'))
  error(['variable not found : ',check_var_list{i}])
 end
end

%% Start building the pieces needed for calibration
% each of the following sections are triggered to build the file if that
% specific calibration file doesn't exist. If it exists, it will load that
% file, and you can get as specific as you want, adding version control to
% outputs, for example. This is just a basic main that is meant to be
% adapted by the user.

% global flag to rerun all processing - use with great caution!
rerun_all=0;

%% run transform for cameras to collocate points
TFORMFILN=fullfile(TFORM_SAVDIR,['TFORM_case',num2str(casenum),'.mat']);
run_tform=~exist(fullfile(TFORMFILN),'file');

if(run_tform || rerun_all)
 grid_size=9; % mm of the grid used
 nimg=1;      % number of images
 
 [tform,scale]=collocate_twocam(TF_PIV,TF_PLIF,grid_size,nimg,TFORMFILN);
 pause();
 close all;

end

load(TFORMFILN);
disp(['loaded file : ',TFORMFILN]);
PLIFCAL.tform=tform; PLIFCAL.scale=scale;
clear tform scale;

%% build D and FF images for PLIF calibration
DFILN=fullfile(D_SAVDIR,['D_case',num2str(casenum),'.mat']);
FFFILN=fullfile(FF_SAVDIR,['FF_case',num2str(casenum),'.mat']);
run_D_FF=(~exist(DFILN,'file') || ~exist(FFFILN,'file'));

if(run_D_FF || rerun_all)
 build_D_FF(DDIR,DFILN,FFDIR,FFFILN,FF_FLG)
end

D=load(DFILN);
disp(['loaded file : ',DFILN])
FF=load(FFFILN);
disp(['loaded file : ',FFFILN])

PLIFCAL.D=D.IMG;   PLIFCAL.FF=FF.IMG;

%% find PLIF correction constants
PLIFCORFILN=fullfile(PLIFCOR_SAVDIR,['PLIFCORE_case',num2str(casenum),'.mat']);
run_PLIFCORE=~exist(PLIFCORFILN,'file');

if(run_PLIFCORE || rerun_all)
 calibrate_PLIF(C,BDIR,PLIF_FLG,PLIFCAL,PLIFCORFILN);
end

PC=load(PLIFCORFILN);
disp(['loaded file : ',PLIFCORFILN])

PLIFCAL.c=PC.PLIFC;

%% run PIV
clearvars -except PLIFCAL;
run('inputs');
run('output_paths');
% the PIV processing may take a few hours, so the code usually run with a
% check to avoid overwriting an existing file

% determine starting and ending images, this is usually consistent through
% all experiments to make post-processing easy
st=696;
ed=696;

[~,infiln,~]=fileparts(PIVDIR);
fext=sprintf('%06d_%06d',st,ed);
PIVFILN=[PIVOUTDIR,'/PIV_',infiln,'__',fext,'_3p_',date,'.mat'];

runPIV=~exist(PIVFILN);

if(runPIV) 
 run('run_mat_PIV_3p');
end

%% combine with PLIF to get density and velocity fields 

clearvars -except PLIFCAL PIVFILN;
run('inputs');
run('output_paths');

load(PIVFILN,'st','ed');
[~,infiln,~]=fileparts(PIVDIR);
fext=sprintf('%06d_%06d',st,ed);
PLIFFILN=[COMBFLUXOUT,'/FLUX_',infiln,'_',fext,date,'.mat'];

runComb=~exist(PLIFFILN);
if(runComb)
 run('comb_PIV_PLIF')
end

return


