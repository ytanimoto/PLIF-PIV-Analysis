function build_D_FF(DDIR,DFILN,FFDIR,FFFILN,FF_FLG)
% function build_D_FF(DDIR,DFILN,FFDIR,FFFILN,FF_FLG)
%
% Function will generate a dark response and flat field image needed for
% PLIF calibration, these are the D and FF images referred to in Crimaldi &
% Koseff, 2001.
%
% Inputs:
% DDIR  : directory with dark response images
% DFILN : filename to save dark response image as
% FFDIR : directory with flat-field  images
% FFFILN: filename to save flat field images
% FF_FLG: either 'odd' or 'even' is expected, and denotes whether PLIF
% images are odd (images 1,3,5) or even. 
%
% Outputs: none

% build the dark response image, D
filnlist=fullfile(DDIR,{dir([DDIR,'/C*tif']).name});
D=get_mean_img(filnlist,DFILN);

% build flat-field image
filnlist=fullfile(FFDIR,{dir([FFDIR,'/C*tif']).name});
nfiln=numel(filnlist);

switch FF_FLG
 case('odd')
  offset=0;
 case('even')
   offset=1;
end
FF=get_mean_img({filnlist{(offset+1):2:nfiln}},FFFILN);



end