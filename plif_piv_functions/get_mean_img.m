function bg=BackgroundImage(inputnames,outputname)
%% function bg=BackgroundImage(inputnames,outputname)
%
% Function will generate a mean image based on inputs, and will save if
% desired.
%
% Inputs:
% inputnames : list of file paths of images to be averaged
% outputname : optional argument to save the averaged image as a matlab
% file
%
% Outputs:
% bg : averaged image, if desired to be used instead of saving/loading from file\

FF0=double(imread(inputnames{1})).*0;
NF=numel(inputnames);

parfor i=1:NF
 FF0=FF0+double(imread(inputnames{i}));
end

IMG=FF0./NF;
if(nargin>1)
 save(outputname,'IMG');
 disp(['Saved file : ',outputname]);
end

bg=IMG;
