function corenum=get_corenum(filn)
%% function corenum=get_corenum(filn)
%
% Function to return the core event number from an IO industries Core device
%
% Inputs
% filn : string filename, such as 'CoreView_3_bobcat_0_0696
%
% Outputs
% corenum : event core number, 3 for the input example

 c=strsplit(filn,'_');
 corenum=str2double(char(c{2}));
end