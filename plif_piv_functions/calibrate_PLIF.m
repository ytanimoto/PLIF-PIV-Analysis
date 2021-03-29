function calibrate_PLIF(C,PLIFDIR,PLIF_FLG,PLIFCAL,PLIFCORFILN)
%% function calibrate_PLIF(C,PLIFDIR,PLIF_FLG,PLIFCAL,PLIFCORFILN)
%
% Function will do a PLIF calibration using the top and bottom layers of a
% stratified interface, see Troy & Koseff 2005 for details
%
% Inputs:
% C       : struct of concentrations, ctop and cbot of the salinity of the
% stratification
% PLIDFIR : directory of PLIF images
% PLIF_FLG:either 'odd' or 'even' is expected, and denotes whether PLIF
% images are odd (images 1,3,5) or even. 
% PLIFCAL : struct with D, FF. See build_D_FF.m for details 
% PLIFCORFILN : filename where the correction constants will be saved


%%
nimg=50;
filnlist=fullfile(PLIFDIR,{dir([PLIFDIR,'/C*tif']).name});
nimg=min([nimg,numel(filnlist)]);

% parse directory to read file structure
[~,infiln,~]=fileparts(filnlist{1});
corenum=str2double(char(infiln(2)));

FF=PLIFCAL.FF;
D=PLIFCAL.D;


%% build calibration image, B

switch PLIF_FLG
 case('odd')
  offset=0;
 case('even')
  offset=1;
end

% FF and D CORRECT
B=get_mean_img({filnlist{(offset+1):2:nimg}});
B=(B-PLIFCAL.D)./(PLIFCAL.FF-PLIFCAL.D);

%%

%% run calibration

xmid=1000;
BPROF=nanmean(B(:,xmid-25:xmid+25),2);
ztop=500; zbot=1500;
ZS=[ztop zbot];
BBOT=nanmean(BPROF(abs((1:numel(BPROF))-ZS(2))<50));

x0=[1 0];
fun=@(const) findcorrerr(const,C,BPROF,ZS);
xfind=fminsearch(fun,x0);


clear PLIFC;
PLIFC=[C.cbot/BBOT xfind(1) xfind(2)];

save(PLIFCORFILN,'PLIFC');
disp(['file saved in ',PLIFCORFILN]);

% debug help
if(0)
 newC=BPROF.*PLIFC(1).*PLIFC(2)+PLIFC(3);
 z=1:numel(BPROF);
 CTOP=nanmean(newC(abs(z-ZS(1))<50));
 CBOT=nanmean(newC(abs(z-ZS(2))<50));
 
 figure;
 imshow(B.*PLIFC.c1.*PLIFC.c2+PLIFC.c3);
 caxis([C.ctop C.cbot]); colormap(cmocean('haline'));
end

end % end main function

function ERR=findcorrerr(const,C,BPROF,ZS)
ztop=ZS(1);
zbot=ZS(2);

[s1,~]=size(BPROF);
z=1:s1;

BBOT=nanmean(BPROF(abs(z-zbot)<50));

CPROFC=BPROF.*(C.cbot./BBOT).*const(1)+const(2);
CTOP=nanmean(CPROFC(abs(z-ztop)<50));
CBOT=nanmean(CPROFC(abs(z-zbot)<50));

ERR=sqrt(sqrt((CTOP-C.ctop).^2)+sqrt((CBOT-C.cbot).^2));

%disp([C.cbot./BBOT const(1) const(2)]);

end
