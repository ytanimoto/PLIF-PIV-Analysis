
%-------------------------------------------------------------------------%
% Script will take PIV outputs, and  collocate them with the corrected PLIF
% images to get a simultaneous density and velocity measurements. The code
% is also designed to provide the velocity and density field (x,y) in two
% coordinate systems- one that is regular and one that is alog and normal
% to a slope of 6 degrees. 
%-------------------------------------------------------------------------%

%%

showplot=1;

switch PLIF_FLG
 case('odd')
  offset=0;
 case('even')
  offset=1;
end

%% read PIV
% copy over variables from PIV output and dischard rest

PIV=load(PIVFILN);
 
PIVR.u_save=PIV.u_save;
PIVR.v_save=PIV.v_save;
PIVR.c_save=PIV.c_save;

PIVR.rmask=PIV.rmask;
PIVR.cmask=PIV.cmask;
PIVR.st=PIV.st;
PIVR.ed=PIV.ed;

X=PIV.x3;
Y=PIV.y3;

corenum=PIV.corenum;
clear PIV;


%% read image sizes 

%load one image for size
PIV0=double(imread(fullfile(PIVDIR,sprintf('CoreView_%01d_bobcat_0_%04d.tif',corenum,PIVR.st+1))));
PLIF0=double(imread(fullfile(PLIFDIR,sprintf('CoreView_%01d_redlake_0_%04d.tif',corenum,PIVR.st+1))));

pivrotsize=size(imrotate(PIV0,-6));
pivsize=size(PIV0);
plifsize=size(PLIF0);

clear PIV0 PLIF0;

%% Constants for PLIF correction

k_s=1.0850e-04;% 	 (cm g/L)^-1
k_e=1.1625e-05;% 	 (cm g/L)^-1
k_r=3.7200e-04;% 	 (cm ug/L)^-1
kt=k_r+k_s;

dens=@(sal) 9.980*10^(-1) + sal.*(6.980200*10^(-4));

%% Pre-allocate FLUX VARIBALES
imax=floor((PIVR.ed-PIVR.st)/2)+1;


FLUX.dbot=dens(C.cbot);
FLUX.dtop=dens(C.ctop);

% velocity thresholds and percentiles to calculate
FLUX.uthresh=[0:0.5:5];
FLUX.upctl=0:10:100;
FLUX.nuthresh=numel(FLUX.uthresh);
FLUX.nupctl=numel(FLUX.upctl);

% vertical coodinates
FLUX.reg_uz=cell(1,imax);
FLUX.rot_uz=cell(1,imax);

% std dev and counts
FLUX.rot_rho_std=cell(1,imax);
FLUX.rot_u_std=cell(1,imax);
FLUX.reg_rho_std=cell(1,imax);
FLUX.reg_u_std=cell(1,imax);

FLUX.rot_rho_nc=cell(1,imax);
FLUX.rot_u_nc=cell(1,imax);
FLUX.reg_rho_nc=cell(1,imax);
FLUX.reg_u_nc=cell(1,imax);

% velocity statistics
FLUX.reg_uprof_stat=cell(imax,FLUX.nuthresh,FLUX.nupctl+4);
FLUX.rot_uprof_stat=cell(imax,FLUX.nuthresh,FLUX.nupctl+4);

% density measurements
FLUX.rot_rhozall=cell(1,imax);
FLUX.rot_rhoall=cell(1,imax);

FLUX.reg_rhozall=cell(1,imax);
FLUX.reg_rhoall=cell(1,imax);

%%  Time March Loop
if(showplot)
 h=figure;
 set(h,'Position',[680         155        1095         943]);
 set(h,'Visible','on');
 set(h,'color','w');
end

for i=1:1:imax

 dt=1/7.5;
 time(i)=i*dt;
 
 if(showplot)
  clf;
 end
 
 % Read PLIF image and go through the correction process (see Crimaldi and
 % Koseff, (2001), and Troy and Koseff (2005));
 A=double(imread(fullfile(PLIFDIR,sprintf('CoreView_%01d_redlake_0_%04d.tif',corenum,PIVR.st+i*2-1+offset))));
 Cplif=(A-PLIFCAL.D)./(PLIFCAL.FF-PLIFCAL.D).*PLIFCAL.c(1).*PLIFCAL.c(2)+PLIFCAL.c(3);
 Cplif_tf=imwarp(Cplif,PLIFCAL.tform,'OutputView',imref2d(pivsize));
 
 % Due to attenuation of the laser from the dye, need to correct
 % (see Koochesfahani 1984, Tian & Roberts 2003, Odier et al. 2014)
 [s1,s2]=size(Cplif_tf);
 acorr=zeros(size(Cplif_tf));
 acorr=cumsum(Cplif_tf,1);
 acorr(2:end,:)=acorr(1:end-1,:);
 acorr(1,:)=0;
 Ccor=Cplif_tf./exp(-acorr*PLIFCAL.scale/10*kt);
 
 SAL=Ccor;
 SAL(SAL>50)=nan;
 
 % run a median filter on the PLIF outputs in case there are any particles
 SAL=medfilt2(SAL,[3 3]);
 [a,b]=size(SAL);
 xc=[1 b].*PLIFCAL.scale;
 yc=[1 a].*PLIFCAL.scale;
 
 imy=(1:a).*PLIFCAL.scale;
 imx=(1:b).*PLIFCAL.scale;
 
 if(showplot)
  cla;
  imagesc(xc,yc,SAL);hold on;
  xlim([min(PIVR.cmask) max(PIVR.cmask)].*PLIFCAL.scale);
  colormap(cmocean('haline'));
  caxis([0 30]);
  drawnow;
 end
 
 % get velocity field (x,y) for this time step
 u=squeeze(PIVR.u_save(i,:,:));
 v=squeeze(PIVR.v_save(i,:,:));
 c=squeeze(PIVR.c_save(i,:,:));

 %FLIP SIGN OF V HERE
 U=real(u); V=real(v)*-1; % 
 
 U=U.*PLIFCAL.scale./(PIV_DT*1e-3)/10; % cms
 V=V.*PLIFCAL.scale./(PIV_DT*1e-3)/10; % cms
 
 % APPLY MASK HERE TO SPEED UP PROCESSINGF AND KEEP ONLY MIDDLE COLUMNS
 ncolp=50;
 nrowp=0;
 [nrx,ncx]=size(X);
 
 p.x=X(1:nrx-nrowp,floor(ncx/2)-ncolp:floor(ncx/2)+ncolp);
 p.y=Y(1:nrx-nrowp,floor(ncx/2)-ncolp:floor(ncx/2)+ncolp);
 p.u=U(1:nrx-nrowp,floor(ncx/2)-ncolp:floor(ncx/2)+ncolp);
 p.v=V(1:nrx-nrowp,floor(ncx/2)-ncolp:floor(ncx/2)+ncolp);
 p.c=c(1:nrx-nrowp,(floor(ncx/2)-ncolp):floor(ncx/2)+ncolp);
 
 p.u=inpaint_nans(p.u,1);
 p.v=inpaint_nans(p.v,1);

 % cut off the velocities according to the correlation constant
 [p.su,p.sv]=snrfilt(p.x,p.y,p.u,p.v,p.c,0.5);
 
 % run a local filter and universal filter -  see MatPIV1.7 by Sveen
 [p.sun,p.svn]=localfilt(p.x,p.y,p.su,p.sv,5,'median');
 [p.sun2,p.svn2,info1]=unifilter(p.x,p.y,p.sun,p.svn,2);
 
 p.ufilt=PIVnaninterp2(p.sun2,0.5);
 p.vfilt=PIVnaninterp2(p.svn2,0.5);
 
 % the PIV was processed on a rotated, cropped image, so undo the crop first
 p.xp=(p.x+min(PIVR.cmask)); 
 p.yp=(p.y+min(PIVR.rmask)); 
 
 % imrotate rotates images about the center, but rotational matrices rotate 
 % about the origin (0,0), so need to rotate by removing the centroid then
 % add it back in (see lines adding in pivsize(1 and 2))
 p.xp=p.xp-ceil((pivrotsize(1)+1)/2);
 p.yp=p.yp-ceil((pivrotsize(2)+1)/2);
 
 % undo rotation with rotational matrix, vectorized
 [p.nrow,p.ncol]=size(p.xp);
 theta=-6;
 R=[cosd(theta) -sind(theta); sind(theta) cosd(theta)];
 
 p.rot_xy=R*[p.xp(:)' ; p.yp(:)'];
 p.xpr=reshape(p.rot_xy(1,:)',[p.nrow,p.ncol]);
 p.ypr=reshape(p.rot_xy(2,:)',[p.nrow,p.ncol]);
 
 p.rot_uv=R*[p.ufilt(:)' ; p.vfilt(:)'];
 p.ufiltr=reshape(p.rot_uv(1,:)',[p.nrow,p.ncol]);
 p.vfiltr=reshape(p.rot_uv(2,:)',[p.nrow,p.ncol]);

 % Add cetroid back in
 p.xpr=p.xpr+ceil((pivsize(1)+1)/2);
 p.ypr=p.ypr+ceil((pivsize(2)+1)/2);

 % need to save this grid sot hat we can make a new density grid from it
 p.xp_gc=p.xp;
 p.yp_gc=p.yp;
 
 p.xpr=p.xpr.*PLIFCAL.scale;
 p.ypr=p.ypr.*PLIFCAL.scale;

 % define new grid for interflow and map densities
 % keep rough number of columns as original velcity
 ncol_piv=numel(unique(p.xp_gc));
 nrow_piv=numel(unique(p.yp_gc));
 newx=linspace(floor(min(unique(p.xpr))),max(unique(p.xpr)),ncol_piv);
 newy=linspace(floor(min(unique(p.ypr))),max(unique(p.ypr)),nrow_piv);
 [RX,RY]=meshgrid(newx,newy);

 % griddata will nan out things outside of Xq and Yq
 RU=griddata(p.xpr,p.ypr,p.ufiltr,RX,RY,'linear'); 
 RV=griddata(p.xpr,p.ypr,p.vfiltr,RX,RY,'linear'); 
 
 % along+normal to slope corindates
 if(~exist('IMGSX','var') ||~exist('IMGSY','var') )
  [IMGX,IMGY]=meshgrid(imx,imy);
  % take row of RX, oversample and expand 8x, take col of RY, do same
  [orow,ocol]=size(p.xpr);
  IMGSX=nan(orow*8,ocol*8);
  IMGSY=nan(orow*8,ocol*8);
  
  % use linspace to make the new coordites for mapping
  for ir=1:orow
   IMGSX((ir-1)*8+1,:)=linspace(min(p.xpr(ir,:)),max(p.xpr(ir,:)),numel(p.xpr(ir,:))*8);
  end
  for ic=1:ocol
   IMGSY(:,(ic-1)*8+1)=linspace(min(p.ypr(:,ic)),max(p.ypr(:,ic)),numel(p.ypr(:,ic))*8);
  end
 end
 
 % interpolate salinity to rotated grid
 SAL_S=interp2(IMGX,IMGY,SAL,IMGSX,IMGSY,'linear',-99);
 SAL_S(SAL_S==-99)=nan;
 
 %% horizontal average
 % determine columns to average for velocity and density, on a 
 % rotated and regular grid
 xloc=195; 
 xavg=10;
 
 % rotated velocity
 [~,xp_rot_low]=min(abs((unique(p.xpr(ceil((1+end)/2),:)))-(xloc-xavg)));
 [~,xp_rot_high]=min(abs((unique(p.xpr(ceil((1+end)/2),:)))-(xloc+xavg)));

 %rotated density
 [~,im_rot_xlow]=min(abs((unique(IMGSX(ceil((1+end)/2),:))-(xloc-xavg))));
 [~,im_rot_xhigh]=min(abs((unique(IMGSX(ceil((1+end)/2),:))-(xloc+xavg))));

 % regular density 
 [~,im_reg_xlow]=min(abs(unique(imx(ceil((1+end)/2),:))-(xloc-xavg)));
 [~,im_reg_xhigh]=min(abs(unique(imx(ceil((1+end)/2),:))-(xloc+xavg)));
 
 % regular velocity
 [~,xu_reg_low]=min(abs(unique(RX(ceil((1+end)/2),:))-(xloc-xavg)));
 [~,xu_reg_high]=min(abs(unique(RX(ceil((1+end)/2),:))-(xloc+xavg)));


 if(showplot)
  % downsample using qf
  qf=4;
  
  quiver(p.xpr(1:qf:end,1:qf:end),p.ypr(1:qf:end,1:qf:end),...
    p.ufiltr(1:qf:end,1:qf:end),p.vfiltr(1:qf:end,1:qf:end),2,'Color','w');
  ylim([0 a].*PLIFCAL.scale);
  ylabel('z (mm)');
  xlabel('x (mm)');

 end
 
 %% calculate verlocity statistics
 
 % RX AND RY ARE REGULAR VELOCIES, ALONG WITH RU AND RV
 % USE WITH SAL, IMX,IMY
 
 % IMGSX AND IMGSY GO WITH SAL_S
 % ALONG WITH P.XPR P.YPR AND P.UFILTR P.VFILTR
 
 % rotated velocity 
 FLUX.rot_u_std{i}=std(p.ufiltr(:,xp_rot_low:xp_rot_high),0,2,'omitnan');
 FLUX.rot_u_nc{i}=nansum(~isnan(p.ufiltr(:,xp_rot_low:xp_rot_high)),2);
 
 rot_uprof=nanmean(p.ufiltr(:,xp_rot_low:xp_rot_high),2);
 rot_vprof=nanmean(p.vfiltr(:,xp_rot_low:xp_rot_low),2);
 z_velo_rot=nanmean(p.ypr(:,xp_rot_low:xp_rot_low),2);
 
 % regular velocity 
 FLUX.reg_u_std{i}=std(RU(:,xu_reg_low:xu_reg_high),0,2,'omitnan');
 FLUX.reg_u_nc{i}=nansum(~isnan(RU(:,xu_reg_low:xu_reg_high)),2);
 
 reg_uprof=nanmean(RU(:,xu_reg_low:xu_reg_high),2);
 reg_vprof=nanmean(RV(:,xu_reg_low:xu_reg_high),2);
 z_velo_reg=nanmean(RY(:,xu_reg_low:xu_reg_high),2);

 FLUX.reg_uz{i}=z_velo_reg;
 FLUX.rot_uz{i}=z_velo_rot;
 
 % run velocity threshold for both profiles (rotated and regular)
 % this is a quick way to look at statistics of velocites with a certain
 % lower bound threshold applied
 for ui=1:FLUX.nuthresh
  uthresh=FLUX.uthresh(ui);
  
  % rotated velocities
  uhight_rot=p.ufiltr(:,xp_rot_low:xp_rot_high);
  uhight_rot(abs(uhight_rot)<uthresh)=nan;

  FLUX.rot_uprof_stat{i,ui,1}=nanmean(uhight_rot,2);
  FLUX.rot_uprof_stat{i,ui,2}=nanmedian(uhight_rot,2);
  FLUX.rot_uprof_stat{i,ui,3}=nanmax(uhight_rot,[],2);
  FLUX.rot_uprof_stat{i,ui,4}=nanmin(uhight_rot,[],2);
  
  % regular velocities
  uhight_reg=RU(:,xu_reg_low:xu_reg_high);
  uhight_reg(abs(uhight_reg)<uthresh)=nan;

  FLUX.reg_uprof_stat{i,ui,1}=nanmean(uhight_reg,2);
  FLUX.reg_uprof_stat{i,ui,2}=nanmedian(uhight_reg,2);
  FLUX.reg_uprof_stat{i,ui,3}=nanmax(uhight_reg,[],2);
  FLUX.reg_uprof_stat{i,ui,4}=nanmin(uhight_reg,[],2);
  
  for pi=1:FLUX.nupctl
   FLUX.rot_uprof_stat{i,ui,pi+4}=prctile(uhight_rot,FLUX.upctl(pi),2);
   FLUX.reg_uprof_stat{i,ui,pi+4}=prctile(uhight_reg,FLUX.upctl(pi),2);
  end

 end
 clear uhight_rot;
 clear uhight_reg;
 
 %% calculate density statistics
  
 FLUX.reg_rhozall{i}=nanmean(IMGY(:,im_reg_xlow:im_reg_xhigh),2);
 FLUX.reg_rhoall{i}=nanmean(SAL(:,im_reg_xlow:im_reg_xhigh),2);
   
 FLUX.reg_rho_std{i}=std(SAL(:,im_reg_xlow:im_reg_xhigh),0,2,'omitnan');
 FLUX.reg_rho_nc{i}=nansum(~isnan(SAL(:,im_reg_xlow:im_reg_xhigh)),2);
 
 FLUX.rot_rhozall{i}=nanmean(IMGSY(:,im_rot_xlow:im_rot_xhigh),2);
 FLUX.rot_rhoall{i}=nanmean(SAL_S(:,im_rot_xlow:im_rot_xlow),2);
 
 FLUX.rot_rho_std{i}=std(SAL_S(:,im_rot_xlow:im_rot_xlow),0,2,'omitnan');
 FLUX.rot_rho_nc{i}=nansum(~isnan(SAL_S(:,im_rot_xlow:im_rot_xlow)),2);
 
%% 
 clear A SAL;
 clear p;
 clear u v snr pkh;
 
 if(showplot)
  T=title([infiln,sprintf('\\\\ time=%4.2f s',i*dt)]);
  set(T,'interpreter','none');
  
 end
  
 if(mod(i,100)==0)
  fprintf('%4.f percent done...',i/imax*100);
  disp(datetime(now,'ConvertFrom','datenum'))
 end
 
end



FLUX.time=time;
FLUX.scale=PLIFCAL.scale;

save(PLIFFILN,'FLUX','-v7.3');
disp(['flux file saved in :',PLIFFILN]);

return
