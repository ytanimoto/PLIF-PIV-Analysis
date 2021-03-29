%-------------------------------------------------------------------------%
% Script will run a 3 pass PIV using a sub-pixel cross-correlation peak
% locating PIV algorithm developed by Liao & Cowen (2005) and based on 
% Cowen & Monismith (1997). The images are cropped and rotated, and the
% outputs of previous PIV passes are smoothed, and used as guesses for the
% next. The code is written wit matlab PARFOR to take advantages of
% computers with multiple cores to cut down on processing time.

% Note: There is an option to remove the background,
% which has been turned off in this example because there is only one image
% provided, and removing the background would remove the image.

%-------------------------------------------------------------------------%

addpath matpiv

%% parse directory to find the core number, and the experiment name


filnlist=fullfile({dir([PIVDIR,'/C*tif']).name});
corenum=get_corenum(filnlist{1});

bkdir=fullfile(mainout,'PIV_BKGD');
if(~exist(bkdir,'dir'))
 mkdir(bkdir)
end

switch PIV_FLG
 case('odd')
  offset=0;
 case('even')
  offset=1;
end

%% build background image for PIV

bkremove=0; 

if(bkremove)
 bkfiln1=[infiln,'_A.mat'];
 bkfiln2=[infiln,'_B.mat'];
 
 bk{1}=fullfile(bkdir,bkfiln1);
 bk{2}=fullfile(bkdir,bkfiln2);
 
 filnlist=fullfile(PIVDIR,{dir([PIVDIR,'/C*tif']).name});
 nimg=numel(filnlist);
 
 if(exist(bk{1},'file')~=2 ||exist(bk{2},'file')~=2  ) % if it doesn't exist, build.
  BK{1}=get_mean_img({filnlist{1:2:nimg}},bk{1});
  BK{2}=get_mean_img({filnlist{2:2:nimg}},bk{2});
 else
  for i=1:2
   tmp=load(bk{i});
   BK{i}=tmp.IMG; 
   clear tmp;
  end
 end
 
 % keep track of which background image goes with which
 bkorder=[2 1];
 if(offset)
  bkorder=sort(bkorder);
 end
  
end 

% Define region to perform PIV 

% for full profile
rmask=190:2100;
cmask=795:2010;

%% Define PIV parameters outside of the loop- this will depend on the 
% PIV code used

% PIV parameters for the first pass
pp.ny=64;   pp.nx=64;   pp.max_pass=3;
pp.method=0;
pp.Vel_Range=[-100 100 -40 40];
pp.Min_Corr_Coef=0.15;
pp.Sub_Pixel_method = 0;

% Second pass
pp2.ny=32;   pp2.nx=32;   pp2.max_pass=3;
pp2.method=0;
pp2.Vel_Range=[-100 100 -40 40];
pp2.Min_Corr_Coef=0.1;
pp2.Sub_Pixel_method = 0;

pp3.ny=32;   pp3.nx=32;   pp3.max_pass=0;
pp3.method=0;
pp3.Vel_Range=[-100 100 -40 40];
pp3.Min_Corr_Coef=0.1;
pp3.Sub_Pixel_method = 1;

% Meshgrid for the different passes
[x,y]=meshgrid(1:32:1217,1:32:1857); 
[x2,y2]=meshgrid(1:16:1217,1:16:1857);  
[x3,y3]=meshgrid(1:8:1217,1:8:1857); 

[m,n]=size(x);
[m2,n2]=size(x2);

% Estimated velocity field for the first pass
ue=zeros(m,n);
ve=zeros(m,n);

%% Preallocate arrays for PIV outputs, and parallel processing

nvec=ceil((ed-st)/2);
u_save=zeros([nvec size(x3)]);
v_save=zeros([nvec size(x3)]);
c_save=zeros([nvec size(x3)]);

x_save=x3;
y_save=y3;

%% run PIV with matlab's parallel processing.
% because images need to be processed in pairs, but the order for the
% processing of the pairs doesn't matter, so this can run in parallel

system('touch start.process');
disp(['start at ',datestr(datetime(now,'ConvertFrom','datenum'))]);
parfor i=1:1:(st-ed)+1
 
 if(mod(i,2)==1)

  A0=double(imread(fullfile(PIVDIR,sprintf('CoreView_%01d_bobcat_0_%04d.tif',corenum,i+st-1+offset))));
  B0=double(imread(fullfile(PIVDIR,sprintf('CoreView_%01d_bobcat_0_%04d.tif',corenum,i+st+offset))));
  
  if(bkremove)  
   A=A0-bk(bkorder(1));
   B=B0-bk(bkorder(2));
   end

  
  A=A0; B=B0;
  A=A.*(A>0);
  B=B.*(B>0);
 
  theta=-6;
  A=imrotate(A,theta);
  B=imrotate(B,theta);
  
  As=A(rmask,cmask);
  Bs=B(rmask,cmask);
   
  %in case anything weird came out of subtracting the background
  As=As.*(As>2); 
  Bs=Bs.*(Bs>2);
  
  % Intensity capping of Shavit et al (2007)
  ncon=2;
  medA=median(As(As~=0));
  stdA=std(As(As~=0));
  As(As>(medA+ncon*stdA))=medA+2*stdA;
  
  medB=median(Bs(Bs~=0));
  stdB=std(Bs(Bs~=0));
  Bs(Bs>(medB+ncon*stdB))=medB+2*stdB;
  
  
  
  %% MATPIV
  % this section of the code should be swapped out with whatever 
  % PIV code is available for the user
  
  [u,v,c]=mat_piv(As,Bs,x,y,ue,ve,pp);
  
  um=medfilt2(u,[4 4],'symmetric');
  vm=medfilt2(v,[4 4],'symmetric');
  %do a spatial smoothing filter to remove erroneous vectors instead of
  %feeding noise into the second refined pass
  
  % Map the result to the second pass 
  ue2=inpaint_nans(interp2(x,y,um,x2,y2));
  ve2=inpaint_nans(interp2(x,y,vm,x2,y2));
  [u2,v2,c2]=mat_piv(As,Bs,x2,y2,ue2,ve2,pp2);
  
  um2=medfilt2(u2,[4 4],'symmetric');
  vm2=medfilt2(v2,[4 4],'symmetric');
  
  ue3=inpaint_nans(interp2(x2,y2,vm2,x3,y3));
  ve3=inpaint_nans(interp2(x2,y2,vm2,x3,y3));
  
  [u3,v3,c3]=mat_piv(As,Bs,x3,y3,ue3,ve3,pp3);
  
  u_save(i,:,:)=u3;
  v_save(i,:,:)=v3;
  c_save(i,:,:)=c3;
  %clear um2 vm2 ue3 ve3 ue2 ve2 um vm u v c;

  
 end
 
end
system('touch end.process');

%% save file 

disp(['finished at ',datestr(datetime(now,'ConvertFrom','datenum'))]);

save(PIVFILN);
disp(['file saved in :',PIVFILN]);


