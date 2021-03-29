function [tform,scale]=collocate_twocam(TF_PIV,TF_PLIF,grid_size,nimg,TFORMFILN)
%% function [tform,scale]=collocate_twocam(TF_PIV,TF_PLIF,grid_size,nimg,TFORMFILN)
% will take an image of a grid from two very closely aligned cameras
% and find the transform to match the two. TF_PLIF will be matched to
% TF_PIV. The code uses "detectCheckerboardPoints" which is a part of the
% Computer Vvision Toolbox. 
%
% Inputs:
% TF_PIV : directory of PIV camera images
% TF_PLIF : directory of PLIF camera images
% grid_size: grid size in mm of the grid in images
% nimg : number of images to be used
%
% Outputs: 
% tform : the translation and scaling matrices
% scale : spatial resolution of the images (in mm/px)

%% code will match D1 to D2
D1=TF_PLIF;
D2=TF_PIV;


%% begin auto rectify

cam1.files=dir(fullfile([D1,'/Core*tif']))';
cam1.filename = fullfile(D1, {cam1.files.name});

cam2.files=dir(fullfile([D2,'/Core*tif']));
cam2.filename = fullfile(D2, {cam2.files.name});


%% read images for cam 1 and 2, if less than nimg, use that number

cam1.nimg=min([nimg numel(cam1.filename)]);
cam2.nimg=min([nimg numel(cam2.filename)]);

% Detect checkerboards in images
[imagePoints, boardSize, imagesUsed] = detectCheckerboardPoints({cam1.filename{1:cam1.nimg}},'MinCornerMetric',0.15);
cam1.imagePoints=imagePoints;
cam1.boardSize=boardSize;
cam1.imagesUsed=imagesUsed;
cam1.imageFileNames={cam1.filename{1:cam1.nimg}};
clear imagePoints boardSize imagesUsed imageFileNames

% Detect checkerboards in images
[imagePoints, boardSize, imagesUsed] = detectCheckerboardPoints({cam2.filename{1:cam2.nimg}},'MinCornerMetric',0.15);
cam2.imagePoints=imagePoints;
cam2.boardSize=boardSize;
cam2.imagesUsed=imagesUsed;
cam2.imageFileNames={cam2.filename{1:cam2.nimg}};

clear imagePoints boardSize imagesUsed imageFileNames

cam2.worldPoints = generateCheckerboardPoints(cam2.boardSize, grid_size);
scale=mean(pdist(cam2.worldPoints)./pdist(cam2.imagePoints(:,:,1)));

% grab the points from the first image
c1pt=cam1.imagePoints(:,:,1);
c2pt=cam2.imagePoints(:,:,1);

% detectCheckerboardPoints is unable to tell the direction of the grid,
% so this part reshuffles the points to be in order
n=sqrt(length(c1pt));
count=1;
for i=1:n
 for j=1:n
 % disp([i j j*n-(i-1)])
  c2pt_new(count,:)=c2pt(j*n-(i-1),:);
  count=count+1;
 end
end
c2pt=flipud(c2pt_new);

% debug animation, the points should appear together in the same order
if(0)
 figure; clf;
 pause();
 for i=1:length(c1pt)
  plot(c1pt(i,1),c1pt(i,2),'xb');
  hold on;
  plot(c2pt(i,1),c2pt(i,2),'xr');
  pause(0.01);
  xlim([0 2500]); ylim([0 2500]);
 end
 title('red = c1pt, blue=c2pt');
end


%% show full effects
c1img=imread(cam1.filename{1});
c2img=imread(cam2.filename{1});
%match bobcat to photron, fixed=photron
tform = fitgeotrans(c1pt,c2pt,'NonreflectiveSimilarity');
c1new = imwarp(c1img,tform,'OutputView',imref2d(size(c2img)));

figure;
imshowpair(c1new,c2img,'Scaling','joint')
%imshowpair(c1new,c2img,'diff');

save(TFORMFILN,'tform','scale');
disp(['saved TFORM file to :',TFORMFILN]);

return
