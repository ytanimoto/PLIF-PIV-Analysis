function [uout]=PIVnaninterp2(u,trsh)
%%function [uout]=PIVnaninterp1(u,trsh)
% Function will interpolate if trsh (0.5=half) of available neighbors are not NaNs.
% this is not iterative, so that the original neighbor count is honored.
% this can be run iteratively, but too many is not recommended.
% based on naninterp1 from MATPIV 1.7, rewritten YT 2020Oct19
%
% MATPIV1.7 copyright is with J.K Sveen (below)
% Copyright 1999-2001 by J.K.Sveen (jks@math.uio.no)
% Dept. of Mathematics, Mechanics Division, University of Oslo, Norway
%
% Inputs:
% u  : velocity field, 2d matrix
% trsh: fraction of neighboring points that should be present for an
% interpolation to occur
%
% Outputs:
% uout : velocity field with interpolations done

uout=u;
if nargin==2
  
end  
[py,px]=find(isnan(u)==1);
numm=size(py); [dy,dx]=size(u);

  % check number of neighbors
  for i=1:length(py)
    %correction if vector is on edge of matrix
    corx1=0; corx2=0; cory1=0; cory2=0;
    if py(i)==1, cory1=1; cory2=0;
    elseif py(i)==dy, cory1=0; cory2=-1; end
    if px(i)==1, corx1=1; corx2=0;
    elseif px(i)==dx, corx1=0; corx2=-1; end
      
    ma = u( py(i)-1+cory1:py(i)+1+cory2,...
	    px(i)-1+corx1:px(i)+1+corx2 );
    nei(i,1)=sum(~isnan(ma(:)));
    nei(i,2)=px(i); nei(i,3)=py(i);
  end
  % now sort the rows of NEI to interpolate the vectors with the
  % fewest spurious neighbors.
  nei=flipud(sortrows(nei,1));

  %locate only the NaNs with the most "true" neighbors
  ind=find(nei(:,1)>=floor(8*trsh)); 

  % only interpolate these few vectors first.
   py2=nei(ind,3); px2=nei(ind,2);

   teller=0;
   
  % main interpolation loop
  for j=1:size(py2,1)
    corx1=0; corx2=0; cory1=0; cory2=0;
    if py2(j)==1
      cory1=1; cory2=0;
    elseif py2(j)==dy
      cory1=0; cory2=-1;
    end
    if px2(j)==1
      corx1=1; corx2=0;
    elseif px2(j)==dx
      corx1=0; corx2=-1;
    end
    uout(py2(j),px2(j))=nanmean(nanmean(u(py2(j)-1+cory1:py2(j)+1+cory2,...
				     px2(j)-1+corx1:px2(j)+1+corx2)));
    teller=teller+1;
  end 
%disp(teller)
