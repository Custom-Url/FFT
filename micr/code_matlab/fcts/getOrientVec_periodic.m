function [iN,iR] = getOrientVec_periodic(x0,y0,iX,iY,r0,UC)
%determine the local orientation vector of a set of voxels describing the
%edge of a circle
% --> vectors pointing from inclusion center to the voxel center
%
%   

% dimensions of the unit cell
xlim0 = UC.xlim0;
xlim1 = UC.xlim1;
ylim0 = UC.ylim0;
ylim1 = UC.ylim1;
vxsiz = UC.vxsiz;
Lx = xlim1-xlim0;
Ly = ylim1-ylim0;

ncv = length(iX);

% absolute distance
% d = sqrt( (x(iPos)-x0).^2 + (y(iPos)-y0).^2 ); 
d = sqrt( (iX-x0).^2 + (iY-y0).^2 ); 

% outliers (the voxels separated from the centre due to the periodicity) 
id_outliers = d-r0 > vxsiz*1; %sqrt(3)/2 --> 1 to be sure

% initialisation of orientation vectors
iN = zeros(ncv,2);
iN(~id_outliers,1) = iX(~id_outliers) - x0;
iN(~id_outliers,2) = iY(~id_outliers) - y0;

% special treatments for the outliers
iN(id_outliers,1) = iX(id_outliers) - x0;
iN(id_outliers,2) = iY(id_outliers) - y0;

%if the circle centre is close to a side & far from corners
if x0-xlim0<=r0 && (y0-ylim0>r0 && ylim1-y0>r0) 
    iN(id_outliers,1) = iX(id_outliers) - (x0+Lx);
elseif y0-ylim0<=r0 && (x0-xlim0>r0 && xlim1-x0>r0)
    iN(id_outliers,2) = iY(id_outliers) - (y0+Ly);
elseif xlim1-x0<=r0 && (y0-ylim0>r0 && ylim1-y0>r0)
    iN(id_outliers,1) = iX(id_outliers) - (x0-Lx);
elseif ylim1-y0<=r0 && (x0-xlim0>r0 && xlim1-x0>r0)
    iN(id_outliers,2) = iY(id_outliers) - (y0-Ly);
    
%if the circle centre is close to a corner (more complicated)
else 
    % create vectors saving the outliers' coords and orientation vectors
    tmpX = iX(id_outliers); %X,Y coordinates of outliers
    tmpY = iY(id_outliers);
    tmpN = iN(id_outliers,:);
    
    % distinguish the outliers according to their positions
    id1 = tmpX<(xlim0+xlim1)/2 & tmpY<(ylim0+ylim1)/2;
    id2 = tmpX>(xlim0+xlim1)/2 & tmpY<(ylim0+ylim1)/2;
    id3 = tmpX<(xlim0+xlim1)/2 & tmpY>(ylim0+ylim1)/2;
    id4 = tmpX>(xlim0+xlim1)/2 & tmpY>(ylim0+ylim1)/2;
    
    % different treatments on the orientation accoring to the location of
    % circle centre
    if x0-xlim0<=r0 && y0-ylim0<=r0
        tmpN(id2,1) = tmpX(id2) - (x0+Lx);
        tmpN(id3,2) = tmpY(id3) - (y0+Ly);
        tmpN(id4,1) = tmpX(id4) - (x0+Lx);
        tmpN(id4,2) = tmpY(id4) - (y0+Ly);
    elseif xlim1-x0<=r0 && ylim1-y0<=r0
        tmpN(id1,1) = tmpX(id1) - (x0-Lx);
        tmpN(id1,2) = tmpY(id1) - (y0-Ly);
        tmpN(id2,2) = tmpY(id2) - (y0-Ly);
        tmpN(id3,1) = tmpX(id3) - (x0-Lx);
    elseif x0-xlim0<=r0 && ylim1-y0<=r0
        tmpN(id1,2) = tmpY(id1) - (y0-Ly);
        tmpN(id2,1) = tmpX(id2) - (x0+Lx);
        tmpN(id2,2) = tmpY(id2) - (y0-Ly);
        tmpN(id4,1) = tmpX(id4) - (x0+Lx);
    elseif xlim1-x0<=r0 && y0-ylim0<=r0
        tmpN(id1,1) = tmpX(id1) - (x0-Lx);
        tmpN(id3,1) = tmpX(id3) - (x0-Lx);
        tmpN(id3,2) = tmpY(id3) - (y0+Ly);
        tmpN(id4,2) = tmpY(id4) - (y0+Ly);
    end
    
    % give the treated orientation vectors back to the global vector
    iN(id_outliers,1) = tmpN(:,1);
    iN(id_outliers,2) = tmpN(:,2);
end

% normalise the vectors
iR = sqrt(sum(iN.^2,2)); %distance between inclusion center and voxel center
iN = iN ./ iR;

