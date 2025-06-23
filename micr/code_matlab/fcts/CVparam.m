function [T,N,fv1,fv2,fv3,vxlst,zID] = CVparam(V,theta,fvI,matID)
%compute the local parameters for composite voxels
%   [T,N,fv1,fv2,fv3,vxlst,zID] = CVparam(V,theta)
%   ----------------------------------------------
%
% Yang CHEN 2019.07.02

%recompute the distance map
[nR,nC] = size(V);
x0 = nC/2;
y0 = nR/2;
[y,x] = ndgrid(1:nR,1:nC);
d = ( tan(theta).*(x-x0) - (y-y0) ) ./ sqrt(1+tan(theta)^2);

%local orientation vectors
T = [ones(nnz(V==matID),1).*cos(theta), ...
     ones(nnz(V==matID),1).*sin(theta), ...
     zeros(nnz(V==matID),1)];
N = [ones(nnz(V==matID),1).*cos(theta+pi/2), ...
     ones(nnz(V==matID),1).*sin(theta+pi/2), ...
     zeros(nnz(V==matID),1)];

%distance criterion (the max. distance between the voxel centre and the voxel conner)
e = sqrt(2)/2*sin(pi/4+theta);

%mass centre in the voxel local coordinates
G = []

%local surface area
nbins = length(N);
xV  = ones(nbins,1)*[-0.5 0.5 0.5 -0.5 -0.5 0.5 0.5 -0.5];
yV  = ones(nbins,1)*[-0.5 -0.5 0.5 0.5 -0.5 -0.5 0.5 0.5];
zV  = ones(nbins,1)*[-0.5 -0.5 -0.5 -0.5 0.5 0.5 0.5 0.5];
S = area_element_2(N,G,xV,yV,zV);


%local volume fractions
id = d<0 & d>=e;
fv2 = (-d+e)./2*e;
fv2(V~=matID) = NaN;
fv1 = 1-fv2;
    %use constant volume fractions
    fv1 = fv1.*0 + 0.5;
    fv2 = fv2.*0 + 0.5;

fv1 = fv1(~isnan(fv1));
fv2 = fv2(~isnan(fv2));
fv3 = ones(length(fv1),1).*fvI;
fv1 = fv1.*(1-fv3);
fv2 = fv2.*(1-fv3);

%linear position of interphase voxels
vxlst = find(V==matID);

%zone ids for each phase (?) --> what's the point of these files
zID = ones(length(fv1),1);




%% visu
figure;imagesc(V,[0 5]);axis equal;axis tight
[figy,figx] = find(V==matID);
hold on;quiver(figx,figy,T(:,1),T(:,2))
[figy,figx] = find(V==matID);
hold on;quiver(figx,figy,N(:,1),N(:,2))
figure;imagesc(fv1);axis equal;axis tight;title('fv1')
figure;imagesc(fv2);axis equal;axis tight;title('fv2')