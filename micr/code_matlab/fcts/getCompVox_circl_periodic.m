function [id] = getCompVox_circl_periodic(x0,y0,x,y,r0,UC)
%identify the composite voxels of circles defined by centres and radii

% unit cell
vxsiz = UC.vxsiz;
xlim0 = UC.xlim0;
xlim1 = UC.xlim1;
ylim0 = UC.ylim0;
ylim1 = UC.ylim1;
siz = size(x);
npts = siz(1)*siz(2);

% create vectors saving the info about therelative position of each pixel corner
id_sign = zeros(npts,4,'int8');
id0 = false(npts,4); %when d==r0 (the plane exactly passing through the vertice)

% loop for each corner of a pixel
for i=1:4
    switch i
        case 1
            xV = -0.5;  yV = -0.5; 
        case 2
            xV = 0.5;  yV = -0.5;
        case 3
            xV = 0.5;  yV = 0.5;
        case 4
            xV = -0.5;  yV = 0.5;
    end
    xV = xV*vxsiz;
    yV = yV*vxsiz;
    %d = sqrt( (x(:)+xV-x0).^2 + (y(:)+yV-y0).^2 );
    d = dist2D_periodic_v3([x0,y0],[x(:)+xV,y(:)+yV],[xlim0,xlim1;ylim0,ylim1]);
    id = d>=r0;    id_sign(id,i) = 1; %d>=r0 modified by YC 22/03/2020, to take into account the "ghost" interface voxels
    id = d<r0;    id_sign(id,i) = -1;
    id0(:,i) = d==r0;
end

% identify the composite voxels
id = abs(sum(id_sign,2)) ./ sum(int8(~id0),2) < 1;