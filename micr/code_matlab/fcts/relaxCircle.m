function [x0_new,y0_new] = relaxCircle(x0,y0,Cij,x1,y1,ind,ROI)
%calculate the new position of a circle centre by relaxation (repulsive
%force)
%
%   

% unit-cell dimensions
xlim0 = ROI(1,1);
xlim1 = ROI(1,2);
ylim0 = ROI(2,1);
ylim1 = ROI(2,2);
Lx = xlim1 - xlim0;
Ly = ylim1 - ylim0;

%calculate the repulsive force exerting on the considered fibre
ni = length(x1);
Nx = (x0-x1) ./ (1-Cij) .* Cij;
Ny = (y0-y1) ./ (1-Cij) .* Cij;
% Nx = (1/ni) * sum(x1+(x0-x1)./(1-Cij));
% Ny = (1/ni) * sum(y1+(y0-y1)./(1-Cij));


%if periodicity used, then the vector should be opposite
id = ind==1 | ind==2 | ind>=5;
Nx(id) = Nx(id) * (-1);
id = ind==3 | ind==4 | ind>=5;
Ny(id) = Ny(id) * (-1);

%
x0_new = x0 + sum(Nx)/3;
y0_new = y0 + sum(Ny)/3;
% x0_new = x0 + sum(Nx);
% y0_new = y0 + sum(Ny);

%check if the new positions are inside the domain (periodicity considered)
if (x0_new<xlim0); x0_new = x0_new + Lx; end
if (y0_new<ylim0); y0_new = y0_new + Ly; end
if (x0_new>xlim1); x0_new = x0_new - Lx; end
if (y0_new>ylim1); y0_new = y0_new - Ly; end