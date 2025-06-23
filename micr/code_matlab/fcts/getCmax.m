function [Cij] = getCmax(x0_all,y0_all,r0_all,eps,ROI)
%calculate the overlaying parameters of a set of circles
%
%   Cmax = getCmax(x0_all,y0_all,r0_all,eps,ROI)
%   --------------------------------------------
%
%   Inputs:
%       > x0_all,y0_all: the x,y coordinates of a set of circle centres
%       > r0_all : the radii
%       > eps : minimum accessible distance between two circles
%       > ROI : dimension of the region of interest
%               ROI = [xlim0,xlim1;
%                      ylim0,ylim1]
%
%   Outputs:
%       > Cij : the overlaying parameter between each two circles
%
% Yang CHEN, 2020.07.10



xlim0 = ROI(1,1);
xlim1 = ROI(1,2);
ylim0 = ROI(2,1);
ylim1 = ROI(2,2);
nn = length(r0_all);
Cij = zeros(nn,nn-1);
for i=1:nn
    x0 = x0_all(i);
    y0 = y0_all(i);
    r0 = r0_all(i);
    %check list
    chklst = [1:1:nn];  chklst(i) = [];
    %admissible distance
    d_ok = r0 + r0_all(chklst) + eps;
    %periodic distance (2D)
    d = dist2D_periodic_v3([x0,y0],...
                           [x0_all(chklst),y0_all(chklst)],...
                           [xlim0,xlim1;ylim0,ylim1]);
    %find the overlaying fibres
%     Cij(i,:) = (d_ok-d) ./ d;
    Cij(i,:) = (d_ok-d) ./ d_ok;
end


