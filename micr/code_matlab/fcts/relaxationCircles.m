function [x0_all,y0_all,lstRelaxed] = relaxationCircles(x0_all,y0_all,r0_all,id,eps,ROI)
%re-arrange the overlying circles by "relaxation"
%
%   [x0_all,y0_all,lstRelaxed] = relaxationCircles(x0_all,y0_all,r0_all,eps,ROI)
%   ----------------------------------------------------------------------------
%
%   Inputs:
%       > x0_all,y0_all : initial positions of the circle centres
%       > r0_all : radii of the circles
%       > id : indices indicating whether the ith fibre is unstable (overlaying)
%       > eps : minimum admissible distance between two circles
%       > ROI : dimension of the unit cell
%               ROI = [xlim0,xlim1;
%                      ylim0,ylim1]
%
%   Ouputs:
%       > x0_all,y0_all : relaxed positions of the circle centres
%       > lstRelaxed : indicate which circles have been relaxed
%
%  Yang CHEN, 2020.07.10


xlim0 = ROI(1,1);
xlim1 = ROI(1,2);
ylim0 = ROI(2,1);
ylim1 = ROI(2,2);
nn = length(r0_all);

lstRelaxed = false(nn,1);

% calculate the maximum overlaying parameter before relaxation
Cij = getCmax(x0_all,y0_all,r0_all,eps,ROI);
Ci = max(Cij,[],2);
Cmax = max(Ci);

% check for each circle and relaxe it if neccessary
for i=1:nn
    % if stable fibres --> skip 
%     if id(i)==0; continue;end
    
    %
    x0 = x0_all(i);
    y0 = y0_all(i);
    r0 = r0_all(i);
    chklst = [1:1:nn];  chklst(i) = []; %check list

    %admissible distance
    d_ok = r0 + r0_all(chklst) + eps;
    
    %periodic distance (2D)
    [d,ind] = dist2D_periodic_v3([x0,y0],...
                                 [x0_all(chklst),y0_all(chklst)],...
                                 [xlim0,xlim1;ylim0,ylim1]);
                             
    %calculate the overlaying parameter
    Cij = (d_ok-d) ./ d_ok;
    Ci = max( Cij );
    
    %relaxation of the overlaying fibres
    if Ci>0
        %relaxation through the repulsive force
        Iol = Cij>0;
        x1 = x0_all(chklst(Iol));
        y1 = y0_all(chklst(Iol));
        [x0_new,y0_new] = relaxCircle(x0,y0,Cij(Iol),x1,y1,ind(Iol),ROI);
        
        %calculate the overlaying parameter
        d = dist2D_periodic_v3([x0_new,y0_new],...
                               [x0_all(chklst),y0_all(chklst)],...
                               [xlim0,xlim1;ylim0,ylim1]);
        Ci = max( (d_ok-d) ./ d_ok );
        
        %save the relaxed position if C<Cmax
        if Ci<Cmax
            x0_all(i) = x0_new;
            y0_all(i) = y0_new;
            lstRelaxed(i) = 1;
        else
            continue;
        end
    else
        continue;
    end
end

