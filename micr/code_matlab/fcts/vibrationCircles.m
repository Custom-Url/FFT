function [x0_all,y0_all] = vibrationCircles(x0_all,y0_all,r0_all,eps,ROI,lstRelaxed,ur_sc)
%Vibration - arbitrarily move the relaxed by still overlying circles
%


xlim0 = ROI(1,1);
xlim1 = ROI(1,2);
ylim0 = ROI(2,1);
ylim1 = ROI(2,2);
Lx = xlim1 - xlim0;
Ly = ylim1 - ylim0;
nn = length(r0_all);

% calculate the overlaying parameters
Cij = getCmax(x0_all,y0_all,r0_all,eps,ROI);
Ci = max(Cij,[],2);
Cij_relaxed = Cij(lstRelaxed,:);
NBov = sum(Cij_relaxed>0,2);


% list of the relaxed circles
lstVIB = find(lstRelaxed==1);
% lstVIB = find( lstRelaxed==1 & Ci>0);
% lstVIB = find( Ci<=0);


% vibration if no more than 3 circles are overlaying
for j=1:length(lstVIB)
%     if NBov(j)>3; continue; end
    i = lstVIB(j);
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
    Cij = (d_ok-d) ./ d_ok;
    if max(Cij>0)
        %vibrate the overlaying fibres
        ur = rand(1).*ur_sc;
        ut = rand(1).*(2*pi);
        x0_all(i) = x0_all(i) + ur.*cos(ut);
        y0_all(i) = y0_all(i) + ur.*sin(ut);
        %check if the new positions are inside the domain
        if (x0_all(i)<xlim0); x0_all(i) = x0_all(i) + Lx; end
        if (y0_all(i)<ylim0); y0_all(i) = y0_all(i) + Ly; end
        if (x0_all(i)>xlim1); x0_all(i) = x0_all(i) - Lx; end
        if (y0_all(i)>ylim1); y0_all(i) = y0_all(i) - Ly; end
    else
        continue;
    end
end
