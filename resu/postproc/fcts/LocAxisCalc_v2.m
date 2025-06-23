function [er,et] = LocAxisCalc_v2(Ot,ez,M)
%Compute the local axes for each micro pore's CoG
%
%   [er,et] = LocAxisCalc(Ot,ez,M)
%
%  Inputs:
%       > Ot: the origin of the tube (vector of 3x1)
%       > ez: the Z-axis of the tube (vector of 3x1)
%       > M: the coodinates of a point (matrix of (nptsx3))
%

Mp = sum([(M(:,1) - Ot(1))*ez(1), (M(:,2) - Ot(2))*ez(2), (M(:,3) - Ot(3))*ez(3)],2)*ez;
Mp = [Mp(:,1)+Ot(1), Mp(:,2)+Ot(2), Mp(:,3)+Ot(3)];
er = M - Mp;
l_er = sqrt(sum(er.*er,2));
er = [er(:,1)./l_er er(:,2)./l_er er(:,3)./l_er];
et = [ez(:,2).*er(:,3)-ez(:,3).*er(:,2) ez(:,3).*er(:,1)-ez(:,1).*er(:,3) ez(:,1).*er(:,2)-ez(:,2).*er(:,1)];

end