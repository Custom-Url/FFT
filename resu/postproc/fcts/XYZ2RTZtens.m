function [sigrtz] = XYZ2RTZtens(sigxyz,posXYZ,Ot,ez)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Transform 3x3 tensor(s) from XYZ coordinates to RTZ coodinates   %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       [sigrtz] = XYZ2RTZtens(sigxyz,posXYZ,Ot,ez)
%      ---------------------------------------------
%   Inputs:
%       - sigxyz: (3x3xn matrix) the set of 3x3 tensors to be transformed
%       - posXYZ: (nx3 matrix) the postions of points associated to tensors
%       - Ot: the origin of cylindrical coordinates
%       - ez: the tube axis
%   Output:
%       - sigrtz: the set of tensors transformed into RTZ coordinates
%
% Y.Chen 01/02/2016  31/05/2016

% comput the tube-local coodinates for each voxel of tube region
[er,et]=LocAxisCalc_v2(Ot,ez,[posXYZ(:,1),posXYZ(:,2),posXYZ(:,3)]);
    nvx = length(posXYZ(:,1));
Q = cat(3,er,et,ones(nvx,1)*ez);
Q = permute(Q,[2 3 1]);

% multiplication of transfomation matrices
%tic
%sigrtz = zeros(3,3,nvx,'single');
%for i=1:nvx
%    sigrtz(:,:,i) = Q(:,:,i)'*sigxyz(:,:,i)*Q(:,:,i);
%end
%toc

% multiplication of transfomation matrices
sigrtz=sigxyz.*0;
for m=1:3
    for n=1:3
        tmp=zeros(1,1,nvx,'single');
        for i=1:3
            for j=1:3
                tmp = tmp + Q(i,m,:).*sigxyz(i,j,:).*Q(j,n,:);
            end
        end
        sigrtz(m,n,:) = tmp;
    end
end


% % ----test : the test should be equal to [1 0 0] for all points, if not
% % the transformation matrix is wrongly defined and you must correct the script!!! 
% test=zeros(nvx,3);
% for i=1:nvx
%    test(i,:) = Q(:,:,i)*er(i,:)';
% end
