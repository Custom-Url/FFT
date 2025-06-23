function [dV] = volume_element(N,P,xV,yV,zV,eI)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Compute the volumes of the two parts of a cuboid cut by a "thick" plane surface
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	
%	[dS] = volume_element_2(N,P,xV,yV,zV,eI)
%	----------------------------------------
%
%    Inputs:
%    -------
%	>          N: (matrix nx3) the normal vector of the plane
%	>          P: (matrix nx3) the point that the median plane is passing
%	              through 
%	> [xV yV zV]: (matrix nx8) the 8 vertices of the cuboid
%   >         eI: thickness of the plane surface [unit: voxel]
%
%    Output:
%    -------
%	>         dV: (nx2 matrix) volumes of [phase1,phase2]
%             note-1: normal vector is pointing from phase1 to phase2
%             note-2: volume of the thick plane can be computed by
%                     subtraction
%
%    written by Yang CHEN 2019.07.02
%    (inspired from "area_element_2.m")
%

nbox = size(N,1);

%% find the intersection vertices
% for face 1 of the plane ------------------
P_face1 = P;
P_face1(:,1) = P(:,1) - eI.*N(:,1);
P_face1(:,2) = P(:,2) - eI.*N(:,2);
P_face1(:,3) = P(:,3) - eI.*N(:,3);
[xI_face1,yI_face1,zI_face1] = findInterSectVertices(P_face1,N,xV,yV,zV);

% for face 2 of the plane ------------------
P_face2 = P;
P_face2(:,1) = P(:,1) + eI.*N(:,1);
P_face2(:,2) = P(:,2) + eI.*N(:,2);
P_face2(:,3) = P(:,3) + eI.*N(:,3);
[xI_face2,yI_face2,zI_face2] = findInterSectVertices(P_face2,N,xV,yV,zV);


%% identify on which side the cuboid vertices are standing
% id = N(:,1).*(xV-P(:,1)) + N(:,2).*(yV-P(:,2)) + N(:,3).*(zV-P(:,3));
% id = id > 0;  %if id=0, face 1; if id=1, face 2
id = zeros(nbox,8);
for i=1:8
    id(:,i) = N(:,1).*(xV(:,i)-P(:,1)) + N(:,2).*(yV(:,i)-P(:,2)) + N(:,3).*(zV(:,i)-P(:,3));
end
id = id > 0;  %if id=0, face 1; if id=1, face 2

%% compute the volume of each phase (phase 1, 2 and interphase)
dV=zeros(nbox,2);
for i=1:nbox
    %phase 1 --------------------------------------------------------------
    i0 = ~isnan(xI_face1(i,:));
    if nnz(i0)<3
        dV(i,1)=0; continue;
    else
        vertices = [];
        vertices(:,1) = xI_face1(i,i0);
        vertices(:,2) = yI_face1(i,i0);
        vertices(:,3) = zI_face1(i,i0);
        %together with the cuboid vertices on the same side
        vertices = [ vertices; [xV(i,~id(i,:))',yV(i,~id(i,:))',zV(i,~id(i,:))'] ];
        %compute the volume
        try %added by YC 2020.05.07
            [junk,tmp] = convhulln(vertices);
        catch
            warning(['Error in convhulln for point:',num2str(P(i,:))])
            tmp = NaN;
        end
        dV(i,1) = tmp;
    end
    %phase 2 --------------------------------------------------------------
    i0 = ~isnan(xI_face2(i,:));
    if nnz(i0)<3
        dV(i,2)=0; continue;
    else
        vertices = [];
        vertices(:,1) = xI_face2(i,i0);
        vertices(:,2) = yI_face2(i,i0);
        vertices(:,3) = zI_face2(i,i0);
        %together with the cuboid vertices on the same side
        vertices = [ vertices; [xV(i,id(i,:))',yV(i,id(i,:))',zV(i,id(i,:))'] ];
        %compute the volume
        try %added by YC 2020.05.07
            [junk,tmp] = convhulln(vertices);
        catch
            warning(['Error in convhulln for point:',num2str(P(i,:))])
            tmp = NaN;
        end
        dV(i,2) = tmp;
    end
end

