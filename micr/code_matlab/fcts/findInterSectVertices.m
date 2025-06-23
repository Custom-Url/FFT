function [xI,yI,zI] = findInterSectVertices(P,N,xV,yV,zV)
%find the intersection vertices of a plane and a cuboid
%
%   [xI,yI,zI] = findInterSectVertices(P,N,xV,yV,zV)
%   ------------------------------------------------
%
%   Inputs
%       - P : coordinates of the point passed through by the plane
%       - N : normal vector of the plane
%       - xV,yV,zV: x y z coordinates of the cuboid vertices
%
%   Outputs
%       - xI,yI,zi: x y z coordinates of the intersection vertices
%
%   written by Yang CHEN 2019.07.02
%       copied from "area_element_2.m"
%

nbox = size(N,1);
xI=zeros(nbox,12);
yI=xI;
zI=xI;
for i=1:12
    [n1,n2] = num2vertices(i); %each side line is associated to two vertices
    V1 = [xV(:,n1), yV(:,n1), zV(:,n1)];
    V2 = [xV(:,n2), yV(:,n2), zV(:,n2)];
    a = V2(:,1) - V1(:,1);
    b = V2(:,2) - V1(:,2);
    c = V2(:,3) - V1(:,3);
    tmp1 = ( N(:,1).*(P(:,1)-V1(:,1)) + ...
             N(:,2).*(P(:,2)-V1(:,2)) + ...
             N(:,3).*(P(:,3)-V1(:,3)) );
    tmp2 = ( a.*N(:,1) + b.*N(:,2) + c.*N(:,3) );
    
    xI(:,i) = V1(:,1) + a .* tmp1 ./ tmp2;
    yI(:,i) = V1(:,2) + b .* tmp1 ./ tmp2;
    zI(:,i) = V1(:,3) + c .* tmp1 ./ tmp2;
    
        % note: in the case that the S passes through the L (tmp2==0 <->
        % S//L;  tmp1==1 <-> V1 is in S), there will be no intersection
        % points for this being-passed L, however, the two vertices of the
        % being-passed L will be double-counted by the L's four
        % intersecting lines. Since the function "polyarea3D.m" is able to
        % provide the accurate area value even the polygon's vetices are
        % repeately presented, so this special case (S passing through L)
        % can be also correctly considered in this code. 
    
    % intersection points should be contained within the cuboid sides
    i0 = ( xI(:,i)-V1(:,1)>=0 & xI(:,i)-V2(:,1)<=0 ) & ...
         ( yI(:,i)-V1(:,2)>=0 & yI(:,i)-V2(:,2)<=0 ) & ...
         ( zI(:,i)-V1(:,3)>=0 & zI(:,i)-V2(:,3)<=0 );
    i0 = ~i0;
    xI(i0,i) = NaN;
    yI(i0,i) = NaN;
    zI(i0,i) = NaN;
    i0 = xI==Inf;   %plan//cuboid side&cuboid side does not belong to plan
    xI(i0,i) = NaN;
    yI(i0,i) = NaN;
    zI(i0,i) = NaN;
end


%%     Each side line of the cuboid is associated to two vertices
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [n1,n2] = num2vertices(n)
switch n
    case 1
        n1 = 1;    n2 = 2;
    case 2
        n1 = 2;    n2 = 3;
    case 3
        n1 = 4;    n2 = 3;
    case 4
        n1 = 1;    n2 = 4;
    case 5
        n1 = 5;    n2 = 6;
    case 6
        n1 = 6;    n2 = 7;
    case 7
        n1 = 8;    n2 = 7;
    case 8
        n1 = 5;    n2 = 8;
    case 9
        n1 = 1;    n2 = 5;
    case 10
        n1 = 2;    n2 = 6;
    case 11
        n1 = 3;    n2 = 7;
    case 12
        n1 = 4;    n2 = 8;
end

