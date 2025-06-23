function [varargout] = dist2D_periodic_v3(A,B,box)
%compute the distance between two points A and B (2D) considering the
%periodicity
%
%   d = dist2D_periodic(A,B,box)
%   [d,ind] = dist2D_periodic(A,B,box)
%   ----------------------------------
%
%   Inputs
%   ------
%       >   A : [x,y], point A
%       >   B : [x(:),y(:)], point B
%               B may be an array of points
%       > box : [x1,x2;
%                y1,y2], defining the 2D periodic box
%
%   Output
%   ------
%       >   d : distance between the two points A and B
%       > ind : integer (0,1,2,3,4) indicating the formula for periodic
%               distance
%               for ind=0, Xi and Xj away from the borders
%               for ind=1, Xi~left border and Xj~right border
%               for ind=2, Xi~right border and Xj~left border
%               for ind=3, Yi~left border and Yj~right border
%               for ind=4, Yi~right border and Yj~left border
%                   with Xi the point to which periodicity is applied if
%                   required
%
% written by Yang CHEN, 2020.04.01
%

xi = A(1);
yi = A(2);
xj = B(:,1);
yj = B(:,2);
x1 = box(1,1);  x2 = box(1,2);
y1 = box(2,1);  y2 = box(2,2);
npts = size(B,1); %number of points in B
Lx = x2-x1;
Ly = y2-y1;

%ordinary distance (without periodicity)
d0 = sqrt( ( xi - xj ).^2 + ( yi - yj ).^2 ); 

% %periodic distance (2D) --> this was not exhausitive (2020.07.11)
% d1 = Inf.*ones(npts,1);
% d2 = d1;
% d3 = d1;
% d4 = d1;
% d2x1 = xi - x1; %distance to left border (x1)
% d2x2 = x2 - xi; %distance to right border (x2)
% d2y1 = yi - y1; %distance to top border (y1)
% d2y2 = y2 - yi; %distance to bottom border (y2)
%                %Note: these d2? should be positive, otherwise --> error
% if d2x1<=d2x2
%     d1 = sqrt( ( x2+d2x1 - xj ).^2 + ( yi -yj ).^2 );
% else
%     d2 = sqrt( ( x1-d2x2 - xj ).^2 + ( yi - yj ).^2 );
% end
% if d2y1<=d2y2
%     d3 = sqrt( ( xi - xj ).^2 + ( y2+d2y1 - yj ).^2 );
% else
%     d4 = sqrt( ( xi - xj ).^2 + ( y1-d2y2 - yj).^2 );
% end

% periodic distance (2D) modified on 2020.07.11
d1 = sqrt( ( xi+Lx - xj ).^2 + ( yi - yj ).^2 ); 
d2 = sqrt( ( xi-Lx - xj ).^2 + ( yi - yj ).^2 );
d3 = sqrt( ( xi - xj ).^2 + ( yi+Ly - yj ).^2 );
d4 = sqrt( ( xi - xj ).^2 + ( yi-Ly - yj).^2 );

d5 = sqrt( ( xi+Lx - xj ).^2 + ( yi+Ly -yj ).^2 ); 
d6 = sqrt( ( xi+Lx - xj ).^2 + ( yi-Ly -yj ).^2 ); 
d7 = sqrt( ( xi-Lx - xj ).^2 + ( yi-Ly -yj ).^2 ); 
d8 = sqrt( ( xi-Lx - xj ).^2 + ( yi+Ly -yj ).^2 ); 

%with periodicity, the distance is defined as the minimum of d0...4
[d,Ind] = min([d0,d1,d2,d3,d4,d5,d6,d7,d8],[],2);

%output
varargout{1} = d;
if nargout==2
    varargout{2} = Ind-1;
end

