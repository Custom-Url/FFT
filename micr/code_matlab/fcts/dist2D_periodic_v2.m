function d = dist2D_periodic_v2(A,B,box)
%compute the distance between two points A and B (2D) considering the
%periodicity
%
%   d = dist2D_periodic(A,B,box)
%   ----------------------------
%
%   Inputs
%   ------
%       >   A : [x,y], point A
%       >   B : [x,y], point B
%       > box : [x1,x2;
%                y1,y2], defining the 2D periodic box
%
%   Output
%   ------
%       >   d : distance between the two points A and B
%
% written by Yang CHEN, 2020.04.01
%

xi = A(1);
yi = A(2);
xj = B(1);
yj = B(2);
x1 = box(1,1);  x2 = box(1,2);
y1 = box(2,1);  y2 = box(2,2);

%ordinary distance (without periodicity)
d0 = sqrt( ( xi - xj ).^2 + ( yi - yj ).^2 ); 

%periodic distance (2D)
d1 = Inf;
d2 = Inf;
d3 = Inf;
d4 = Inf;
d2x1 = xi - x1; %distance to left border (x1)
d2x2 = x2 - xi; %distance to right border (x2)
d2y1 = yi - y1; %distance to top border (y1)
d2y2 = y2 - yi; %distance to bottom border (y2)
               %Note: these d2? should be positive, otherwise --> error
if d2x1<=d2x2
    d1 = sqrt( ( x2+d2x1 - xj ).^2 + ( yi -yj ).^2 );
else
    d2 = sqrt( ( x1-d2x2 - xj ).^2 + ( yi - yj ).^2 );
end
if d2y1<=d2y2
    d3 = sqrt( ( xi - xj ).^2 + ( y2+d2y1 - yj ).^2 );
else
    d4 = sqrt( ( xi - xj ).^2 + ( y1-d2y2 - yj).^2 );
end

%with periodicity, the distance is defined as the minimum of d0...4
d = min([d0,d1,d2,d3,d4]);

