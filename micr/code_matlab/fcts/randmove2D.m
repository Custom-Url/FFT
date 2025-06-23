
function [x1, y1] = randmove2D(x0, y0, L, ur_sc, varargin)
%random movement of points in a 2D rectangle region
% --> periodicity is considered
%
%   [x1, y1] = randmove2D(x0, y0, L, ur_sc)
%   [x1, y1] = randmove2D(x0, y0, L, ur_sc, pad)
%   --------------------------------------------
%
%   Inputs
%       > x0, y0 : coordinates of the points to be moved
%       >      L : Size of the 2D rectangle region (origin at 0)
%       >  ur_sc : Scale of the movement step
%       >    pad : (optional), if absent, pad='periodic'
%                   'periodic' - periodic padding
%                   'wall' - region outside [0,L]x[0,L] not allwed
%
%   Outputs
%       > x1, y1 : coordinates of new point positions
%
% written by Yang Chen, University of Bath
%

pad = 'periodic';
if nargin>4
    pad = varargin{1};
end

%
[nfy, nfx] = size(x0);

if strcmp(pad,'periodic')
    % random movement of fibres
    ur = rand(nfy,nfx) .* ur_sc;
    ut = rand(nfy,nfx) .* (2*pi);
    x1 = x0 + ur.*cos(ut);
    y1 = y0 + ur.*sin(ut);

    % periodicity
    id = x1<0;  x1(id) = x1(id) + L;
    id = y1<0;  y1(id) = y1(id) + L;
    id = x1>=L;  x1(id) = x1(id) - L;
    id = y1>=L;  y1(id) = y1(id) - L;
    
elseif strcmp(pad,'wall')
    notOkay = true;
    while notOkay
        % random movement of fibres
        ur = rand(nfy,nfx) .* ur_sc;
        ut = rand(nfy,nfx) .* (2*pi);
        x1 = x0 + ur.*cos(ut);
        y1 = y0 + ur.*sin(ut);

        % check the wall-condition
        if nnz( x1<0 | y1<0 | x1>=L| y1>=L ) ==0
            notOkay = false;
        end
    end
else
    error('padding needs to be either "periodic" or "wall"')
end

end    
