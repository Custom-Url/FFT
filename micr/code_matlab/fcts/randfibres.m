function [x0, y0, r0] = randfibres(L, vf, rParam, algo, varargin)
%randomly generate a spatial distribution of fibres in a 2D periodic square
%
%   [x0, y0] = randfibres(L, vf, rParam, algo)
%   [x0, y0] = randfibres(L, vf, rParam, algo, visu)
%   ------------------------------------------------
%
%   Inputs
%       >     L  : side length of the 2D unit-cell square
%       >    vf  : volume fraction of fibres
%                  The true vf will be slightly different to this value
%       > rParam : parameters for statistical distribution of fibre radii
%                  .AV  - average fibre radius
%                  .SIG - std of fibre radii
%                  .MIN - min radius allowed
%                  .MAX - max radius allowed
%       >   algo : parameters for the iterative algorithm
%                  .eps_fibres  - minimum distance between fibre edges
%                  .itMax_ovlap - max nb of iterations in case of overlap 
%                   This parameter is useful for non-constant fibre radii, 
%                   for which initial overlap exists
%                  .num_repeats - number of repeats to reach the randomness
%                                 of fibre distribution 
%
%   Outputs
%       > x0, y0 : 2D coordinates of fibre centres
%       >     r0 : list of fibre radii (index correspoinding to x0, y0)
%
% written by Yang Chen, University of Bath
% 2022.01.21

if nargin>4
    visu = varargin{1};
else
    visu = 'no';
end
fprintf('visu=%s\n', visu)

K = algo.num_repeats;
eps_fibres = algo.eps_fibres;
itMax_ovlap = algo.itMax_ovlap;

% number of fibres in each direction
N = round( sqrt(vf/pi) * L/rParam.AV );

% compact packing
[y0,x0] = ndgrid([0:2*rParam.AV:2*rParam.AV*(N-1)], ...
                 [0:2*rParam.AV:2*rParam.AV*(N-1)]); %N-1 for periodicity

[nfy, nfx] = size(x0);
nf = nfx * nfy;

x0 = x0(:);
y0 = y0(:);

% stretch to the desired vf
x0 = x0 .* L/(2*rParam.AV*N);
y0 = y0 .* L/(2*rParam.AV*N);

% variation of fibre radii
if rParam.SIG == 0
    r0 = ones(nf,1) .* rParam.AV;
else
    rParam.nn = nf;
    r0 = GenDataWithGaussianDistrib(rParam);
end

% true fibre volume fraction
vf_true = sum(pi.*(r0.^2)) / (L*L);
disp(['True vf: ', num2str(vf_true)]);


% random labels of fibres
labels = randperm(nf);


% loop
ovlap = 'true';
Kcounter = 0;
while strcmp(ovlap,'true') || (Kcounter<=K)
    
    if strcmpi(visu,'yes')
        hold off;
        plot(x0, y0, '*', 'LineWidth',4); axis equal
        for i=1:nf
            hold on;rectangle('Position', [x0(i)-r0(i),y0(i)-r0(i),...
                                           2*r0(i),2*r0(i)], ...
                              'Curvature',[1 1])
            text(x0(i), y0(i), num2str(i))
        end
        hold on;rectangle('Position', [0 0 L L], 'EdgeColor', 'r')
        pause(0.1)
    end
    
    % loop over fibres
    tag_ovlap = 0;
    for i=1:nf

        %identify this fibre and the other fibres
        id = [1:1:nf];  id(i) = [];
        xi=x0(labels(i)); yi=y0(labels(i)); ri=r0(labels(i)); %this fibre
        x1=x0(labels(id)); y1=y0(labels(id)); r1=r0(labels(id)); %the other fibres
        
        % pre-determine the movement amplitude by the closest neighbors
        di = dist2D_periodic([xi,yi], [x1,y1], [0 L; 0 L]);
        di = sort(di);
        ur_sc = mean(di(1:3));
        
        % randomly move the fibre, until no overlap
        dimin = 0;
        counter=0;
        n_ovlap0 = Inf;

        while (dimin < eps_fibres)
            % random movement of the fibre
            [xinew, yinew] = randmove2D(xi, yi, L, ur_sc);

            % distance between this fibre and the other fibres
            di = dist2D_periodic([xinew,yinew], [x1,y1], [0 L; 0 L]);
            di = di - ri - r1;

            % the closest neighbor
            dimin = min(di);

            % memorise the "least worst" option
            n_ovlap1 = nnz(di<0);
            if (n_ovlap1<n_ovlap0)
                xi_leastworst = xinew;
                yi_leastworst = yinew;
                n_ovlap0 = n_ovlap1;
            end
            
            % locking
            counter = counter+1;
            if counter > itMax_ovlap
                %use the "least worst" option
                xinew = xi_leastworst;
                yinew = yi_leastworst;
                %tag this
                tag_ovlap = tag_ovlap+1;
                break
            end
            
        end

        % update the new position
        x0(labels(i)) = xinew;
        y0(labels(i)) = yinew;

   end
   
   if tag_ovlap>0
       ovlap = 'true';
   else
       ovlap = 'false';
   end
   fprintf('repeats: %d;  Overlap?: %s \n', Kcounter, ovlap);
   Kcounter = Kcounter + 1;

end
