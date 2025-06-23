function V = generate_ply_noInterface(x0_all, y0_all, r0_all, L, vxsiz0)
%generate one ply, from a 2D fibre distribution map

%%

% discretisation
nR = ceil(L(2)/vxsiz0);  if (mod(nR,2)==0); nR=nR+1; end
nC = ceil(L(1)/vxsiz0);  if (mod(nC,2)==0); nC=nC+1; end
V = ones(nR,nC,1,'uint8').*2;
siz = size(V);

UC.xlim0 = 0;
UC.xlim1 = L(1);
UC.ylim0 = 0;
UC.ylim1 = L(2);
UC.vxsiz = vxsiz0;

% [y,x] = ndgrid(0.5:nR-0.5,0.5:nC-0.5);
[y,x] = ndgrid(1:nR,1:nC);
x = x.*vxsiz0;
y = y.*vxsiz0;

nn = length(x0_all);
for icirc = 1:nn
    
    x0 = x0_all(icirc);
    y0 = y0_all(icirc);
    r0 = r0_all(icirc);
    
    % fibre region
    d = dist2D_periodic_v3([x0,y0],[x(:),y(:)],[0,L(1);0,L(2)]);
    id = d <= r0;
    V(id) = 1; %phase#1 -- fibre
%         hold on;imshow(V,[0 2]);
end

%% extrusion in Z direction
Lz = max(L);
nB = ceil(Lz/vxsiz0);  if (mod(nB,2)==0); nB=nB+1; end
V = repmat(V,1,1,nB);



