addpath 'fcts'

% clear all; close all


% path to vtk files (strain/stress tensors)
prefix = ['../FM_vf0.6/h0.0001_def']; %This needs to be adjusted to your directory structure
% prefix = ['../FM_vf0.6/h0.001_sig']; %This needs to be adjusted to your directory structure
loadStep = '1';

% file name where you saved your fibre centre locations
fname0 = '../../micr/FM_vf0.6/iUC1_ply0.mat';

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% load the vtk files (tensor components)
T = cell(6,1);
for icomp=1:6
    fname = [prefix, num2str(icomp), '_', loadStep, '.vtk'];
    T{icomp} = readVolvtk_amitex(fname);
end


% transform XYZ-tensor (Cartesian) to rthetaz-tensor (cylindrical)
A = load(fname0);
L = A.L;
xc = A.x0_all;  yc = A.y0_all;
r0 = A.r0_all;

npts = length(xc); %number of fibres

% 
nvx = size(T{1},1);
h = L(1)/nvx;
[y,x] = ndgrid(h:h:L, h:h:L);
siz = size(T{1});
T2 = T;
for i=1:6; T2{i}(:)=NaN; end
srr_av = zeros(npts,1);
stt_av = zeros(npts,1);
szz_av = zeros(npts,1);
srt_av = zeros(npts,1);
srz_av = zeros(npts,1);
stz_av = zeros(npts,1);
for i=1:npts
    
    %xyz coordinates
    id = sqrt( (x-xc(i)).^2 + (y-yc(i)).^2 ) <= r0(i);
    xi = x(id);
    yi = y(id);
    zi = xi.*0;
    
    %cylinder origin & axis
    Ot = [xc(i), yc(i), 0];
    ez = [0 0 1];
    
    %coordinate change
    s1 = T{1}(id);
    s2 = T{2}(id);
    s3 = T{3}(id);
    s4 = T{4}(id);
    s5 = T{5}(id);
    s6 = T{6}(id);
    [srr,stt,szz,srt,srz,stz] = tensLst_XYZ2RTZ(s1,s2,s3,s4,s5,s6,...
                                                [xi,yi,zi],Ot,ez);
    
    % assign the transformed tensor back to the 3D mech (for visualisation)
    T2{1}(id) = srr;
    T2{2}(id) = stt;
    T2{3}(id) = szz;
    T2{4}(id) = srt;
    T2{5}(id) = srz;
    T2{6}(id) = stz;
    
    % calculate average of each fibre
    srr_av(i) = mean(srr);
    stt_av(i) = mean(stt);
    szz_av(i) = mean(szz);
    srt_av(i) = mean(srt);
    srz_av(i) = mean(srz);
    stz_av(i) = mean(stz);
end
    
    %some plots (adjust these according to your needs)
    icomp=1;
    figure;imshow(T2{icomp}(:,:,1),[min(T2{icomp}(:)) max(T2{icomp}(:))])
    colormap 'jet'
    figure;imshow(T{icomp}(:,:,1),[min(T{icomp}(:)) max(T{icomp}(:))])
    colormap 'jet'
    
    figure;
    hold on;plot(srr_av,'<','LineWidth',2)
    hold on;plot(stt_av,'o','LineWidth',2)
    hold on;plot(szz_av,'*','LineWidth',2)
    hold on;plot(srt_av,'^','LineWidth',2)
    hold on;plot(srz_av,'s','LineWidth',2)
    hold on;plot(stz_av,'d','LineWidth',2)
    xlabel('fibre number')
    ylabel('value of the component')
    legend('\sigma_{rr}','\sigma_{\theta\theta}','\sigma_{zz}',...
           '\sigma_{r\theta}','\sigma_{rz}','\sigma_{\thetaz}')
       
% overall average values
sigmaRTZ = zeros(6,1);
sigmaRTZ(1) = mean(srr_av);
sigmaRTZ(2) = mean(stt_av);
sigmaRTZ(3) = mean(szz_av);
sigmaRTZ(4) = mean(srt_av);
sigmaRTZ(5) = mean(srz_av);
sigmaRTZ(6) = mean(stz_av);

disp(sigmaRTZ)
%
