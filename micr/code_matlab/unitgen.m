addpath ~/headQuarter_matlab/0_fct_basic/
addpath ~/headQuarter_matlab/0_fct_basic/netcodes/
addpath ~/headQuarter_matlab/0_fct_basic/netcodes/imshow3D
addpath ~/headQuarter_matlab/1a_unitGen/fcts
addpath fcts


clear all

% Parameters --------------------------------------------------------------
% basics
vf = 0.6; %fibre volume fraction
L = [0.1, 0.1]; %mm, side length of the RVE cross-section square
H = 0.1; %mm, length of the RVE

% fibre diameter
rParam.AV = 3.5e-3; %mm, average fibre radius
rParam.MIN = 1.5e-3; %mm, min fibre radius
rParam.MAX = 5.5e-3; %mm, max fibre radius
rParam.SIG = 0.5e-3; %std of radii distribution (Gaussian)

% algorithm parameters
algo.eps_fibres = 3e-4; %minimum distance between fibres
algo.itMax_ovlap = 5e3; %max number of iterations in case of overlap
algo.num_repeats = 100; %number of repeats to reach the randomness of fibre distribution
% -------------------------------------------------------------------------

UC.xlim0 = 0;
UC.xlim1 = L(1);
UC.ylim0 = 0;
UC.ylim1 = L(2);

for iUC=1:1
    % generate random fibre positions -------------------------------------
    figure;
    [x0_all, y0_all, r0_all] = randfibres(L, vf, rParam, algo, 'yes');
    vf_true = sum(pi.*(r0_all.^2)) / (L(1)*L(2));
    disp(['True vf: ', num2str(vf_true)]);
    
    % plot the circles
    hold off; plotPeriodicCircles(x0_all,y0_all,r0_all,UC);
    
    % save the generated data
    dirname = ['../FM_vf',num2str(vf)];
    if ~exist(dirname, 'dir')
        mkdir(dirname);
    end
    
    if mod(iUC,2)==0
        fnameS = ['../FM_vf',num2str(vf),'/iUC',num2str(iUC-1),'_ply90.mat'];
    else
        fnameS = ['../FM_vf',num2str(vf),'/iUC',num2str(iUC),'_ply0.mat'];
    end
    save(fnameS,...
         'x0_all','y0_all','r0_all','vf_true','L')
    disp(['--> data saved in ',fnameS])
end



%% ########################################################################
%% ########################################################################
%% 3D laminate
%% ===========

vxsiz0 = 0.1e-3; %mm
% vxsiz0 = 0.04; %mm

dirname0 = '../FM_vf0.6';

%ply - 0
load([dirname0,'/iUC1_ply0.mat'])
V_0 = generate_ply_noInterface(x0_all, y0_all, r0_all, L, vxsiz0);

% %ply - 90
% load([dirname0,'/iUC1_ply90.mat'])
% V_90 = generate_ply_noInterface(x0_all, y0_all, r0_all, L, vxsiz0);
% 
% %assembly
% siz = size(V_0);
% 
% %rotate w.r.t. y-axis by 90 degree
% V_0(V_0==2) = 3;
% V_90_rot = permute(V_90, [3 2 1]);
% V_90_rot(V_90_rot==2) = 3;
% V_90_rot(V_90_rot==1) = 2;
% V_a = cat(2, V_0, V_90_rot);
% 
% 
%     isl = 11;
%     SL = V_a(:,:,isl);
%     figure;imshow(SL, [0 3])
 


%% save

%vtk
fname = [dirname0,'/h',num2str(vxsiz0),'.vtk'];
saveVolvtk_amitex(V_0(:,:,1),fname,'uint8','cell',ones(3,1).*vxsiz0);
% saveVolvtk_amitex(V_a,fname,'uint8','cell',ones(3,1).*vxsiz0);



