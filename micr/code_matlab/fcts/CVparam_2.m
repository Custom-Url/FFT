function [dS,dV] = CVparam_2(P,N,eI)
%compute the local parameters for composite voxels
%   [dS,dV] = CVparam(P,N,eI)
%   -------------------------
%
%   Inputs
%   ------
%       -  P : (nx3 matrix) projection position in the voxel local coordinates
%       -  N : (nx3 matrix) normal vector of the interphase
%       - eI : (scalar) thickness of the interphase [unit: voxel]
%
%   Outputs
%   -------
%       - dS : (nx1 matrix) elementary surface of the interface (intersection area
%              between the plane and the voxel)
%       - dV : (nx2 matrix) volume fractions of each phase(phase1, phase2)
%              note-1: normal vector is pointing from phase1 to phase 2
%              note-2: volume of the interphase can be computed by
%                      subtraction
%
% Yang CHEN 2019.07.02
%

%local surface area
nbins = size(P,1);
xV  = ones(nbins,1)*[-0.5 0.5 0.5 -0.5 -0.5 0.5 0.5 -0.5];
yV  = ones(nbins,1)*[-0.5 -0.5 0.5 0.5 -0.5 -0.5 0.5 0.5];
zV  = ones(nbins,1)*[-0.5 -0.5 -0.5 -0.5 0.5 0.5 0.5 0.5];
dS = area_element_2(N,P,xV,yV,zV);

%local volumes : phase1, phase2, interphase
% (normal vector is pointing from phase1 to phase2)
dV = volume_element(N,P,xV,yV,zV,eI);

