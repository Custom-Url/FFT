function [srr,stt,szz,srt,srz,stz] = tensLst_XYZ2RTZ(s1,s2,s3,s4,s5,s6,XYZ,Ot,ez)
%% ========================================================================
%% transform the stress/strain tensor from Cartesian to cylindrical
%% ========================================================================
%
%  [srr,stt,szz,srt,srz,stz] = tensLst_XYZ2RTZ(s1,s2,s3,s4,s5,s6,XYZ,Ot,ez)
%  ------------------------------------------------------------------------
%
%   Inputs
%       > s1,...,s6 : components of the stress/strain tensor in the global
%                     system
%                     tensor = [s1 s4 s5
%                               s4 s2 s6
%                               s5 s6 s3]
%       >       XYZ : coordinates of each point
%       >        Ot : (3x1) vector tube centre coordinates
%       >        ez : tube axis
%
%   Outputs
%       > srr,...stz : components of the stress/strain tensor in the local
%                      coordinate system
%                      tensor = [srr srt srz
%                                srt stt stz
%                                srz stz szz]
%
%  28/09/2016 written by Yang Chen
% comments added on 2019.12.01
%


%% regroup into tensor form
sxyz = cat(3,[s1(:) s4(:) s5(:)],...
             [s4(:) s2(:) s6(:)],...
             [s5(:) s6(:) s3(:)]);
sxyz = permute(sxyz,[2 3 1]);

%% coordinate change
srtz = XYZ2RTZtens(sxyz,XYZ,Ot,ez);

%% regroup into point list form

srr = squeeze(srtz(1,1,:));
stt = squeeze(srtz(2,2,:));
szz = squeeze(srtz(3,3,:));
srt = squeeze(srtz(1,2,:));
srz = squeeze(srtz(1,3,:));
stz = squeeze(srtz(2,3,:));


