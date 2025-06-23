function [x0_all,y0_all,Cmax] = generationFibres(r0_all,param,ROI)
%generate a set of fibres with random locations and given radii
%qqa
%   [x0_all,y0_all] = generationFibres(r0_all,param,ROI)
%   ----------------------------------------------------
%
%   Inputs:
%       > r0_all : vector of radii of the fibres
%       >  param : parameters for the generation procedure
%                  - param.eps : admissible min. distance between two fibres
%                  - param.trialsMAX : max. number of trials
%                  - param.ur_sc : vibration magnitude
%       >    ROI : dimensions of the region of interest
%                  ROI = [xlim0,xlim1;
%                         ylim0,ylim1]
%
%   Outputs:
%       > x0_all,y0_all : x,y coordinates of the generated fibre centres
%
% written by Yang CHEN, 2020.07.10
% inspired by C.Chateau's PhD thesis, chapter 1.2
%

% unit cell dimensions
xlim0 = ROI(1,1);
xlim1 = ROI(1,2);
ylim0 = ROI(2,1);
ylim1 = ROI(2,2);
Lx = xlim1 - xlim0;
Ly = ylim1 - ylim0;

% number of fibres
nn = length(r0_all);

% parameters for the generation procedure
eps = param.eps;
trialsMAX = param.trialsMAX;
ur_sc = param.ur_sc;

% random initial fibre center positions
x0_all = rand(nn,1).*Lx;
y0_all = rand(nn,1).*Ly;

% determine the stable fibres (the non-overlaying ones)
Ci = getCmax(x0_all,y0_all,r0_all,eps,ROI);
Ci = max(Ci,[],2);
id = Ci>0;

% reTry=1; %re-try with a new set of initial positions
% while reTry==1
%     % randomly re-distribute the ovelying fibres
%     x0_all(id) = rand(nn2,1).*Lx;
%     y0_all(id) = rand(nn2,1).*Ly;
%     
    %the initial maximum overlaying parameter
    Cmax = getCmax(x0_all,y0_all,r0_all,eps,ROI);
    Cmax = max(Cmax,[],2);
        id = Cmax>0 & id;
        disp(['initially overlapping fibres: ',int2str(transpose(find(id)))])
    Cmax = max(Cmax);
        
    %loop
    trials = 0;
    while (Cmax>0 && trials<=trialsMAX)
        
        trials = trials + 1;
        
        %------------------------------------------relaxe the fibre centres
        [x0_all,y0_all,lstRelaxed] = ...
         relaxationCircles(x0_all,y0_all,r0_all,id,eps,ROI);
                      
        %---------------------vibrate the moved but still overlaying fibres
        [x0_all,y0_all] = ...
           vibrationCircles(x0_all,y0_all,r0_all,eps,ROI,lstRelaxed,ur_sc);
       
        %-------------------------------------identify the overlying fibres
        Ci = getCmax(x0_all,y0_all,r0_all,eps,ROI);
        Ci = max(Ci,[],2);
%         id = Ci>0 & id;
        id = Ci>0;
        
        %----------------------the overlaying parameter after rearrangement
        Cmax = max(Ci);
        
        if ~mod(trials,10)
            disp([int2str(trials),' trials, overlaping fibres: '])
            disp(['      ',int2str(transpose(find(id)))])
%             disp([' x: ',num2str(transpose(x0_all(id)))])
%             disp([' y: ',num2str(transpose(y0_all(id)))])
%             disp([' Ci: ',num2str(transpose(Ci(id)))])
%             disp([' moved: ',num2str(transpose(find(lstRelaxed)))])
        end
    end

    %------if too many trials, re-start with a new set of initial positions
    if trials>trialsMAX
        disp(['after ',num2str(trials), ...
             ' trials, Cmax=',num2str(Cmax),...
             ' --> try again with a new set of initial positions'])
        reTry=1;
    else
        reTry=0;
    end

% end

