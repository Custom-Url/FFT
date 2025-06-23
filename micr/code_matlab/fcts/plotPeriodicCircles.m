function [] = plotPeriodicCircles(x0_all,y0_all,r0_all,UC,varargin)
%plot circles described by centres and radii over a periodic region
%
%   plotPeriodicCircles(x0_all,y0_all,r0_all,UC,offsetIM)
%   -----------------------------------------------------
%
%   Inputs:
%       > x0_all,y0_all: circle centres
%       > r0_all : circle radii
%       > UC : parameters defining the unit-cell region
%              UC.xlim0
%              UC.xlim1
%              UC.ylim0
%              UC.ylim1
%       > offsetIM : (optional) offset for the periodicity to be correctly
%                    considered in an image-type unit-cell region
%                    offsetIM='offsetIM' if plotting over an image
%                    by default, offsetIM absent or any other characters
%
%  Yang CHEN, 2020.07.12
%  chenyangpic@gmail.com

% offset for plotting over an image
offset = 0;
if nargin==5
    if strcmpi(varargin{1},'offsetIM')
        offset = 1;
    end
end


% unit cell
xlim0 = UC.xlim0;
xlim1 = UC.xlim1;
ylim0 = UC.ylim0;
ylim1 = UC.ylim1;

nn = length(x0_all);
nbs = strsplit(int2str([1:1:nn]));

% circle centres
plot(x0_all,y0_all,'xg');axis equal;axis tight

% labels
hold on;for i=1:nn; text(x0_all(i),y0_all(i),nbs{i}); end

% periodicity
hold on;
for ifig=1:nn
    h=plot_circle(x0_all(ifig),y0_all(ifig),r0_all(ifig),'b');
    if x0_all(ifig)-r0_all(ifig) < xlim0
        h=plot_circle(x0_all(ifig)+xlim1-offset,y0_all(ifig),r0_all(ifig),'b');
    end
    if x0_all(ifig)+r0_all(ifig) > xlim1
        h=plot_circle(x0_all(ifig)-xlim1+offset,y0_all(ifig),r0_all(ifig),'b');
    end
    if y0_all(ifig)-r0_all(ifig) < ylim0
        h=plot_circle(x0_all(ifig),y0_all(ifig)+ylim1-offset,r0_all(ifig),'b');
    end
    if y0_all(ifig)+r0_all(ifig) > ylim1
        h=plot_circle(x0_all(ifig),y0_all(ifig)-ylim1+offset,r0_all(ifig),'b');
    end

    if x0_all(ifig)-r0_all(ifig)<xlim0 && y0_all(ifig)-r0_all(ifig)<ylim0
        h=plot_circle(x0_all(ifig)+xlim1-offset,y0_all(ifig)+ylim1-offset,r0_all(ifig),'b');
    end
    if x0_all(ifig)-r0_all(ifig)<xlim0 && y0_all(ifig)+r0_all(ifig)>ylim1
        h=plot_circle(x0_all(ifig)+xlim1-offset,y0_all(ifig)-ylim1+offset,r0_all(ifig),'b');
    end
    if x0_all(ifig)+r0_all(ifig)>xlim1 && y0_all(ifig)+r0_all(ifig)>ylim1
        h=plot_circle(x0_all(ifig)-xlim1+offset,y0_all(ifig)-ylim1+offset,r0_all(ifig),'b');
    end
    if x0_all(ifig)+r0_all(ifig)>xlim1 && y0_all(ifig)-r0_all(ifig)<ylim0
        h=plot_circle(x0_all(ifig)-xlim1+offset,y0_all(ifig)+ylim1-offset,r0_all(ifig),'b');
    end
end

    set(gca,'XLim',[xlim0,xlim1],'YLim',[ylim0,ylim1])

hold off