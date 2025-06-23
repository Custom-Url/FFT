function [r0_all] = GenDataWithGaussianDistrib(r)
%generate a set of data with a given Gaussian distribution
%
%   [r0_all] = GenDataWithGaussianDistrib(r,nn)
%   -------------------------------------------
%
%   Inputs:
%       >  r : struct
%              r.nn - number of data points to be generated
%              r.AV - average
%              r.MIN - min
%              r.MAX - max
%              r.SIG - std
%
%   Outputs:
%       > r0_all : the generated data set following the Gaussian distrib.
%
% Yang CHEN, 2020.07.10
%

nn = r.nn;

x0 = linspace(r.MIN,r.MAX,1000); %mm
r0 = normpdf(x0,r.AV,r.SIG);  %mm

cdf = cumsum(r0);
cdf = cdf ./ max(cdf);
x = linspace(0,1,length(cdf));

r0_all = zeros(nn,1);
for i=1:nn
    tmp = rand(1);
    r0_all(i) = interp1(cdf,x,tmp);
end
r0_all = r0_all.*(max(x0)-min(x0))+min(x0);
