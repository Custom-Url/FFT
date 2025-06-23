function [sig,d]=CDM_interphase_fct(lbda,mu,gc,lc,def)
%compute stress from a CDM model

% %material parameters
% lbda = 121.15; %[GPa]
% mu = 80.77; %[GPa]
% gc = 2.7e-3; %[kN/mm]
% lc = 5e-2; %[mm]

%strain decomposition
defP = ( def + abs(def) ) ./ 2;
defN = ( def - abs(def) ) ./ 2;
Tr_def = def(1) + def(2) + def(3);
Tr_def_P = ( Tr_def + abs(Tr_def) ) / 2;
Tr_def_N = ( Tr_def - abs(Tr_def) ) / 2;


%strain energy
psiP = 0.5*lbda * ( Tr_def_P )^2 ...
           + mu * ( sum(defP(1:3).^2) + sum(defP(4:6).^2)/2 );
       
%damage variable
d = 1-gc/(gc+2*lc*psiP);

sig = zeros(6,1);
%stress
sig(1:3) = ( (1-d)^2 )*( lbda*Tr_def_P + 2*mu*defP(1:3) ) ...
                      +( lbda*Tr_def_N + 2*mu*defN(1:3) );
sig(4:6) = ( (1-d)^2 )*( mu*defP(4:6) ) ...
                      +( mu*defN(4:6) );
                  
