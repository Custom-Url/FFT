function [C] = GetStiffIso(E,NU)

C = 0;
G = E/2/(1+NU);

% Inverse the compliance tensor
DET = (1 + NU) * ((1 - NU)/E - 2*NU^2/E) / E / E;
A = (1 - NU^2)/E/E;
B = NU * (1 + NU)/E/E;
H = (NU/E)^2 + NU/E/E;

C1 = A/DET; %C11=C22=C33
C2 = B/DET; %C12=C21=C13=C31=C23=C32
C3 = G;     %C44=C55=C66

C(1) = C1;
C(2) = C2;
C(3) = C2;
C(4) = C1;
C(5) = C2;
C(6) = C1;
C(7) = C3;
C(8) = C3;
C(9) = C3;
