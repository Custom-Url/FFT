function sig = GetStress(C,STRANT)

sig(1) = C(1)*STRANT(1) + C(2)*STRANT(2) + C(3)*STRANT(3);
sig(2) = C(2)*STRANT(1) + C(4)*STRANT(2) + C(5)*STRANT(3);
sig(3) = C(3)*STRANT(1) + C(5)*STRANT(2) + C(6)*STRANT(3);
sig(4) = C(7)*STRANT(4);
sig(5) = C(8)*STRANT(5);
sig(6) = C(9)*STRANT(6);