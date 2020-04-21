clc; clearvars;
syms deltat kp kv ka Vprimeh z
a0 = [1 -deltat; 0 1];
asig = [-kp*Vprimeh*deltat^2/2, (kp+kv)*deltat^2/2; ...
        kp*Vprimeh*deltat, -(kp+kv)*deltat];
bsig1 = [ka*deltat/2; -ka];
bsig = [-(kv*deltat+ka)*deltat/2, (kv*deltat+ka)];
b0 = [deltat/2; 0];
b1 = [deltat/2; 0];
A = sym(zeros(10));
A(1:2,1:2) = a0;
A(1:2,7:8) = asig;
A(3:end,1:end-2) = eye(8);
Bsig1 = sym(zeros(10,1)); Bsig1(1:2) = bsig1;
Bsig = sym(zeros(10,1)); Bsig(1:2) = bsig;
B0 = sym(zeros(10,1)); B0(1:2) = b0;
B1 = sym(zeros(10,1)); B1(1:2) = b1;
C = sym(zeros(1,10)); C(2) = 1;
T = simplify(C*inv(z*eye(10)-A)*(Bsig1*z^-4+Bsig*z^-3+B0+B1*z))