clc
clearvars
syms deltat Vprime kp kv ka z

a0 = [1 -deltat; 0 1];
asigma1 = [-deltat^2*kp*Vprime/2, deltat^2*(kp+kv+ka/deltat); deltat*kp*Vprime, -deltat*(kp+kv+ka/deltat)];
asigma2 = [0 -deltat*ka/2; 0 ka];

b0 = [deltat/2; 0];
b1 = [deltat/2; 0];
bsigma1 = [-deltat^2*(kv+ka/deltat)/2; deltat*(kv+ka/deltat)];
bsigma2 = [deltat*ka/2; -ka];

c = [0 1];

A = sym(zeros(10));
A(3:end,1:end-2) = eye(8);
A(1:2,1:2) = a0; A(1:2,7:8) = asigma1; A(1:2,end-1:end) = asigma2;
Bsigma1 = sym(zeros(10,1)); Bsigma1(1:2) = bsigma1;
Bsigma2 = sym(zeros(10,1)); Bsigma2(1:2) = bsigma2;
B0 = sym(zeros(10,1)); B0(1:2) = b0;
B1 = sym(zeros(10,1)); B1(1:2) = b1;
C = sym(zeros(1,10)); C(1:2) = c;
disp(A)
disp(Bsigma1)
disp(Bsigma2)
disp(B0)
disp(B1)
disp(C)

TF = simplify(C*inv(z*eye(10)-A)*(Bsigma2*z^(-4)+Bsigma1*z^(-3)+B0+B1*z))