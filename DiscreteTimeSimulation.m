clc
clearvars

kv = 0.25; kp = 0.4; ka = 0.1;
amin = 8; amax = 3; kappa = 0.6; vmax = 30; hst = 5; hgo = hst+vmax/kappa;
deltat = 0.1; vstar = 15; hstar = 30;
t = 0-deltat*3:deltat:200;
vL = 15+2*cos(0.3*t); vL(1:3) = 15;
hdot = zeros(length(t),1); vdot = zeros(length(t),1);
h = zeros(length(t),1); v = zeros(length(t),1);

V=@(h)vmax*(hgo<=h)+(kappa*(h-hst))*(hst<h & h<hgo);
W=@(vL)vmax*(vmax<=vL)+vL.*(vL<vmax);
sat=@(u)(u<-amin).*(-amin)+(-amin<=u & u<=amax).*u+(amax<u).*amax;


for i = 0:length(t)-1
    if i < 4
        v(i+1) = vstar;
        h(i+1) = hstar;
    elseif i == 4
        utemp = kp*(V(h(i-3))-v(i-3))+kv*(W(vL(i-3))-v(i-3));
        v(i+1) = v(i)+sat(utemp)*deltat; 
        h(i+1) = h(i)-v(i)*deltat-0.5*sat(utemp)*(deltat^2)+(vL(i)+vL(i+1))*deltat/2;
    else
        utemp = kp*(V(h(i-3))-v(i-3))+kv*(W(vL(i-3))-v(i-3))+ka*(vL(i-3)-vL(i-4))/deltat;
        v(i+1) = v(i)+sat(utemp)*deltat; 
        h(i+1) = h(i)-v(i)*deltat-0.5*sat(utemp)*(deltat^2)+(vL(i)+vL(i+1))*deltat/2;
    end
end
figure(1)
subplot(1,2,1)
plot(t,h)
hold on
xlabel('Time (s)')
ylabel('Headway (m)')
title('Headway vs Time')
xlim([0-3*deltat 200])

subplot(1,2,2)
plot(t,v,t,vL)
legend('Follower Velocity', 'Lead Velocity')
hold on
xlabel('Time (s)')
ylabel('Velocity (m/s)')
title('Velocity vs Time')
xlim([0-3*deltat 200])