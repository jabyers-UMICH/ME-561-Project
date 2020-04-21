clearvars
clc
%% Parameters
% number of vehicles
veh_num=1;
% distribution of ACC and human driven vehicles
driver=cell(1,veh_num); driver(:)={'CCC'};  
amin = 8; amax = 3;
sat=@(u)(u<-amin).*(-amin)+(-amin<=u & u<=amax).*u+(amax<u).*amax;
% ACC parameters
% control gains
kv=0.5; kp=0.4; ka=0.4;
% delay
sigma=0.6;
% range policy parameters
kappa=0.6; vmax=35; hst=5; hgo=55;
% equilibrium velocity and headway
vstar=15; hstar=10;
% range policy
V=@(h)vmax*(hgo<=h)+kappa*(h-hst).*(hst<h & h<hgo);
W=@(vL)vmax*(vmax<=vL)+vL.*(vL<vmax);
% eq's for CCC vehicle
ccc=@(t,x,xdelay,xdotdel,vL,vLdelay,vLdot,sat,kp,kv,ka,V,W)...
    [vL(t)-x(2);...
    sat(kp*(V(xdelay(1))-xdelay(2))+kv*(W(vLdelay(t))-xdelay(2))+ka*(vLdot(t-sigma)))];
% initial condition - history for human driven vehicle
v0 = 15; h0 = 10; cccinit=@(t)[h0;v0];
% Lead Car Input
vLead=@(t) (v0-10*t).*(t>=0 & t<1.5)+(t<0)*v0;
vLdot =@(tdelay) (tdelay>=0 & tdelay<=1.5)*(-10)+0.*tdelay;
% regions for plotting
hminplot=-5; hmaxplot=20; vminplot=-5; vmaxplot=30;
% time
t0=0; tend=10; deltat=0.001;
%% Simulation
% simulation time
hist=t0-sigma:deltat:t0;	% time of history (initial conditions)
time=t0:deltat:tend;	% simulation time
% simulate the motion of each vehicle one by one
hFollow=zeros(veh_num,length(time));
vFollow=zeros(size(hFollow));
for kveh=1:veh_num
    disp(['Simulation of vehicle #',num2str(kveh)]);
    % pick initial condition and model for ACC or human driven vehicle
    delay=sigma;
    DDErhs=@(t,x,xdelay,xdotdel,vL,vLdelay)ccc(t,x,xdelay,xdotdel,vL,vLdelay,vLdot,sat,kp,kv,ka,V,W);
    DDEic=@(t)cccinit(t);
    % pick leader velocity for simulations
    if kveh==1  % first vehicle with prescribed velocity
        vL=@(t)vLead(t);
        vLdelay=@(t)vLead(t-delay);
    else        % preceding vehicle
        vL=@(t)interp1(time,vFollow(kveh-1,:),t,'linear','extrap');
        select=@(u,ku)u(ku);
        vLdelay=@(t)(t-delay<=t0).*select(DDEic(t-delay),2)+...
           (t-delay>t0).*interp1(time,vFollow(kveh-1,:),t-delay,'linear','extrap');
    end
    % perform simulation
    sol=ddensd(@(t,x,xdelay,xdotdel)DDErhs(t,x,xdelay,xdotdel,vL,vLdelay),@(t,x)t-sigma,@(t,x)t-sigma,DDEic,[time(1) time(end)]);
    [x, xp] = deval(sol,time); 
    % extract headways and velocities from the results
    hFollow(kveh,:)=x(1,:);
    vFollow(kveh,:)=x(2,:);
    aFollow(kveh,:)=xp(2,:);
end
%% Plot of solution
% headways
subplot(3,1,1); hold on; box on;
plot(time,hFollow);
plot([t0,tend],[hstar,hstar],'k--');
axis([t0,tend,-2,12]);
xlabel('t');
ylabel('h(t)');
title(['Headways vs Time - vLdot and vdot Feedback',10,...
       'Parameters: ','sigma=',num2str(sigma,'%3.2f'),' [s]   ',...
       'k_p=',num2str(kp,'%3.2f'),'   ',...
       'k_v=',num2str(kv,'%3.2f'),'   ',...
       'k_a=',num2str(ka,'%3.2f'),'   ']);
driverlegend=cellfun(@(x,y)['#',num2str(x),': ',y],...
             num2cell(1:veh_num),driver,'UniformOutput',false);
legend(horzcat(driverlegend,{'h*'}),'Location','northeastoutside');

% velocities
subplot(3,1,2); hold on; box on;
plot(time,vLead(time),'k');
plot(time,vFollow);
plot([t0,tend],[vstar,vstar],'k--');
axis([t0,tend,vminplot,vmaxplot]);
xlabel('t');
ylabel('v(t)');
title('Velocity vs Time - vLdot and vdot Feedback');
driverlegend=cellfun(@(x,y)['#',num2str(x),': ',y],...
             num2cell(1:veh_num),driver,'UniformOutput',false);
legend(horzcat({'#0: leader'},driverlegend,{'v*'}),'Location','northeastoutside');

% acceleration
subplot(3,1,3); hold on; box on;
plot(time,vLdot(time),time,aFollow);
plot([t0,tend],[0,0],'k--');
axis([t0,tend,-10,4]);
xlabel('t');
ylabel('a(t)');
title('Acceleration vs Time - vLdot and vdot Feedback');
driverlegend=cellfun(@(x,y)['#',num2str(x),': ',y],...
             num2cell(1:veh_num),driver,'UniformOutput',false);
legend(horzcat({'#0: leader'},driverlegend,{'a*'}),'Location','northeastoutside');