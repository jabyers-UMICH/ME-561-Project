%% Create stability charts for CCC with ZOH
clc
clear; %close all; clc;
clearvars
%% Parameters
kappa=0.4619; deltat=0.1; ka = 0.4;
% range of parameters
kvmin=-2; kvmax=3; kpmin=-1; kpmax=4; kvst=100; kpst=100;
% range of frequencies
ommin=0; ommax=2*pi/deltat;
dom=2*pi/deltat/2000; om_v=ommin:dom:ommax;   

% grid in plane of parameters
kv_v=linspace(kvmin,kvmax,kvst+1);kp_v=linspace(kpmin,kpmax,kpst+1);    
kv_m = zeros(length(kv_v),length(kp_v));kp_m = zeros(size(kv_m)); 
% exponential of frequencies
expom=exp(1i*om_v*deltat);
% list of parameters to put on figures
parlist=['\kappa=',num2str(kappa,'%3.2f'),' [1/s]   ',...
         '\Deltat=',num2str(deltat,'%3.2f'),' [s]'];

%% Stability analysis\
% coefficient matrices
a0=[1,-deltat;0,1]; b0=[deltat/2;0]; b1=[deltat/2;0]; c=[0,1];
% construction of parameter-independent parts of system matrices
A=zeros(10);
B0=zeros(10,1); B1=zeros(10,1); Bsig=zeros(10,1); Bsig1=zeros(10,1); C=zeros(1,10); II=eye(10);
A(1:2,1:2)=a0; A(3:end,1:end-2)=eye(8); B0(1:2)=b0; B1(1:2)=b1; C(1:2)=c;

% matrix of eigenvalues
eigall = zeros(length(kv_v),length(kp_v),size(A,1));
% matrix of transfer function values
T = zeros(length(kv_v),length(kp_v),length(om_v));

% computation of eigenvalues and frequency response
for kkv=1:length(kv_v)
    kv=kv_v(kkv);
    for kkp=1:length(kp_v)
        kp=kp_v(kkp);
        % construction of remaining parts of system matrices
        asig=[-kp*kappa*deltat^2/2,	 (kp+kv)*deltat^2/2;
               kp*kappa*deltat,     -(kp+kv)*deltat];
        bsig=[-(kv*deltat+ka)*deltat/2;kv*deltat+ka];
        bsig1=[ka*deltat/2; -ka];
        A(1:2,end-3:end-2)=asig;
        Bsig(1:2)=bsig;
        Bsig1(1:2)=bsig1;
        % computation of eigenvalues
        eigall(kkv,kkp,:) = eig(A);
        % computation of transfer function values
        warning('off');
        T(kkv,kkp,:) = arrayfun(@(z)C/(z*II-A)*(Bsig1*z^(-4)+Bsig*z^(-3)+B0+B1*z),expom);
        % store parameter values
        kv_m(kkv,kkp) = kv;
        kp_m(kkv,kkp) = kp;
    end
    % show progress of computations
    disp(['Computation of row #',num2str(kkv),'/',...
            num2str(length(kv_v)),' of the stability chart.']);
end

% matrix of critical eigenvalues
eigcr = max(abs(eigall),[],3);
% matrix of maximum transfer function magnitudes
Mmax = max(abs(T(:,:,0<om_v & om_v<2*pi/deltat)),[],3);

%% Plot stability chart
figure(1); clf; hold on; box on;
% plant stability
contourf(kv_m,kp_m,-eigcr,-[1 1],'r','Linewidth',1.5);
% string stability
contourf(kv_m,kp_m,-Mmax,-[1 1],'b','Linewidth',1.5);
colormap gray
axis([kvmin kvmax kpmin kpmax]);
pbaspect([1 1 1]);
xlabel('k_V [1/s]');
ylabel('k_P [1/s]');
title(['Stability chart of CCC with ZOH',10,parlist]);
legend('plant stability','string stability','Location','southeast');

%% Plot eigenvalues and frequency response
% select parameters of interest
figure(1); [kv0,kp0]=ginput(1);    % by clicking on stability chart
% kv0=0.4; kp0=0.2;                % or by giving manually

% get closest parameter values to the selected ones
kv0=interp1(kv_v,kv_v,kv0,'nearest');
kp0=interp1(kp_v,kp_v,kp0,'nearest');

% get corresponding eigenvalues and frequency response
eig_v=squeeze(eigall(kv_v==kv0,kp_v==kp0,:));
M_v=squeeze(abs(T(kv_v==kv0,kp_v==kp0,:)));

% plot eigenvalues
figure(2); clf; hold on; box on;
plot(real(eig_v),imag(eig_v),'rx','Linewidth',2);
plot(cos(0:2*pi/100:2*pi),sin(0:2*pi/100:2*pi),'k--');    % unit circle
pbaspect([1 1 1]);
xlabel('Re');
ylabel('Im');
title(['Eigenvalues plot of CCC with ZOH and time delay',10,parlist,...
       '   k_V=',num2str(kv0,'%3.2f'),' [1/s]',...
       '   k_P=',num2str(kp0,'%3.2f'),' [1/s]']);

% plot frequency response
figure(3); clf; hold on; box on;
plot(om_v,M_v,'b');
plot([ommin ommax],[1 1],'k--');    % plot unit gain
xlim([ommin ommax]);
xlabel('\omega [rad/s]');
ylabel('|T(e^{j\omega\Deltat})|');
title(['Frequency response of CCC with ZOH and time delay',10,parlist,...
       '   k_V=',num2str(kv0,'%3.2f'),' [1/s]',...
       '   k_P=',num2str(kp0,'%3.2f'),' [1/s]']);