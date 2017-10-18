clear; %close all; clc;
global gamma

%% Parameters
CFL     = 0.75;	% CFL number
tFinal	= 0.80;	% Final time
nEx      = 100;  % Number of cells in x
nEy      = 100;  % Number of cells in y
gamma   = 1.4;  % Ratio of specific heats for ideal di-atomic gas
IC      = 1;	% 10 IC cases are available
plot_fig= 1;

% Discretize spatial domain
a=0; b=1; dx=(b-a)/nEx; nx=nEx+1; x=linspace(a,b,nx);
c=0; d=1; dy=(d-c)/nEy; ny=nEy+1; y=linspace(c,d,ny);
[Y,X] = meshgrid(y,x);


% Set IC
[rho0,u0,v0,p0,tFinal,CFL] = Euler_IC2d(X,Y,IC);
E0 = p0./((gamma-1)*rho0)+0.5*(u0.^2+v0.^2);  % Total Energy density
a0 = sqrt(gamma*p0./rho0);            % Speed of sound
q0=reshape([rho0 rho0.*u0 rho0.*v0 rho0.*E0],nx,ny,4);        % vec. of conserved properties


% Discretize time domain
lambda0=max(max(abs(u0)+abs(v0)+a0)); dt0=CFL*dx/lambda0;  % using the system's largest eigenvalue

%% Solver Loop
% Load initial condition
q=q0; it=0; dt=dt0; t=0; lambda=lambda0;

%% Solver Loop
figure;
% h=createButton;
while t<tFinal

    % RK Initial step
    qo = q;

    % 1st stage
    dFx=WENO5LF2d(lambda,q,dx,1);
    dFy=WENO5LF2d(lambda,q,dy,2);
    dF = dFx+dFy;
    q = qo-dt*dF;
    q(1,:,:)=q0(1,:,:); q(end,:,:)=q0(end,:,:); q(:,1,:)=q0(:,1,:); q(:,end,:)=q0(:,end,:);
    q(1,:,:)=q(2,:,:); q(end,:,:)=q(end-1,:,:); q(:,1,:)=q(:,2,:); q(:,end,:)=q(:,end-1,:);


    % 2nd Stage
    dFx=WENO5LF2d(lambda,q,dx,1);
    dFy=WENO5LF2d(lambda,q,dy,2);
    dF = dFx+dFy;
    q = 0.75*qo+0.25*(q-dt*dF);
    q(1,:,:)=q0(1,:,:); q(end,:,:)=q0(end,:,:); q(:,1,:)=q0(:,1,:); q(:,end,:)=q0(:,end,:);
    q(1,:,:)=q(2,:,:); q(end,:,:)=q(end-1,:,:); q(:,1,:)=q(:,2,:); q(:,end,:)=q(:,end-1,:);

    % 3rd stage
    dFx=WENO5LF2d(lambda,q,dx,1);
    dFy=WENO5LF2d(lambda,q,dy,2);
    dF = dFx+dFy;
    q = (qo+2*(q-dt*dF))/3;
    q(1,:,:)=q0(1,:,:); q(end,:,:)=q0(end,:,:); q(:,1,:)=q0(:,1,:); q(:,end,:)=q0(:,end,:);
    q(1,:,:)=q(2,:,:); q(end,:,:)=q(end-1,:,:); q(:,1,:)=q(:,2,:); q(:,end,:)=q(:,end-1,:);

    % compute primary properties
    rho=q(:,:,1); u=q(:,:,2)./rho; v=q(:,:,3)./rho; E=q(:,:,4)./rho; p=(gamma-1)*rho.*(E-0.5*(u.^2+v.^2));
    a=sqrt(gamma*p./rho); if min(p)<0; error('negative pressure found!'); end

    % Update dt and time
    lambda=max(max(abs(u)+abs(v)+a)); dt=CFL*dx/lambda; if t+dt>tFinal; dt=tFinal-t; end

    % Update time and iteration counter
	t=t+dt; it=it+1;


%     cla;plot(x,u); grid on; hold on; scatter(x,u);
    cla;fdisplay(X,Y,rho);title(['t=' num2str(t)]);
    pause(0.001)

end

% % Calculation of flow parameters
% a = sqrt(gamma*p./rho); M = u./a; % Mach number [-]
% p_ref = 101325;             % Reference air pressure (N/m^2)
% rho_ref= 1.225;             % Reference air density (kg/m^3)
% s = 1/(gamma-1)*(log(p/p_ref)+gamma*log(rho_ref./rho));
%                             % Entropy w.r.t reference condition
% ss = log(p./rho.^gamma);    % Dimensionless Entropy
% Q = rho.*u;                 % Mass Flow rate per unit area
% e = p./((gamma-1)*rho);     % internal Energy
