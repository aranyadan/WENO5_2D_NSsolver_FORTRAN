function [dFx,dFy] = WENOcaller(lambda,q,dx,dy)
%% Calls WENO

% Setting up ghost cells
[nx,ny,eqns] = size(q);
qnew = zeros(nx+5,ny+5,eqns);
qnew(3:2+nx,3:2+ny,:) = q(:,:,:);
qnew(1:2,:,:) = repmat(qnew(3,:,:),[2,1,1]); qnew(nx+3:nx+5,:,:) = repmat(qnew(nx+2,:,:),[3,1,1]);
qnew(:,1:2,:) = repmat(qnew(:,3,:),[1,2,1]); qnew(:,ny+3:ny+5,:) = repmat(qnew(:,ny+2,:),[1,3,1]);
w = qnew;

% Constructing fluxes to send
[G1,G2] = F(qnew,dx,dy);

% Sending for derivative along x
turn=[1 0 0];
% G = F(qnew,1);
dFx=WENO5LF1d(lambda,qnew,G1,dx,turn);
dFx = dFx(3:2+nx,3:2+ny,:);

% Sending for derivative along y
turn=[0 1 0];
% G = F(qnew,2);
dFy=WENO5LF1d(lambda,qnew,G2,dy,turn);
dFy = dFy(3:2+nx,3:2+ny,:);



end


%% Compute flux vector
function [G1,G2] = F(q,dx,dy)
    global gamma
    global R_gas_const
    % primary properties
    rho=q(:,:,1); u=q(:,:,2)./rho; v=q(:,:,3)./rho; E=q(:,:,4)./rho;
    p=(gamma-1)*rho.*(E-0.5*(u.^2 + v.^2));

    temp = p./(rho.*R_gas_const);
    [tauxx,tauyy,tauxy] = visc_stress(u,v,temp,dx,dy);


    % flux vector of conserved properties

    G1=reshape([rho.*u (rho.*u.^2+p-tauxx) (rho.*u.*v-tauxy)    (u.*(rho.*E+p-tauxx)-v.*tauxy)], size(q));

    G2=reshape([rho.*v (rho.*u.*v-tauxy)   (rho.*v.^2+p-tauyy)  (v.*(rho.*E+p-tauyy)-u.*tauxy)], size(q));

end


%% Compute viscous stress terms
function [tauxx,tauyy,tauxy] = visc_stress(u,v,T,dx,dy)
    global Suth_const
    global Re

    [u_x,u_y] = first_der(u,dx,dy);
    [v_x,v_y] = first_der(v,dx,dy);

    meu = 1;%(T.^(3/2)) .* ((1 + Suth_const)./(T + Suth_const));

    lambda = (-2/3) .* meu;

    tauxx = (lambda.*(u_x+v_y) + 2.*meu.*u_x)./Re;
    tauyy = (lambda.*(u_x+v_y) + 2.*meu.*v_y)./Re;
    tauxy = (meu.*(u_y+v_x))./Re;
end

%% Compute first derivatives
function [u_x,u_y] = first_der(u,dx,dy)
u_x = [ (1/20).*circshift(u,[2,0]) + (-0.5).*circshift(u,[1,0]) + (-1/3).*u + ...
		circshift(u,[-1,0]) + (-0.25).*circshift(u,[-2,0]) + (1/30).*circshift(u,[-3,0]) ]./dx;
u_y = [ (1/20).*circshift(u,[0,2]) + (-0.5).*circshift(u,[0,1]) + (-1/3).*u + ...
        circshift(u,[0,-1]) + (-0.25).*circshift(u,[0,-2]) + (1/30).*circshift(u,[0,-3]) ]./dy;
end
