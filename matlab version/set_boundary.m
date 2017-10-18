%% Function to set BC
function q0 = set_boundary(q)
    global gamma
    rho=q(:,:,1); u=q(:,:,2)./rho; v=q(:,:,3)./rho; E=q(:,:,4)./rho; p=(gamma-1)*rho.*(E-0.5*(u.^2+v.^2));

    u(:,end) = 1;
    v(:,end) = 0;

    u(:,1) = 0;
    v(:,1) = 0;

    u(1,:) = 0;
    v(1,:) = 0;

    u(end,:) = 0;
    v(end,:) = 0;

    q0=reshape([rho rho.*u rho.*v rho.*E],size(q));
end
