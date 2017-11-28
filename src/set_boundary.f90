subroutine set_boundary(q,x,y,n_x,n_y,Cv,case_id)
  real, dimension(n_x,n_y,4) :: q
  real,dimension(n_x,n_y) :: rho,u,v,E,p,a
  real, parameter :: PI = 4.0*ATAN(1.0)
  real, dimension(n_x) :: x
  real, dimension(n_y) :: y
  real :: Cv, gamma=1.4
  integer :: n_x,n_y,case_id, temp

  call primitives(q,n_x,n_y,rho,u,v,E,p,a)

  select case (case_id)
  case(1)
    u(:,n_y) = (1.0 - (2.0*x(:)-1.0)**18)**2
    v(:,n_y) = 0
    E(:,n_y) = 1.0/(gamma-1.0)

    u(:,1) = -1.0*u(:,2)
    v(:,1) = -1.0*v(:,2)
    E(:,1) = 1.0/(gamma-1.0)

    u(1,:) = 0
    v(1,:) = 0
    E(1,:) = 1.0/(gamma-1.0)

    u(n_x,:) = 0
    v(n_x,:) = 0
    E(n_x,:) = 1.0/(gamma-1.0)


  case(2)

  case(3)
    q(1,:,:)=q(2,:,:)
    q(n_x,:,:)=q(n_x-1,:,:)
    q(:,1,:)=q(:,2,:)
    q(:,n_y,:)=q(:,n_y-1,:)

  case(4)
    u(:,n_y) = 1
    v(:,n_y) = 0
    E(:,n_y) = 1.0/(gamma-1.0)

    u(:,1) = -1.0*u(:,2)
    v(:,1) = -1.0*v(:,2)
    E(:,1) = 1.0/(gamma-1.0)

  case(5)
    ! Inflow
    u(1,:) = (COS(PI*y(:)))**2
    v(1,:) = 0
    E(1,:) = 1.0/(gamma-1.0) + 0.5 * u(1,:) * u(1,:)

    ! Walls
    u(:,n_y) = 0
    v(:,n_y) = 0
    E(:,n_y) = 1.0/(gamma-1.0)

    u(:,1) = 0
    v(:,1) = 0
    E(:,1) = 1.0/(gamma-1.0)

    ! Outflow
    E(n_x,:) = 1.0/((gamma-1.0)*rho(n_x,:)) + 0.5*(u(n_x,:)**2 + v(n_x,:)**2)

  case(6)
    ! Inflow
    rho(1,:) = 1.0
    u(1,:) = 1.0
    v(1,:) = 0
    E(1,:) = 1.0/(gamma-1.0) + 0.5 * u(1,:) * u(1,:)

    ! Upper Wall
    rho(:,n_y) = 1.0
    u(:,n_y) = 1.0
    v(:,n_y) = 0
    E(:,n_y) = 1.0/(gamma-1.0) + 0.5 * u(1,:) * u(1,:)
    ! u(:,n_y) = 0
    ! v(:,n_y) = 0
    ! E(:,n_y) = 1.0/(gamma-1.0)

    ! Lower wall
    temp = NINT(0.3*n_x)
    u(1:temp,1) = u(1:temp,2)
    ! E(temp:n_x,2) = 1.0/(gamma-1.0)
    rho(temp:n_x,2) = 1.0
    u(temp:n_x,1) = -1.0 * u(temp:n_x,2)
    v(temp:n_x,1) = -1.0 * v(temp:n_x,2)
    E(temp:n_x,1) = E(temp:n_x,2)
    rho(temp:n_x,1) = rho(temp:n_x,2)

    ! Outflow
    ! E(n_x,:) = 1.0/((gamma-1.0)*rho(n_x,:)) + 0.5*(u(n_x,:)**2 + v(n_x,:)**2)



  end select

  q(:,:,1) = rho(:,:)
  q(:,:,2) = rho(:,:)*u(:,:)
  q(:,:,3) = rho(:,:)*v(:,:)
  q(:,:,4) = rho(:,:)*E(:,:)




end subroutine set_boundary
