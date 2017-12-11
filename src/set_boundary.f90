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
  case(1)              ! Lid driven cavity
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


  case(2)                  ! Taylor Green Vortex

  case(3)              ! Reimann Problem
    q(1,:,:)=q(2,:,:)
    q(n_x,:,:)=q(n_x-1,:,:)
    q(:,1,:)=q(:,2,:)
    q(:,n_y,:)=q(:,n_y-1,:)

  case(4)              ! Couette Flow
    u(:,n_y) = 1
    v(:,n_y) = 0
    E(:,n_y) = 1.0/(gamma-1.0)

    u(:,1) = -1.0*u(:,2)
    v(:,1) = -1.0*v(:,2)
    E(:,1) = 1.0/(gamma-1.0)

  case(5)              ! Poissuelli Flow
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

  case(6)               ! Flat plat boundary layer
    ! Inflow
    rho(1,:) = 1.0
    u(1,:) = 1.0
    v(1,:) = 0
    p(1,:) = 1.0

    ! Upper Wall
    rho(:,n_y) = 1.0
    ! u(:,n_y) = 1.0
    ! v(:,n_y) = 0
    p(:,n_y) = 1.0

    ! Lower wall

    temp = NINT(0.1*n_x)
    u(1:temp,1) = u(1:temp,2)
    u(temp+1:n_x,1) = -1.0 * u(temp+1:n_x,2)
    v(:,1) = -1.0 * v(:,2)
    p(:,1) = p(:,2)
    rho(:,1) = (rho(:,2)/p(:,2))*p(:,1)

    ! Outflow
    u(n_x,:) = u(n_x-1,:)
    v(n_x,:) = v(n_x-1,:)
    rho(n_x,:) = (rho(n_x-1,:)/p(n_x-1,:))*p(n_x,:)
    p(n_x,:) = p(n_x-1,:)
    ! rho(n_x,:) = rho(n_x-1,:)

  case(7)                 ! 2D Viscous shock tube SWBLI
    ! Symmetry
    u(:,n_y) = u(:,n_y-1)
    v(:,n_y) = v(:,n_y-1)
    p(:,n_y) = p(:,n_y-1)
    rho(:,n_y) = rho(:,n_y-1)

    ! Noslip
    u(1,:) = 0
    v(1,:) = 0
    u(n_x,:) = 0
    v(n_x,:) = 0
    u(:,1) = 0
    v(:,1) = 0

    !Adiabatic
    p(1,:) = (p(2,:)/rho(2,:))*rho(1,:)
    p(n_x,:) = (p(n_x-1,:)/rho(n_x-1,:))*rho(n_x,:)
    p(:,1) = (p(:,2)/rho(:,2))*rho(:,1)


  end select

  E = p/((gamma-1.0)*rho) + 0.5*(u**2+v**2)
  q(:,:,1) = rho(:,:)
  q(:,:,2) = rho(:,:)*u(:,:)
  q(:,:,3) = rho(:,:)*v(:,:)
  q(:,:,4) = rho(:,:)*E(:,:)




end subroutine set_boundary
