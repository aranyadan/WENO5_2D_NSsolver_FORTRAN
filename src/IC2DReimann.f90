subroutine IC2DReimann(Prim,q,n_x,n_y,x,y,case_id,tend,Re,Pr,Suth,Cv)
  implicit none
  integer :: n_x,n_y,mid_x,mid_y,i,j, case_id
  real :: delx,del_y,gamma=1.4,tend,cfl,Re,Suth,p_ref,rho_ref,&
          T_ref,R_gas_const=287.0,Cp,Cv,Pr
  real, parameter :: PI = 4.0*ATAN(1.0)
  real, dimension(4) :: p,u,v,rho
  real,dimension(n_x,n_y,4) :: Prim,q
  real, dimension(n_x) :: x
  real, dimension(n_y) :: y
  real,dimension(n_x,n_y) :: E
  select case (case_id)
  case (1)              ! Lid driven cavity
    tend = 0.56
    cfl = 0.6
    p_ref = 101325             ! Reference air pressure (N/m^2)
    rho_ref= 1.225             ! Reference air density (kg/m^3)
    T_ref = p_ref / (rho_ref * R_gas_const);
    Cp = gamma * R_gas_const / (gamma-1);
    Cv = Cp - gamma;
    Re = 100.0;
    Suth = 110.4/T_ref;
    Pr = 0.7
    Prim(:,:,1) = 0.0
    Prim(:,:,2) = 0.0
    Prim(:,:,3) = 1.0
    Prim(:,:,4) = 1.0

  case (2)                  ! Taylor Green Vortex
    tend = 6
    cfl = 0.475
    p_ref = 101325             ! Reference air pressure (N/m^2)
    rho_ref= 1.225             ! Reference air density (kg/m^3)
    T_ref = p_ref / (rho_ref * R_gas_const);
    Cp = gamma * R_gas_const / (gamma-1);
    Cv = Cp - gamma;
    Re = 100;
    Suth = 110.4/T_ref;
    Pr = 0.7
    do j=1,n_y
      do i=1,n_x
        Prim(i,j,1) = SIN(PI*x(i)) * COS(PI*y(j));
        Prim(i,j,2) = -1.0 * COS(PI*x(i)) * SIN(PI*y(j))
        Prim(i,j,3) = 1.0 + (1.0/8.0) * (COS(2.0*PI*x(i)) + COS(2.0*PI*y(j)))
        Prim(i,j,4) = Prim(i,j,3)
      end do
    end do

  case (3)              ! Reimann Problem
    p =   (/ 1.5, 0.3,    0.029, 0.3    /)
    rho = (/ 1.5, 0.5323, 0.138, 0.5323 /)
    u =   (/ 0.0, 1.206,  1.206, 0.0    /)
    v =   (/ 0.0, 0.0,    1.206, 1.206  /)
    tend = 0.3
    mid_x =(n_x+1)/2
    mid_y =(n_y+1)/2
    cfl = 0.475
    p_ref = 101325             ! Reference air pressure (N/m^2)
    rho_ref= 1.225             ! Reference air density (kg/m^3)
    T_ref = p_ref / (rho_ref * R_gas_const);
    Cp = gamma * R_gas_const / (gamma-1);
    Cv = Cp - gamma;
    Re = 10000;
    Suth = 110.4/T_ref;
    Pr = 0.7
    ! Quadrant 1
    Prim(mid_x+1:n_x, mid_y+1:n_y, 1) = u(1)
    Prim(mid_x+1:n_x, mid_y+1:n_y, 2) = v(1)
    Prim(mid_x+1:n_x, mid_y+1:n_y, 3) = p(1)
    Prim(mid_x+1:n_x, mid_y+1:n_y, 4) = rho(1)
    ! Quadrant 2
    Prim(1:mid_x, mid_y+1:n_y, 1) = u(2)
    Prim(1:mid_x, mid_y+1:n_y, 2) = v(2)
    Prim(1:mid_x, mid_y+1:n_y, 3) = p(2)
    Prim(1:mid_x, mid_y+1:n_y, 4) = rho(2)
    ! Quadrant 3
    Prim(1:mid_x, 1:mid_y, 1) = u(3)
    Prim(1:mid_x, 1:mid_y, 2) = v(3)
    Prim(1:mid_x, 1:mid_y, 3) = p(3)
    Prim(1:mid_x, 1:mid_y, 4) = rho(3)
    ! Quadrant 4
    Prim(mid_x+1:n_x, 1:mid_y, 1) = u(4)
    Prim(mid_x+1:n_x, 1:mid_y, 2) = v(4)
    Prim(mid_x+1:n_x, 1:mid_y, 3) = p(4)
    Prim(mid_x+1:n_x, 1:mid_y, 4) = rho(4)

  case (4)              ! Couette Flow
    tend = 100
    cfl = 0.6
    p_ref = 101325             ! Reference air pressure (N/m^2)
    rho_ref= 1.225             ! Reference air density (kg/m^3)
    T_ref = p_ref / (rho_ref * R_gas_const);
    Cp = gamma * R_gas_const / (gamma-1);
    Cv = Cp - gamma;
    Re = 100.0;
    Suth = 110.4/T_ref;
    Pr = 0.7
    Prim(:,:,1) = 0.0
    Prim(:,:,2) = 0.0
    Prim(:,:,3) = 1.0
    Prim(:,:,4) = 1.0


  case (5)              ! Poissuelli Flow
    tend = 10
    cfl = 0.6
    p_ref = 101325             ! Reference air pressure (N/m^2)
    rho_ref= 1.225             ! Reference air density (kg/m^3)
    T_ref = p_ref / (rho_ref * R_gas_const);
    Cp = gamma * R_gas_const / (gamma-1);
    Cv = Cp - gamma;
    Re = 100.0;
    Suth = 110.4/T_ref;
    Pr = 0.7
    ! Set the primitive variables
    do i=1,n_x
      Prim(i,:,1) = (COS(PI*y(:)))**2
    end do
    Prim(:,:,2) = 0.0
    Prim(:,:,3) = 1.0
    Prim(:,:,4) = 1.0


  case(6)               ! Flat plat boundary layer
    tend = 5
    cfl = 0.6
    p_ref = 101325             ! Reference air pressure (N/m^2)
    rho_ref= 1.225             ! Reference air density (kg/m^3)
    T_ref = p_ref / (rho_ref * R_gas_const);
    Cp = gamma * R_gas_const / (gamma-1);
    Cv = Cp - gamma;
    Re = 10000.0;
    Suth = 110.4/T_ref;
    Pr = 0.7
    ! Set the primitive variables
    Prim(:,:,1) = 1.0
    Prim(:,:,2) = 0.0
    Prim(:,:,3) = 1.0
    Prim(:,:,4) = 1.0

  case(7)                 ! 2D Viscous shock tube SWBLI
    tend = 1
    cfl = 1.1
    p_ref = 101325             ! Reference air pressure (N/m^2)
    rho_ref= 1.225             ! Reference air density (kg/m^3)
    T_ref = p_ref / (rho_ref * R_gas_const);
    Cp = gamma * R_gas_const / (gamma-1);
    Cv = Cp - gamma;
    Re = 200.0;
    Suth = 110.4/T_ref;
    Pr = 0.72

    p =   (/ 120.0/gamma, 1.2/gamma,0.0,0.0   /)
    rho = (/ 120.0, 1.2,0.0,0.0  /)
    u =   (/ 0.0, 0.0,0.0,0.0    /)
    v =   (/ 0.0, 0.0,0.0,0.0  /)


    ! Set the primitive variables
    ! First half
    mid_x =(n_x+1)/2

    Prim(1:mid_x, :, 1) = u(1)
    Prim(1:mid_x, :, 2) = v(1)
    Prim(1:mid_x, :, 3) = p(1)
    Prim(1:mid_x, :, 4) = rho(1)
    ! 2nd Half
    Prim(mid_x+1:n_x, :, 1) = u(2)
    Prim(mid_x+1:n_x, :, 2) = v(2)
    Prim(mid_x+1:n_x, :, 3) = p(2)
    Prim(mid_x+1:n_x, :, 4) = rho(2)


  end select

  do j = 1,n_y
    do i = 1,n_x
      q(i,j,1) = Prim(i,j,4)
      q(i,j,2) = Prim(i,j,1) * Prim(i,j,4)
      q(i,j,3) = Prim(i,j,2) * Prim(i,j,4)
      E(i,j) = Prim(i,j,3) / ( (gamma-1)*Prim(i,j,4) ) &
               + 0.5 * (Prim(i,j,1)**2 + Prim(i,j,2)**2)
      q(i,j,4) = E(i,j) * Prim(i,j,4)
    end do
  end do
  ! Check IC setup
  ! open(unit=48,file='data')
  ! do i=1,n_x
  !   write(48,*)x(i),F(i,3)
  ! end do
  ! close(48)
  ! call system('gnuplot -p plotter.plt')
  return
end subroutine IC2DReimann
