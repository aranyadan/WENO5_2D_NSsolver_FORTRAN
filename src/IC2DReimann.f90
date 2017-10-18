subroutine IC2DReimann(Prim,q,n_x,n_y,x,y,case_id,tend,Re,Suth)
  implicit none
  integer :: n_x,n_y,mid_x,mid_y,i,j, case_id
  real :: delx,del_y,gamma=1.4,tend,cfl,Re,Suth,p_ref,rho_ref,&
          T_ref,R_gas_const=287.0,Cp,Cv
  real, dimension(4) :: p,u,v,rho
  real,dimension(n_x,n_y,4) :: Prim,q
  real, dimension(n_x) :: x
  real, dimension(n_y) :: y
  real,dimension(n_x,n_y) :: E

  tend = 5.0
  cfl = 0.6
  p_ref = 101325             ! Reference air pressure (N/m^2)
  rho_ref= 1.225             ! Reference air density (kg/m^3)
  T_ref = p_ref / (rho_ref * R_gas_const);
  Cp = gamma * R_gas_const / (gamma-1);
  Cv = Cp - gamma;
  Re = 200;
  Suth = 110.4/T_ref;




  ! Set the primitive variables

  Prim(:,:,1) = 0.0
  Prim(:,:,2) = 0.0
  Prim(:,:,3) = 1.0
  Prim(:,:,4) = 1.0



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
