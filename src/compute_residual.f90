subroutine compute_residual(q,qo,residual,n_x,n_y)
  real, dimension(n_x,n_y,4) :: q,qo
  real, dimension(n_x,n_y) :: u,v,E,p,a,u2,rho
  real, dimension(4) :: res
  real :: residual
  integer :: dim

  call primitives(q,n_x,n_y,rho,u,v,E,p,a)
  call primitives(qo,n_x,n_y,rho,u2,v,E,p,a)
  do dim=1,1
    res(dim) = MAXVAL(MAXVAL(ABS(u-u2),1));
  end do
  residual = res(1);


end subroutine compute_residual
