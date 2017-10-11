subroutine primitives(q,n_x,n_y,rho,u,v,E,p,a)
  integer :: n_x,n_y
  real,dimension(n_x,n_y,4) :: q
  real :: gamma = 1.4
  real,dimension(n_x,n_y) :: rho,u,v,E,p,a

  rho = q(:,:,1)
  u = q(:,:,2)/rho
  v = q(:,:,3)/rho
  E = q(:,:,4)/rho
  p = (gamma-1)*rho*(E-0.5*(u*u+v*v))
  a = SQRT(gamma*p/rho)

end subroutine primitives
