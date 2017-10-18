subroutine set_boundary(q,n_x,n_y)
  real, dimension(n_x,n_y,4) :: q
  real,dimension(n_x,n_y) :: rho,u,v,E,p,a
  integer :: n_x,n_y

  call primitives(q,n_x,n_y,rho,u,v,E,p,a)

  u(:,n_y) = 1
  v(:,n_y) = 0

  u(:,1) = 0
  v(:,1) = 0

  u(1,:) = 0
  v(1,:) = 0

  u(n_x,:) = 0
  v(n_x,:) = 0

  q(:,:,1) = rho(:,:)
  q(:,:,2) = rho(:,:)*u(:,:)
  q(:,:,3) = rho(:,:)*v(:,:)
  q(:,:,4) = rho(:,:)*E(:,:)


end subroutine set_boundary
