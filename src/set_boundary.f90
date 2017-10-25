subroutine set_boundary(q,x,n_x,n_y,Cv)
  real, dimension(n_x,n_y,4) :: q
  real,dimension(n_x,n_y) :: rho,u,v,E,p,a
  real, dimension(n_x) :: x
  real :: Cv, gamma=1.4
  integer :: n_x,n_y

  ! call primitives(q,n_x,n_y,rho,u,v,E,p,a)
  !
  ! u(:,n_y) = (1.0 - (2.0*x(:)-1.0)**18)**2
  ! v(:,n_y) = 0
  ! E(:,n_y) = 1.0/(gamma-1.0)
  !
  ! u(:,1) = 0
  ! v(:,1) = 0
  ! E(:,1) = 1.0/(gamma-1.0)
  !
  ! u(1,:) = 0
  ! v(1,:) = 0
  ! E(1,:) = 1.0/(gamma-1.0)
  !
  ! u(n_x,:) = 0
  ! v(n_x,:) = 0
  ! E(n_x,:) = 1.0/(gamma-1.0)
  !
  ! q(:,:,1) = rho(:,:)
  ! q(:,:,2) = rho(:,:)*u(:,:)
  ! q(:,:,3) = rho(:,:)*v(:,:)
  ! q(:,:,4) = rho(:,:)*E(:,:)



  q(1,:,:)=q(2,:,:)
  q(n_x,:,:)=q(n_x-1,:,:)
  q(:,1,:)=q(:,2,:)
  q(:,n_y,:)=q(:,n_y-1,:)

end subroutine set_boundary
