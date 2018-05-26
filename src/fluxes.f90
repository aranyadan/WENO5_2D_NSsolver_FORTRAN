subroutine build_flux(q,n_x,n_y,delx,dely,Re,Pr,Suth,F,G,VISCOUS)
  integer :: n_x,n_y,VISCOUS
  real :: gamma=1.4,Re,Suth,Pr
  real, dimension(n_x,n_y,4) :: q,F, G
  real, dimension(n_x,n_y) :: u,v,p,rho,E,a,tauxx,tauxy,tauyy,q_x,q_y,T
  real, dimension(n_x,n_y,3) :: vel

  call primitives(q,n_x,n_y,rho,u,v,E,p,a)
  T = p/rho

  if(VISCOUS==1) then
    vel(:,:,1) = u(:,:)
    vel(:,:,2) = v(:,:)
    vel(:,:,3) = T(:,:)
    call viscous_fluxes(vel,p,rho,n_x,n_y,delx,dely,Suth,Re,Pr,tauxx,tauyy,tauxy,q_x,q_y)
  endif

  F(:,:,1) = rho*u
  F(:,:,2) = rho*u*u + p + (- tauxx)*VISCOUS
  F(:,:,3) = rho*u*v + (- tauxy)*VISCOUS
  F(:,:,4) = rho*u*E + p*u + (- u*tauxx - v*tauxy - q_x)*VISCOUS

  G(:,:,1) = rho*v
  G(:,:,2) = rho*u*v + (- tauxy)*VISCOUS
  G(:,:,3) = rho*v*v + p + (-tauyy)*VISCOUS
  G(:,:,4) = rho*v*E + p*v + (- u*tauxy - v*tauyy - q_y)*VISCOUS

  ! write(*,*)MAXVAL(MAXVAL(tauxx(:,:),1))
  ! write(*,*)MINVAL(MINVAL(rho(:,:),1))


end subroutine build_flux
