subroutine viscous_fluxes(vel,p,rho,n_x,n_y,delx,dely,Suth,Re,tauxx,tauyy,tauxy)
  implicit none
  real, dimension(n_x,n_y,2) :: vel     ! 1-xvel, 2-yvel
  real, dimension(n_x,n_y) :: p,rho,tauxx,tauyy,tauxy
  real, dimension(-1:n_x+3,-1:n_y+3,2) :: velnew
  real :: delx,dely,Suth,R_gas_const = 287.0, T,meu,lambda,Re
  real :: u_x,u_y,v_x,v_y
  real, dimension(1,2) :: temp_der,temp_der2
  real, dimension(6,2) :: temp,temp2
  integer :: n_x,n_y,i,j

  call pad_fluxes(vel,velnew,n_x,n_y,2)
  ! Calculate the derivatives in x and y directions
  do j=1,n_y
    do i =1,n_x

      temp(1:6,:) = velnew(i-2:i+3,j,:)
      call get_first_deriv(temp,delx,temp_der)
      u_x = temp_der(1,1); v_x = temp_der(1,2)

      temp(1:6,:) = velnew(i,j-2:j+3,:)
      call get_first_deriv(temp,dely,temp_der)
      u_y = temp_der(1,1); v_y = temp_der(1,2)

      T = p(i,j)/(rho(i,j)*R_gas_const)
      meu = (T**(3.0/2.0)) * ( (1+Suth)/(T+Suth) )
      lambda = (-2.0/3.0)*meu

      tauxx(i,j) = (lambda*(u_x+v_y) + 2*meu*u_x)/Re
      tauyy(i,j) = (lambda*(u_x+v_y) + 2*meu*v_y)/Re
      tauxy(i,j) = (meu*(u_y+v_x))/Re

    end do
  end do

end subroutine viscous_fluxes


subroutine get_first_deriv(vel,delta,deriv)
  implicit none
  real, dimension(6,2) :: vel
  real :: delta
  real, dimension(1,2) :: deriv

  deriv(1,:) = ( (1.0/20.0)*vel(1,:) + (-0.5)*vel(2,:) + (-1.0/3.0)*vel(3,:) &
               +	vel(4,:) + (-0.25)*vel(5,:) + (1.0/30.0)*vel(6,:) )/delta

end subroutine get_first_deriv
