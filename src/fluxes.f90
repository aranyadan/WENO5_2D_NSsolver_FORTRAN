module flux
contains
  function turn(F,n_x,n_y,times,direction)
    integer :: n_x,times,direction
    real,dimension(n_x,n_y,4) :: F,F_new,turn

    if(times == 0) then
      turn = F
    else
      if(direction==1) then
        do i = 1,n_x
          if (mod(i+times,n_x)>0) then
            F_new(mod(i+times,n_x),:,:) = F(i,:,:)
          else
            F_new(n_x,:,:) = F(i,:,:)
          end if
        end do
        turn = F_new
      else
        do i = 1,n_y
          if (mod(i+times,n_y)>0) then
            F_new(:,mod(i+times,n_y),:) = F(:,i,:)
          else
            F_new(:,n_y,:) = F(:,i,:)
          end if
        end do
        turn = F_new
      end if
    end if

    return
  end function turn

  function turnx(F,n_x,times)
    integer :: n_x,times
    real,dimension(n_x,4) :: F,F_new,turn
    if(times == 0) then
      turn = F
    else
      do i = 1,n_x
        if (mod(i+times,n_x)>0) then
          F_new(mod(i+times,n_x),:) = F(i,:)
        else
          F_new(n_x,:) = F(i,:)
        end if
      end do
      turn = F_new
    end if
    return
  end function turnx

  subroutine build_flux(q,n_x,n_y,delx,dely,Re,Pr,Suth,F,G)
    integer :: n_x,n_y
    real :: gamma=1.4,Re,Suth,Pr
    real, dimension(n_x,n_y,4) :: q,F, G
    real, dimension(n_x,n_y) :: u,v,p,rho,E,a,tauxx,tauxy,tauyy,q_x,q_y,T
    real, dimension(n_x,n_y,3) :: vel

    call primitives(q,n_x,n_y,rho,u,v,E,p,a)
    T = p/rho

    vel(:,:,1) = u(:,:)
    vel(:,:,2) = v(:,:)
    vel(:,:,3) = T(:,:)
    call viscous_fluxes(vel,p,rho,n_x,n_y,delx,dely,Suth,Re,Pr,tauxx,tauyy,tauxy,q_x,q_y)

    F(:,:,1) = rho*u
    F(:,:,2) = rho*u*u + p - tauxx
    F(:,:,3) = rho*u*v - tauxy
    F(:,:,4) = rho*u*E + p*u - u*tauxx - v*tauxy - q_x

    G(:,:,1) = rho*v
    G(:,:,2) = rho*u*v - tauxy
    G(:,:,3) = rho*v*v + p -tauyy
    G(:,:,4) = rho*v*E + p*v - u*tauxy - v*tauyy - q_y

    ! write(*,*)MAXVAL(MAXVAL(tauxx(:,:),1))
    ! write(*,*)MINVAL(MINVAL(rho(:,:),1))


  end subroutine build_flux

  ! function get_deriv(lambda,hp,hn,q,n_x,dx)
  !   integer :: n_x
  !   real :: lambda,dx
  !   real, dimension(n_x,3) :: hp,hn,res1,res2,res,get_deriv,q
  !   res1 = 0.5 * ( build_flux(hp,n_x) + build_flux(hn,n_x) - lambda*(hn - hp) )
  !   res2 = 0.5 * ( build_flux((turn(hp,n_x,1)),n_x) + build_flux((turn(hn,n_x,1)),n_x) - lambda*( turn(hn,n_x,1) - turn(hp,n_x,1)) )
  !   res = (res1-res2)/dx
  !
  !   get_deriv= res
  !
  ! end function get_deriv
end module flux
