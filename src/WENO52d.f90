subroutine WENO52d(lambda,F,q,n_x,n_y,hp,hn,dir)
  use flux
  use transform
  implicit none
  integer :: n_x,i,j,dir,k,n_y
  real,dimension(n_x,n_y,4) :: F,q,hp,hn,Ftemp,qtemp
  real :: lambda
  real,dimension(1) :: rho,u,v,E,p,a
  real,dimension(6,4) :: F_i,q_i,g_i,v_i
  real,dimension(4,4) :: R_i,Rinv_i
  real,dimension(1,4) :: q_h,hpr,hnr


  do j=1,n_y
    do i=1,n_x
      do k = 6,1,-1
        Ftemp = turn(F,n_x,n_y,3-k,dir)
        qtemp = turn(q,n_x,n_y,3-k,dir)
        F_i(k,:) = Ftemp(i,j,:)
        q_i(k,:) = qtemp(i,j,:)
      end do

      if(dir==1) then
        if(i == n_x) then
          q_h(1,:) = (q(i,j,:) + q(1,j,:))/2.0
        else
          q_h(1,:) = (q(i,j,:) + q(i+1,j,:))/2.0
        end if
      else
        if(j == n_y) then
          q_h(1,:) = (q(i,j,:) + q(i,1,:))/2.0
        else
          q_h(1,:) = (q(i,j,:) + q(i,j+1,:))/2.0
        end if
      end if

      call primitives(q_h,1,1,rho,u,v,E,p,a)
      R_i = Rcalc(u(1),a(1))
      Rinv_i = Rinv(u(1),a(1))
      print*, 'Sending ',i,',',j
      call WENO(lambda,F_i,q_i,R_i,Rinv_i,hpr,hnr)
      hp(i,j,:) = hpr(1,:)
      hn(i,j,:) = hnr(1,:)
    end do
  end do


  return
end subroutine WENO52d
