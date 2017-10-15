subroutine WENO52d(lambda,F,q,n_x,n_y,hp,hn,dir)
  use flux
  use transform
  implicit none
  integer :: n_x,i,j,dir,k,n_y
  real,dimension(n_x,n_y,4) :: F,q,hp,hn,Ftemp,qtemp
  real,dimension(-1:n_x+3,-1:n_y+3,4) :: Fnew, qnew
  real :: lambda
  real,dimension(1) :: rho,u,v,E,p,a
  real,dimension(6,4) :: F_i,q_i,g_i,v_i
  real,dimension(4,4) :: R_i,Rinv_i
  real,dimension(1,4) :: q_h,hpr,hnr

  Fnew(1:n_x,1:n_y,:) = F(1:n_x,1:n_y,:)
  qnew(1:n_x,1:n_y,:) = q(1:n_x,1:n_y,:)

  qnew(-1,:,:) = qnew(1,:,:); qnew(0,:,:) = qnew(1,:,:)
  qnew(:,-1,:) = qnew(:,1,:); qnew(:,0,:) = qnew(:,1,:)
  qnew(n_x+1,:,:) = qnew(n_x,:,:); qnew(n_x+2,:,:) = qnew(n_x,:,:); qnew(n_x+3,:,:) = qnew(n_x,:,:);
  qnew(:,n_y+1,:) = qnew(:,n_y,:); qnew(:,n_y+2,:) = qnew(:,n_y,:); qnew(:,n_y+3,:) = qnew(:,n_y,:);

  Fnew(-1,:,:) = Fnew(1,:,:); Fnew(0,:,:) = Fnew(1,:,:)
  Fnew(:,-1,:) = Fnew(:,1,:); Fnew(:,0,:) = Fnew(:,1,:)
  Fnew(n_x+1,:,:) = Fnew(n_x,:,:); Fnew(n_x+2,:,:) = Fnew(n_x,:,:); Fnew(n_x+3,:,:) = Fnew(n_x,:,:);
  Fnew(:,n_y+1,:) = Fnew(:,n_y,:); Fnew(:,n_y+2,:) = Fnew(:,n_y,:); Fnew(:,n_y+3,:) = Fnew(:,n_y,:);


  if (dir==1) then
    do j=1,n_y
      do i=1,n_x
        F_i(1:6,:) = Fnew(i-2:i+3,j,:)
        q_i(1:6,:) = qnew(i-2:i+3,j,:)

        q_h(1,:) = (qnew(i,j,:) + qnew(i+1,j,:))/2.0

        call primitives(q_h,1,1,rho,u,v,E,p,a)
        R_i = Rcalc(u(1),a(1))
        Rinv_i = Rinv(u(1),a(1))
        call WENO(lambda,F_i,q_i,R_i,Rinv_i,hpr,hnr)
        hp(i,j,:) = hpr(1,:)
        hn(i,j,:) = hnr(1,:)
      end do
    end do
  else
    do j=1,n_y
      do i=1,n_x
        F_i(1:6,:) = Fnew(i,j-2:j+3,:)
        q_i(1:6,:) = qnew(i,j-2:j+3,:)

        q_h(1,:) = (qnew(i,j,:) + qnew(i,j+1,:))/2.0

        call primitives(q_h,1,1,rho,u,v,E,p,a)
        R_i = Rcalc(u(1),a(1))
        Rinv_i = Rinv(u(1),a(1))
        call WENO(lambda,F_i,q_i,R_i,Rinv_i,hpr,hnr)
        hp(i,j,:) = hpr(1,:)
        hn(i,j,:) = hnr(1,:)
      end do
    end do
  end if





  return
end subroutine WENO52d
