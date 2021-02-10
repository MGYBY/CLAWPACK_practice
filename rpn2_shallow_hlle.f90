! =====================================================
subroutine rpn2(ixy,maxm,meqn,mwaves,maux,mbc,mx,ql,qr,auxl,auxr,wave,s,amdq,apdq)
! =====================================================

! HLLC solver for the 2D shallow water equations.
! reconstruction of the middle state which has two states

! waves: 3
! equations: 3

! Conserved quantities:
!       1 depth
!       2 x-momentum
!       3 y-momentum

    implicit none

    integer, intent(in) :: ixy, maxm, meqn, mwaves, maux, mbc, mx
    double precision, dimension(meqn,1-mbc:maxm+mbc), intent(in) :: ql, qr
    double precision, dimension(maux,1-mbc:maxm+mbc), intent(in) :: auxl, auxr
    double precision, dimension(meqn, mwaves, 1-mbc:maxm+mbc), intent(out) :: wave
    double precision, dimension(mwaves, 1-mbc:maxm+mbc), intent(out) :: s
    double precision, dimension(meqn, 1-mbc:maxm+mbc), intent(out) :: amdq, apdq

    double precision :: u_l, u_r, c_l, c_r, u_hat, c_hat, v_l, v_r, v_hat, h_hat
    double precision :: h_l, h_r, hsqrt_l, hsqrt_r, hsq2
    double precision :: grav
    double precision :: h_m, hu_m, hv_m, s1, s2
    double precision :: h_ml, hu_ml, hv_ml, h_mr, hu_mr, hv_mr
    double precision :: qqleft, qqright, hstar, sstar, s_l, s_r
    integer :: depth, mu, mv
    integer :: i, m, mw

    common /cparam/  grav

!     # set mu to point to  the component of the system that corresponds
!     # to momentum in the direction of this slice, mv to the orthogonal
!     # momentum:
!

    depth = 1
    if (ixy.eq.1) then
        mu = 2
        mv = 3
    else
        mu = 3
        mv = 2
    endif

    do i=2-mbc,mx+mbc
        h_l = qr(depth,i-1)
        h_r = ql(depth,i)
        ! Velocity
        u_l = qr(mu,i-1) / qr(depth,i-1)
        u_r = ql(mu,i  ) / ql(depth,i  )
        v_l = qr(mv,i-1) / qr(depth,i-1)
        v_r = ql(mv,i  ) / ql(depth,i  )
        ! Sound speed
        c_l = dsqrt(grav*h_l)
        c_r = dsqrt(grav*h_r)

        hstar = 1.0/2.0*(h_l+h_r)-(1.0/4.0*(u_r-u_l)*(h_l+h_r))/(c_l+c_r)
        if (hstar>h_l) then
            qqleft = dsqrt(1.0/2.0*((hstar+h_l)*hstar/h_l**2))
        else
            qqleft = 1.0
        endif
        ! speed of non-shear wave left
        s_l = u_l-c_l*qqleft

        if (hstar>h_r) then
            qqright = dsqrt(1.0/2.0*((hstar+h_r)*hstar/h_r**2))
        else
            qqright = 1.0
        endif
        ! speed of non-shear wave right
        s_r = u_r+c_r*qqright

        ! speed of approximated shear wave sstar
        sstar = (s_l*h_r*(u_r-s_r)-s_r*h_l*(u_l-s_l))/(h_r*(u_r-s_r)-h_l*(u_l-s_l))
        ! Another choice of sstar, not advised
        ! sstar = 1.0/2.0*(u_l+u_r)- ((h_r-h_l)*(c_l+c_r))/(h_l+h_r)

        ! Middle state
        h_ml = h_l*((s_l-u_l)/(s_l-sstar))*1.0
        hu_ml = h_l*((s_l-u_l)/(s_l-sstar))*(sstar)
        hv_ml = h_l*((s_l-u_l)/(s_l-sstar))*(v_l)
        h_mr = h_r*((s_r-u_r)/(s_r-sstar))*1.0
        hu_mr = h_r*((s_r-u_r)/(s_r-sstar))*(sstar)   
        hv_mr = h_r*((s_r-u_r)/(s_r-sstar))*(v_r)
        
        wave(1,1,i) = h_ml - h_l
        wave(2,1,i) = hu_ml - qr(2,i-1)
        wave(3,1,i) = hv_ml - qr(3,i-1)
        s(1,i) = s_l

        wave(1,2,i) = h_mr - h_ml
        wave(2,2,i) = hu_mr - hu_ml
        wave(3,2,i) = hv_mr - hv_ml
        s(2,i) = sstar
    
        wave(1,3,i) = h_r - h_mr
        wave(2,3,i) = ql(2,i) - hu_mr
        wave(3,3,i) = ql(3,i) - hv_mr
        s(3,i) = s_r

    end do


    do m=1, meqn
        do i=2-mbc, mx+mbc
            amdq(m,i) = 0.d0
            apdq(m,i) = 0.d0
            do mw=1, mwaves
                if (s(mw,i) .lt. 0.d0) then
                    amdq(m,i) = amdq(m,i) + s(mw,i)*wave(m,mw,i)
                 else
                   apdq(m,i) = apdq(m,i) + s(mw,i)*wave(m,mw,i)
                 endif
            end do
        end do
    end do

end subroutine rpn2
