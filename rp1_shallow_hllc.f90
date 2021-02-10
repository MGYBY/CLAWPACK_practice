! =========================================================
subroutine rp1(maxmx,meqn,mwaves,maux,mbc,mx,ql,qr,auxl,auxr,wave,s,amdq,apdq)
! =========================================================

! solve Riemann problem for the 1D shallow water equations using the HLLC
! approximate Riemann solver.

! waves: 3
! equations: 2
! reconstruction of the middle state which has two states

! Conserved quantities:
!       1 depth (h)
!       2 momentum (hu)

    implicit none

    integer, intent(in) :: maxmx, meqn, mwaves, mbc, mx, maux
    double precision, dimension(meqn,1-mbc:maxmx+mbc), intent(in) :: ql, qr
    double precision, dimension(maux,1-mbc:maxmx+mbc), intent(in) :: auxl, auxr
    double precision, dimension(meqn, mwaves, 1-mbc:maxmx+mbc), intent(out) :: wave
    double precision, dimension(meqn, 1-mbc:maxmx+mbc), intent(out) :: amdq, apdq
    double precision, dimension(mwaves, 1-mbc:maxmx+mbc), intent(out) :: s

    double precision :: u_l, u_r, h_l, h_r, c_l, c_r
    double precision :: hsqrt_l, hsqrt_r, u_hat, h_hat, c_hat, grav
    double precision :: h_m, hu_m, h_ml, hu_ml, h_mr, hu_mr
    ! double precision :: rho_m, rhou_m, E_m, s1, s2
    double precision :: qqleft, qqright, hstar, sstar, s_l, s_r

    integer :: m, i, mw

    common /cparam/  grav

    do i=2-mbc,mx+mbc
        h_l = qr(1,i-1)
        h_r = ql(1,i)
        u_l = qr(2,i-1) / qr(1,i-1)
        u_r = ql(2,i  ) / ql(1,i  )
        c_l = dsqrt(grav*h_l)
        c_r = dsqrt(grav*h_r)

        hstar = 1.0/2.0*(h_l+h_r)-(1.0/4.0*(u_r-u_l)*(h_l+h_r))/(c_l+c_r)
        if (hstar>h_l) then
            qqleft = dsqrt(1.0/2.0*((hstar+h_l)*hstar/h_l**2))
        else
            qqleft = 1.0
        endif
        s_l = u_l-c_l*qqleft

        if (hstar>h_r) then
            qqright = dsqrt(1.0/2.0*((hstar+h_r)*hstar/h_r**2))
        else
            qqright = 1.0
        endif
        s_r = u_r+c_r*qqright

        sstar = (s_l*h_r*(u_r-s_r)-s_r*h_l*(u_l-s_l))/(h_r*(u_r-s_r)-h_l*(u_l-s_l))
        ! Another choice of sstar, not advised
        ! sstar = 1.0/2.0*(u_l+u_r)- ((h_r-h_l)*(c_l+c_r))/(h_l+h_r)
        
        ! Middle state
        h_ml = h_l*((s_l-u_l)/(s_l-sstar))*1.0
        hu_ml = h_l*((s_l-u_l)/(s_l-sstar))*(sstar)
        h_mr = h_r*((s_r-u_r)/(s_r-sstar))*1.0
        hu_mr = h_r*((s_r-u_r)/(s_r-sstar))*(sstar)

        wave(1,1,i) = h_ml - h_l
        wave(2,1,i) = hu_ml - qr(2,i-1)
        s(1,i) = s_l

        wave(1,2,i) = h_mr - h_ml
        wave(2,2,i) = hu_mr - hu_ml
        s(2,i) = sstar
    
        wave(1,3,i) = h_r - h_mr
        wave(2,3,i) = ql(2,i) - hu_mr
        s(3,i) = s_r

    end do 

    do m=1, meqn
        do i=2-mbc, mx+mbc
            amdq(m,i) = 0.d0
            apdq(m,i) = 0.d0
            do mw=1, mwaves
                if (s(mw,i) < 0.d0) then
                    amdq(m,i) = amdq(m,i) + s(mw,i)*wave(m,mw,i)
                else
                    apdq(m,i) = apdq(m,i) + s(mw,i)*wave(m,mw,i)
                endif
            end do
        end do
    end do

end subroutine rp1
