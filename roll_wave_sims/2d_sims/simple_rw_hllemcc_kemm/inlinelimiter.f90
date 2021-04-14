!
! Inlined version of limiter code
!
!     =====================================================
!     # Apply a limiter to the waves.
!     # The limiter is computed by comparing the 2-norm of each wave with
!     # the projection of the wave from the interface to the left or
!     # right onto the current wave.  For a linear system this would
!     # correspond to comparing the norms of the two waves.  For a
!     # nonlinear problem the eigenvectors are not colinear and so the
!     # projection is needed to provide more limiting in the case where the
!     # neighboring wave has large norm but points in a different direction
!     # in phase space.
!
!     # The specific limiter used in each family is determined by the
!     # value of the corresponding element of the array mthlim, as used in
!     # the function philim.
!     # Note that a different limiter may be used in each wave family.
!
!     # dotl and dotr denote the inner product of wave with the wave to
!     # the left or right.  The norm of the projections onto the wave are then
!     # given by dotl/wnorm2 and dotr/wnorm2, where wnorm2 is the 2-norm
!     # of wave.

subroutine limiter(maxm,num_eqn,num_waves,num_ghost,mx,wave,s,mthlim,dtdx)

    ! implicit none
    implicit double precision (a-h,o-z)

    ! Arguments
    integer, intent(in) :: maxm, num_eqn, num_waves, num_ghost, mx
    real(kind=8), intent(in out) :: wave(num_eqn, num_waves, 1-num_ghost:maxm+num_ghost)
    real(kind=8), intent(in) :: s(num_waves, 1-num_ghost:maxm+num_ghost)
    integer, intent(in) :: mthlim(num_waves)

    ! Local storage
    integer :: m, mw, i
    real(kind=8) :: r, c, wlimiter, wnorm2, dotl
    real(kind=8), dimension(num_waves) :: dotr
    real(kind=8) :: caut, cflloc, s1, s2, phimax, cfmod1, cfmod2, left, middle, rmin, ultra
    dimension dtdx(1-num_ghost:maxm+num_ghost)

    dotr = 0.d0

    x_loop: do i = 0, mx+1

        wave_loop: do mw=1,num_waves
            if (mthlim(mw) == 0) then
                cycle wave_loop
            endif

            ! Construct dot products
            wnorm2 = 0.d0
            dotl = dotr(mw)
            dotr(mw) = 0.d0
            do m=1,num_eqn
                wnorm2 = wnorm2 + wave(m,mw,i)**2
                dotr(mw) = dotr(mw) + wave(m,mw,i)*wave(m,mw,i+1)
            end do

            ! Skip this loop if it's on the boundary or the size of the wave is
            ! zero (but still want dot products to be initialized above)
            if (i == 0) cycle wave_loop
            if (wnorm2 == 0.d0) cycle wave_loop
        
            ! Compute ratio of this wave's strength to upwind wave's strength
            if (s(mw,i) > 0.d0) then
                r = dotl / wnorm2
            else
                r = dotr(mw) / wnorm2
            endif

            ! Compute value of limiter function
            select case(mthlim(mw))
                
                ! Minmod
                case(1)
                    wlimiter = max(0.d0, min(1.d0, r))

                ! Superbee
                case(2)
                    wlimiter = max(0.d0, min(1.d0, 2.d0*r), min(2.d0, r))

                ! Van Leer
                case(3)
                    wlimiter = (r + abs(r)) / (1.d0 + abs(r))

                ! Monotonized - Centered
                case(4)
                    c = (1.d0 + r)/2.d0
                    wlimiter = max(0.d0, min(c, 2.d0, 2.d0*r))

                ! Beam Warming
                case(5)
                    wlimiter = r

                ! Arora-Row
                case(6)
                    caut = 0.99d0
!     
!                   # For limiters depending on the local CFL-number:
                    cflloc = dabs(s(i,mw) * dtdx(i))
!     
                    s1     = caut * 2.d0/cflloc
                    s2     = (1.d0 + cflloc)/3.d0
                    phimax = caut * 2.d0/(1.d0 - cflloc)
                    wlimiter = dmax1(0.d0, dmin1(s1*r, 1.d0 + s2*(r-1.d0), phimax)) 

                ! Arora-Roe theta=3/4
                case(7)
!     ------------------------------
!     # Theta Limiter, theta=0.75
!     ------------------------------
!     
!     # For limiters depending on the local CFL-number:
                    cflloc = dabs(s(i,mw) * dtdx(i))
                    !     Avoiding division by zero
                                cfmod1 = dmax1(0.001,cflloc)
                                cfmod2 = dmin1(0.999,cflloc)
                                s1     = 2.d0/cfmod1
                                s2     = (1.d0 + cflloc)/3.d0
                                phimax = 2.d0/(1.d0 - cfmod2)
                    
                    
                                left = dmax1(-0.25d0*s1, 1.0d0 + s2*(r-1.0d0))
                                middle = dmax1(-0.25d0*phimax*r,0.75d0*s1*r)
                                
                                wlimiter = dmin1(left,middle, 0.75d0*phimax)

                ! Superbee beta=2/3
                case(8)
!     ------------------------------
!     # beta=2/3 limiter
!     ------------------------------
!     
!     # For limiters depending on the local CFL-number:
            cflloc = dabs(s(i,mw) * dtdx(i))
!     Avoiding division by zero
            cfmod1 = dmax1(0.001,cflloc)
            cfmod2 = dmin1(0.999,cflloc)
            s1     = 2.d0/cfmod1
            s2     = (1.d0 + cflloc)/3.d0
            phimax = 2.d0/(1.d0 - cfmod2)
!     
            rmin = r-1.0d0
            ultra = dmax1(0.0d0,dmin1(s1*r,phimax))
            wlimiter = dmin1(ultra,dmax1(1.0d0 + (s2-0.333333333d0)*rmin, 1.0d0 + (s2+0.333333333d0)*rmin))
            wlimiter = dmax1(0.0d0,wlimiter)

                case default
                    stop "Invalid limiter method."

            end select

            ! Apply resulting limit
            wave(:,mw,i) = wlimiter * wave(:,mw,i)
        
        end do wave_loop
    end do x_loop

end subroutine limiter
