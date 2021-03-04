subroutine setaux(mbc,mx,my,xlower,ylower,dx,dy,maux,aux)

    ! Called at start of computation before calling qinit, and
    ! when AMR is used, also called every time a new grid patch is created.
    ! Use to set auxiliary arrays aux(1:maux, 1-mbc:mx+mbc, 1-mbc:my+mbc).
    ! Note that ghost cell values may need to be set if the aux arrays
    ! are used by the Riemann solver(s).
    !
    ! This default version does nothing. 
 
    implicit none
    integer, intent(in) :: mbc,mx,my,maux
    real(kind=8), intent(in) :: xlower,ylower,dx,dy
    real(kind=8), intent(out) ::  aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
    integer :: i, j, m
    real(kind=8) :: sigma, x0, bPrime, gaussianHump

    real(kind=8) :: hn, xcell, ycell, So

    So = 0.05011
    hn = 0.00798
    sigma = 1.0d0
    x0 = 20.0d0
    bPrime = 0.05d0*hn
    
    do j=1-mbc, my+mbc
        do i=1-mbc, mx+mbc
            xcell = xlower + (i-0.5d0)*dx
            ycell = ylower + (j-0.5d0)*dy
                do m=1,maux
                    ! Gaussian hump on the channel
                    aux(m,i,j) = -So*xcell + &
                    gaussianHump(sigma, x0, bPrime, xcell)
                end do
        enddo
    enddo

end subroutine setaux

! -----------------------------------------------------------
! ------    friction term and slope term functions
! -----------------------------------------------------------
real(kind=8) function gaussianHump(sigma, x0, bPrime, xCoord)
    real(kind=8), intent(in) :: sigma, x0, bPrime, xCoord
    gaussianHump = (1.0d0)*bprime*exp(-((xCoord-x0)**2.0d0)/(sigma**2.0d0))
end
