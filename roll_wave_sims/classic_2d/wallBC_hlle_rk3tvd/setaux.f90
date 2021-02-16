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

    real(kind=8) :: x1, y1, x2, y2, xcell, ycell

    x1 = 40.0d0
    y1 = -0.200d0
    x2 = 40.40d0
    y2 = 0.20d0
    
    do j=1-mbc, my+mbc
        do i=1-mbc, mx+mbc
            xcell = xlower + (i-0.5d0)*dx
            ycell = ylower + (j-0.5d0)*dy

            ! internal wall spatial BC
            if ((xcell .ge. x1) .and. (xcell .le. x2) .and. (ycell .ge. y1) .and. (ycell .le. y2)) then
                ! mask the region inside the internal walls
                do m=1,maux
                    ! inside the block, zero
                    aux(m,i,j) = 0.0d0
                end do
            else
                ! the region outside the internal walls
                do m=1,maux
                    ! outside the block, unity
                    aux(m,i,j) = 1.0d0
                end do
            end if
        enddo
    enddo

end subroutine setaux
