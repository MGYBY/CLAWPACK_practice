subroutine qinit(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux)

    ! Set initial conditions for the q array.
    ! This default version prints an error message since it should
    ! not be used directly.  Copy this to an application directory and
    ! loop over all grid cells to set values of q(1:meqn, 1:mx, 1:my).

    implicit none
    
    integer, intent(in) :: meqn,mbc,mx,my,maux
    real(kind=8), intent(in) :: xlower,ylower,dx,dy
    real(kind=8), intent(in) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
    real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
    integer :: i,j
    real(kind=8) :: xcell, ycell

    ! Parameters for problem
    real(kind=8) :: slope,hn,un, x1, y1, x2, y2

    slope = 0.05011d0
    hn = 0.00798d0
    un = 1.03774d0
    x1 = 40.0d0
    y1 = -0.10d0
    x2 = 40.20d0
    y2 = 0.10d0

    ! do j=1,my
    !     do i=1,mx
    !         xcell = xlower + (i-0.5d0)*dx
    !         ycell = ylower + (j-0.5d0)*dy

    !         ! mask the region in the internal walls
    !         if ((xcell>=x1) .and. (xcell<=x2) .and. (ycell>=y1) .and. (ycell<=y2)) then
    !             q(1,i,j) = 0.0d0
    !             q(2,i,j) = 0.0d0
    !             q(3,i,j) = 0.0d0
    !         else
    !             q(1,i,j) = hn
    !             q(2,i,j) = hn*un
    !             q(3,i,j) = 0.0d0
    !         end if
    !     enddo
    ! enddo

    do j=1,my
        do i=1,mx
            ! xcell = xlower + (i-0.5d0)*dx
            ! ycell = ylower + (j-0.5d0)*dy
            q(1,i,j) = hn
            q(2,i,j) = hn*un
            q(3,i,j) = 0.0d0
        enddo
    enddo

end subroutine qinit
