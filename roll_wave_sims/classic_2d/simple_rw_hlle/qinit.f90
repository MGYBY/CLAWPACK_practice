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

    ! Parameters for problem
    real(kind=8) :: slope,hn,un

    slope = 0.05011d0
    hn = 0.00798d0
    un = 1.03774d0

    do j=1,my
        do i=1,mx
            ! xcell = xlower + (i-0.5d0)*dx
            q(1,i,j) = hn
            q(2,i,j) = hn*un
            q(3,i,j) = 0.0d0
        enddo
    enddo

end subroutine qinit
