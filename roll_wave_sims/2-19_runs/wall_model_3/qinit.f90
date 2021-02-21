! qinit routine for parabolic bowl problem, only single layer
subroutine qinit(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux)

    use geoclaw_module, only: grav

    implicit none

    ! Subroutine arguments
    integer, intent(in) :: meqn,mbc,mx,my,maux
    real(kind=8), intent(in) :: xlower,ylower,dx,dy
    real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
    real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

    ! ! Parameters for problem
    integer :: i,j
    real(kind=8) :: x,y,slope,hn,un

    slope = 0.05011d0
    hn = 0.00798d0
    un = 1.03774d0
    
    do i=1-mbc,mx+mbc
        x = xlower + (i - 0.5d0)*dx
        do j=1-mbc,my+mbc
            y = ylower + (j - 0.5d0) * dy
            ! eta = sigma * h0 / a**2 * (2.d0 * x - sigma)
            q(1,i,j) = hn
            q(2,i,j) = un*hn
            q(3,i,j) = 0.0d0
        enddo
    enddo
    
end subroutine qinit
