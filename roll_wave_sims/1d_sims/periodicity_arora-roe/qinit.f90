subroutine qinit(meqn,mbc,mx,xlower,dx,q,maux,aux)

    ! Set initial conditions for the q array.
    ! This default version prints an error message since it should
    ! not be used directly.  Copy this to an application directory and
    ! loop over all grid cells to set values of q(1:meqn, 1:mx).

    implicit none
    
    integer, intent(in) :: meqn,mbc,mx,maux
    real(kind=8), intent(in) :: xlower,dx
    real(kind=8), intent(in) :: aux(maux,1-mbc:mx+mbc)
    real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc)

    integer :: i
    real(kind=8) :: xcell
    ! Parameters for problem
    real(kind=8) :: slope,hn,un, amp, domainX, piVal

    slope = 0.05011d0
    hn = 0.00798d0
    un = 1.03774d0
    piVal = 3.1416d0
    amp = 0.00d0
    domainX = 2.0d0
 
 
      do i=1,mx+mbc
         xcell = xlower + (i-0.5d0)*dx
         q(1,i) = hn*(1.0d0+amp*sin(2*piVal*xcell/domainX))
         q(2,i) = hn*un
      enddo

end subroutine qinit

