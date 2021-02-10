subroutine src2(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux,t,dt)

    ! Called to update q by solving source term equation 
    ! $q_t = \psi(q)$ over time dt starting at time t.
    !
    ! This default version does nothing. 
 
    implicit none
    integer, intent(in) :: mbc,mx,my,meqn,maux
    real(kind=8), intent(in) :: xlower,ylower,dx,dy,t,dt
    real(kind=8), intent(in) ::  aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
    real(kind=8), intent(inout) ::  q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
    integer :: i, j
    real(kind=8) :: frictionTermX, frictionTermY, slopeTerm
    ! problem-related params
    real(kind=8) :: distPeriod, hn, un, piVal, gravity, So, cf

    distPeriod = 0.933d0
    hn = 0.00798d0
    un = 1.03774d0
    piVal = 3.1415926d0
    gravity = 9.81d0
    So = 0.05011d0
    cf = gravity * So * 2 * hn / (un**2)    

    do j=1-mbc,my+mbc
        do i=1-mbc,mx+mbc
            frictionTermX = -dt*cf/(2.0d0)*(q(2, i, j))*sqrt(q(2, i, j)**2+q(3, i, j)**2)/((q(1, i, j))**2)
            frictionTermY = -dt*cf/(2.0d0)*(q(3, i, j))*sqrt(q(2, i, j)**2+q(3, i, j)**2)/((q(1, i, j))**2)
            slopeTerm = -dt*gravity*(-So)*q(1, i, j)

            q(2, i, j) = q(2, i, j) + frictionTermX + slopeTerm
            q(3, i, j) = q(3, i, j) + frictionTermY
        enddo
    enddo


end subroutine src2
