subroutine src1(meqn,mbc,mx,xlower,dx,q,maux,aux,t,dt)

    ! Called to update q by solving source term equation 
    ! $q_t = \psi(q)$ over time dt starting at time t.
    !
    ! This default version does nothing. 
 
    implicit none
    integer, intent(in) :: mbc,mx,meqn,maux
    real(kind=8), intent(in) :: xlower,dx,t,dt
    real(kind=8), intent(in) ::  aux(maux,1-mbc:mx+mbc)
    real(kind=8), intent(inout) ::  q(meqn,1-mbc:mx+mbc)

    real(kind=8) :: frictionTerm, slopeTerm

    ! problem-related params
    real(kind=8) :: distPeriod, hn, un, piVal, grav, So, cf
    integer :: i, j

    distPeriod = 0.933d0
    hn = 0.00798d0
    un = 1.03774d0
    piVal = 3.1415926d0
    grav = 9.81d0
    So = 0.05011d0
    cf = grav * So * 2 * hn / (un**2)

    do i = 1, mx+mbc
        frictionTerm = -dt*cf/(2.0d0)*(q(2, i))**2/((q(1, i))**2)
        slopeTerm = -dt*grav*(-So)*q(1,i)

        q(2, i) = q(2, i) + frictionTerm + slopeTerm
    end do

end subroutine src1
