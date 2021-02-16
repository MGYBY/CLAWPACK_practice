subroutine src2(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux,t,dt)
    ! added feature: RK3TVD time integration

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
    real(kind=8) :: calcFrictionTermX, calcFrictionTermY, calcSlope
    real(kind=8) ::    qmed2
    real(kind=8) ::    qmed3

    distPeriod = 0.933d0
    hn = 0.00798d0
    un = 1.03774d0
    piVal = 3.1415926d0
    gravity = 9.81d0
    So = 0.05011d0
    cf = gravity * So * 2 * hn / (un**2)    

    do j=1-mbc,my+mbc
        do i=1-mbc,mx+mbc
            frictionTermX = calcFrictionTermX(dt, cf, q(1, i, j), q(2, i, j), q(3, i, j))
            frictionTermY = calcFrictionTermY(dt, cf, q(1, i, j), q(2, i, j), q(3, i, j))
            slopeTerm = calcSlope(dt, gravity, So, q(1, i, j))
            qmed2 = q(2, i, j) + frictionTermX + slopeTerm
            qmed3 = q(3, i, j) + frictionTermY

            frictionTermX = calcFrictionTermX(dt/(4.0d0), cf, q(1, i, j), qmed2, qmed3)
            frictionTermY = calcFrictionTermY(dt/(4.0d0), cf, q(1, i, j), qmed2, qmed3)
            slopeTerm = calcSlope(dt/(4.0d0), gravity, So, q(1, i, j))
            qmed2 = (3.0d0/4.0d0)*q(2, i, j) + (1.0d0/4.0d0)*qmed2 + frictionTermX + slopeTerm
            qmed3 = (3.0d0/4.0d0)*q(3, i, j) + (1.0d0/4.0d0)*qmed3 + frictionTermY

            frictionTermX = calcFrictionTermX(dt*(2.0d0/3.0d0), cf, q(1, i, j), qmed2, qmed3)
            frictionTermY = calcFrictionTermY(dt*(2.0d0/3.0d0), cf, q(1, i, j), qmed2, qmed3)
            slopeTerm = calcSlope(dt*(2.0d0/3.0d0), gravity, So, q(1, i, j))
            q(2, i, j) = (1.0d0/3.0d0)*q(2, i, j) + (2.0d0/3.0d0)*qmed2 + frictionTermX + slopeTerm
            q(3, i, j) = (1.0d0/3.0d0)*q(3, i, j) + (2.0d0/3.0d0)*qmed3 + frictionTermY
        enddo
    enddo


end subroutine src2

! -----------------------------------------------------------
! ------    friction term and slope term functions
! -----------------------------------------------------------
real(kind=8) function calcFrictionTermX(dt, cf, q1, q2, q3)
    real(kind=8), intent(in) :: dt, cf, q1, q2, q3
    calcFrictionTermX = (-1.0d0)*dt*cf/(2.0d0)*(q2)*sqrt(q2**2+q3**2)/((q1)**2)
end

real(kind=8) function calcFrictionTermY(dt, cf, q1, q2, q3)
    real(kind=8), intent(in) :: dt, cf, q1, q2, q3
    calcFrictionTermY = (-1.0d0)*dt*cf/(2.0d0)*(q3)*sqrt(q2**2+q3**2)/((q1)**2)
end

real(kind=8) function calcSlope(dt, g, S, q1)
    real(kind=8), intent(in) :: dt, g, S, q1
    calcSlope = (-1.0d0)*dt*g*((-1.0d0)*S)*q1
end