subroutine b4step1(mbc,mx,meqn,q,xlower,dx,t,dt,maux,aux)

    ! Called before each call to step1.
    ! Use to set time-dependent aux arrays or perform other tasks.
    !
    ! This default version does nothing. 
 
    implicit none
    integer, intent(in) :: mbc,mx,meqn,maux
    real(kind=8), intent(in) :: xlower,dx,t,dt
    real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc)
    real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc)
    
    real(kind=8) :: outputTime
    character( len = 2 ) :: cTemp
    integer :: i , j, tint
    
    outputTime = 1.0d0
    
    if (mod(t, outputTime) .lt. (dt)) then
    tint = int(t)
    write( cTemp,'(i2)' ) tint
    open (1, file = 'out-' // trim(adjustl( cTemp )) // '.txt',status='unknown',form='formatted')
    do i = 1, mx
    write(1,'(*(e15.7))') xlower+(i-1.0/2.0), q(1,i)
    enddo
!     610 format(e15.7,e15.7,e15.7)
    close( 1 )
    end if
    
    

end subroutine b4step1

