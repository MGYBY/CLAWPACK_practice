
subroutine b4step2(mbc,mx,my,meqn,q,xlower,ylower,dx,dy,t,dt,maux,aux)

    ! Called before each call to step2.
    ! Use to set time-dependent aux arrays or perform other tasks.
    !
    ! This default version does nothing. 
 
    implicit none
    integer, intent(in) :: mbc,mx,my,meqn,maux
    real(kind=8), intent(in) :: xlower,ylower,dx,dy,t,dt
    real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
    real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
    
    real(kind=8) :: xcell, ycell, x1, x2, y1, y2
    integer :: i , j, i1Solid=0, i2Solid=0, j1Solid=0, j2Solid=0, jCenter=0
    logical :: foundFile
    character*15 fname

    x1 = 10.0d0
    y1 = -0.200d0
    x2 = 10.40d0
    y2 = 0.20d0

    fname = 'centerline'

    do j=1-mbc, my+mbc
        do i=1-mbc, mx+mbc
            xcell = xlower + (i-0.5d0)*dx
            ycell = ylower + (j-0.5d0)*dy

            ! in the solid
            if ((xcell .ge. x1) .and. (xcell .le. x2) .and. (ycell .ge. y1) .and. (ycell .le. y2)) then
                if (i .ge. i2Solid) then
                    i2Solid = i
                end if
                if ((i .le. i1Solid) .or. (i1Solid .le. 0)) then
                    i1Solid = i
                end if
                if ((j .le. j1Solid) .or. (j1Solid .le. 0)) then
                    j1Solid = j
                end if
                if (j .ge. j2Solid) then
                    j2Solid = j
                end if
            end if
        enddo
    enddo
    jCenter = (j1Solid+j2Solid)/2

    ! output centerline gauges
    ! outside solid
    inquire(file=fname, exist=foundFile)
    if (.not. foundFile) then
        open(61, file = fname,status='new')
        write(61,610) t, q(1,i1Solid-1,jCenter), q(1,i2Solid+1,jCenter), i1Solid, i2Solid, jCenter
    610 format(e26.16,e26.16,e26.16,I6,I6,I6)
    else
        open(61, file = fname,status='old', form='formatted')
        write(61,620) t, q(1,i1Solid-1,jCenter), q(1,i2Solid+1,jCenter), i1Solid, i2Solid, jCenter
    620 format(e26.16,e26.16,e26.16,I6,I6,I6)
    end if

end subroutine b4step2
