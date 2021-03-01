
subroutine b4step2(mbc,mx,my,meqn,q,xlower,ylower,dx,dy,t,dt,maux,aux)

    ! Called before each call to step2.
    ! Use to set time-dependent aux arrays or perform other tasks.
    !
    ! This default version does nothing. 

    ! implementation of the zero normal BC of the block
 
    implicit none
    integer, intent(in) :: mbc,mx,my,meqn,maux
    real(kind=8), intent(in) :: xlower,ylower,dx,dy,t,dt
    real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
    real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

    ! problem-related parameters
    real(kind=8) :: x1, x2, y1, y2, bLength, bWidth, xcell, ycell, xCenter, yCenter
    real(kind=8) :: outputTime
    integer :: i , j, i1Solid=0, i2Solid=0, j1Solid=0, j2Solid=0, jCenter=0
    ! integer :: frontGauge1, frontGauge2, frontGauge3, frontGauge4, frontGauge5
    ! integer :: rearGauge1, rearGauge2, rearGauge3, rearGauge4, rearGauge5
    logical :: foundFile
    character*20 fname, fname1, fname2, fname3

    x1 = 10.00d0
    y1 = -0.200d0
    x2 = 10.40d0
    y2 = 0.20d0
    bLength = x2-x1
    bWidth = y2-y1
    xCenter = (x1+x2)/(2.0d0)
    yCenter = (y1+y2)/(2.0d0)
    outputTime = 0.10d0
    
    fname = 'centerline'
    fname1 = 'centerline_front'
    fname2 = 'centerline_rear'
    fname3 = 'centerline_time'

    do j=1,my
        do i=1,mx
            xcell = xlower + (i-0.5d0)*dx
            ycell = ylower + (j-0.5d0)*dy

            ! inside the solid
            ! set velocity normal to similarity rectangles to zero
            if ((xcell>=x1) .and. (xcell<=x2) .and. (ycell>=y1) .and. (ycell<=y2)) then
                ! right side
                if (((xcell-xCenter) .ge. 0.0d0) .and. (abs(ycell-yCenter)/abs(xcell-xCenter) .le. (bLength/bWidth))) then
                    q(2, i, j) = 0.0d0
                ! left side
                else if (((xcell-xCenter) .le. 0.0d0) .and. (abs(ycell-yCenter)/abs(xcell-xCenter) .le. (bLength/bWidth))) then
                    q(2, i, j) = 0.0d0
                ! top side
                else if (((ycell-yCenter) .ge. 0.0d0) .and. (abs(ycell-yCenter)/abs(xcell-xCenter) .gt. (bLength/bWidth))) then
                    q(3, i, j) = 0.0d0
                ! bottom side
                else if (((ycell-yCenter) .le. 0.0d0) .and. (abs(ycell-yCenter)/abs(xcell-xCenter) .gt. (bLength/bWidth))) then
                    q(3, i, j) = 0.0d0
                end if
            end if
        enddo
    enddo
    
    ! do j=1-mbc, my+mbc
    !         do i=1-mbc, mx+mbc
    !             xcell = xlower + (i-0.5d0)*dx
    !             ycell = ylower + (j-0.5d0)*dy

    !             ! in the solid
    !             if ((xcell .ge. x1) .and. (xcell .le. x2) .and. (ycell .ge. y1) .and. (ycell .le. y2)) then
    !                 if (i .ge. i2Solid) then
    !                     i2Solid = i
    !                 end if
    !                 if ((i .le. i1Solid) .or. (i1Solid .le. 0)) then
    !                     i1Solid = i
    !                 end if
    !                 if ((j .le. j1Solid) .or. (j1Solid .le. 0)) then
    !                     j1Solid = j
    !                 end if
    !                 if (j .ge. j2Solid) then
    !                     j2Solid = j
    !                 end if
    !             end if
    !         enddo
    ! enddo
    !     jCenter = int((j1Solid+j2Solid)/2)

    !     ! output centerline gauges
    !     ! outside solid
    !     inquire(file=fname, exist=foundFile)
    !     if (.not. foundFile) then
    !         open(61, file = fname,status='new')
    !         write(61,610) t, q(1,i1Solid-1,jCenter), q(1,i2Solid+1,jCenter)
    !     610 format(e15.7,e15.7,e15.7)
    !     else
    !         open(61, file = fname,status='old', form='formatted')
    !         write(61,620) t, q(1,i1Solid-1,jCenter), q(1,i2Solid+1,jCenter)
    !     620 format(e15.7,e15.7,e15.7)
    !     end if

        ! not appropriate for AMR
        ! ! if (mod(t, outputTime) .lt. (dt)) then
        !     inquire(file=fname1, exist=foundFile)
        !     if (.not. foundFile) then
        !         open(62, file = fname1,status='new')
        !         write(62,'(*(e15.7))') (q(1,i1Solid-1,j),j=j1Solid,j2Solid)
        !     else
        !         open(62, file = fname1,status='old', form='formatted')
        !         write(62,'(*(e15.7))') (q(1,i1Solid-1,j),j=j1Solid,j2Solid)
        !     end if

        !     inquire(file=fname2, exist=foundFile)
        !     if (.not. foundFile) then
        !         open(63, file = fname2,status='new')
        !         write(63,'(*(e15.7))') (q(1,i2Solid+1,j),j=j1Solid,j2Solid)
        !     else
        !         open(63, file = fname2,status='old', form='formatted')
        !         write(63,'(*(e15.7))') (q(1,i2Solid+1,j),j=j1Solid,j2Solid)
        !     end if
            
        !     inquire(file=fname3, exist=foundFile)
        !     if (.not. foundFile) then
        !         open(64, file = fname3,status='new')
        !         write(64,'(*(e15.7))') t
        !     else
        !         open(64, file = fname3,status='old', form='formatted')
        !         write(64,'(*(e15.7))') t
        !     end if

        ! ! end if


end subroutine b4step2
