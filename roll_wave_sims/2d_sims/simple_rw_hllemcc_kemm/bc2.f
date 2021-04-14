c
c
c     =====================================================
      subroutine bc2(meqn,mbc,mx,my,xlower,ylower,
     &               dx,dy,q,maux,aux,t,dt,mthbc)
c     =====================================================
c
c     # Standard boundary condition choices for claw2
c
c     # At each boundary  k = 1 (left),  2 (right),  3 (top), 4 (bottom):
c     #   mthbc(k) =  0  for user-supplied BC's (must be inserted!)
c     #            =  1  for zero-order extrapolation
c     #            =  2  for periodic boundary coniditions
c     #            =  3  for solid walls, assuming this can be implemented
c     #                  by reflecting the data about the boundary and then
c     #                  negating the 2'nd (for k=1,2) or 3'rd (for k=3,4)
c     #                  component of q.
c     ------------------------------------------------
c
c     # Extend the data from the interior cells (1:mx, 1:my)
c     # to the ghost cells outside the region:
c     #   (i, 1-jbc)   for jbc = 1,mbc,  i = 1-mbc, mx+mbc
c     #   (i, my+jbc)  for jbc = 1,mbc,  i = 1-mbc, mx+mbc
c     #   (1-ibc, j)   for ibc = 1,mbc,  j = 1-mbc, my+mbc
c     #   (mx+ibc, j)  for ibc = 1,mbc,  j = 1-mbc, my+mbc
c
      implicit double precision (a-h,o-z)
      dimension q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      dimension aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
      dimension mthbc(4)
      
      double precision hn,un,amp,ghostPhase,piVal,gravity
      double precision period, gp, xCenter
      hn = 0.00798d0
      un = 1.03774d0
      amp = 0.05d0
      piVal = 3.1415926d0
      gravity = 9.81d0
      period = 0.933d0      
      gp = sqrt(gravity*hn)

c
c
c-------------------------------------------------------
c     # left boundary:
c-------------------------------------------------------
      go to (100,110,120,130) mthbc(1)+1
c
  100 continue
c     # user-specified boundary conditions go here in place of error output
      do j = 1-mbc, my+mbc
         do ibc=1,mbc
            xCenter = (0.50+(ibc-1))*dx
            ghostPhase=(2.0d0*piVal*(t+(xCenter/gp)))/period
            q(1,1-ibc,j) = hn*(1.d0+amp*sin(ghostPhase))
            q(2,1-ibc,j) = hn*un
            q(3,1-ibc,j) = 0.0d0
         end do
      end do
      go to 199
c
  110 continue
c     # zero-order extrapolation:
      do j = 1-mbc, my+mbc
         do ibc=1,mbc
            do m=1,meqn
               q(m,1-ibc,j) = q(m,1,j)
            end do
         end do 
      end do
      go to 199

  120 continue
c     # periodic:  
      do j = 1-mbc, my+mbc
         do ibc=1,mbc
            do m=1,meqn
               q(m,1-ibc,j) = q(m,mx+1-ibc,j)
            end do
         end do
      end do
      go to 199

  130 continue
c     # solid wall (assumes 2'nd component is velocity or momentum in x):
      do j = 1-mbc, my+mbc
         do ibc=1,mbc
            do m=1,meqn
               q(m,1-ibc,j) = q(m,ibc,j)
            end do
         end do
      end do
c     # negate the normal velocity:
      do j = 1-mbc, my+mbc
         do ibc=1,mbc
            q(2,1-ibc,j) = -q(2,ibc,j)
         end do
      end do
      go to 199

  199 continue
c
c-------------------------------------------------------
c     # right boundary:
c-------------------------------------------------------
      go to (200,210,220,230) mthbc(2)+1
c
  200 continue
c     # user-specified boundary conditions go here in place of error output
      write(6,*) '*** ERROR *** mthbc(2)=0 and no BCs specified in bc2'
      stop
      go to 299

  210 continue
c     # zero-order extrapolation:
      do j = 1-mbc, my+mbc
         do ibc=1,mbc
            do m=1,meqn
               q(m,mx+ibc,j) = q(m,mx,j)
            end do
         end do
      end do
      go to 299

  220 continue
c     # periodic:  
      do j = 1-mbc, my+mbc
         do ibc=1,mbc
            do m=1,meqn
               q(m,mx+ibc,j) = q(m,ibc,j)
            end do
         end do
      end do
      go to 299

  230 continue
c     # solid wall (assumes 2'nd component is velocity or momentum in x):
      do j = 1-mbc, my+mbc
         do ibc=1,mbc
            do m=1,meqn
               q(m,mx+ibc,j) = q(m,mx+1-ibc,j)
            end do
         end do
      end do
c     # negate the normal velocity:
      do j = 1-mbc, my+mbc
         do ibc=1,mbc
            q(2,mx+ibc,j) = -q(2,mx+1-ibc,j)
         end do 
      end do
      go to 299

  299 continue
c
c-------------------------------------------------------
c     # bottom boundary:
c-------------------------------------------------------
      go to (300,310,320,330) mthbc(3)+1
c
  300 continue
c     # user-specified boundary conditions go here in place of error output
      write(6,*) '*** ERROR *** mthbc(3)=0 and no BCs specified in bc2'
      stop
      go to 399
c
  310 continue
c     # zero-order extrapolation:
      do jbc=1,mbc
         do i = 1-mbc, mx+mbc
            do m=1,meqn
               q(m,i,1-jbc) = q(m,i,1)
            end do
         end do
      end do
      go to 399

  320 continue
c     # periodic:  
      do jbc=1,mbc
         do i = 1-mbc, mx+mbc
            do m=1,meqn
               q(m,i,1-jbc) = q(m,i,my+1-jbc)
            end do
         end do
      end do
      go to 399

  330 continue
c     # solid wall (assumes 3'rd component is velocity or momentum in y):
      do jbc=1,mbc
         do i = 1-mbc, mx+mbc
            do m=1,meqn
               q(m,i,1-jbc) = q(m,i,jbc)
            end do
         end do
      end do
c     # negate the normal velocity:
      do jbc=1,mbc
         do i = 1-mbc, mx+mbc
            q(3,i,1-jbc) = -q(3,i,jbc)
         end do
      end do
      go to 399

  399 continue
c
c-------------------------------------------------------
c     # top boundary:
c-------------------------------------------------------
      go to (400,410,420,430) mthbc(4)+1
c
  400 continue
c     # user-specified boundary conditions go here in place of error output
      write(6,*) '*** ERROR *** mthbc(4)=0 and no BCs specified in bc2'
      stop
      go to 499

  410 continue
c     # zero-order extrapolation:
      do jbc=1,mbc
         do i = 1-mbc, mx+mbc
            do  m=1,meqn
               q(m,i,my+jbc) = q(m,i,my)
            end do
         end do
      end do
      go to 499

  420 continue
c     # periodic:  
      do jbc=1,mbc
         do i = 1-mbc, mx+mbc
            do m=1,meqn
               q(m,i,my+jbc) = q(m,i,jbc)
            end do
         end do
      end do
      go to 499

  430 continue
c     # solid wall (assumes 3'rd component is velocity or momentum in y):
      do jbc=1,mbc
         do i = 1-mbc, mx+mbc
            do m=1,meqn
               q(m,i,my+jbc) = q(m,i,my+1-jbc)
            end do
         end do
      end do
c     # negate the normal velocity:
      do jbc=1,mbc
         do i = 1-mbc, mx+mbc
            q(3,i,my+jbc) = -q(3,i,my+1-jbc)
         end do
      end do
      go to 499

  499 continue

      return
      end

