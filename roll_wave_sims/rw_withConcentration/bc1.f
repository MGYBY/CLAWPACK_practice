
c
c
c     =================================================================
      subroutine bc1(meqn,mbc,mx,xlower,dx,q,maux,aux,t,dt,mthbc)
c     =================================================================
c
c     # Standard boundary condition choices for claw2
c
c     # At each boundary  k = 1 (left),  2 (right):
c     #   mthbc(k) =  0  for user-supplied BC's (must be inserted!)
c     #            =  1  for zero-order extrapolation
c     #            =  2  for periodic boundary coniditions
c     #            =  3  for solid walls, assuming this can be implemented
c     #                  by reflecting the data about the boundary and then
c     #                  negating the 2'nd component of q.
c     ------------------------------------------------
c
c     # Extend the data from the computational region
c     #      i = 1, 2, ..., mx2
c     # to the virtual cells outside the region, with
c     #      i = 1-ibc  and   i = mx+ibc   for ibc=1,...,mbc
c
      implicit double precision (a-h,o-z)
      dimension q(meqn,1-mbc:mx+mbc)
      dimension aux(maux,1-mbc:mx+mbc)

      dimension mthbc(2)
      double precision hn, un, amp, ghostPhase, piVal, gravity
      double precision period, gp, xCenter, conc
      hn = 0.00798d0
      un = 1.03774d0
      amp = 0.05d0
      piVal = 3.1415926d0
      gravity = 9.81d0
      period = 0.933d0
      conc = 1.0d0

c
c
c-------------------------------------------------------
c     # left boundary:
c-------------------------------------------------------
      go to (100,110,120,130) mthbc(1)+1
c
  100 continue
c     # user-specified boundary conditions go here in place of error output
c     # inflow wave maker BC 
c     # extrap tracer
      do ibc=1,mbc
         xCenter = (0.50+(ibc-1))*dx
         gp = sqrt(gravity*hn)
         ghostPhase = (2.0d0 * piVal*(t+(xCenter/gp)))/period
         q(1, 1-ibc) = hn*(1.d0 + amp*sin(ghostPhase))
         q(2, 1-ibc) = hn*un
         q(3, 1-ibc) = hn*conc
      end do
      go to 199
c
  110 continue
c     # zero-order extrapolation:
      do ibc=1,mbc
         do m=1,meqn
            q(m,1-ibc) = q(m,1)
         end do
      end do
      go to 199

  120 continue
c     # periodic:  
      do ibc=1,mbc
         do m=1,meqn
            q(m,1-ibc) = q(m,mx+1-ibc)
         end do
      end do
      go to 199

  130 continue
c     # solid wall (assumes 2'nd component is velocity or momentum in x):
      do ibc=1,mbc
         do m=1,meqn
            q(m,1-ibc) = q(m,ibc)
         end do
         q(2,1-ibc) = -q(2,ibc)
      end do
c     # negate the normal velocity:
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
      do ibc=1,mbc
         do m=1,meqn
            q(m,mx+ibc) = q(m,mx)
         end do
      end do
      go to 299

  220 continue
c     # periodic:  
      do ibc=1,mbc
         do m=1,meqn
            q(m,mx+ibc) = q(m,ibc)
         end do
      end do
      go to 299

  230 continue
c     # solid wall (assumes 2'nd component is velocity or momentum in x):
      do ibc=1,mbc
         do m=1,meqn
            q(m,mx+ibc) = q(m,mx+1-ibc)
         end do
         q(2,mx+ibc) = -q(2,mx+1-ibc)
      end do
      go to 299

  299 continue
c
      return
      end

