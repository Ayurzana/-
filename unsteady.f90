! Lattice Boltzmann sample, written in Fortran 90
    ! http://www.palabos.org/forum/read.php?3,7698,7698
    ! Copyright (C) 2006 Orestis Malaspinas
    ! Address: EPFL STI ISE LIN, ME A2 398, 1015 Lausanne, Switzerland
    ! E-mail: orestis.malaspinas@epfl.ch
    !
    ! This program is free software; you can redistribute it and/or
    ! modify it under the terms of the GNU General Public License
    ! as published by the Free Software Foundation; either version 2
    ! of the License, or (at your option) any later version.
    !
    ! This program is distributed in the hope that it will be useful,
    ! but WITHOUT ANY WARRANTY; without even the implied warranty of
    ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    ! GNU General Public License for more details.
    !
    ! You should have received a copy of the GNU General Public 
    ! License along with this program; if not, write to the Free 
    ! Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
    ! Boston, MA  02110-1301, USA.
    ! This example examines an unsteady flow past a cylinder placed in a channel.
    ! The cylinder is offset somewhat from the center of the flow to make the
    ! steady-state symmetrical flow unstable. At the inlet and outlet, a Poiseuille
    ! profile is imposed on the velocity. At Reynolds numbers around 100,
    ! an unstable periodic pattern is created, the Karman vortex street.
    ! Note that with the implemented Zou/He boundary condition, you must
    ! increase the resolution to keep the simulation stable if you increase
    ! the Reynolds number.
    !	 ========================================================
    !	 Constants that identify different cell-types according
    !        to the dynamics they implement
    !	 ========================================================
      MODULE cellConst
      integer, parameter:: fluid=0, wall=100, inlet=50, outlet=80
      END MODULE cellConst
 
    !	 ========================================================
    !	 Lattice constants for the D2Q9 lattice
    !	 ========================================================
      MODULE D2Q9Const
    !	 D2Q9 Weights
      double precision,parameter:: t(0:8) = (/4.0d0/9.0d0,1.0d0/9.0d0&
      ,1.0d0/9.0d0,1.0d0/9.0d0,1.0d0/9.0d0,1.0d0/36.0d0,1.0d0/36.0d0&
      ,1.0d0/36.0d0,1.0d0/36.0d0/)
    !	D2Q9 Directions
      integer:: v(0:8,0:1)&
           = (/(/0,1,0,-1,0,1,-1,-1,1/),(/0,0,1,0,-1,1,1,-1,-1/)/)
 
      integer, parameter:: opposite(0:8) = (/0,3,4,1,2,7,8,5,6/)
      END MODULE D2Q9Const
 
 
    !	 ========================================================
    !	 Constants for simulation setup
    !	 ========================================================
      MODULE simParam
      integer, parameter:: xDim = 150*5 
      integer, parameter:: yDim = 20*5  
      integer, parameter:: obstX = xDim/5
      integer, parameter:: obstY = yDim/2
      integer, parameter:: obstR = min(yDim/10+0.5,10.5)
 
      integer, parameter:: tMax = 20000
     
      double precision, parameter:: uMax = 0.1d0
    !  double precision, parameter:: Re = 40.0d0
      END MODULE simParam
 
 
    !	 ========================================================
    !	 The main program, implementing a flow past a cylinder
    !	 ========================================================
 
      PROGRAM unsteady
      USE simParam, ONLY: xDim, yDim, tMax
      use cellConst
      implicit none
 
      double precision:: omega, time1, time2, timeTot
      double precision, dimension(:,:,:), allocatable:: f,f_previous&
                                                     , fEq, u,u0
      double precision, dimension(:,:), allocatable:: rho,uSqr,lift,drag
      integer, dimension(:,:), allocatable:: image
      integer:: tStep
      integer:: i, j, k ,x, y
      double precision :: Err
      double precision:: Re
 
      allocate(f(yDim,xDim,0:8))
      allocate(f_previous(yDim,xDim,0:8))
      allocate(fEq(yDim,xDim,0:8))
      allocate(u(yDim,xDim,0:1))
      allocate(u0(yDim,xDim,0:1))
      allocate(uSqr(yDim,xDim))
      allocate(rho(yDim,xDim))
      allocate(lift(yDim,xDim))
      allocate(drag(yDim,xDim))
      allocate(image(yDim,xDim))
 
      CALL constructImage(image)
      CALL computeOmega(omega)
	!  Re1=Re
      CALL writeInput(omega)
      CALL initMacro(rho,u,u0,uSqr,image)
      CALL computeFeq(fEq,rho,u,uSqr)
 
      f = fEq
 
      timeTot = 0.0d0
      do tStep = 1, tMax
 
      CALL CPU_TIME(time1)
      CALL collide(f,fEq,omega,image)
      CALL stream(f)
      CALL inletOutlet(f,rho,u,image)  
      CALL boundaries(f,image,u)
      CALL computeMacros(f,rho,u,uSqr)
      CALL computeFeq(fEq,rho,u,uSqr)
 
    !f_previous(1:yDim,1:xDim,0:8)=f(1:yDim,1:xDim,0:8)
 
  !  call computeForce(f,rho,u,uSqr,f_previous,lift,drag)
 
      CALL CPU_TIME(time2)
      timeTot = timeTot + (time2-time1)
 
 
     ! if( mod(tStep,500) == 0) then
       ! call output(tStep,image,u,rho,lift,drag)
       ! write(*,*) "tStep ",tStep," done "," err ",Err(u,u0,image)
 
      !end if
      if( Err(u,u0,image) < 1.0d-6 .or. isnan(Err(u,u0,image))) then
        exit
      end if
      end do
      print*,"TCPU=",timeTot
      CALL writeImage(image)
      CALL Output(tStep,image,u,rho,lift,drag)
      write(*,*) dble(tMax)*(dble(yDim*xDim))/timeTot,'cells per second'
 
      deallocate(f)
      deallocate(fEq)
      deallocate(u)
      deallocate(uSqr)
      deallocate(rho)
      deallocate(image)
      END PROGRAM unsteady
 
! program finished. functions and subroutines are started.
! *********************************************************************************
      double precision function Err(u,u0,image)
      use cellConst
      USE simParam
      double precision, INTENT(INOUT):: u(yDim,xDim,0:1)
      double precision, INTENT(INOUT):: u0(yDim,xDim,0:1)
      integer, INTENT(INOUT)::image(yDim, xDim)
 
      integer:: i, j, k ,x, y
      double precision ::  e1,e2
      e1=0.0
      e2=0.0
      do i=1,yDim    
      do j=1,xDim
      if(image(i,j) /= wall) then
 
      e1=e1+ sqrt((u(i,j,0)-u0(i,j,0))**2+(u(i,j,1)-u0(i,j,1))**2)
 
            !e2+=sqrt( ux[j][i]*ux[j][i]+uy[j][i]*uy[j][i]);
      e2=e2+sqrt(u(i,j,0)**2+u(i,j,1)**2)
 
            !u0[j][i]=ux[j][i];v0[j][i]=uy[j][i];
 
      end if          
      end do
      end do
      u0=u
      e1=sqrt(e1)
      e1=sqrt(e1)
 
      Err =e1 /( e2+1d-30) 
      end function Err
! *********************************************************************************
      subroutine output(tStep,image,u,rho,lift,drag)
      use cellConst
      USE simParam, ONLY: xDIm, yDim
      double precision, INTENT(INOUT):: u(yDim,xDim,0:1),rho(yDim, xDim)&
                                     ,lift(yDim, xDim),drag(yDim, xDim)
      integer, INTENT(INOUT)::image(yDim, xDim)
      integer:: i,j
      integer:: tStep
 
      Character( Len = 1200 ) :: filename
 
      write(filename,"(I7,A4)") tStep,".plt"
 
 
      open(1,file="filename.dat")
      write(1,*)"VARIABLES =X, Y, U, V, RHO, IMAGE"!, Vel
      write(1,*)"ZONE ","I=",xDim,"J=",yDim,",","F=POINT"
 
      do i=1,yDim    
        do j=1,xDim
            if(image(i,j)==wall) then
                u(i,j,0)=0
                u(i,j,1)=0
                rho(i,j)=wall
            end if
 
            write(1,*) j, i,u(i,j,0), u(i,j,1), rho(i,j),image(i,j)!,sqrt(u(i,j,0)**2+u(i,j,1)**2)
        end do
      end do
 
      end subroutine output
! *********************************************************************************
    !	 ========================================================
    !	 Compute the relaxation parameter from the Reynolds number
    !	 ========================================================
      SUBROUTINE computeOmega(omega)
      USE simParam, ONLY: uMax,obstR !,Re
 
      implicit none
 
      double precision, INTENT(INOUT):: omega
      double precision:: nu
      double precision:: Re
      nu    =0.04d0 !  uMax * 2.0d0 * dble(obstR) / Re
      Re=uMax * 2.0d0 * dble(obstR) / nu
      print*,"Re1=",Re
      omega = 1.0d0 / (3.0d0*nu+0.5d0)
      END SUBROUTINE computeOmega
! *********************************************************************************
    !	 ========================================================
    !	 Construct an array the defines the flow geometry
    !	 ========================================================
      SUBROUTINE constructImage(image)
      USE cellConst
      USE simParam, ONLY: xDim, yDim, obstX, obstY, obstR
      USE D2Q9Const, ONLY: v
 
      implicit none
 
      integer, INTENT(INOUT):: image(yDim,xDim)
      integer:: x,y
 
     ! v(0:8,0) = (/0,1,0,-1,0,1,-1,-1,1/)
     ! v(0:8,1) = (/0,0,1,0,-1,1,1,-1,-1/)
 
      image          = fluid
      image(:,1)     = inlet
      image(:,xDim)  = outlet
      image(1,:)     = wall
      image(yDim,:)  = wall
      do x = 1, xDim
        do y = 1, yDim
      if (((x-obstX)**2+(y-obstY)**2)<=(obstR**2)) image(y,x)=wall
        end do
      end do
 
      END SUBROUTINE constructImage
! *********************************************************************************
    !	 ========================================================
    !	 Initialize the simulation to Poiseuille profile at
    !        an equilibrium distribution
    !	 ========================================================
      SUBROUTINE initMacro(rho,u,u0,uSqr,image)
 
      USE simParam!, ONLY: xDim, yDim
      use cellConst
 
      implicit none
      integer, INTENT(INOUT):: image(yDim,xDim)
 
      double precision, INTENT(INOUT):: rho(yDim,xDim), u(yDim,xDim,0:1)&
                                    , u0(yDim,xDim,0:1),uSqr(yDim,xDim)
      double precision:: uProf
      integer:: x,y
 
      do y = 1, yDim
        u(y,:,0) = uProf(y)
        u(y,:,1) = 0.0d0
        u0(y,:,0) = 0.0d0
        u0(y,:,1) = 0.0d0
      end do
 
 
      rho  = 1.0d0
      uSqr = u(:,:,0) * u(:,:,0) + u(:,:,1) * u(:,:,1)
      END SUBROUTINE initMacro
! *********************************************************************************
    !	 ========================================================
    !	 Compute equilibrium distribution
    !	 ========================================================
      SUBROUTINE computeFeq(fEq,rho,u,uSqr)
      USE D2Q9COnst, ONLY: t, v
      USE simParam, ONLY: xDim, yDim
      implicit none
 
      double precision, INTENT(IN):: rho(yDim,xDim), uSqr(yDim,xDim)&
                                  , u(yDim,xDim,0:1)
      double precision, INTENT(INOUT):: fEq(yDim,xDim,0:8)
      integer:: i, x, y
      double precision:: uxy
 
      do i = 0, 8
        do x = 1, xDim
            do y = 1, yDim
                uxy = u(y,x,0) * v(i,0) + u(y,x,1) * v(i,1)
                fEq(y,x,i)=t(i)*rho(y,x)*(1.0d0+3.0d0*uxy+4.5d0&
                          *uxy*uxy-1.5d0*uSqr(y,x))
            end do
        end do
      end do
      END SUBROUTINE computeFeq
! *********************************************************************************
    !	 ========================================================
    !	 Compute drag and lift from distribution functions
    !	 ========================================================
      SUBROUTINE computeForce(f,rho,u,uSqr,f_previous,lift,drag)
      USE simParam, ONLY: xDIm, yDim
      implicit none
 
      double precision, INTENT(IN):: f(yDim,xDim,0:8)&
                          ,f_previous(yDim,xDim,0:8)
      double precision, INTENT(INOUT)::u(yDim,xDim,0:1), rho(yDim, xDim)&
                  , uSqr(yDim, xDim), lift(yDim, xDim),drag(yDim, xDim)
      integer:: x,y
 
      do x = 1, xDim-1
        do y = 1, yDim-1
            !rho(y,x)  = f(y,x,0) + f(y,x,1) + f(y,x,2) + f(y,x,3) + f(y,x,4) + f(y,x,5) + f(y,x,6) + f(y,x,7) + f(y,x,8)
            !u(y,x,0)  = (f(y,x,1) - f(y,x,3) + f(y,x,5) - f(y,x,6) - f(y,x,7) + f(y,x,8)) / rho(y,x)
            !u(y,x,1)  = (f(y,x,2) - f(y,x,4) + f(y,x,5) + f(y,x,6) - f(y,x,7) - f(y,x,8)) / rho(y,x)
            !uSqr(y,x) = u(y,x,0) * u(y,x,0) + u(y,x,1) * u(y,x,1)
      drag(y,x)=(f(y,x+1,1)-f(y,x+1,3)+f(y,x+1,5)-f(y,x+1,6)-f(y,x+1,7)&
               +f(y,x+1,8))+(f_previous(y,x,1)-f_previous(y,x,3)&
               +f_previous(y,x,5)-f_previous(y,x,6)-f_previous(y,x,7)&
               +f_previous(y,x,8))
      lift(y,x)=(f(y+1,x,2)-f(y+1,x,4)+f(y+1,x,5)+f(y+1,x,6)-f(y+1,x,7)&
               -f(y+1,x,8))+(f_previous(y,x,2)-f_previous(y,x,4)&
               +f_previous(y,x,5)+f_previous(y,x,6)-f_previous(y,x,7)&
               -f_previous(y,x,8))
        end do
      end do
      END SUBROUTINE computeForce
! *********************************************************************************
    !	 ========================================================
    !	 Compute density and velocity from distribution functions
    !	 ========================================================
      SUBROUTINE computeMacros(f,rho,u,uSqr)
      USE simParam, ONLY: xDIm, yDim
      implicit none
 
      double precision, INTENT(IN):: f(yDim,xDim,0:8)
      double precision, INTENT(INOUT)::u(yDim,xDim,0:1),rho(yDim,xDim)&
                                   ,uSqr(yDim,xDim)
      integer:: x,y
 
      do x = 1, xDim
        do y = 1, yDim
            !if (image(y,x) /= wall) then
            rho(y,x)=f(y,x,0)+f(y,x,1)+f(y,x,2)+f(y,x,3)+f(y,x,4)&
                    +f(y,x,5)+f(y,x,6)+f(y,x,7)+f(y,x,8)
            u(y,x,0)=(f(y,x,1)-f(y,x,3)+f(y,x,5)-f(y,x,6)-f(y,x,7)&
                    +f(y,x,8))/rho(y,x)
            u(y,x,1)=(f(y,x,2)-f(y,x,4)+f(y,x,5)+f(y,x,6)-f(y,x,7)&
                    -f(y,x,8))/rho(y,x)
            uSqr(y,x) = u(y,x,0) * u(y,x,0) + u(y,x,1) * u(y,x,1)
            !end if
        end do
      end do
      END SUBROUTINE computeMacros
! *********************************************************************************
    !	 ========================================================
    !	 Implement Bounce-back on upper/lower boundaries
    !	 ========================================================
      SUBROUTINE boundaries(f,image,u)
      USE D2Q9Const, ONLY: opposite
      USE cellConst, ONLY: wall
      USE simParam, ONLY: xDim, yDim
      implicit none
 
      integer, INTENT(IN):: image(yDim,xDim)
      double precision, INTENT(INOUT):: f(yDim,xDim,0:8)&
                                      , u(yDim,xDim,0:1)
   ! double precision, INTENT(INOUT):: rho(yDim,xDim)
      double precision:: fTmp(0:8)
      integer:: i, x, y
 
      do x = 1, xDim
        do y = 1, yDim
            if (image(y,x) == wall) then
 
            u(y,x,0)=0.0
            u(y,x,1)=0.0
            do i = 0, 8
                fTmp(i) = f(y,x,opposite(i))
            end do
            do i = 0, 8
                f(y,x,i) = fTmp(i)
            end do
 
           ! rho(y,x)=wall
 
 
            end if
        end do
      end do
      END SUBROUTINE boundaries
! *********************************************************************************
    !	 ========================================================
    !	 Use Zou/He boundary condition to implement Dirichlet
    !        boundaries on inlet/outlet
    !	 ========================================================
      SUBROUTINE inletOutlet(f,rho,u,image)
      USE cellConst, ONLY: inlet, outlet
      USE simParam
 
      implicit none
 
      double precision, INTENT(INOUT):: f(yDim,xDim,0:8)&
                      , u(yDim,xDim,0:1), rho(yDim,xDim)
      integer, INTENT(IN):: image(yDim,xDim)
 
      double precision:: uProf
      integer:: x, y
 
      do x = 1, xDim
        do y = 1, yDim
            if (image(y,x) == inlet) then
                u(y,x,0) = uProf(y)
                u(y,x,1) = 0.0d0
                CALL inletZou(f(y,x,:),u(y,x,:),rho(y,x))
            else if (image(y,x) == outlet) then
                u(y,x,0) = uProf(y)
                u(y,x,1) = 0.0d0
                CALL outletZou(f(y,x,:),u(y,x,:),rho(y,x))
            end if
        end do
      end do
 
!      CONTAINS
      END SUBROUTINE
! *********************************************************************************
    !	 ========================================================
    !	 Zou/He boundary on inlet
    !	 ========================================================
      SUBROUTINE inletZou(f,u,rho)
      implicit none
 
      double precision, INTENT(INOUT):: f(0:8),rho
      double precision, INTENT(IN):: u(0:1)
      double precision:: fInt, fInt2
 
      fInt   = f(0) + f(2) + f(4)
      fInt2  = f(3) + f(6) + f(7)
      rho    = (fInt + 2.0d0 * fInt2) / (1.0d0 - u(0))
      CALL zouWestWall(f,rho,u)
      END SUBROUTINE inletZou
! *********************************************************************************
      SUBROUTINE zouWestWall(f,rho,u)
      implicit none
 
      double precision, INTENT(INOUT):: f(0:8)
      double precision, INTENT(IN):: rho, u(0:1)
      double precision:: fDiff, rhoUx, rhoUy
 
      fDiff = 0.5d0 * (f(2) - f(4))
      rhoUx = rho * u(0) / 6.0d0
      rhoUy = 0.5d0 * rho * u(1)
 
      f(1) = f(3) + 4.0d0 * rhoUx
      f(5) = f(7) - fDiff + rhoUx + rhoUy
      f(8) = f(6) + fDiff + rhoUx - rhoUy
      END SUBROUTINE zouWestWall
! *********************************************************************************
    !	 ========================================================
    !	 Zou/He boundary on outlet
    !	 ========================================================
      SUBROUTINE outletZou(f,u,rho)
      implicit none
 
      double precision, INTENT(INOUT):: f(0:8),rho,u(0:1)
      double precision:: fInt, fInt2
 
      fInt  = f(0) + f(2) + f(4)
      fInt2 = f(1) + f(8) + f(5)
      rho   = (fInt + 2.0d0 * fInt2) / (1.0d0 + u(0))
      CALL zouEastWall(f,rho,u)
      END SUBROUTINE outletZou
! *********************************************************************************
      SUBROUTINE zouEastWall(f,rho,u)
      implicit none
 
      double precision, INTENT(INOUT):: f(0:8)
      double precision, INTENT(IN):: rho, u(0:1)
      double precision:: fDiff, rhoUx, rhoUy
 
      fDiff = 0.5d0 * (f(2) - f(4))
      rhoUx = rho * u(0) / 6.0d0
      rhoUy = 0.5d0 * rho * u(1)
 
      f(3) = f(1) - 4.0d0 * rhoUx
      f(7) = f(5) + fDiff - rhoUx - rhoUy
      f(6) = f(8) - fDiff - rhoUx + rhoUy
      END SUBROUTINE zouEastWall
! *********************************************************************************
  !    END SUBROUTINE inletOutlet
    !	 ========================================================
    !	 Computation of Poiseuille profile for the inlet/outlet
    !	 ========================================================
      FUNCTION uProf(y)
      USE simParam, ONLY: yDIm, uMax
      implicit none
 
      integer, INTENT(IN):: y
      double precision:: radius, uProf
 
      radius = dble(yDim-1) * 0.5d0
      uProf  = -uMax * ((abs(1 - dble(y-1) / radius))**2 - 1.0d0)
      END FUNCTION uProf
! *********************************************************************************
    !	 ========================================================
    !	 Streaming step: the population functions are shifted
    !        one site along their corresponding lattice direction
    !        (no temporary memory is needed)
!	 ========================================================
      SUBROUTINE stream(f)
      USE simParam
      implicit none
 
      double precision, INTENT(INOUT):: f(yDim,xDim,0:8)
      double precision:: periodicHor(yDim), periodicVert(xDim)
 
!	 -------------------------------------
!	 right direction
      periodicHor   = f(:,xDim,1)
      f(:,2:xDim,1) = f(:,1:xDim-1,1)
      f(:,1,1)      = periodicHor
!	 -------------------------------------
!	 up direction
      periodicVert    = f(yDim,:,2)
      f(2:yDim,:,2) = f(1:yDim-1,:,2)
      f(1,:,2)     = periodicVert
!	 -------------------------------------
!	 left direction
      periodicHor     = f(:,1,3)
      f(:,1:xDim-1,3) = f(:,2:xDim,3)
      f(:,xDim,3)     = periodicHor
!	 -------------------------------------
!	 down direction
      periodicVert     = f(1,:,4)
      f(1:yDim-1,:,4)  = f(2:yDim,:,4)
      f(yDim,:,4)      = periodicVert
!	 -------------------------------------
!	 up-right direction
      periodicVert = f(yDim,:,5)
      periodicHor  = f(:,xDim,5)
      f(2:yDim,2:xDim,5)   = f(1:yDim-1,1:xDim-1,5)
      f(1,2:xDim,5)        = periodicVert(1:xDim-1)
      f(1,1,5)             = periodicVert(xDim)
      f(2:yDim,1,5)        = periodicHor(1:yDim-1)
!	 -------------------------------------
!	 up-left direction
      periodicVert = f(yDim,:,6)
      periodicHor  = f(:,1,6)
      f(2:yDim,1:xDim-1,6)   = f(1:yDim-1,2:xDim,6)
      f(1,1:xDim-1,6)        = periodicVert(2:xDim)
      f(1,xDim,6)            = periodicVert(1)
      f(2:yDim,xDim,6)       = periodicHor(1:yDim-1)
!	 -------------------------------------
!	 down-left direction
      periodicVert = f(1,:,7)
      periodicHor  = f(:,1,7)
      f(1:yDim-1,1:xDim-1,7) = f(2:yDim,2:xDim,7)
      f(yDim,1:xDim-1,7)     = periodicVert(2:xDim)
      f(yDim,xDim,7)         = periodicVert(1)
      f(1:yDim-1,xDim,7)     = periodicHor(2:yDim)
!	 -------------------------------------
!	 down-right direction
      periodicVert = f(1,:,8)
      periodicHor  = f(:,xDim,8)
      f(1:yDim-1,2:xDim,8)  = f(2:yDim,1:xDim-1,8)
      f(yDim,2:xDim,8)      = periodicVert(1:xDim-1)
      f(yDim,1,8)           = periodicVert(xDim)
      f(1:yDim-1,1,8)         = periodicHor(2:yDim)
      END SUBROUTINE stream
! *********************************************************************************
    !	 ========================================================
    !	 LBGK collision step
    !	 ========================================================
      SUBROUTINE collide(f,fEq,omega,image)
      USE simParam, ONLY: xDim, yDim
      USE cellConst, ONLY: wall
      implicit none
 
      integer, INTENT(IN):: image(yDim,xDim)
      double precision, INTENT(IN):: fEq(yDim,xDim,0:8), omega
      double precision, INTENT(INOUT):: f(yDim,xDim,0:8)
 
      integer:: x,y,i
 
      do i = 0, 8
        do x = 1, xDim
            do y = 1, yDim
                !if (image(y,x) /= wall) 
      f(y,x,i)= (1.0d0-omega)*f(y,x,i)+omega*feq(y,x,i)
            end do
        end do
      end do
      END SUBROUTINE collide
! *********************************************************************************
    !	 ========================================================
    !	 Write the flow geometry to a file
    !	 ========================================================
      SUBROUTINE writeImage(image)
      USE simParam, ONLY: xDim, yDim
      implicit none
 
      integer, INTENT(IN):: image(yDim,xDim)
 
      integer:: x,y
 
      open(13,file='outputImage.dat')
	  write(13,*)"variables=x, y, image"
	  write(13,*)"zone","x=",xDim,"y=",yDim, "F=Point"
      do x=1, xDim
        do y=1, yDim
            write(13,102)x,y, image(y,x)
        end do
      end do
  102 format(3i10)
      close(1)
      END SUBROUTINE writeImage
! *********************************************************************************
    !	 ========================================================
    !	 Print out simulation parameters to screen
    !	 ========================================================
      SUBROUTINE writeInput(omega)
      USE simParam
      implicit none
 
      double precision, INTENT(IN):: omega
       double precision:: Re
      write(*,*) 'xDim                 = ', xDim
      write(*,*) 'yDim                 = ', yDim
      write(*,*) 'Obstacle X           = ', obstX
      write(*,*) 'Obstacle Y           = ', obstY
      write(*,*) 'Obstacle Radius      = ', obstR
      write(*,*) 'tMax                 = ', tMax
      write(*,*) 'uMax                 = ', uMax
      write(*,*) 'Re                   = ', Re
      write(*,*) 'omega                = ', omega
      END SUBROUTINE writeInput
