  module diff_module
!------------------------------------------------------------------------
! computes dudy, dvdx, dudx, dvdy, for vorticity and div/conv calculations.
!
! program log:
!   February, 2020    Jesse Meng   Initial code
!------------------------------------------------------------------------
!  use masks, only: dx, dy
!
  implicit none
  
  real dx1, dy1
  real u1, u2
  real v1, v2
  real dudx1
  real dudy1
  real dvdx1
  real dvdy1
! 
!-------------------------------------------------------------------------------------
!
  contains
!
!-------------------------------------------------------------------------------------
  subroutine dvdxdudy(u1,u2,v1,v2,dx1,dy1,dvdx1,dudy1)
!
    implicit none

    real u1, u2, v1, v2, dx1, dy1
    real r2dx, r2dy
    real dvdx1, dudy1
!
!-- local variables
!--
              R2DX   = 1./(2.*DX1)
              R2DY   = 1./(2.*DY1)
              DVDX1  = (V2-V1)*R2DX
              DUDY1  = (U2-U1)*R2DY
!
!--
!
  end subroutine dvdxdudy
!-------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
!
  subroutine dudxdvdy(u1,u2,v1,v2,dx1,dy1,dudx1,dvdy1)
!
    implicit none

    real u1, u2, v1, v2, dx1, dy1
    real r2dx, r2dy
    real dudx1, dvdy1
!   
!-- local variables
!--
              R2DX    = 1./(2.*DX1)
              R2DY    = 1./(2.*DY1)
              DUDX1   = (U2-U1)*R2DX
              DVDY1   = (V2-V1)*R2DY
!
!--
!
  end subroutine dudxdvdy
!
!--
!
  end module diff_module

