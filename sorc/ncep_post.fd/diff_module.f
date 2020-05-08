  module diff_module
!------------------------------------------------------------------------
! computes dudy, dvdx, dudx, dvdy, for vorticity and div/conv calculations.
!
! program log:
!   February, 2020    Jesse Meng   Initial code
!------------------------------------------------------------------------
  use masks,        only: dx, dy
  use ctlblk_mod,   only: im, jsta_2l, jend_2u, jsta_m, jend_m, spval
  use gridspec_mod, only: gridtype
!
  implicit none

  REAL, ALLOCATABLE :: DDVDX(:,:)
  REAL, ALLOCATABLE :: DDUDY(:,:)
  REAL, ALLOCATABLE :: UUAVG(:,:)

! 
!-------------------------------------------------------------------------------------
!
  contains
!
!-------------------------------------------------------------------------------------
  subroutine dvdxdudy(uwnd,vwnd)
!
      implicit none

      REAL, dimension(im,jsta_2l:jend_2u), intent(in)    :: UWND, VWND
!
!-- local variables
!--
      integer i, j
      real r2dx, r2dy
      INTEGER, allocatable ::  IHE(:),IHW(:)
!
      IF(GRIDTYPE == 'A')THEN
!$omp parallel do  private(i,j,r2dx,r2dy)
        DO J=JSTA_M,JEND_M
        DO I=2,IM-1
           IF(VWND(I+1,J).LT.SPVAL.AND.VWND(I-1,J).LT.SPVAL.AND.              &
     &        UWND(I,J+1).LT.SPVAL.AND.UWND(I,J-1).LT.SPVAL) THEN
              R2DX   = 1./(2.*DX(I,J))
              R2DY   = 1./(2.*DY(I,J))
              DDVDX(I,J)   = (VWND(I+1,J)-VWND(I-1,J))*R2DX
              DDUDY(I,J)   = (UWND(I,J+1)-UWND(I,J-1))*R2DY
              UUAVG(I,J)   = UWND(I,J)
           END IF
        END DO
        END DO
      ELSE IF (GRIDTYPE == 'E')THEN
        allocate(ihw(JSTA_2L:JEND_2U), IHE(JSTA_2L:JEND_2U))
!$omp  parallel do private(j)
        DO J=JSTA_2L,JEND_2U
          IHW(J) = -MOD(J,2)
          IHE(J) = IHW(J)+1
        ENDDO
!$omp parallel do  private(i,j,r2dx,r2dy)
        DO J=JSTA_M,JEND_M
          DO I=2,IM-1
            IF(VWND(I+IHE(J),J) < SPVAL.AND.VWND(I+IHW(J),J) < SPVAL .AND.   &
     &         UWND(I,J+1) < SPVAL     .AND.UWND(I,J-1) < SPVAL) THEN
               R2DX   = 1./(2.*DX(I,J))
               R2DY   = 1./(2.*DY(I,J))
               DDVDX(I,J)   = (VWND(I+IHE(J),J)-VWND(I+IHW(J),J))*R2DX
               DDUDY(I,J)   = (UWND(I,J+1)-UWND(I,J-1))*R2DY
               UUAVG(I,J)   = 0.25*(UWND(I+IHE(J),J)+UWND(I+IHW(J),J)               &
     &                      +       UWND(I,J+1)+UWND(I,J-1))
            END IF
          END DO
        END DO
        deallocate(ihw, IHE)
      ELSE IF (GRIDTYPE == 'B')THEN
        CALL EXCH_F(VWND)
        DO J=JSTA_M,JEND_M
          DO I=2,IM-1
            R2DX = 1./DX(I,J)
            R2DY = 1./DY(I,J)
            if(VWND(I,  J)==SPVAL .or. VWND(I,  J-1)==SPVAL .or. &
               VWND(I-1,J)==SPVAL .or. VWND(I-1,J-1)==SPVAL .or. &
               UWND(I,  J)==SPVAL .or. UWND(I-1,J  )==SPVAL .or. &
               UWND(I,J-1)==SPVAL .or. UWND(I-1,J-1)==SPVAL) cycle
            DDVDX(I,J) = (0.5*(VWND(I,J)+VWND(I,J-1))-0.5*(VWND(I-1,J) &
     &           +       VWND(I-1,J-1)))*R2DX
            DDUDY(I,J) = (0.5*(UWND(I,J)+UWND(I-1,J))-0.5*(UWND(I,J-1) &
     &           +       UWND(I-1,J-1)))*R2DY
            UUAVG(I,J) = 0.25*(UWND(I-1,J-1)+UWND(I-1,J)               &
     &           +       UWND(I,  J-1)+UWND(I,  J))
          END DO
        END DO
      END IF
!
!--
!
  end subroutine dvdxdudy
!
!--
!
  end module diff_module

