! *****************************************************************************
MODULE ELMFIRE_SUPPRESSION
! *****************************************************************************

USE ELMFIRE_VARS
USE SORT

IMPLICIT NONE

CONTAINS

! *****************************************************************************
SUBROUTINE CENTROID(IT)
! *****************************************************************************
! Computes the centroid (IXCEN, IYCEN cell indices) of the active fire perimeter
! for suppression resource IT by averaging the positions of tagged cells that are
! not yet burned or suppressed. Result is stored in SUPP(IT).

INTEGER, INTENT(IN) :: IT
INTEGER :: I, IXCEN, IYCEN, COUNT
TYPE(NODE), POINTER :: C

IXCEN=0
IYCEN=0
COUNT=0
C => LIST_TAGGED%HEAD

DO I = 1, LIST_TAGGED%NUM_NODES
   IF (C%BURNED) THEN
      C => C%NEXT
      CYCLE
   ENDIF
   IF (C%TIME_SUPPRESSED .GT. 0.) THEN
      C => C%NEXT
      CYCLE
   ENDIF

   COUNT = COUNT + 1
   IXCEN = IXCEN + C%IX
   IYCEN = IYCEN + C%IY
   C => C%NEXT
ENDDO

IF (COUNT .EQ. 0) COUNT = 1
IXCEN=NINT(REAL(IXCEN)/REAL(COUNT))
IYCEN=NINT(REAL(IYCEN)/REAL(COUNT))

SUPP(IT)%IXCEN = IXCEN
SUPP(IT)%IYCEN = IYCEN

! *****************************************************************************
END SUBROUTINE CENTROID
! *****************************************************************************

! *****************************************************************************
SUBROUTINE CONTAINMENT(IT,T)
! *****************************************************************************
! Bins tagged/suppressed fireline cells into 360 angular sectors about the centroid
! and computes current containment from per-bin suppressed and fireline fractions.
! If below the target containment, selects the lowest-spread-velocity sectors and
! marks their tagged cells as suppressed at time T (sets TIME_SUPPRESSED).

INTEGER, INTENT(IN) :: IT
REAL(8), INTENT(IN) :: T
INTEGER :: I, IDEG, J, COUNT, N_NONZERO
REAL :: DX, DY, CURRENT_CONTAINMENT, VELOCITY_SUM
REAL, DIMENSION (0:359) :: SUPPRESSED_FRACTION, FIRELINE_FRACTION, DEG, &
                           VELOCITY_SMOOTHED, VELOCITY0_SMOOTHED
LOGICAL, DIMENSION (0:359) :: SELECTED
INTEGER, PARAMETER :: DEG_SMOOTH_WIDTH=22
TYPE (DLL), POINTER :: L
TYPE(NODE), POINTER :: C

! Start by zeroing arrays
SUPP(IT)%NCELLS(:)=0
SUPP(IT)%VELOCITY(:)=0.
SUPP(IT)%SUPPRESSED_FRACTION(:)=0.

! Now loop over tagged / suppressed cells and determine which degree bin they fall in 
COUNT=0

DO J = 1, 2
   IF (J .EQ. 1) L => LIST_TAGGED
   IF (J .EQ. 2) L => LIST_SUPPRESSED
   C => L%HEAD

   DO I = 1, L%NUM_NODES
      DX = REAL(C%IX - SUPP(IT)%IXCEN) 
      DY = REAL(C%IY - SUPP(IT)%IYCEN)
!      IDEG = NINT(ATAN2D(DY,DX) - 90.0)
      IDEG = NINT(ATAN2(DY,DX)/PIO180 - 90.0)   ! PIO180 = PI/180 (full-precision), was 3.14159
      IF (IDEG .LT.   0) IDEG = IDEG + 360
      IF (IDEG .EQ. 360) IDEG = 0
      C%SUPPRESSION_IDEG = IDEG

      COUNT = COUNT + 1
      SUPP(IT)%NCELLS(IDEG) = SUPP(IT)%NCELLS(IDEG) + 1

      IF (C%TIME_SUPPRESSED .GT. 0.) THEN
         SUPP(IT)%SUPPRESSED_FRACTION(IDEG) = SUPP(IT)%SUPPRESSED_FRACTION(IDEG) + 1.
      ELSE
         SUPP(IT)%VELOCITY(IDEG) = SUPP(IT)%VELOCITY(IDEG) + C%VELOCITY
      ENDIF
      C => C%NEXT
   ENDDO
ENDDO

CONTINUE

! Determine current containment
CURRENT_CONTAINMENT = 0.
DO IDEG = 0, 359
   IF (SUPP(IT)%NCELLS(IDEG) .EQ. 0) THEN
      SUPP(IT)%VELOCITY(IDEG) = 999999.
      SUPP(IT)%SUPPRESSED_FRACTION(IDEG) = 999999.
   ELSE
      SUPP(IT)%VELOCITY(IDEG) = SUPP(IT)%VELOCITY(IDEG) / REAL(SUPP(IT)%NCELLS(IDEG))
      SUPP(IT)%SUPPRESSED_FRACTION(IDEG) = SUPP(IT)%SUPPRESSED_FRACTION(IDEG) / REAL(SUPP(IT)%NCELLS(IDEG))
   ENDIF
   SUPP(IT)%FIRELINE_FRACTION(IDEG) = REAL(SUPP(IT)%NCELLS(IDEG)) / REAL(COUNT)
   CURRENT_CONTAINMENT = CURRENT_CONTAINMENT + SUPP(IT)%SUPPRESSED_FRACTION(IDEG) * SUPP(IT)%FIRELINE_FRACTION(IDEG) 
ENDDO

! Smooth out velocity
DO IDEG = 0, 359
   N_NONZERO=0
   VELOCITY_SUM=0.
   DO I = -DEG_SMOOTH_WIDTH, DEG_SMOOTH_WIDTH
      J = IDEG + I
      IF (J .LT. 0  ) J = J + 360
      IF (J .GT. 359) J = J - 360
      IF (SUPP(IT)%NCELLS(J) .GT. 0) THEN
         N_NONZERO = N_NONZERO + 1
         VELOCITY_SUM = VELOCITY_SUM + SUPP(IT)%VELOCITY(J)
      ENDIF
   ENDDO
   IF (N_NONZERO .GT. 0) THEN
      VELOCITY_SUM = VELOCITY_SUM / REAL(N_NONZERO)
   ELSE
      VELOCITY_SUM = 999999.
   ENDIF
   SUPP(IT)%VELOCITY_SMOOTHED(IDEG) = VELOCITY_SUM
ENDDO

IF (SUPP(IT)%TARGET_CONTAINMENT .GT. CURRENT_CONTAINMENT) THEN
   DO IDEG = 0, 359
      DEG(IDEG) = REAL(IDEG)
   ENDDO
   SUPPRESSED_FRACTION(:) = SUPP(IT)%SUPPRESSED_FRACTION(:)
   FIRELINE_FRACTION  (:) = SUPP(IT)%FIRELINE_FRACTION(:)
   VELOCITY_SMOOTHED  (:) = SUPP(IT)%VELOCITY_SMOOTHED(:)
   VELOCITY0_SMOOTHED (:) = SUPP(IT)%VELOCITY_SMOOTHED(:)
   CALL DSORT(VELOCITY_SMOOTHED(0:), SUPPRESSED_FRACTION(0:), 360, 2); VELOCITY_SMOOTHED(:) = VELOCITY0_SMOOTHED(:)
   CALL DSORT(VELOCITY_SMOOTHED(0:), FIRELINE_FRACTION  (0:), 360, 2); VELOCITY_SMOOTHED(:) = VELOCITY0_SMOOTHED(:)
   CALL DSORT(VELOCITY_SMOOTHED(0:), DEG                (0:), 360, 2)

   ! Select degree bins (in increasing-velocity order) until the containment target is met.
   ! Flag the chosen bins here, then suppress their tagged cells in a single pass afterwards.
   ! This avoids re-walking the entire tagged list once per selected bin (was O(bins x nodes)).
   SELECTED(:) = .FALSE.
   I=-1
   DO WHILE (CURRENT_CONTAINMENT .LT. SUPP(IT)%TARGET_CONTAINMENT .AND. I .LT. 359)
      I = I + 1
      IDEG = INT(DEG(I))
      IF (SUPP(IT)%NCELLS(IDEG) .EQ. 0) CYCLE

      SELECTED(IDEG) = .TRUE.

      CURRENT_CONTAINMENT = CURRENT_CONTAINMENT + (1. - SUPP(IT)%SUPPRESSED_FRACTION(IDEG)) * SUPP(IT)%FIRELINE_FRACTION(IDEG)
      SUPP(IT)%SUPPRESSED_FRACTION(IDEG) = 1.0

   ENDDO

   ! Single pass over the tagged list: suppress every not-yet-suppressed cell whose bin was selected.
   C => LIST_TAGGED%HEAD
   DO J = 1, LIST_TAGGED%NUM_NODES
      IF (SELECTED(C%SUPPRESSION_IDEG) .AND. C%TIME_SUPPRESSED .LT. 0.) THEN
         C%TIME_SUPPRESSED = T
         C%SUPPRESSION_ADJUSTMENT_FACTOR = 0.0
      ENDIF
      C => C%NEXT
   ENDDO
ENDIF

CONTINUE

! *****************************************************************************
END SUBROUTINE CONTAINMENT
! *****************************************************************************





! FUNCTION AND SUBROUTINE FOR NEW SUPPRESSION MODEL
! *****************************************************************************
LOGICAL FUNCTION IS_CLOSE(a, b, tol1, tol2, tol3, tol4)
! *****************************************************************************
    TYPE(Node), INTENT(IN) :: a, b
    REAL, INTENT(IN) :: tol1, tol2, tol3, tol4

    IS_CLOSE = ABS(a%VELOCITY - b%VELOCITY) <= tol1 .AND. &
               ABS(a%FLAME_LENGTH - b%FLAME_LENGTH) <= tol2 .AND. &
               ABS(a%SDI - b%SDI) <= tol3 .AND. &
               ABS(a%PCL - b%PCL) <= tol4 .AND. &
               (a%TYPE_GROUP .EQ. b%TYPE_GROUP)
! *****************************************************************************  
END FUNCTION IS_CLOSE
! *****************************************************************************

! *****************************************************************************
SUBROUTINE ABSORB_SMALL_SEGMENTS(NODES, GRID, N, NX, NY, IX_MIN, IY_MIN, &
                                 NUM_SEGMENTS, SMALL_LIMIT)
! *****************************************************************************

TYPE(NODE_PTR), INTENT(INOUT) :: NODES(:)
INTEGER, INTENT(IN) :: GRID(:,:)
INTEGER, INTENT(IN) :: N, NX, NY
INTEGER, INTENT(IN) :: IX_MIN, IY_MIN
INTEGER, INTENT(IN) :: NUM_SEGMENTS
INTEGER, INTENT(IN) :: SMALL_LIMIT

INTEGER, ALLOCATABLE :: SEG_SIZE(:)
INTEGER, ALLOCATABLE :: SEG_COUNT(:)
INTEGER :: I, S, DX, DY, GX, GY, IX, IY
INTEGER :: NB, NB_SEG
INTEGER :: BEST_SEG, BEST_COUNT

ALLOCATE(SEG_SIZE(NUM_SEGMENTS))
ALLOCATE(SEG_COUNT(NUM_SEGMENTS))

SEG_SIZE = 0

DO I = 1, N
   S = NODES(I)%P%SEGMENT_GROUP
   IF (S .GE. 1 .AND. S .LE. NUM_SEGMENTS) THEN
      SEG_SIZE(S) = SEG_SIZE(S) + 1
   ENDIF
ENDDO

DO I = 1, N

   S = NODES(I)%P%SEGMENT_GROUP

   IF (S .LT. 1 .OR. S .GT. NUM_SEGMENTS) CYCLE
   IF (SEG_SIZE(S) > SMALL_LIMIT) CYCLE

   SEG_COUNT = 0

   DO DX = -1, 1
      DO DY = -1, 1

         IF (DX .EQ. 0 .AND. DY .EQ. 0) CYCLE

         IX = NODES(I)%P%IX + DX
         IY = NODES(I)%P%IY + DY

         GX = IX - IX_MIN + 1
         GY = IY - IY_MIN + 1

         IF (GX .LT. 1 .OR. GX .GT. NX) CYCLE
         IF (GY .LT. 1 .OR. GY .GT. NY) CYCLE

         NB = GRID(GX, GY)

         IF (NB .EQ. 0) CYCLE

         NB_SEG = NODES(NB)%P%SEGMENT_GROUP

         IF (NB_SEG .LT. 1 .OR. NB_SEG .GT. NUM_SEGMENTS) CYCLE
         IF (NB_SEG .EQ. S) CYCLE
         IF (SEG_SIZE(NB_SEG) .LE. SMALL_LIMIT) CYCLE

         SEG_COUNT(NB_SEG) = SEG_COUNT(NB_SEG) + 1

      ENDDO
   ENDDO

   BEST_SEG = S
   BEST_COUNT = 0

   DO NB_SEG = 1, NUM_SEGMENTS
      IF (SEG_COUNT(NB_SEG) .GT. BEST_COUNT) THEN
         BEST_COUNT = SEG_COUNT(NB_SEG)
         BEST_SEG = NB_SEG
      ENDIF
   ENDDO

   IF (BEST_SEG .NE. S) THEN
      NODES(I)%P%SEGMENT_GROUP = BEST_SEG
   ENDIF

ENDDO

DEALLOCATE(SEG_SIZE)
DEALLOCATE(SEG_COUNT)

! *****************************************************************************
END SUBROUTINE ABSORB_SMALL_SEGMENTS
! *****************************************************************************


! *****************************************************************************
SUBROUTINE SEGMENT_FIRELINE
! *****************************************************************************
! REAL(8), INTENT(IN) :: T
! REAL,    INTENT(IN) :: DT
TYPE(NODE), POINTER :: C, C2
REAL :: C_ELLIPSE, X_PNT, COSANG
INTEGER  :: IX_MIN, IX_MAX, IY_MIN, IY_MAX
INTEGER :: I, I2, SEGMENT_ID, N, NX, NY, GX, GY, DX, DY, IX, IY, HEAD, TAIL, CURRENT, NB
INTEGER :: IX1, IY1, IX2, IY2, IX3, IY3, IX4, IY4, IX5, IY5, IX6, IY6, IX7, IY7, IX8, IY8, TEST_INDX

TYPE(NODE_PTR), ALLOCATABLE :: NODES(:)
INTEGER, ALLOCATABLE :: GRID(:,:)
INTEGER, ALLOCATABLE :: QUEUE(:)
INTEGER, ALLOCATABLE :: MAP_SEG(:)
INTEGER :: OLD_SEG

! detect fireline
C => LIST_TAGGED%HEAD
DO I = 1, LIST_TAGGED%NUM_NODES
   IF (C%TIME_OF_ARRIVAL .EQ. -1) THEN
      DO DX = -FIRE_LINE_THICKNESS, FIRE_LINE_THICKNESS
         DO DY = -FIRE_LINE_THICKNESS, FIRE_LINE_THICKNESS
            IF (DX .EQ. 0 .AND. DY .EQ. 0) CYCLE
            IX = C%IX + DX
            IY = C%IY + DY
            IF (TIME_OF_ARRIVAL(IX,IY) .GT. 0) C%FIRE_LINE = .TRUE.
         ENDDO
      ENDDO
   ENDIF
   
   C => C%NEXT
ENDDO

! remove isolated fireline
C => LIST_TAGGED%HEAD
DO I = 1, LIST_TAGGED%NUM_NODES
   IF (C%FIRE_LINE) THEN
      IX1 = C%IX - 1
      IY1 = C%IY
      IX2 = C%IX - 1
      IY2 = C%IY + 1
      IX3 = C%IX
      IY3 = C%IY + 1
      IX4 = C%IX + 1
      IY4 = C%IY + 1
      IX5 = C%IX + 1
      IY5 = C%IY
      IX6 = C%IX + 1
      IY6 = C%IY - 1
      IX7 = C%IX
      IY7 = C%IY - 1
      IX8 = C%IX - 1
      IY8 = C%IY - 1   
      
      C2 => LIST_TAGGED%HEAD
      TEST_INDX = 0
      DO I2 = 1, LIST_TAGGED%NUM_NODES
         IF (C2%IX .EQ. IX1 .AND. C2%IY .EQ. IY1 .AND. C2%FIRE_LINE) TEST_INDX = TEST_INDX + 1
         IF (C2%IX .EQ. IX2 .AND. C2%IY .EQ. IY2 .AND. C2%FIRE_LINE) TEST_INDX = TEST_INDX + 1
         IF (C2%IX .EQ. IX3 .AND. C2%IY .EQ. IY3 .AND. C2%FIRE_LINE) TEST_INDX = TEST_INDX + 1
         IF (C2%IX .EQ. IX4 .AND. C2%IY .EQ. IY4 .AND. C2%FIRE_LINE) TEST_INDX = TEST_INDX + 1
         IF (C2%IX .EQ. IX5 .AND. C2%IY .EQ. IY5 .AND. C2%FIRE_LINE) TEST_INDX = TEST_INDX + 1
         IF (C2%IX .EQ. IX6 .AND. C2%IY .EQ. IY6 .AND. C2%FIRE_LINE) TEST_INDX = TEST_INDX + 1
         IF (C2%IX .EQ. IX7 .AND. C2%IY .EQ. IY7 .AND. C2%FIRE_LINE) TEST_INDX = TEST_INDX + 1
         IF (C2%IX .EQ. IX8 .AND. C2%IY .EQ. IY8 .AND. C2%FIRE_LINE) TEST_INDX = TEST_INDX + 1
         C2 => C2%NEXT
      ENDDO

      IF (TEST_INDX .EQ. 0) C%FIRE_LINE = .FALSE.

   ENDIF
   
   C => C%NEXT
ENDDO


! C => LIST_TAGGED%HEAD
! DO I = 1, LIST_TAGGED%NUM_NODES
!    IF (C%TIME_OF_ARRIVAL .GT. T-DT .AND. C%TIME_OF_ARRIVAL .LT. T+DT) THEN
!        C%FIRE_LINE = .TRUE.
!    ENDIF

!    C => C%NEXT
! ENDDO


N = LIST_TAGGED%NUM_NODES
IF (N .LE. 0) THEN
   LIST_TAGGED%NUM_SEGMENTS = 0
   RETURN
END IF


! determine fireline type.
C => LIST_TAGGED%HEAD
DO I = 1, LIST_TAGGED%NUM_NODES

   IF (.NOT. C%FIRE_LINE) THEN
      C => C%NEXT
      CYCLE
   ENDIF

   ! We can get sin(theta - dms) and cos(theta - dms) directly:
   COSANG   = C%NORMVECTORY*C%NORMVECTORY_DMS + C%NORMVECTORX*C%NORMVECTORX_DMS

   C_ELLIPSE = (C%VELOCITY_DMS - C%VBACK)*0.5
   X_PNT = C%VELOCITY*COSANG

   ! TYPE_GROUP = 2:head; 1:flank; 0:back
   IF (X_PNT .GE. C_ELLIPSE) THEN
      C%TYPE_GROUP = 2
   ELSE IF (X_PNT .GE. 0) THEN
      C%TYPE_GROUP = 1
   ELSE
      C%TYPE_GROUP = 0
   ENDIF

   C => C%NEXT
ENDDO

! find IX_MIN, IX_MAX, IY_MIN, IY_MAX
IX_MIN = 9999; IX_MAX = -9999
IY_MIN = 9999; IY_MAX = -9999

C => LIST_TAGGED%HEAD
DO I = 1, LIST_TAGGED%NUM_NODES
   IF (C%IX .GT. IX_MAX) IX_MAX = C%IX
   IF (C%IX .LT. IX_MIN) IX_MIN = C%IX
   IF (C%IY .GT. IY_MAX) IY_MAX = C%IY
   IF (C%IY .LT. IY_MIN) IY_MIN = C%IY
   C => C%NEXT
ENDDO

! segment the fireline

NX = IX_MAX - IX_MIN + 1
NY = IY_MAX - IY_MIN + 1

ALLOCATE(NODES(N))
ALLOCATE(GRID(NX, NY))
ALLOCATE(QUEUE(N))

GRID = 0

C => LIST_TAGGED%HEAD

DO I = 1, N

   NODES(I)%P => C

   IF (.NOT. C%FIRE_LINE) THEN
      C => C%NEXT
      CYCLE
   ENDIF


   GX = C%IX - IX_MIN + 1
   GY = C%IY - IY_MIN + 1

   IF (GX .LT. 1 .OR. GX .GT. NX .OR. GY .LT. 1 .OR. GY .GT. NY) THEN
      WRITE(*,*) 'SUPPRESSION MODULE ERROR: NODE OUTSIDE GRID BOUNDS'
      STOP
   END IF

   GRID(GX, GY) = I
   C%SEGMENT_GROUP = 0

   C => C%NEXT
ENDDO

SEGMENT_ID = 0

DO I = 1, N

   IF (NODES(I)%P%SEGMENT_GROUP .NE. 0) CYCLE

   SEGMENT_ID = SEGMENT_ID + 1

   HEAD = 1
   TAIL = 1
   QUEUE(TAIL) = I
   NODES(I)%P%SEGMENT_GROUP = SEGMENT_ID

   DO WHILE (HEAD .LE. TAIL)
      CURRENT = QUEUE(HEAD)
      HEAD = HEAD + 1

      DO DX = -1, 1
         DO DY = -1, 1
            IF (DX .EQ. 0 .AND. DY .EQ. 0) CYCLE

            IX = NODES(CURRENT)%P%IX + DX
            IY = NODES(CURRENT)%P%IY + DY

            GX = IX - IX_MIN + 1
            GY = IY - IY_MIN + 1

            IF (GX .LT. 1 .OR. GX .GT. NX) CYCLE
            IF (GY .LT. 1 .OR. GY .GT. NY) CYCLE

            NB = GRID(GX, GY)

            IF (NB .EQ. 0) CYCLE

            IF (NODES(NB)%P%SEGMENT_GROUP .NE. 0) CYCLE

            IF (IS_CLOSE(NODES(CURRENT)%P, NODES(NB)%P, DELTA_ROS, DELTA_FL, DELTA_SDI, DELTA_PCL)) THEN
               TAIL = TAIL + 1
               QUEUE(TAIL) = NB
               NODES(NB)%P%SEGMENT_GROUP = SEGMENT_ID
            ENDIF
         ENDDO
      ENDDO
   ENDDO
ENDDO

CALL ABSORB_SMALL_SEGMENTS(NODES, GRID, N, NX, NY, IX_MIN, IY_MIN, SEGMENT_ID, 3)
LIST_TAGGED%NUM_SEGMENTS = SEGMENT_ID
DEALLOCATE(NODES)
DEALLOCATE(GRID)
DEALLOCATE(QUEUE)

!correcting segment_id and numbers
ALLOCATE(MAP_SEG(LIST_TAGGED%NUM_SEGMENTS))
MAP_SEG = 0
C => LIST_TAGGED%HEAD
SEGMENT_ID = 0

DO I = 1, LIST_TAGGED%NUM_NODES

   IF (.NOT. C%FIRE_LINE) THEN
      C => C%NEXT
      CYCLE
   ENDIF

   OLD_SEG = C%SEGMENT_GROUP

   IF (MAP_SEG(OLD_SEG) .EQ. 0) THEN
      SEGMENT_ID = SEGMENT_ID + 1
      MAP_SEG(OLD_SEG) = SEGMENT_ID
   ENDIF

   C%SEGMENT_GROUP = MAP_SEG(OLD_SEG)

   C => C%NEXT

ENDDO


LIST_TAGGED%NUM_SEGMENTS = SEGMENT_ID
DEALLOCATE(MAP_SEG)

! *****************************************************************************
END SUBROUTINE SEGMENT_FIRELINE
! *****************************************************************************

! *****************************************************************************
SUBROUTINE CALCULATE_STS
! *****************************************************************************

! STS: Suppression Type Score

TYPE(NODE), POINTER :: C
REAL, ALLOCATABLE ::  SDI_AVG(:), FL_AVG(:), N_CELLS(:)
INTEGER :: I
REAL :: FL_NI, SDI_NI

ALLOCATE(SUPPRESSION_TYPE_SCORE(LIST_TAGGED%NUM_SEGMENTS))
ALLOCATE(FL_AVG(LIST_TAGGED%NUM_SEGMENTS))
ALLOCATE(SDI_AVG(LIST_TAGGED%NUM_SEGMENTS))
ALLOCATE(N_CELLS(LIST_TAGGED%NUM_SEGMENTS))

SDI_AVG = 0.0
FL_AVG  = 0.0
N_CELLS = 0.0


! calculate max of FL and SDI over fireline
! calculate mean of SDI and FL over each segment
C => LIST_TAGGED%HEAD
DO I = 1, LIST_TAGGED%NUM_NODES

   IF (.NOT. C%FIRE_LINE) THEN
      C => C%NEXT
      CYCLE
   ENDIF

   SDI_AVG(C%SEGMENT_GROUP) = SDI_AVG(C%SEGMENT_GROUP) + C%SDI
   FL_AVG(C%SEGMENT_GROUP)  = FL_AVG(C%SEGMENT_GROUP)  + C%FLAME_LENGTH
   N_CELLS(C%SEGMENT_GROUP) = N_CELLS(C%SEGMENT_GROUP) + 1.0

   C => C%NEXT
ENDDO

! calculate STS for seg
DO I = 1, LIST_TAGGED%NUM_SEGMENTS

   IF (N_CELLS(I) .EQ. 0) CYCLE

   FL_AVG(I) = FL_AVG(I)/N_CELLS(I)
   FL_NI = (FL_AVG(I)/FL_MAX_DIRECT_ATTACK)

   SDI_AVG(I) = SDI_AVG(I)/N_CELLS(I)
   SDI_NI = (SDI_AVG(I)/SDI_MAX_DIRECT_ATTACK)

   SUPPRESSION_TYPE_SCORE(I) = (0.7*FL_NI) + (0.3*SDI_NI)
ENDDO


! assign sts to each node
C => LIST_TAGGED%HEAD
DO I = 1, LIST_TAGGED%NUM_NODES

   IF (.NOT. C%FIRE_LINE) THEN
      C => C%NEXT
      CYCLE
   ENDIF

   C%STS = SUPPRESSION_TYPE_SCORE(C%SEGMENT_GROUP)

   C => C%NEXT
ENDDO


DEALLOCATE(FL_AVG)
DEALLOCATE(SDI_AVG)
DEALLOCATE(N_CELLS)


! *****************************************************************************
END SUBROUTINE CALCULATE_STS
! *****************************************************************************

! *****************************************************************************
SUBROUTINE SORT_STS
! *****************************************************************************
INTEGER :: I, J, TEMP_INDX
ALLOCATE(SUPPRESSION_TYPE_SCORE_RANK(LIST_TAGGED%NUM_SEGMENTS))

! initiate the rank array
DO I = 1, LIST_TAGGED%NUM_SEGMENTS
   SUPPRESSION_TYPE_SCORE_RANK(I) = i
ENDDO

! sorting index
DO I = 1, LIST_TAGGED%NUM_SEGMENTS-1
   DO J = I+1, LIST_TAGGED%NUM_SEGMENTS
      IF (SUPPRESSION_TYPE_SCORE(SUPPRESSION_TYPE_SCORE_RANK(J)) .LT. SUPPRESSION_TYPE_SCORE(SUPPRESSION_TYPE_SCORE_RANK(I))) THEN
         TEMP_INDX = SUPPRESSION_TYPE_SCORE_RANK(I)
         SUPPRESSION_TYPE_SCORE_RANK(I) = SUPPRESSION_TYPE_SCORE_RANK(J)
         SUPPRESSION_TYPE_SCORE_RANK(J) = TEMP_INDX
      ENDIF
   ENDDO
ENDDO

! *****************************************************************************
END SUBROUTINE SORT_STS
! *****************************************************************************


! *****************************************************************************
SUBROUTINE DIRECT_ATTACK(T, IT, rank_finished, DT, ICASE, TSTOP)
! *****************************************************************************
REAL(8), INTENT(IN) :: T
INTEGER, INTENT(IN) :: IT
INTEGER, INTENT(IN) :: ICASE
REAL, INTENT(INOUT)    :: DT, TSTOP
INTEGER, INTENT(INOUT) :: rank_finished
TYPE(NODE), POINTER :: C
INTEGER :: I, J, J2, J3, N_CELL_AVAIL
REAL :: CAP_0, L_CAP, L_AVAIL
REAL, ALLOCATABLE :: L_SEG(:), L_REQ(:), FL_MIN(:)
LOGICAL :: IS_MIN

IF (IT .NE. 1) THEN
   SUPP(IT)%FIRE_LINE_LENGTH = SUPP(IT-1)%FIRE_LINE_LENGTH
   SUPP(IT)%SUPPRESSED_FIRELINE_LENGTH = SUPP(IT-1)%SUPPRESSED_FIRELINE_LENGTH
ENDIF


CAP_0 = FIRE_LINE_THICKNESS * AVAILABLE_SUPPRESSION_CAPACITY * (1/3600.0)
L_CAP = CAP_0*DT_EXTENDED_ATTACK
ALLOCATE(L_SEG(LIST_TAGGED%NUM_SEGMENTS))
ALLOCATE(L_REQ(LIST_TAGGED%NUM_SEGMENTS))

L_SEG = 0
L_REQ = 0

! Calculate fireline length in each segement
C => LIST_TAGGED%HEAD
DO I = 1, LIST_TAGGED%NUM_NODES

   IF (.NOT. C%FIRE_LINE) THEN
      C => C%NEXT
      CYCLE
   ENDIF

   SUPP(IT)%FIRE_LINE_LENGTH = SUPP(IT)%FIRE_LINE_LENGTH + ANALYSIS_CELLSIZE

   L_SEG(C%SEGMENT_GROUP) = L_SEG(C%SEGMENT_GROUP) + ANALYSIS_CELLSIZE

   C => C%NEXT
ENDDO

! consider effect of sts on length in each segment
DO I = 1, LIST_TAGGED%NUM_SEGMENTS
   IF (1-SUPPRESSION_TYPE_SCORE(I) .GT. 0) THEN
      L_REQ(I) = L_SEG(I)/(1-SUPPRESSION_TYPE_SCORE(I))
   ELSEIF (1-SUPPRESSION_TYPE_SCORE(I) .EQ. 0) THEN
      L_REQ(I) = L_SEG(I)
   ELSE
      L_REQ(I) = 0
   ENDIF
ENDDO

! Apply direct attack
I = 0
L_CAP = MIN(L_CAP, SUPP(IT)%FIRE_LINE_LENGTH)
DO WHILE (L_CAP .GT. 0 .AND. I .LT. LIST_TAGGED%NUM_SEGMENTS)

   I = I + 1

   IF (L_REQ(SUPPRESSION_TYPE_SCORE_RANK(I)) .EQ. 0) CYCLE

   IF (L_CAP .GT. L_REQ(SUPPRESSION_TYPE_SCORE_RANK(I))) THEN

      ! SUPP(IT)%SUPPRESSED_FIRELINE_LENGTH = SUPP(IT)%SUPPRESSED_FIRELINE_LENGTH + L_REQ(SUPPRESSION_TYPE_SCORE_RANK(I))

      C => LIST_TAGGED%HEAD
      DO J = 1, LIST_TAGGED%NUM_NODES

         IF (.NOT. C%FIRE_LINE .OR. C%SEGMENT_GROUP .NE. SUPPRESSION_TYPE_SCORE_RANK(I)) THEN
            C => C%NEXT
            CYCLE
         ENDIF

         C%TIME_SUPPRESSED = T
         C%SUPPRESSION_ADJUSTMENT_FACTOR = 0.0
         SUPP(IT)%SUPPRESSED_FIRELINE_LENGTH = SUPP(IT)%SUPPRESSED_FIRELINE_LENGTH + ANALYSIS_CELLSIZE

         C => C%NEXT

      ENDDO

      L_CAP = L_CAP - L_REQ(SUPPRESSION_TYPE_SCORE_RANK(I))

   ELSE

      L_AVAIL = L_CAP * (1-SUPPRESSION_TYPE_SCORE(SUPPRESSION_TYPE_SCORE_RANK(I)))

      ! SUPP(IT)%SUPPRESSED_FIRELINE_LENGTH = SUPP(IT)%SUPPRESSED_FIRELINE_LENGTH + L_AVAIL

      N_CELL_AVAIL = CEILING(L_AVAIL/ANALYSIS_CELLSIZE)
      ALLOCATE(FL_MIN(N_CELL_AVAIL))
      DO J = 1, N_CELL_AVAIL
         FL_MIN(J) = 999999.0
      ENDDO

      C => LIST_TAGGED%HEAD
      DO J = 1, LIST_TAGGED%NUM_NODES
         IF (.NOT. C%FIRE_LINE .OR. C%SEGMENT_GROUP .NE. SUPPRESSION_TYPE_SCORE_RANK(I)) THEN
            C => C%NEXT
            CYCLE
         ENDIF


         DO J2 = 1, N_CELL_AVAIL

            IF (C%FLAME_LENGTH .LT. FL_MIN(J2)) THEN

               DO J3 = N_CELL_AVAIL, J2 + 1, -1
                  FL_MIN(J3) = FL_MIN(J3 - 1)
               ENDDO

               FL_MIN(J2) = C%FLAME_LENGTH
               EXIT

            ENDIF
         ENDDO

         C => C%NEXT
      ENDDO

      C => LIST_TAGGED%HEAD
      DO J = 1, LIST_TAGGED%NUM_NODES

         IF (.NOT. C%FIRE_LINE .OR. C%SEGMENT_GROUP .NE. SUPPRESSION_TYPE_SCORE_RANK(I) .OR. L_AVAIL .LE. 0) THEN
            C => C%NEXT
            CYCLE
         ENDIF

         IS_MIN = .FALSE.

         DO J2 = 1, N_CELL_AVAIL
            IF (C%FLAME_LENGTH .EQ. FL_MIN(J2)) IS_MIN = .TRUE.
         ENDDO

         IF (.NOT. IS_MIN) THEN
            C => C%NEXT
            CYCLE
         ENDIF

         C%TIME_SUPPRESSED = T
         C%SUPPRESSION_ADJUSTMENT_FACTOR = 0.0
         SUPP(IT)%SUPPRESSED_FIRELINE_LENGTH = SUPP(IT)%SUPPRESSED_FIRELINE_LENGTH + ANALYSIS_CELLSIZE

         L_AVAIL = L_AVAIL - ANALYSIS_CELLSIZE

         C => C%NEXT

      ENDDO

      L_CAP = -1.0
      DEALLOCATE(FL_MIN)

   ENDIF


ENDDO

SUPP(IT)%CONTAINMENT = (SUPP(IT)%SUPPRESSED_FIRELINE_LENGTH / SUPP(IT)%FIRE_LINE_LENGTH)

WRITE(*,*) SUPP(IT)%T, SUPP(IT)%SUPPRESSED_FIRELINE_LENGTH, SUPP(IT)%FIRE_LINE_LENGTH, SUPP(IT)%CONTAINMENT

IF (SUPP(IT)%CONTAINMENT .GE. 0.99) THEN
   rank_finished = 1
   DT = DT_METEOROLOGY
   STATS_FINAL_CONTAINMENT_FRAC(ICASE) = 1.0
   STATS_SIMULATION_TSTOP_HOURS(ICASE) = T / 3600.0
   TSTOP = T
ENDIF

! *****************************************************************************
END SUBROUTINE DIRECT_ATTACK
! *****************************************************************************


! *****************************************************************************
END MODULE ELMFIRE_SUPPRESSION
! *****************************************************************************
