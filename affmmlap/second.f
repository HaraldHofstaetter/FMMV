      FUNCTION SECOND()
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c  purpose:
c
c    a subroutine to get the cpu time for timing purposes.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
      IMPLICIT NONE
      REAL *8 SECOND
      REAL *4 TIMER
c
      CALL CPU_TIME(TIMER)
      SECOND=TIMER
c
      RETURN
      END
