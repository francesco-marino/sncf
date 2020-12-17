C     Information                                 p.cusdin@qub.ac.uk
C     ==============================================================
C     This code is  available for use free of charge.  QUB School of
C     Aeronautical Engineering  is not liable for the correctness or
C     accuracy of data that results from the use of this routine nor
C     is it accountable for damages that may be casued by its use.
C     The author will be pleased to hear how  the code is being used
C     in any research or commercial enterprise.
C     ==============================================================
      
C*********************************************************************
C Timer programme to create a log file which conatins the events
C of a script as called from the main routine
C     
C     USEAGE
C For correct timing, the log file must be initialised at the start
C of your programme using:
C     CALL TIMELOG(.true.,msg,time1,sum)
C      
C Subsequent entries are applied using the following:
C     msg = 'My Message'
C     CALL TIMELOG(.false.,msg,time1,sum)
C      
C You will also need to declare the following variables:
C     double precision time1
C     double precision sum
C     character msg*20
C     
C Compile and link in the normal way TIMELOG cleans the atime.log
C file each time the primal  program starts: to keep old data you
C must rename that file at the end of each program
C *******************
C P.Cusdin 05.12.2002
C*********************************************************************

      SUBROUTINE TIMELOG(first,msg,etime1,last,elap)

      implicit none

c     external time

      double precision last
      double precision this
      double precision etime1
      double precision elap
      double precision sum

      character datum*24
      character msg*20

      integer tfile
      logical first

      parameter (tFile = 9)
C
C     Set log file format
C     -------------------
 101  FORMAT (A20,2X,2(F9.2,2X))
 102  FORMAT (A20,2X,A9,2X,A9,2X)
      
      IF(first) THEN
C
C       Open (and wipe) time log file
C       -----------------------------
        OPEN (tFile,FILE='atime.log',STATUS='UNKNOWN')
        REWIND tFile
        
c       call cpu_time(etime1)
        WRITE (tFile,*) 'Courtesy of QUB Sc. of Aeronautical Eng.'
        WRITE (tFile,*) '***      http://www.ea.qub.ac.uk     ***'
        WRITE (tFile,*) '----------------------------------------'
        WRITE (tFile,102) 'NOTE','CALL','TOTAL'
        CLOSE(tFile)

      ELSE
C
C       Open (and append) time log file
C       -------------------------------
        OPEN (tFile,FILE='atime.log',STATUS='UNKNOWN',ACCESS='append')
C
C       Get CPU data
C       ------------
c       call cpu_time(this)
        elap = (this - last)
        sum = (this - etime1)
        WRITE (tFile,101) msg,elap,sum
C
C       Reset last time
C       --------------
        last = this
        CLOSE(tFile)
        
      END IF

      RETURN
      END

