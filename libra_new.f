C****************************************************************************
C
C     BENT LAURITZEN LIBRARY (received May 1993)
C     (Modified for Alpha in Jan 1996)
C     (Modified for Digital Unix in Nov 1998, G. Colo`) 
C****************************************************************************


      FUNCTION WUE(IA,L)
      PI = 3.14159
      XL = L
      A = IA
      AEL = (1.2*(A**(1./3.)))**(2*XL)
      Y = 3./(3.+XL)
      WUE = (1./(4.*PI))*Y*Y*AEL
      RETURN
      END

C*******************************************************************************

      FUNCTION WUM(IA,L)
      PI = 3.14159
      XL = L
      A = IA
      AML = (1.2*(A**(1./3.)))**(2*XL-2.)
      Y = 3./(3.+XL)
      WUM = (10./(PI))*Y*Y*AML
      RETURN
      END

C*******************************************************************************

      FUNCTION OSCRAD(I1,XR,BHBW)
      COMMON/BINDEX/XL(3,200)
      F(X,Y) = SQRT(DFAC(2.*X+2.*Y+1.)/FAC(X))/DFAC(2.*Y+1.)
      G(X,Y,Z) = PARI(Z)*BINOM(X,Z)*DFAC(2.*Y+1.)/DFAC(2.*Y+2.*Z+1.)
      B2 = BHBW*BHBW
CBAB      V = 41.465/HBW
      V = 1./B2
      XN = XL(1,I1)-1.
      XLL = XL(2,I1)
      CC = V**(.75+XLL/2.)
      XX = F(XN,XLL)
      C = 1.50225*(2.**((XLL/2.)-(XN/2.)))
      IXN = XN+1.
      R = XR
      S = 0.
      DO 21 J = 1,IXN
      XK = J-1
      IF(R.EQ.0..AND.XK.EQ.0.) RRR=1.
      IF(R.EQ.0..AND.XK.NE.0.) RRR=0.
      IF(R.NE.0.) RRR=R**(2.*XK)
21    S = S + (V**XK)*G(XN,XLL,XK)*(RRR)*(2.**XK)
      YY = V*R*R/2.
      IF(R.EQ.0..AND.XLL.EQ.0.) RRR=1.
      IF(R.EQ.0..AND.XLL.NE.0.) RRR=0.
      IF(R.NE.0.) RRR=R**XLL
      SS = XX*EXP(-YY)*S*(RRR)
      SS = SS*C
      OSCRAD = SS*CC
      RETURN
      END

C*******************************************************************************

      FUNCTION XLSJJ(X1,X2,X3,X4,X5,X6,X7,X8,X9)
      X = (2.*X3+1.)*(2.*X6+1.)*(2.*X7+1.)*(2.*X8+1.)
      XLSJJ = SQRT(X)*QJ(X1,X2,X3,X4,X5,X6,X7,X8,X9)
      RETURN
      END

C*******************************************************************************

      FUNCTION QJ(X1,X2,X3,X4,X5,X6,X7,X8,X9)
      QJ = 0.
      J1 = 2.*X1
      J2 = 2.*X2
      J3 = 2.*X3
      J4 = 2.*X4
      J5 = 2.*X5
      J6 = 2.*X6
      J7 = 2.*X7
      J8 = 2.*X8
      J9 = 2.*X9
      KMIN=MAX0(IABS(J1-J9),IABS(J2-J6),IABS(J4-J8))
      KMAX=MIN0(J1+J9,J2+J6,J4+J8)
      IF (KMIN.GT.KMAX) RETURN
      DO 100 K=KMIN,KMAX,2
      XK = K/2.
      A=SJ(X1,X4,X7,X8,X9,XK)
      B=SJ(X2,X5,X8,X4,XK,X6)
      C=SJ(X3,X6,X9,XK,X1,X2)
      QJ=QJ+PARI(2.*XK)*(2.*XK+1.)*A*B*C
100   CONTINUE
      RETURN
      END

C*******************************************************************************

      FUNCTION DELTA(X)     ! kronecker delta function (real argument)
      IF(X) 10,20,10
20    DELTA = 1.0
      RETURN
10    DELTA = 0.0
      RETURN
      END

C*******************************************************************************

      FUNCTION IDELTA(I)     ! kronecker delta function (integer argument)
      IF(I) 10,20,10
20    IDELTA = 1
      RETURN
10    IDELTA = 0
      RETURN
      END

C*******************************************************************************

      FUNCTION TR(A1,A2,A3)
      IF(A1.LT.0.) GO TO 2
      IF(A2.LT.0.) GO TO 2
      IF(A3.LT.0.) GO TO 2
      AMIN = ABS(A1-A2)
      AMAX = A1+A2
      XA = AMIN
10    CONTINUE
      IF(A3.EQ.XA) GO TO 20
      XA = XA+1.
      IF(XA.GT.AMAX) GO TO 2
      GO TO 10
20    TR=1.
      RETURN
2     TR=-1.
      RETURN
      END

C*******************************************************************************

      FUNCTION ITR(IA1,IA2,IA3)
      IF(IA1.LT.0.) GO TO 2
      IF(IA2.LT.0.) GO TO 2
      IF(IA3.LT.0.) GO TO 2
      IAMIN = IABS(IA1-IA2)
      IAMAX = IA1+IA2
      IXA = IAMIN
10    CONTINUE
      IF(IA3.EQ.IXA) GO TO 20
      IXA = IXA+1
      IF(IXA.GT.IAMAX) GO TO 2
      GO TO 10
20    ITR=1
      RETURN
2     ITR=-1
      RETURN
      END

C*******************************************************************************

      FUNCTION DFACN(A)
      DFACN = DFAC(A)
      RETURN
      END

C*******************************************************************************

      FUNCTION FACS(X)
      DIMENSION FACS_SAVE(100)
      IS = X+1
      IF(X) 3,4,4
3     FACS = 0.0
      RETURN
4     FACS = FACS_SAVE(IS)
      IF(FACS.NE.0.) RETURN
      S = 1.
      DO 1 I = 1,1000
      Y = I
      IF(X-Y) 2,1,1
1     S = S*SQRT(Y)
2     FACS = S
      FACS_SAVE(IS) = FACS
      RETURN
      END

C*******************************************************************************

      FUNCTION FAC(X)
      DIMENSION FAC_SAVE(100)
      IS = X+1
      IF(X) 3,4,4
3     FAC = 0.0
      RETURN
4     FAC = FAC_SAVE(IS)
      IF(FAC.NE.0.) RETURN
      S = 1.
      DO 1 I = 1,1000
      Y = I
      IF(X-Y) 2,1,1
1     S = S*Y
2     FAC = S
      FAC_SAVE(IS) = FAC
      RETURN
      END

C*******************************************************************************

      FUNCTION FACI(IX)
      X = IX
      FACI = FAC(X)
      RETURN
      END

C*******************************************************************************

      FUNCTION FACIS(IX)
      X = IX
      FACIS = FACS(X)
      RETURN
      END

C*******************************************************************************

        FUNCTION BINOM(A,B)
        IF(A.EQ.B)GO TO 12
        IF(B.EQ.0.)GO TO 12
        IF(A-B.GT.B)GO TO 9
        M=A-B
        GO TO 11
9       M=B
11      BINOM=1.
        X=A-FLOAT(M)
        DO 10 I=1,M
        AI=I
10      BINOM=BINOM*(X+AI)/AI
        RETURN
12      BINOM=1.
        RETURN
        END

C*******************************************************************************

      FUNCTION DFAC(A)
      DIMENSION DFAC_SAVE(100)
      IS = A
      IF(A)6,7,8
6     DFAC=0.
      RETURN
7     DFAC=1.
      RETURN
8     N=A
      IF(IFIRST.EQ.0.AND.PARI(A).EQ.1.) WRITE(*,7777) A
7777  FORMAT(1X,'WARNING IN DFAC: INPUT = ',F10.5,
     1/,1X,'THE INPUT IS USUALLY ODD',
     1/,1X,'THIS MESSAGE WILL NOW BE SUPRESSED')
      IF(IFIRST.EQ.0) IFIRST = 999
CBAB      IF(PARI(A).EQ.1.) STOP
      DFAC = DFAC_SAVE(IS)
      IF(DFAC.NE.0.) RETURN
      DFAC=1.
CBAB      DO 9 I=1,N,2
      DO 9 I=N,1,-2
      AI=I
9     DFAC=AI*DFAC
      DFAC_SAVE(IS) = DFAC
      RETURN
      END

C*******************************************************************************

      FUNCTION PARI(A)
      I = A
      I = I/2
      P = A/2.0
      Q = I
      IF(P-Q) 1,2,1
1     PARI = -1.0
      RETURN
2     PARI = 1.0
      RETURN
      END

C*******************************************************************************

      FUNCTION IPOT(I)     ! sign factor : (-1)^I
      A = I
      I1 = I/2
      P = A/2.0
      Q = I1
      IF(P-Q) 1,2,1
1     IPOT = -1.0
      RETURN
2     IPOT = 1.0
      RETURN
      END

C*******************************************************************************

      FUNCTION IPARI(I)     ! sign factor : (-1)^I
      A = I
      I1 = I/2
      P = A/2.0
      Q = I1
      IF(P-Q) 1,2,1
1     IPARI = -1.0
      RETURN
2     IPARI = 1.0
      RETURN
      END

C*******************************************************************************

      FUNCTION FINDEX(XN,XL,XJ)
CBAB XN STARTS AT 1
      FINDEX = 0
      IF(XN.EQ.0.) RETURN
      IF(TR(XL,XJ,0.5).EQ.-1.) RETURN
      N = XN-1.
      L = XL
      J2 = 2.*XJ
      FINDEX = ((2*N+L)*(2*N+L+3)-(J2-1)+2)/2
      RETURN
      END

C*******************************************************************************

      SUBROUTINE INDEX
      CHARACTER*2 ZLABEL 
      CHARACTER*1 LABEL
      COMMON/CLABEL/LABEL(20)
      COMMON/ZLABEL/ZLABEL(200)
      COMMON/BINDEX/XL(3,200)
      I = 0
      DO 100 I1 = 1,1000
      N = I1-1
      LN = 0
      L = N-2*LN
200   IF(L.LT.0) GO TO 100
      I = I + 1
      IF(I.GT.200) GO TO 101
      XJ2 = 2*L+1
      XL(1,I) = LN+1
      XL(2,I) = L
      XL(3,I) = XJ2/2.
      IF(L.EQ.0) GO TO 201
      I = I + 1
      IF(I.GT.200) GO TO 101
      XJ2 = 2*L-1
      XL(1,I) = LN+1
      XL(2,I) = L
      XL(3,I) = XJ2/2.
201   L = L-2
      LN = LN + 1
      GO TO 200
100   CONTINUE
101   CONTINUE
      LABEL(1) = 's'
      LABEL(2) = 'p'
      LABEL(3) = 'd'
      LABEL(4) = 'f'
      LABEL(5) = 'g'
      LABEL(6) = 'h'
      LABEL(7) = 'i'
      LABEL(8) = 'j'
      LABEL(9) = 'k'
      LABEL(10) = 'l'
      LABEL(11) = 'm'
      LABEL(12) = 'n'
      LABEL(13) = 'o'
      LABEL(14) = 'p'
      LABEL(15) = 'q'
      LABEL(16) = 'r'
      LABEL(17) = 's'
      LABEL(18) = 't'
      LABEL(19) = 'u'
      LABEL(20) = 'v'
      ZLABEL(1) = 'H'
      ZLABEL(2) = 'He'
      ZLABEL(3) = 'Li'
      ZLABEL(4) = 'Be'
      ZLABEL(5) = 'B'
      ZLABEL(6) = 'C'
      ZLABEL(7) = 'N'
      ZLABEL(8) = 'O'
      ZLABEL(9) = 'F'
      ZLABEL(10) = 'Ne'
      ZLABEL(11) = 'Na'
      ZLABEL(12) = 'Mg'
      ZLABEL(13) = 'Al'
      ZLABEL(14) = 'Si'
      ZLABEL(15) = 'P'
      ZLABEL(16) = 'S'
      ZLABEL(17) = 'Cl'
      ZLABEL(18) = 'Ar'
      ZLABEL(19) = 'K'
      ZLABEL(20) = 'Ca'
      ZLABEL(21) = 'Sc'
      ZLABEL(22) = 'Ti'
      ZLABEL(23) = 'V'
      ZLABEL(24) = 'Cr'
      ZLABEL(25) = 'Mn'
      ZLABEL(26) = 'Fe'
      ZLABEL(27) = 'Co'
      ZLABEL(28) = 'Ni'
      ZLABEL(29) = 'Cu'
      ZLABEL(30) = 'Zn'
      ZLABEL(31) = 'Ga'
      ZLABEL(32) = 'Ge'
      ZLABEL(33) = 'As'
      ZLABEL(34) = 'Se'
      ZLABEL(35) = 'Br'
      ZLABEL(36) = 'Kr'
      ZLABEL(37) = 'Rb'
      ZLABEL(38) = 'Sr'
      ZLABEL(39) = 'Y'
      ZLABEL(40) = 'Zr'
      ZLABEL(41) = 'Nb'
      ZLABEL(42) = 'Mo'
      ZLABEL(43) = 'Tc'
      ZLABEL(44) = 'Ru'
      ZLABEL(45) = 'Rh'
      ZLABEL(46) = 'Pd'
      ZLABEL(47) = 'Ag'
      ZLABEL(48) = 'Cd'
      ZLABEL(49) = 'In'
      ZLABEL(50) = 'Sn'
      ZLABEL(51) = 'Sb'
      ZLABEL(52) = 'Te'
      ZLABEL(53) = 'I'
      ZLABEL(54) = 'Xe'
      ZLABEL(55) = 'Cs'
      ZLABEL(56) = 'Ba'
      ZLABEL(57) = 'La'
      ZLABEL(58) = 'Ce'
      ZLABEL(59) = 'Pr'
      ZLABEL(60) = 'Nd'
      ZLABEL(61) = 'Pm'
      ZLABEL(62) = 'Sm'
      ZLABEL(63) = 'Eu'
      ZLABEL(64) = 'Gd'
      ZLABEL(65) = 'Tb'
      ZLABEL(66) = 'Dy'
      ZLABEL(67) = 'Ho'
      ZLABEL(68) = 'Er'
      ZLABEL(69) = 'Tm'
      ZLABEL(70) = 'Yb'
      ZLABEL(71) = 'Lu'
      ZLABEL(72) = 'Hf'
      ZLABEL(73) = 'Ta'
      ZLABEL(74) = 'W'
      ZLABEL(75) = 'Re'
      ZLABEL(76) = 'Os'
      ZLABEL(77) = 'Ir'
      ZLABEL(78) = 'Pt'
      ZLABEL(79) = 'Au'
      ZLABEL(80) = 'Hg'
      ZLABEL(81) = 'Tl'
      ZLABEL(82) = 'Pb'
      ZLABEL(83) = 'Bi'
      ZLABEL(84) = 'Po'
      ZLABEL(85) = 'At'
      ZLABEL(86) = 'Rn'
      ZLABEL(87) = 'Fr'
      ZLABEL(88) = 'Ra'
      ZLABEL(89) = 'Ac'
      ZLABEL(90) = 'Th'
      ZLABEL(91) = 'Pa'
      ZLABEL(92) = 'U'
      ZLABEL(93) = 'Np'
      ZLABEL(94) = 'Pu'
      ZLABEL(95) = 'Am'
      ZLABEL(96) = 'Cm'
      ZLABEL(97) = 'Bk'
      ZLABEL(98) = 'Cf'
      ZLABEL(99) = 'Es'
      ZLABEL(100) = 'Fm'
      ZLABEL(101) = 'Md'
      ZLABEL(102) = 'No'
      ZLABEL(103) = 'Lr'
      ZLABEL(104) = 'X1'
      ZLABEL(105) = 'X2'
      ZLABEL(106) = 'X3'
      ZLABEL(107) = 'X4'
      ZLABEL(108) = 'X5'
      ZLABEL(109) = 'X6'

      END

C*******************************************************************************

      FUNCTION SJ(A1,A2,B2,B1,A3,B3)
      SJ = PARI(A1+A2+A3+B1)*RACAHR(A1,A2,A3,B1,B2,B3)
      RETURN
      END
CBAB ORIGINAL
CBAB      FUNCTION W(A1,A2,A3,B1,B2,B3)
CBAB      W = PARI(A1+A2+A3+B1)*SJ(A1,A2,B2,B1,A3,B3)
CBAB      RETURN
CBAB      END

C*******************************************************************************

      FUNCTION TJ(A1,A2,A3,B1,B2,B3)
      TJ = PARI(A1-A2-B3)*CLEBR(A1,B1,A2,B2,A3,-B3)/SQRT(2.*A3+1.)
      RETURN
      END
CBAB ORIGINAL
CBAB      FUNCTION CG(A1,B1,A2,B2,A3,B3)
CBAB      CG = PARI(A1-A2+B3)*SQRT(2.*A3+1.)*TJ(A1,A2,A3,B1,B2,-B3)
CBAB      RETURN
CBAB      END

C*******************************************************************************

      FUNCTION CG(X1,X2,X3,X4,X5,X6)
      CG = CLEBR(X1,X2,X3,X4,X5,X6)
      RETURN
      END
      FUNCTION W(X1,X2,X3,X4,X5,X6)
      W = RACAHR(X1,X2,X3,X4,X5,X6)
      RETURN
      END

C*******************************************************************************

      SUBROUTINE TH_FACINIT
C ...  SET UP LOG OF FACTORIALS
      PARAMETER (LFACTC=200)
      LOGICAL FIRST
      COMMON / LOGFAC / FIRST,FACLOG(LFACTC)
      COMMON / LOGFAC1 / LFACT
      DATA FIRST/.TRUE./,LFACT/LFACTC/
C      REAL*8 FACLOG,FN
      FIRST=.FALSE.
      FACLOG(1)=0.0
      FACLOG(2)=0.0
      FN=1.0
      DO 10 I=3,LFACTC
      FN=FN+1.0
      FACLOG(I)=FACLOG(I-1)+LOG(FN)
   10 CONTINUE
      RETURN
      END
C
C*******************************************************************************

C  ***********************  CLEB, CLEBI, CLEB2I ************************
C
      REAL FUNCTION CLEBR(A,B,C,D,E,F)
C
C      ARGUMENTS ARE REAL AND OF TRUE VALUE; J1,M1,J2,M2,J3,M3
C
      PARAMETER (LFACTC=200)
      COMMON / LOGFAC / FIRST,FACLOG(LFACTC)
      COMMON / LOGFAC1 / LFACT
      LOGICAL FIRST
C      REAL*8 FACLOG
CCCCCC      IA=2*J1,ID=2*M1 ETC.(J1 IS OF TRUE VALUE)
      IA=NINT(2.*A)
      IB=NINT(2.*C)
      IC=NINT(2.*E)
      ID=NINT(2.*B)
      IE=NINT(2.*D)
      IF=NINT(2.*F)
      GOTO 7000
C ...............  CLEBI  ...........................
C
C      ARGUMENTS ARE INTEGER AND OF TRUE VALUE
C
      ENTRY CLEBI(LL1,LM1,LL2,LM2,LL3,LM3)
      IA=2*LL1
      IB=2*LL2
      IC=2*LL3
      ID=2*LM1
      IE=2*LM2
      IF=2*LM3
      GOTO 7000
C ..............  CLEB  ..............................
C
C      ARGUMENTS ARE INTEGER AND REPRESENT TWICE THE REAL VALUE
C
      ENTRY CLEB(I2J1,I2M1,I2J2,I2M2,I2J3,I2M3)
      IA=I2J1
      IB=I2J2
      IC=I2J3
      ID=I2M1
      IE=I2M2
      IF=I2M3
 7000 IF (FIRST) CALL TH_FACINIT
      RAC=0.0
      IF(ID+IE-IF) 1000,105,1000
  105 K1=IA+IB+IC
      IF((-1)**K1) 1000,110,110
  110 K1=IA+IB-IC
      K2=IC+IA-IB
      K3=IB+IC-IA
      K4=IA-IABS (IB-IC)
      K5=IB-IABS (IC-IA)
      K6=IC-IABS (IA-IB)
      K7= MIN0 (K1,K2,K3,K4,K5,K6)
      IF(K7) 1000,120,120
  120 IF((-1)**(IA+ID)) 1000,1000,130
  130 IF((-1)**(IB+IE)) 1000,1000,140
  140 IF((-1)**(IC+IF)) 1000,1000,150
  150 IF(IA-IABS (ID)) 1000,152,152
  152 IF(IB-IABS (IE)) 1000,154,154
  154 IF(IC-IABS (IF)) 1000,160,160
  160 SIGNFC=1.0
      IAM=IA
      IBM=IB
      ICM=IC
      IDM=ID
      IEM=IE
      IFM=IF
      IF(IA-IB) 210,220,220
  210 IF(IA-IC) 215,225,225
  215 IT=IA
      IA=IB
      IB=IT
      IT=ID
      ID=IE
      IE=IT
      SIGNFC=(-1.0)**((IA+IB-IC)/2)
      GO TO 235
  220 IF(IC-IB) 225,235,235
  225 IT=IC
      IC=IB
      IB=IT
      IT=IF
      IF=-IE
      IE=-IT
      FIBM=IBM+1
      FICM=ICM+1
      SIGNFC=(-1.)**((IAM-IDM)/2)*SQRT (FICM/FIBM)
  235 IF(IB) 237,236,237
  236 RAC=SIGNFC
      GO TO 900
  237 IF(IE) 250,250,240
  240 SIGNFC=SIGNFC*((-1.0)**((IA+IB-IC)/2))
      ID=-ID
      IE=-IE
      IF=-IF
  250 FC2=IC+1
      IABCP=(IA+IB+IC)/2+1
      IABC=IABCP-IC
      ICAB=IABCP-IB
      IBCA=IABCP-IA
      IAPD=(IA+ID)/2+1
      IAMD=IAPD-ID
      IBPE=(IB+IE)/2+1
      IBME=IBPE-IE
      ICPF=(IC+IF)/2+1
      ICMF=ICPF-IF
      SQFCLG=0.5*(ALOG(FC2)-FACLOG(IABCP+1)
     1      +FACLOG(IABC)+FACLOG(ICAB)+FACLOG(IBCA)
     2      +FACLOG(IAPD)+FACLOG(IAMD)+FACLOG(IBPE)
     3      +FACLOG(IBME)+FACLOG(ICPF)+FACLOG(ICMF))
      NZMIC2=(IB-IC-ID)/2
      NZMIC3=(IA-IC+IE)/2
      NZMI= MAX0 (0,NZMIC2,NZMIC3)+1
      NZMX= MIN0 (IABC,IAMD,IBPE)
      IF(NZMI-NZMX) 310,310,900
  310 SS=0.0
      S1=(-1.0)**(NZMI-1)
      DO 400 NZ=NZMI,NZMX
      NZM1=NZ-1
      NZT1=IABC-NZM1
      NZT2=IAMD-NZM1
      NZT3=IBPE-NZM1
      NZT4=NZ-NZMIC2
      NZT5=NZ-NZMIC3
      TERMLG=SQFCLG-FACLOG(NZ)-FACLOG(NZT1)-FACLOG(NZT2)
     1           -FACLOG(NZT3)-FACLOG(NZT4)-FACLOG(NZT5)
      SSTERM=S1*EXP (TERMLG)
      SS=SS+SSTERM
  400 S1=-S1
      RAC=SIGNFC*SS
  900 IA=IAM
      IB=IBM
      IC=ICM
      ID=IDM
      IE=IEM
      IF=IFM
 1000 CLEB=RAC
      RETURN
      END

C*******************************************************************************

      FUNCTION RACAH(JAD,JBD,JCD,JDD,JED,JFD)
C
C        CALCULATES RACAH COEFFICIENTS
C
C        SOURCE : UNKNOWN
C        MODIFIED : MARCH 1982 , OLAF SCHOLTEN
C              RUN TIME OPTIMIZED FOR VAX780 MACHINE
C
C        ENTRIES : RACAH , RACAHI , RACAHR
C            RACAH  : INTEGER ARGUMENTS = 2*J
C            RACAHI : ARGUMENTS = TRUE INTEGER VALUE
C            RACAHR : ARGUMENTS = TRUE REAL VALUE
C        EXTERNAL : TH_FACINIT , GENERATES FACTORIAL TABLE
C
      PARAMETER (LFACTC=200)
      DIMENSION I(16)
      LOGICAL FIRST
C      REAL*8 G,S
      COMMON / LOGFAC / FIRST,FACLOG(LFACTC)
      COMMON / LOGFAC1 / LFACT
      COMMON / LOGFAC2 / G(1)
      EQUIVALENCE(I(1),I1),(I(2),I2),(I(3),I3),(I(4),I4),(I(5),I5),
     1 (I(6),I6),(I(7),I7),(I(8),I8),(I(9),I9),(I(10),I10),(I(11),I11),
     2 (I(12),I12),(I(13),I13),(I(14),I14),(I(15),I15),(I(16),I16)
C        MAKE USEFULL COMBINATIONS
      K=JAD+JBD-JED+2
      I1=K/2
      IF((2*I1).NE.K) GOTO 300
      K=JCD+JDD-JED+2
      I4=K/2
      IF((2*I4).NE.K) GOTO 300
      K=JAD+JCD-JFD+2
      I7=K/2
      IF((2*I7).NE.K) GOTO 300
      K=JBD+JDD-JFD+2
      I10=K/2
      IF((2*I10).NE.K) GOTO 300
      I13=I1+JED
      I14=I4+JED
      I15=I7+JFD
      I16=I10+JFD
      I2=I13-JAD
      I3=I13-JBD
      I5=I14-JCD
      I6=I14-JDD
      I8=I15-JAD
      I9=I15-JCD
      I11=I16-JBD
      I12=I16-JDD
C       CHECK TRIANGULAR INEQUALITIES,FIND NO. OF TERMS IN SUM
      N=MIN(I1,I2,I3,I4,I5,I6,I7,I8,I9,I10,I11,I12)-1
      IF(N) 300,2,2
C       FIND MINIMUM VALUE OF SUMMATION INDEX
    2 IL=MAX(I13,I14,I15,I16)
      IF(MIN(JAD,JBD,JCD,JDD,JED,JFD)) 300,20,1
C    ..............
      ENTRY RACAHI(JA1,JB1,JC1,JD1,JE1,JF1)
C        MAKE USEFULL COMBINATIONS
      I13=JA1+JB1+JE1+1
      I14=JC1+JD1+JE1+1
      I15=JA1+JC1+JF1+1
      I16=JB1+JD1+JF1+1
      I1=I13-JE1*2
      I2=I13-JA1*2
      I3=I13-JB1*2
      I4=I14-JE1*2
      I5=I14-JC1*2
      I6=I14-JD1*2
      I7=I15-JF1*2
      I8=I15-JA1*2
      I9=I15-JC1*2
      I10=I16-JF1*2
      I11=I16-JB1*2
      I12=I16-JD1*2
C       CHECK TRIANGULAR INEQUALITIES,FIND NO. OF TERMS IN SUM
      N=MIN(I1,I2,I3,I4,I5,I6,I7,I8,I9,I10,I11,I12)-1
      IF(N) 300,4,4
C       FIND MINIMUM VALUE OF SUMMATION INDEX
    4 IL=MAX(I13,I14,I15,I16)
      LMIN=MIN(JA1,JB1,JC1,JD1,JE1,JF1)
      IF(LMIN)300,20,1
C     ............
      ENTRY RACAHR(A,B,C,D,E,F)
C     CONVERT ARGUMENTS TO INTEGER 
      JA=NINT(2.*A)
      JB=NINT(2.*B)
      JC=NINT(2.*C)
      JD=NINT(2.*D)
      JE=NINT(2.*E)
      JF=NINT(2.*F)
C        MAKE USEFULL COMBINATIONS
      K=JA+JB-JE+2
      I1=K/2
      IF((2*I1-K).NE.0) GOTO 300
      K=JC+JD-JE+2
      I4=K/2
      IF((2*I4-K).NE.0) GOTO 300
      K=JA+JC-JF+2
      I7=K/2
      IF((2*I7-K).NE.0) GOTO 300
      K=JB+JD-JF+2
      I10=K/2
      IF((2*I10-K).NE.0) GOTO 300
      I13=I1+JE
      I14=I4+JE
      I15=I7+JF
      I16=I10+JF
      I2=I13-JA
      I3=I13-JB
      I5=I14-JC
      I6=I14-JD
      I8=I15-JA
      I9=I15-JC
      I11=I16-JB
      I12=I16-JD
C       CHECK TRIANGULAR INEQUALITIES,FIND NO. OF TERMS IN SUM
      N=MIN(I1,I2,I3,I4,I5,I6,I7,I8,I9,I10,I11,I12)-1
      IF(N) 300,3,3
C       FIND MINIMUM VALUE OF SUMMATION INDEX
    3 IL=MAX(I13,I14,I15,I16)
      LMIN=MIN(JA,JB,JC,JD,JE,JF)
      IF(LMIN)300,20,1
C      ------------
    1 IF(FIRST) CALL TH_FACINIT
      IF(IL.GE.LFACT) STOP 'RACAH: LENGTH FACTORIAL TABLE INSUFFICIENT'
      J1=IL-I13+1 
      J2=IL-I14+1 
      J3=IL-I15+1 
      J4=IL-I16+1
      J5=I13+I4-IL 
      J6=I15+I5-IL 
      J7=I16+I6-IL
      PH=1.
      IF(2*(J5/2).EQ.J5) PH=-1.
      H=PH*EXP ((G(I1)+G(I2)+G(I3)-G(I13+1)+G(I4)+G(I5)+G(I6)-
     1G(I14+1)+G(I7)+G(I8)+G(I9)-G(I15+1)+G(I10)+G(I11)+G(I12)-G(I16+1))
     2*.5+G(IL+1)-G(J1)-G(J2)-G(J3)-G(J4)-G(J5)-G(J6)-G(J7))
      IF(N)300,110,120
C
  110 RACAH=H 
      RETURN
C
  120 S=1.
      K=N-1
      KL=IL+1
      J5=J5-1
      J6=J6-1
      J7=J7-1
      DO 130 J=1,N   ! K=N-J
      S=1.-((KL+K)*(J5-K)*(J6-K)*(J7-K))*S/((J1+K)*(J2+K)*(J3+K)*(J4+K))
      K=K-1
  130 CONTINUE  
      RACAH=H*S
      RETURN
C
C      ONE OF THE ARGUMENTS =0
   20 IAD=IL
      IBD=IL
      DO 21 J=13,16
      IF(IAD.LT.I(J)) GOTO 22
      IF(IAD.LT.IBD) IBD=IAD
      IAD=I(J)
      GOTO 21
   22 IF(IBD.GT.I(J)) IBD=I(J)
   21 CONTINUE
      J5=I13+I4-IL 
      PH=1.
      IF(2*(J5/2).EQ.J5) PH=-1.
      RACAH=PH/SQRT(FLOAT(IAD*IBD))
      RETURN
C
C      IMPOSSIBLE COMBINATION OF ARGUMENTS
  300 RACAH=0.
      RETURN
      END
C
C*******************************************************************************

C     C O E F 9 J    
C
      FUNCTION COEF9J(J1,J2,J3,J4,J5,J6,J7,J8,J9)
C
CCCC  TAMURAS NUMBERING CONVENTION IS USED HERE.
CCCCC THE ARGUMENTS OF COEF9J, WHEN NUMBERED SEQUENTIALLY 1 THROUGH 9,
CCCC    CORRESPOND TO THE ARRAY
CCCC                               1  2  5
CCCC                               3  4  6
CCCC                               7  8  9
C
      DIMENSION LT(9)
CCCCCC
CCCCCC      ALL L9 MUST BE TWICE AS LARGE AS TRUE ARGUMENTS
CCCCCC
C
C               CHANGED FOR THEORY LIBRARY 4/8/82   HK
C
      U9=0.0
      LT(1)=J1
      LT(2)=J2
      LT(3)=J3
      LT(4)=J4
      LT(5)=J5
      LT(6)=J6
      LT(7)=J7
      LT(8)=J8
      LT(9)=J9
      LMIN=LT(1)
      IMIN=1
      DO 20 I=2,9
      IF(LT(I)-LMIN) 15,20,20
   15 LMIN=LT(I)
      IMIN=I
   20 CONTINUE
      KEX=0
      GO TO (110,110,110,110,150,150,170,170,190),IMIN
  110 MM=(IMIN-1)/2+1
      M1=MM+MM-1
      M2=M1+1
      M3=MM+4
      L1=LT(7)
      LT(7)=LT(M1)
      LT(M1)=L1
      L1=LT(8)
      LT(8)=LT(M2)
      LT(M2)=L1
      L1=LT(9)
      LT(9)=LT(M3)
      LT(M3)=L1
      IMIN=IMIN+(7-M1)
      GO TO 175
  150 KEX=1
      M1=7
      M2=8
      M3=IMIN+IMIN-9
      M4=M3+1
      GO TO 180
  170 KEX=1
  175 M1=5
      M2=6
      M3=IMIN-6
      M4=M3+2
  180 L1=LT(M1)
      L1=LT(M1)
      LT(M1)=LT(M3)
      LT(M3)=L1
      L1=LT(M2)
      LT(M2)=LT(M4)
      LT(M4)=L1
      L1=LT(9)
      LT(9)=LT(IMIN)
      LT(IMIN)=L1
  190 IF(LT(9)) 200,200,300
  200 IF(LT(5)-LT(6)) 1000,210,1000
  210 IF(LT(7)-LT(8)) 1000,220,1000
  220 RT=(LT(5)+1)*(LT(7)+1)
      K=(LT(5)+LT(7)-LT(1)-LT(4))/2
      RAC= RACAH(LT(1),LT(2),LT(3),LT(4),LT(5),LT(7))
      PH=1.
      IF (2*(K/2) .NE. K) PH=-1.
      U9=(RAC/SQRT(RT))*PH
      GO TO 370
  300 K1=IABS(LT(2)-LT(7))
      K2=IABS(LT(3)-LT(5))
      K3=IABS(LT(4)-LT(9))
      NMIN=MAX0(K1,K2,K3)
      K1=LT(2)+LT(7)
      K2=LT(3)+LT(5)
      K3=LT(4)+LT(9)
      NMAX=MIN0(K1,K2,K3)
      IF (NMIN-NMAX) 320, 320, 1000
  320 DO 350 N=NMIN,NMAX,2
      W1=N+1
      RAC= RACAH(LT(2),LT(5),LT(7),LT(3),LT(1),N)
      IF (RAC) 321, 350, 321
  321 W1=W1*RAC
      RAC= RACAH(LT(2),LT(4),LT(7),LT(9),LT(8),N)
      IF (RAC) 322, 350, 322
  322 W1=W1*RAC
      RAC= RACAH(LT(3),LT(4),LT(5),LT(9),LT(6),N)
      IF (RAC) 323, 350, 323
  323 U9=U9+W1*RAC
  350 CONTINUE
  370 IF(KEX) 400,1000,400
  400 KP=0
      DO 410 I=1,9
  410 KP=KP+LT(I)
      K=KP/2
      PH=1.
      IF (2*(K/2) .NE. K) PH=-1.
      U9=U9*PH
 1000 COEF9J=U9
      RETURN
      END

      SUBROUTINE SIXJ(XJ1,XJ2,XJ3,XL1,XL2,XL3,C6J)                      TRO01770
C     IMPLICIT REAL*4 (A-H,O-Z)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION FLOG(301)                                               TRO01780
C                                                                       TRO01790
C     MISE EN DATA DES LOG(FACTORIELLE) ET DE EPS                       TRO01800
C                                                                       TRO01810
      DATA(FLOG(I),I=2,31)/0.D0,.69314718D0,1.7917595D0,3.1780538D0,
     A4.7874917D0,6.
     15792511D0,8.5251613D0,10.604603D0,12.801827D0,15.104413D0,
     B17.502307D0,19.98721
     24D0,22.552163D0,25.191221D0,27.899271D0,30.671860D0,33.505072D0,
     C36.395445D0,39.3
     339884D0,42.335616D0,45.380139D0,48.471180D0,51.606674D0,54.
     D784729D0,58.003604D0,
     461.261702D0,64.557537D0,67.889743D0,71.257038D0,74.658235D0/   
      DATA(FLOG(I),I=32,61)/78.092223D0,81.557959D0,85.054466D0,88.
     A580827D0,92.1
     136175D0,95.719694D0,99.330612D0,102.96820D0,106.63176D0,110.
     B32064D0,114.03421D0,
     2117.77188D0,121.53308D0,125.31727D0,129.12393D0,132.95257D0,
     C136.80272D0,140.67
     3392D0,144.56574D0,148.47776D0,152.40959D0,156.36083D0,160.3311
     D2D0,164.32011D0,16
     48.32744D0,172.35279D0,176.39584D0,180.45629D0,184.53383D0,
     E188.62817D0/        
      DATA(FLOG(I),I=62,91)/192.73904D0,196.86618D0,201.00931D0,
     A205.16820D0,209.
     134258D0,213.53224D0,217.73693D0,221.95644D0,226.19054D0,
     B230.43904D0,234.70172D0,
     2238.97839D0,243.26885D0,247.57291D0,251.89040D0,256.22113D0,
     C260.56494D0,264.92
     3164D0,269.29110D0,273.67312D0,278.06757D0,282.47429D0,286.
     D89313D0,291.32394D0,29
     45.76659D0,300.22094D0,304.68685D0,309.16419D0,313.65283D0,
     E318.15264D0/        
      DATA(FLOG(I),I=92,121)/322.66349D0,327.18529D0,331.71788D0,
     A336.26118D0,340
     1.81505D0,345.37940D0,349.95411D0,354.53908D0,359.13420D0,363.
     B73937D0,368.35449
     2D0,372.97946D0,377.61419D0,382.25859D0,386.91255D0,391.57598D0,
     C396.24881D0,400.9
     33094D0,405.62230D0,410.32277D0,415.03230D0,419.75080D0,424.47819
     DD0,429.21439D0,4
     433.95932D0,438.71291D0,443.47508D0,448.24576D0,453.02489D0,
     E457.81238D0/       
      DATA(FLOG(I),I=122,151)/462.60817D0,467.41220D0,472.22438D0,477.
     A04466D0,48
     11.87298D0,486.70926D0,491.55345D0,496.40547D0,501.26529D0,506.
     B13282D0,511.0080
     22D0,515.89082D0,520.78117D0,525.67901D0,530.58428D0,535.49694D0,
     C540.41692D0,545.
     334417D0,550.27865D0,555.22029D0,560.16905D0,565.12488D0,570.08772
     DD0,575.05753D0,
     4580.03427D0,585.01787D0,590.00830D0,595.00552D0,600.00946D0,
     E605.02010D0/      
      DATA(FLOG(I),I=152,181)/610.03738D0,615.06126D0,620.09170D0,
     A625.12866D0,63
     10.17208D0,635.22193D0,640.27818D0,645.34077D0,650.40968D0,655.
     B48486D0,660.5662
     26D0,665.65385D0,670.74760D0,675.84747D0,680.95341D0,686.06541D0,
     C691.18340D0,696.
     330735D0,701.43726D0,706.57306D0,711.71472D0,716.86221D0,722.
     D01551D0,727.17456D0,
     4732.33934D0,737.50983D0,742.68598D0,747.86776D0,753.05516D0,
     E758.24811D0/      
      DATA(FLOG(I),I=182,211)/763.44661D0,768.65061D0,773.86010D0,
     A779.07503D0,78
     14.29539D0,789.52114D0,794.75224D0,799.98869D0,805.23044D0,810.
     B47747D0,815.7297
     23D0,820.98722D0,826.24991D0,831.51778D0,836.79078D0,842.06890D0,
     C847.35209D0,852.
     364036D0,857.93366D0,863.23199D0,868.53529D0,873.84356D0,879.
     D15676D0,884.47488D0,
     4889.79789D0,895.12577D0,900.45848D0,905.79603D0,911.13836D0,
     E916.48547D0/      
      DATA(FLOG(I),I=212,241)/921.83732D0,927.19391D0,932.55521D0,
     A937.92118D0,94
     13.29181D0,948.66710D0,954.04699D0,959.43148D0,964.82056D0,970.
     B21419D0,975.6123
     25D0,981.01503D0,986.42220D0,991.83385D0,997.24995D0,1002.6705D0,
     C1008.0954D0,1013
     3.5248D0,1018.9585D0,1024.3966D0,1029.8389D0,1035.2857D0,1040.
     D7367D0,1046.1920D0,
     41051.6516D0,1057.1155D0,1062.5836D0,1068.0558D0,1073.5323D0,
     E1079.0129D0/      
      DATA(FLOG(I),I=242,271)/1084.4977D0,1089.9866D0,1095.4797D0,1100.
     A9768D0,11
     106.4781D0,1111.9834D0,1117.4928D0,1123.0063D0,1128.5237D0,1134.
     B0452D0,1139.570
     26D0,1145.1001D0,1150.6335D0,1156.1708D0,1161.7120D0,1167.2573D0,
     C1172.8063D0,1178
     3.3593D0,1183.9161D0,1189.4768D0,1195.0413D0,1200.6097D0,1206.
     D1818D0,1211.7577D0,
     41217.3375D0,1222.9209D0,1228.5082D0,1234.0992D0,1239.6939D0,
     E1245.2924D0/      
      DATA(FLOG(I),I=272,301)/1250.8944D0,1256.5003D0,1262.1097D0,
     A1267.7228D0,12
     173.3396D0,1278.9600D0,1284.5840D0,1290.2117D0,1295.8429D0,
     B1301.4777D0,1307.116
     20D0,1312.7580D0,1318.4034D0,1324.0524D0,1329.7048D0,1335.3609D0,
     C1341.0203D0,1346
     3.6833D0,1352.3497D0,1358.0196D0,1363.6929D0,1369.3697D0,
     D1375.0499D0,1380.7334D0,
     41386.4204D0,1392.1107D0,1397.8045D0,1403.5016D0,1409.2020D0,
     E1414.9058D0/      
      DATA EPS1,EPS2/.1D0,-.2D0/                                            TRO02320
C                                                                       TRO02330
C     CALCUL DES COMBINAISONS J,L                                       TRO02340
C                                                                       TRO02350
      XN=-XJ1+XJ2+XJ3+EPS1                                              TRO02360
      IF (XN.LT.0.D0) XN=XN+EPS2                                          TRO02370
      N1=XN                                                             TRO02380
      XN=-XL1+XL2+XJ3+EPS1                                              TRO02390
      IF (XN.LT.0.D0) XN=XN+EPS2                                          TRO02400
      N2=XN                                                             TRO02410
      XN=-XL1+XJ2+XL3+EPS1                                              TRO02420
      IF (XN.LT.0.D0) XN=XN+EPS2                                          TRO02430
      N3=XN                                                             TRO02440
      XN=-XJ1+XL2+XL3+EPS1                                              TRO02450
      IF (XN.LT.0.D0) XN=XN+EPS2                                          TRO02460
      N4=XN                                                             TRO02470
      XN=XJ1-XJ2+XJ3+EPS1                                               TRO02480
      IF (XN.LT.0.D0) XN=XN+EPS2                                          TRO02490
      N5=XN                                                             TRO02500
      XN=XL1-XL2+XJ3+EPS1                                               TRO02510
      IF (XN.LT.0.D0) XN=XN+EPS2                                          TRO02520
      N6=XN                                                             TRO02530
      XN=XL1-XJ2+XL3+EPS1                                               TRO02540
      IF (XN.LT.0.D0) XN=XN+EPS2                                          TRO02550
      N7=XN                                                             TRO02560
      XN=XJ1-XL2+XL3+EPS1                                               TRO02570
      IF (XN.LT.0.D0) XN=XN+EPS2                                          TRO02580
      N8=XN                                                             TRO02590
      XN=XJ1+XJ2-XJ3+EPS1                                               TRO02600
      IF (XN.LT.0.D0) XN=XN+EPS2                                          TRO02610
      N9=XN                                                             TRO02620
      XN=XL1+XL2-XJ3+EPS1                                               TRO02630
      IF (XN.LT.0.D0) XN=XN+EPS2                                          TRO02640
      N10=XN                                                            TRO02650
      XN=XL1+XJ2-XL3+EPS1                                               TRO02660
      IF (XN.LT.0.D0) XN=XN+EPS2                                          TRO02670
      N11=XN                                                            TRO02680
      XN=XJ1+XL2-XL3+EPS1                                               TRO02690
      IF (XN.LT.0.D0) XN=XN+EPS2                                          TRO02700
      N12=XN                                                            TRO02710
      XN=-XJ1-XL1+XJ3+XL3+EPS1                                          TRO02720
      IF (XN.LT.0.D0) XN=XN+EPS2                                          TRO02730
      N13=XN                                                            TRO02740
      XN=-XJ2-XL2+XJ3+XL3+EPS1                                          TRO02750
      IF (XN.LT.0.D0) XN=XN+EPS2                                          TRO02760
      N14=XN                                                            TRO02770
      XN=XJ1+XL1+XJ2+XL2+EPS1                                           TRO02780
      IF (XN.LT.0.D0) XN=XN+EPS2                                          TRO02790
      N15=XN                                                            TRO02800
      N15=N15+1                                                         TRO02810
      XN=XJ1+XJ2+XJ3+EPS1                                               TRO02820
      IF (XN.LT.0.D0) XN=XN+EPS2                                          TRO02830
      N16=XN                                                            TRO02840
      N16=N16+1                                                         TRO02850
      XN=XL1+XL2+XJ3+EPS1                                               TRO02860
      IF (XN.LT.0.D0) XN=XN+EPS2                                          TRO02870
      N17=XN                                                            TRO02880
      N17=N17+1                                                         TRO02890
      XN=XL1+XJ2+XL3+EPS1                                               TRO02900
      IF (XN.LT.0.D0) XN=XN+EPS2                                          TRO02910
      N18=XN                                                            TRO02920
      N18=N18+1                                                         TRO02930
      XN=XJ1+XL2+XL3+EPS1                                               TRO02940
      IF (XN.LT.0.D0) XN=XN+EPS2                                          TRO02950
      N19=XN                                                            TRO02960
      N19=N19+1                                                         TRO02970
C                                                                       TRO02980
C     TEST SUR LES J ET L                                               TRO02990
C                                                                       TRO03000
      IF(N9.LT.0) GO TO 50                                              TRO03010
      IF(N5.LT.0) GO TO 50                                              TRO03020
      IF(N1.LT.0) GO TO 50                                              TRO03030
      IF(N10.LT.0) GO TO 50                                             TRO03040
      IF(N6.LT.0) GO TO 50                                              TRO03050
      IF(N2.LT.0) GO TO 50                                              TRO03060
      IF(N11.LT.0) GO TO 50                                             TRO03070
      IF(N7.LT.0) GO TO 50                                              TRO03080
      IF(N3.LT.0) GO TO 50                                              TRO03090
      IF(N12.LT.0) GO TO 50                                             TRO03100
      IF(N8.LT.0) GO TO 50                                              TRO03110
      IF(N4.LT.0) GO TO 50                                              TRO03120
      K=N17+N18+N19-3                                                   TRO03130
      IF(K.GT.0) GO TO 54                                               TRO03140
      C6J=1.D0                                                            TRO03150
      RETURN                                                            TRO03160
   50 C6J=0.D0                                                            TRO03170
      RETURN                                                            TRO03180
C                                                                       TRO03190
C      CALCUL DE LA SOMME ALTERNEE                                      TRO03200
C                                                                       TRO03210
   54 K=0                                                               TRO03220
      L=-N13                                                            TRO03230
      IF(L.GT.K) K=L                                                    TRO03240
      L=-N14                                                            TRO03250
      IF(L.GT.K) K=L                                                    TRO03260
      L=N9                                                              TRO03270
      IF(N10.LT.L) L=N10                                                TRO03280
      IF(N11.LT.L) L=N11                                                TRO03290
      IF(N12.LT.L) L=N12                                                TRO03300
      IF(N15.LT.L) L=N15                                                TRO03310
      F=1.D0                                                              TRO03320
      S=1.D0                                                              TRO03330
      I=K+1                                                             TRO03340
   62 IF(I.GT.L) GO TO 80                                               TRO03350
      IM1=I-1                                                           TRO03360
      NN=(N9-IM1)*(N10-IM1)*(N11-IM1)*(N12-IM1)                         TRO03370
      ND=I*(N13+I)*(N14+I)*(N15-IM1)                                    TRO03380
      F=-F*DFLOAT(NN)/DFLOAT(ND)                                          TRO03390
      S=S+F                                                             TRO03400
      I=I+1                                                             TRO03410
      GO TO 62                                                          TRO03420
C                                                                       TRO03430
C     CALCUL DE LA RACINE                                               TRO03440
C                                                                       TRO03450
   80 C2N=FLOG(N1+1)+FLOG(N2+1)+FLOG(N3+1)+FLOG(N4+1)+FLOG(N5+1)+FLOG(N6TRO03460
     1+1)+FLOG(N7+1)+FLOG(N8+1)+FLOG(N9+1)+FLOG(N10+1)+FLOG(N11+1)+FLOG(TRO03470
     2N12+1)                                                            TRO03480
      C2N=.5D0*C2N                                                        TRO03490
      C2D=FLOG(N16+1)+FLOG(N17+1)+FLOG(N18+1)+FLOG(N19+1)               TRO03500
      C2D=.5D0*C2D                                                        TRO03510
      KM1=K-1                                                           TRO03520
      KP1=K+1                                                           TRO03530
      C2N=C2N+FLOG(N15-KM1)                                             TRO03540
      C2D=C2D+FLOG(KP1)+FLOG(N13+KP1)+FLOG(N14+KP1)+FLOG(N9-KM1)+FLOG(N1TRO03550
     10-KM1)+FLOG(N11-KM1)+FLOG(N12-KM1)                                TRO03560
C                                                                       TRO03570
C     CALCUL DU C6J SANS PHASE                                          TRO03580
C                                                                       TRO03590
      F=C2D-C2N                                                         TRO03600
      IF(F.GT.80.D0) GO TO 98                                             TRO03610
      F=C2N/C2D                                                         TRO03620
      IF((F.LT.1.01D0).AND.(F.GT.0.98D0)) GO TO 98                          TRO03630
      C6J=S*DEXP(C2N-C2D)                                                TRO03640
      GO TO 106                                                         TRO03650
   98 IF(S) 100,50,102                                                  TRO03660
  100 S=DLOG(-S)                                                        TRO03670
      C6J=-DEXP(S+C2N-C2D)                                               TRO03680
      GO TO 106                                                         TRO03690
  102 S=DLOG(S)                                                         TRO03700
      C6J=DEXP(S+C2N-C2D)                                                TRO03710
C                                                                       TRO03720
C     CALCUL DE LA PHASE                                                TRO03730
C                                                                       TRO03740
  106 L=N15+KM1                                                         TRO03750
      K=L/2                                                             TRO03760
      K=2*K                                                             TRO03770
      IF(L.NE.K) C6J=-C6J                                               TRO03780
      RETURN                                                            TRO03790
      END                                                               TRO03800

      SUBROUTINE NEUFJ(XJ11,XJ12,XJ13,XJ21,XJ22,XJ23,XJ31,XJ32,XJ33,C9J)TRO03820
C     IMPLICIT REAL*4 (A-H,O-Z)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION FLOG(301)                                               TRO03830
C                                                                       TRO03840
C     MISE EN DATA DES LOG(FACTORIELLE) ET DE EPS                       TRO03850
C                                                                       TRO03860
      DATA(FLOG(I),I=2,31)/0.D0,.69314718D0,1.7917595D0,3.1780538D0,
     A4.7874917D0,6.
     15792511D0,8.5251613D0,10.604603D0,12.801827D0,15.104413D0,
     B17.502307D0,19.98721
     24D0,22.552163D0,25.191221D0,27.899271D0,30.671860D0,33.505072D0,
     C36.395445D0,39.3
     339884D0,42.335616D0,45.380139D0,48.471180D0,51.606674D0,54.
     D784729D0,58.003604D0,
     461.261702D0,64.557537D0,67.889743D0,71.257038D0,74.658235D0/   
      DATA(FLOG(I),I=32,61)/78.092223D0,81.557959D0,85.054466D0,88.
     A580827D0,92.1
     136175D0,95.719694D0,99.330612D0,102.96820D0,106.63176D0,110.
     B32064D0,114.03421D0,
     2117.77188D0,121.53308D0,125.31727D0,129.12393D0,132.95257D0,
     C136.80272D0,140.67
     3392D0,144.56574D0,148.47776D0,152.40959D0,156.36083D0,160.3311
     D2D0,164.32011D0,16
     48.32744D0,172.35279D0,176.39584D0,180.45629D0,184.53383D0,
     E188.62817D0/        
      DATA(FLOG(I),I=62,91)/192.73904D0,196.86618D0,201.00931D0,
     A205.16820D0,209.
     134258D0,213.53224D0,217.73693D0,221.95644D0,226.19054D0,
     B230.43904D0,234.70172D0,
     2238.97839D0,243.26885D0,247.57291D0,251.89040D0,256.22113D0,
     C260.56494D0,264.92
     3164D0,269.29110D0,273.67312D0,278.06757D0,282.47429D0,286.
     D89313D0,291.32394D0,29
     45.76659D0,300.22094D0,304.68685D0,309.16419D0,313.65283D0,
     E318.15264D0/        
      DATA(FLOG(I),I=92,121)/322.66349D0,327.18529D0,331.71788D0,
     A336.26118D0,340
     1.81505D0,345.37940D0,349.95411D0,354.53908D0,359.13420D0,363.
     B73937D0,368.35449
     2D0,372.97946D0,377.61419D0,382.25859D0,386.91255D0,391.57598D0,
     C396.24881D0,400.9
     33094D0,405.62230D0,410.32277D0,415.03230D0,419.75080D0,424.47819
     DD0,429.21439D0,4
     433.95932D0,438.71291D0,443.47508D0,448.24576D0,453.02489D0,
     E457.81238D0/       
      DATA(FLOG(I),I=122,151)/462.60817D0,467.41220D0,472.22438D0,477.
     A04466D0,48
     11.87298D0,486.70926D0,491.55345D0,496.40547D0,501.26529D0,506.
     B13282D0,511.0080
     22D0,515.89082D0,520.78117D0,525.67901D0,530.58428D0,535.49694D0,
     C540.41692D0,545.
     334417D0,550.27865D0,555.22029D0,560.16905D0,565.12488D0,570.08772
     DD0,575.05753D0,
     4580.03427D0,585.01787D0,590.00830D0,595.00552D0,600.00946D0,
     E605.02010D0/      
      DATA(FLOG(I),I=152,181)/610.03738D0,615.06126D0,620.09170D0,
     A625.12866D0,63
     10.17208D0,635.22193D0,640.27818D0,645.34077D0,650.40968D0,655.
     B48486D0,660.5662
     26D0,665.65385D0,670.74760D0,675.84747D0,680.95341D0,686.06541D0,
     C691.18340D0,696.
     330735D0,701.43726D0,706.57306D0,711.71472D0,716.86221D0,722.
     D01551D0,727.17456D0,
     4732.33934D0,737.50983D0,742.68598D0,747.86776D0,753.05516D0,
     E758.24811D0/      
      DATA(FLOG(I),I=182,211)/763.44661D0,768.65061D0,773.86010D0,
     A779.07503D0,78
     14.29539D0,789.52114D0,794.75224D0,799.98869D0,805.23044D0,810.
     B47747D0,815.7297
     23D0,820.98722D0,826.24991D0,831.51778D0,836.79078D0,842.06890D0,
     C847.35209D0,852.
     364036D0,857.93366D0,863.23199D0,868.53529D0,873.84356D0,879.
     D15676D0,884.47488D0,
     4889.79789D0,895.12577D0,900.45848D0,905.79603D0,911.13836D0,
     E916.48547D0/      
      DATA(FLOG(I),I=212,241)/921.83732D0,927.19391D0,932.55521D0,
     A937.92118D0,94
     13.29181D0,948.66710D0,954.04699D0,959.43148D0,964.82056D0,970.
     B21419D0,975.6123
     25D0,981.01503D0,986.42220D0,991.83385D0,997.24995D0,1002.6705D0,
     C1008.0954D0,1013
     3.5248D0,1018.9585D0,1024.3966D0,1029.8389D0,1035.2857D0,1040.
     D7367D0,1046.1920D0,
     41051.6516D0,1057.1155D0,1062.5836D0,1068.0558D0,1073.5323D0,
     E1079.0129D0/      
      DATA(FLOG(I),I=242,271)/1084.4977D0,1089.9866D0,1095.4797D0,1100.
     A9768D0,11
     106.4781D0,1111.9834D0,1117.4928D0,1123.0063D0,1128.5237D0,1134.
     B0452D0,1139.570
     26D0,1145.1001D0,1150.6335D0,1156.1708D0,1161.7120D0,1167.2573D0,
     C1172.8063D0,1178
     3.3593D0,1183.9161D0,1189.4768D0,1195.0413D0,1200.6097D0,1206.
     D1818D0,1211.7577D0,
     41217.3375D0,1222.9209D0,1228.5082D0,1234.0992D0,1239.6939D0,
     E1245.2924D0/      
      DATA(FLOG(I),I=272,301)/1250.8944D0,1256.5003D0,1262.1097D0,
     A1267.7228D0,12
     173.3396D0,1278.9600D0,1284.5840D0,1290.2117D0,1295.8429D0,
     B1301.4777D0,1307.116
     20D0,1312.7580D0,1318.4034D0,1324.0524D0,1329.7048D0,1335.3609D0,
     C1341.0203D0,1346
     3.6833D0,1352.3497D0,1358.0196D0,1363.6929D0,1369.3697D0,
     D1375.0499D0,1380.7334D0,
     41386.4204D0,1392.1107D0,1397.8045D0,1403.5016D0,1409.2020D0,
     E1414.9058D0/      
      DATA EPS1,EPS2/.1D0,-.2D0/                                            TRO04370
C                                                                       TRO04380
C     CALCUL DES COMBINAISONS XJ11,XJ12,XJ13,XJ21,XJ22,XJ23,XJ31,XJ32,XJTRO04390
C                                                                       TRO04400
      XN=-XJ11+XJ21+XJ31+EPS1                                           TRO04410
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO04420
      N1=XN                                                             TRO04430
      XN=-XJ32+XJ33+XJ31+EPS1                                           TRO04440
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO04450
      N2=XN                                                             TRO04460
      XN=XJ11-XJ21+XJ31+EPS1                                            TRO04470
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO04480
      N5=XN                                                             TRO04490
      XN=XJ32-XJ33+XJ31+EPS1                                            TRO04500
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO04510
      N6=XN                                                             TRO04520
      XN=XJ11+XJ21-XJ31+EPS1                                            TRO04530
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO04540
      N9=XN                                                             TRO04550
      XN=XJ32+XJ33-XJ31+EPS1                                            TRO04560
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO04570
      N10=XN                                                            TRO04580
      XN=XJ11+XJ32+XJ21+XJ33+EPS1                                       TRO04590
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO04600
      N15=XN                                                            TRO04610
      N15=N15+1                                                         TRO04620
      XN=XJ11+XJ21+XJ31+EPS1                                            TRO04630
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO04640
      N16=XN                                                            TRO04650
      N16=N16+1                                                         TRO04660
      XN=XJ32+XJ33+XJ31+EPS1                                            TRO04670
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO04680
      N17=XN                                                            TRO04690
      N17=N17+1                                                         TRO04700
      XN=-XJ12+XJ22+XJ32+EPS1                                           TRO04710
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO04720
      N21=XN                                                            TRO04730
      XN=-XJ21+XJ22+XJ23+EPS1                                           TRO04740
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO04750
      N23=XN                                                            TRO04760
      XN=XJ12-XJ22+XJ32+EPS1                                            TRO04770
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO04780
      N25=XN                                                            TRO04790
      XN=XJ21-XJ22+XJ23+EPS1                                            TRO04800
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO04810
      N27=XN                                                            TRO04820
      XN=XJ12+XJ22-XJ32+EPS1                                            TRO04830
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO04840
      N29=XN                                                            TRO04850
      XN=XJ21+XJ22-XJ23+EPS1                                            TRO04860
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO04870
      N31=XN                                                            TRO04880
      XN=-XJ12-XJ21+XJ32+XJ23+EPS1                                      TRO04890
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO04900
      N33=XN                                                            TRO04910
      XN=XJ12+XJ22+XJ32+EPS1                                            TRO04920
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO04930
      N36=XN                                                            TRO04940
      N36=N36+1                                                         TRO04950
      XN=XJ21+XJ22+XJ23+EPS1                                            TRO04960
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO04970
      N38=XN                                                            TRO04980
      N38=N38+1                                                         TRO04990
      XN=-XJ13+XJ23+XJ33+EPS1                                           TRO05000
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO05010
      N41=XN                                                            TRO05020
      XN=-XJ13+XJ11+XJ12+EPS1                                           TRO05030
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO05040
      N44=XN                                                            TRO05050
      XN=XJ13-XJ23+XJ33+EPS1                                            TRO05060
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO05070
      N45=XN                                                            TRO05080
      XN=XJ13-XJ11+XJ12+EPS1                                            TRO05090
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO05100
      N48=XN                                                            TRO05110
      XN=XJ13+XJ23-XJ33+EPS1                                            TRO05120
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO05130
      N49=XN                                                            TRO05140
      XN=XJ13+XJ11-XJ12+EPS1                                            TRO05150
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO05160
      N52=XN                                                            TRO05170
      XN=-XJ23-XJ11+XJ33+XJ12+EPS1                                      TRO05180
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO05190
      N54=XN                                                            TRO05200
      XN=XJ13+XJ23+XJ33+EPS1                                            TRO05210
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO05220
      N56=XN                                                            TRO05230
      N56=N56+1                                                         TRO05240
      XN=XJ13+XJ11+XJ12+EPS1                                            TRO05250
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO05260
      N59=XN                                                            TRO05270
      N59=N59+1                                                         TRO05280
C                                                                       TRO05290
C     TEST SUR LES XJ11,XJ12,XJ13,XJ21,XJ22,XJ23,XJ31,XJ32,XJ33         TRO05300
C                                                                       TRO05310
      IF(N9.LT.0) GO TO 50                                              TRO05320
      IF(N5.LT.0) GO TO 50                                              TRO05330
      IF(N1.LT.0) GO TO 50                                              TRO05340
      IF(N10.LT.0) GO TO 50                                             TRO05350
      IF(N6.LT.0) GO TO 50                                              TRO05360
      IF(N2.LT.0) GO TO 50                                              TRO05370
      IF(N29.LT.0) GO TO 50                                             TRO05380
      IF(N25.LT.0) GO TO 50                                             TRO05390
      IF(N21.LT.0) GO TO 50                                             TRO05400
      IF(N31.LT.0) GO TO 50                                             TRO05410
      IF(N27.LT.0) GO TO 50                                             TRO05420
      IF(N23.LT.0) GO TO 50                                             TRO05430
      IF(N49.LT.0) GO TO 50                                             TRO05440
      IF(N45.LT.0) GO TO 50                                             TRO05450
      IF(N41.LT.0) GO TO 50                                             TRO05460
      IF(N52.LT.0) GO TO 50                                             TRO05470
      IF(N48.LT.0) GO TO 50                                             TRO05480
      IF(N44.LT.0) GO TO 50                                             TRO05490
      K=N1+N2+N5+N6+N9+N10+N21+N23+N25+N27+N29+N31+N41+N44+N45+N48+N49+NTRO05500
     152                                                                TRO05510
      IF(K.GT.0) GO TO 54                                               TRO05520
      C9J=1.D0                                                            TRO05530
      RETURN                                                            TRO05540
  50  C9J=0.D0                                                            TRO05550
      RETURN                                                            TRO05560
C                                                                       TRO05570
C     CALCUL DE LA SOMME SUR XJ                                         TRO05580
C                                                                       TRO05590
   54 XN=2.D0*(XJ21-XJ32)+EPS1                                            TRO05600
      IF(XN.LT.0.D0) XN=-(XN+EPS2)                                        TRO05610
      JMIN=XN                                                           TRO05620
      XN=2.D0*(XJ11-XJ33)+EPS1                                            TRO05630
      IF(XN.LT.0.D0) XN=-(XN+EPS2)                                        TRO05640
      N=XN                                                              TRO05650
      IF(N.GT.JMIN) JMIN=N                                              TRO05660
      XN=2.D0*(XJ12-XJ23)+EPS1                                            TRO05670
      IF(XN.LT.0.D0) XN=-(XN+EPS2)                                        TRO05680
      N=XN                                                              TRO05690
      IF(N.GT.JMIN) JMIN=N                                              TRO05700
      XN=2.D0*(XJ21+XJ32)+EPS1                                            TRO05710
      JMAX=XN                                                           TRO05720
      XN=2.D0*(XJ11+XJ33)+EPS1                                            TRO05730
      N=XN                                                              TRO05740
      IF(N.LT.JMAX) JMAX=N                                              TRO05750
      XN=2.D0*(XJ12+XJ23)+EPS1                                            TRO05760
      N=XN                                                              TRO05770
      IF(N.LT.JMAX) JMAX=N                                              TRO05780
      XJMIN=DFLOAT(JMIN)/2.D0                                              TRO05790
      XJMAX=DFLOAT(JMAX)/2.D0                                              TRO05800
      S=0.D0                                                              TRO05810
      XJ=XJMIN                                                          TRO05820
      XN=-XJ32+XJ21+XJ+EPS1                                             TRO05830
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO05840
      N3=XN                                                             TRO05850
      XN=-XJ11+XJ33+XJ+EPS1                                             TRO05860
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO05870
      N4=XN                                                             TRO05880
      XN=XJ32-XJ21+XJ+EPS1                                              TRO05890
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO05900
      N7=XN                                                             TRO05910
      XN=XJ11-XJ33+XJ+EPS1                                              TRO05920
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO05930
      N8=XN                                                             TRO05940
      XN=XJ32+XJ21-XJ+EPS1                                              TRO05950
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO05960
      N11=XN                                                            TRO05970
      XN=XJ11+XJ33-XJ+EPS1                                              TRO05980
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO05990
      N12=XN                                                            TRO06000
      XN=-XJ11-XJ32+XJ31+XJ+EPS1                                        TRO06010
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO06020
      N13=XN                                                            TRO06030
      XN=-XJ21-XJ33+XJ31+XJ+EPS1                                        TRO06040
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO06050
      N14=XN                                                            TRO06060
      XN=XJ32+XJ21+XJ+EPS1                                              TRO06070
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO06080
      N18=XN                                                            TRO06090
      N18=N18+1                                                         TRO06100
      XN=XJ11+XJ33+XJ+EPS1                                              TRO06110
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO06120
      N19=XN                                                            TRO06130
      N19=N19+1                                                         TRO06140
      XN=-XJ21+XJ+XJ32+EPS1                                             TRO06150
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO06160
      N22=XN                                                            TRO06170
      XN=-XJ12+XJ23+XJ+EPS1                                             TRO06180
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO06190
      N24=XN                                                            TRO06200
      XN=XJ21-XJ+XJ32+EPS1                                              TRO06210
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO06220
      N26=XN                                                            TRO06230
      XN=XJ12-XJ+XJ23+EPS1                                              TRO06240
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO06250
      N28=XN                                                            TRO06260
      XN=XJ21+XJ-XJ32+EPS1                                              TRO06270
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO06280
      N30=XN                                                            TRO06290
      XN=XJ12+XJ-XJ23+EPS1                                              TRO06300
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO06310
      N32=XN                                                            TRO06320
      XN=-XJ22-XJ+XJ32+XJ23+EPS1                                        TRO06330
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO06340
      N34=XN                                                            TRO06350
      XN=XJ12+XJ21+XJ22+XJ+EPS1                                         TRO06360
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO06370
      N35=XN                                                            TRO06380
      N35=N35+1                                                         TRO06390
      XN=XJ21+XJ+XJ32+EPS1                                              TRO06400
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO06410
      N37=XN                                                            TRO06420
      N37=N37+1                                                         TRO06430
      XN=XJ12+XJ+XJ23+EPS1                                              TRO06440
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO06450
      N39=XN                                                            TRO06460
      N39=N39+1                                                         TRO06470
      XN=-XJ+XJ11+XJ33+EPS1                                             TRO06480
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO06490
      N42=XN                                                            TRO06500
      XN=-XJ+XJ23+XJ12+EPS1                                             TRO06510
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO06520
      N43=XN                                                            TRO06530
      XN=XJ-XJ11+XJ33+EPS1                                              TRO06540
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO06550
      N46=XN                                                            TRO06560
      XN=XJ-XJ23+XJ12+EPS1                                              TRO06570
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO06580
      N47=XN                                                            TRO06590
      XN=XJ+XJ11-XJ33+EPS1                                              TRO06600
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO06610
      N50=XN                                                            TRO06620
      XN=XJ+XJ23-XJ12+EPS1                                              TRO06630
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO06640
      N51=XN                                                            TRO06650
      XN=-XJ13-XJ+XJ33+XJ12+EPS1                                        TRO06660
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO06670
      N53=XN                                                            TRO06680
      XN=XJ13+XJ+XJ23+XJ11+EPS1                                         TRO06690
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO06700
      N55=XN                                                            TRO06710
      N55=N55+1                                                         TRO06720
      XN=XJ+XJ11+XJ33+EPS1                                              TRO06730
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO06740
      N57=XN                                                            TRO06750
      N57=N57+1                                                         TRO06760
      XN=XJ+XJ23+XJ12+EPS1                                              TRO06770
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO06780
      N58=XN                                                            TRO06790
      N58=N58+1                                                         TRO06800
      GO TO 10                                                          TRO06810
   52 IF(XJ.GT.XJMAX) GO TO 120                                         TRO06820
      N3=N3+1                                                           TRO06830
      N4=N4+1                                                           TRO06840
      N7=N7+1                                                           TRO06850
      N8=N8+1                                                           TRO06860
      N11=N11-1                                                         TRO06870
      N12=N12-1                                                         TRO06880
      N13=N13+1                                                         TRO06890
      N14=N14+1                                                         TRO06900
      N18=N18+1                                                         TRO06910
      N19=N19+1                                                         TRO06920
      N22=N22+1                                                         TRO06930
      N24=N24+1                                                         TRO06940
      N26=N26-1                                                         TRO06950
      N28=N28-1                                                         TRO06960
      N30=N30+1                                                         TRO06970
      N32=N32+1                                                         TRO06980
      N34=N34-1                                                         TRO06990
      N35=N35+1                                                         TRO07000
      N37=N37+1                                                         TRO07010
      N39=N39+1                                                         TRO07020
      N42=N42-1                                                         TRO07030
      N43=N43-1                                                         TRO07040
      N46=N46+1                                                         TRO07050
      N47=N47+1                                                         TRO07060
      N50=N50+1                                                         TRO07070
      N51=N51+1                                                         TRO07080
      N53=N53-1                                                         TRO07090
      N55=N55+1                                                         TRO07100
      N57=N57+1                                                         TRO07110
      N58=N58+1                                                         TRO07120
C                                                                       TRO07130
C     CALCUL DES SOMMES ALTERNEES S1(K1),S2(K2),S3(K3)                  TRO07140
C                                                                       TRO07150
   10 K1=0.D0                                                             TRO07160
      L1=-N13                                                           TRO07170
      IF(L1.GT.K1) K1=L1                                                TRO07180
      L1=-N14                                                           TRO07190
      IF(L1.GT.K1) K1=L1                                                TRO07200
      L1=N9                                                             TRO07210
      IF(N10.LT.L1) L1=N10                                              TRO07220
      IF(N11.LT.L1) L1=N11                                              TRO07230
      IF(N12.LT.L1) L1=N12                                              TRO07240
      IF(N15.LT.L1) L1=N15                                              TRO07250
      F1=1.D0                                                             TRO07260
      S1=1.D0                                                             TRO07270
      I1=K1+1                                                           TRO07280
  62  IF(I1.GT.L1) GO TO 64                                             TRO07290
      I1M1=I1-1                                                         TRO07300
      NN1=(N9-I1M1)*(N10-I1M1)*(N11-I1M1)*(N12-I1M1)                    TRO07310
      ND1=I1*(N13+I1)*(N14+I1)*(N15-I1M1)                               TRO07320
      F1=-F1*DFLOAT(NN1)/DFLOAT(ND1)                                      TRO07330
      S1=S1+F1                                                          TRO07340
      I1=I1+1                                                           TRO07350
      GO TO 62                                                          TRO07360
   64 K2=0                                                              TRO07370
      L2=-N33                                                           TRO07380
      IF(L2.GT.K2) K2=L2                                                TRO07390
      L2=-N34                                                           TRO07400
      IF(L2.GT.K2) K2=L2                                                TRO07410
      L2=N29                                                            TRO07420
      IF(N30.LT.L2) L2=N30                                              TRO07430
      IF(N31.LT.L2) L2=N31                                              TRO07440
      IF(N32.LT.L2) L2=N32                                              TRO07450
      IF(N35.LT.L2) L2=N35                                              TRO07460
      F2=1.D0                                                             TRO07470
      S2=1.D0                                                             TRO07480
      I2=K2+1                                                           TRO07490
  70  IF(I2.GT.L2) GO TO 80                                             TRO07500
      I2M2=I2-1                                                         TRO07510
      NN2=(N29-I2M2)*(N30-I2M2)*(N31-I2M2)*(N32-I2M2)                   TRO07520
      ND2=I2*(N33+I2)*(N34+I2)*(N35-I2M2)                               TRO07530
      F2=-F2*DFLOAT(NN2)/DFLOAT(ND2)                                      TRO07540
      S2=S2+F2                                                          TRO07550
      I2=I2+1                                                           TRO07560
      GO TO 70                                                          TRO07570
   80 K3=0                                                              TRO07580
      L3=-N53                                                           TRO07590
      IF(L3.GT.K3) K3=L3                                                TRO07600
      L3=-N54                                                           TRO07610
      IF(L3.GT.K3) K3=L3                                                TRO07620
      L3=N49                                                            TRO07630
      IF(N50.LT.L3) L3=N50                                              TRO07640
      IF(N51.LT.L3) L3=N51                                              TRO07650
      IF(N52.LT.L3) L3=N52                                              TRO07660
      IF(N55.LT.L3) L3=N55                                              TRO07670
      F3=1.D0                                                             TRO07680
      S3=1.D0                                                             TRO07690
      I3=K3+1                                                           TRO07700
   84 IF(I3.GT.L3) GO TO 90                                             TRO07710
      I3M3=I3-1                                                         TRO07720
      NN3=(N49-I3M3)*(N50-I3M3)*(N51-I3M3)*(N52-I3M3)                   TRO07730
      ND3=I3*(N53+I3)*(N54+I3)*(N55-I3M3)                               TRO07740
      F3=-F3*DFLOAT(NN3)/DFLOAT(ND3)                                      TRO07750
      S3=S3+F3                                                          TRO07760
      I3=I3+1                                                           TRO07770
      GO TO 84                                                          TRO07780
C                                                                       TRO07790
C     CALCUL DE LA RACINE D'UN TERME DE LA SOMME SUR J                  TRO07800
C                                                                       TRO07810
   90 S2N=FLOG(N3+1)+FLOG(N4+1)+FLOG(N7+1)+FLOG(N8+1)+FLOG(N11+1)+FLOG(NTRO07820
     112+1)+FLOG(N22+1)+FLOG(N24+1)+FLOG(N26+1)+FLOG(N28+1)+FLOG(N30+1)+TRO07830
     2FLOG(N32+1)+FLOG(N42+1)+FLOG(N43+1)+FLOG(N46+1)+FLOG(N47+1)+FLOG(NTRO07840
     350+1)+FLOG(N51+1)                                                 TRO07850
      S2N=.5D0*S2N                                                        TRO07860
      S2D=FLOG(N18+1)+FLOG(N19+1)+FLOG(N37+1)+FLOG(N39+1)+FLOG(N57+1)+FLTRO07870
     1OG(N58+1)                                                         TRO07880
      S2D=.5D0*S2D                                                        TRO07890
      KM1=K1-1                                                          TRO07900
      KP1=K1+1                                                          TRO07910
      KM2=K2-1                                                          TRO07920
      KP2=K2+1                                                          TRO07930
      KM3=K3-1                                                          TRO07940
      KP3=K3+1                                                          TRO07950
      S2N=S2N+FLOG(N15-KM1)+FLOG(N35-KM2)+FLOG(N55-KM3)                 TRO07960
      S2D=S2D+FLOG(KP1)+FLOG(KP2)+FLOG(KP3)+FLOG(N9-KM1)+FLOG(N10-KM1)+FTRO07970
     1LOG(N11-KM1)+FLOG(N12-KM1)+FLOG(N13+KP1)+FLOG(N14+KP1)+FLOG(N29-KMTRO07980
     22)+FLOG(N30-KM2)+FLOG(N31-KM2)+FLOG(N32-KM2)+FLOG(N33+KP2)+FLOG(N3TRO07990
     34+KP2)+FLOG(N49-KM3)+FLOG(N50-KM3)+FLOG(N51-KM3)                  TRO08000
     4+FLOG(N52-KM3)+FLOG(N53+KP3)+FLOG(N54+KP3)                        TRO08010
      F=S2D-S2N                                                         TRO08020
      IF (F.GT.80.D0) GO TO 100                                            TRO08030
      F=S2N/S2D                                                         TRO08040
      IF((F.LT.1.01D0).AND.(F.GT.0.98D0)) GO TO 100                         TRO08050
      F=S1*S2*S3*DEXP(S2N-S2D)                                           TRO08060
      GO TO 110                                                         TRO08070
  100 F=S1*S2*S3                                                        TRO08080
      IF(F) 102,112,104                                                 TRO08090
  102 F=DLOG(-F)                                                        TRO08100
      F=-DEXP(F+S2N-S2D)                                                 TRO08110
      GO TO 110                                                         TRO08120
  104 F=DLOG(F)                                                         TRO08130
      F=DEXP(F+S2N-S2D)                                                  TRO08140
  110 F=F*(2.D0*XJ+1.D0)                                                    TRO08150
C                                                                       TRO08160
C     CALCUL DE LA PHASE D'UN TERME DE LA SOMME SUR J                   TRO08170
C                                                                       TRO08180
      L1=K1+K2+K3                                                       TRO08190
      K1=L1/2                                                           TRO08200
      K1=2*K1                                                           TRO08210
      IF(L1.NE.K1) F=-F                                                 TRO08220
      S=S+F                                                             TRO08230
  112 XJ=XJ+1.D0                                                          TRO08240
      GO TO 52                                                          TRO08250
C                                                                       TRO08260
C     CALCUL DU C9J SANS PHASE                                          TRO08270
C                                                                       TRO08280
  120 C2N=FLOG(N1+1)+FLOG(N2+1)+FLOG(N5+1)+FLOG(N6+1)+FLOG(N9+1)+FLOG(N1TRO08290
     10+1)+FLOG(N21+1)+FLOG(N23+1)+FLOG(N25+1)+FLOG(N27+1)+FLOG(N29+1)+FTRO08300
     2LOG(N31+1)+FLOG(N41+1)+FLOG(N44+1)+FLOG(N45+1)+FLOG(N48+1)+FLOG(N4TRO08310
     39+1)+FLOG(N52+1)                                                  TRO08320
      C2N=.5D0*C2N                                                        TRO08330
      C2D=FLOG(N16+1)+FLOG(N17+1)+FLOG(N36+1)+FLOG(N38+1)+FLOG(N56+1)+FLTRO08340
     1OG(N59+1)                                                         TRO08350
      C2D=.5D0*C2D                                                        TRO08360
      F=C2D-C2N                                                         TRO08370
      IF(F.GT.80.D0) GO TO 122                                            TRO08380
      F=C2N/C2D                                                         TRO08390
      IF((F.LT.1.01D0).AND.(F.GT.0.98D0)) GO TO 122                         TRO08400
      C9J=S*DEXP(C2N-C2D)                                                TRO08410
      GO TO 130                                                         TRO08420
  122 IF(S) 124,50,126                                                  TRO08430
  124 S=DLOG(-S)                                                        TRO08440
      C9J=-DEXP(S+C2N-C2D)                                               TRO08450
      GO TO 130                                                         TRO08460
  126 S=DLOG(S)                                                         TRO08470
      C9J=DEXP(S+C2N-C2D)                                                TRO08480
C                                                                       TRO08490
C     CALCUL DE LA PHASE                                                TRO08500
C                                                                       TRO08510
  130 K=N9+N16+N36+N56-1                                                TRO08520
      L=K/2                                                             TRO08530
      L=2*L                                                             TRO08540
      IF (L.NE.K) C9J=-C9J                                              TRO08550
      RETURN                                                            TRO08560
      END

C***********************************************************************

      SUBROUTINE TRED2(A,NV,NP,D,E)

      implicit real*8 (a-h,o-z)
      dimension a(np,np),d(np),e(np)

      IF(NV.GT.1)THEN
        DO 18 I=NV,2,-1  
          L=I-1
          H=0.d0
          SCALE=0.d0
          IF(L.GT.1)THEN
            DO 11 K=1,L
              SCALE=SCALE+DABS(A(I,K))
11          CONTINUE
            IF(SCALE.EQ.0.d0)THEN
              E(I)=A(I,L)
            ELSE
              DO 12 K=1,L
                A(I,K)=A(I,K)/SCALE
                H=H+A(I,K)**2
12            CONTINUE
              F=A(I,L)
              G=-SIGN(dsqrt(H),F)
              E(I)=SCALE*G
              H=H-F*G
              A(I,L)=F-G
              F=0.d0
              DO 15 J=1,L
                A(J,I)=A(I,J)/H
                G=0.d0
                DO 13 K=1,J
                  G=G+A(J,K)*A(I,K)
13              CONTINUE
                IF(L.GT.J)THEN
                  DO 14 K=J+1,L
                    G=G+A(K,J)*A(I,K)
14                CONTINUE
                ENDIF
                E(J)=G/H
                F=F+E(J)*A(I,J)
15            CONTINUE
              HH=F/(H+H)
              DO 17 J=1,L
                F=A(I,J)
                G=E(J)-HH*F
                E(J)=G
                DO 16 K=1,J
                  A(J,K)=A(J,K)-F*E(K)-G*A(I,K)
16              CONTINUE
17            CONTINUE
            ENDIF
          ELSE
            E(I)=A(I,L)
          ENDIF
          D(I)=H
18      CONTINUE
      ENDIF
      D(1)=0.d0
      E(1)=0.d0
      DO 23 I=1,NV
        L=I-1
        IF(D(I).NE.0.)THEN
          DO 21 J=1,L
            G=0.d0
            DO 19 K=1,L
              G=G+A(I,K)*A(K,J)
19          CONTINUE
            DO 20 K=1,L
              A(K,J)=A(K,J)-G*A(K,I)
20          CONTINUE
21        CONTINUE
        ENDIF
        D(I)=A(I,I)
        A(I,I)=1.d0
        IF(L.GE.1)THEN
          DO 22 J=1,L
            A(I,J)=0.d0
            A(J,I)=0.d0
22        CONTINUE
        ENDIF
23    CONTINUE

      RETURN
      END
C***********************************************************************

      SUBROUTINE TQLI(D,E,NV,NP,Z)

      implicit real*8 (a-h,o-z)
      dimension d(np),e(np),z(np,np)

      IF (NV.GT.1) THEN
        DO 11 I=2,NV
          E(I-1)=E(I)
11      CONTINUE
        E(NV)=0.d0
        DO 15 L=1,NV
          ITER=0
1         DO 12 M=L,NV-1
            DD=DABS(D(M))+DABS(D(M+1))
            IF (DABS(E(M))+DD.EQ.DD) GO TO 2
12        CONTINUE
          M=NV
2         IF(M.NE.L)THEN
            IF(ITER.EQ.30)stop  'too many iterations'
            ITER=ITER+1
            G=(D(L+1)-D(L))/(2.d0*E(L))
            R=dsqrt(G**2+1.d0)
            G=D(M)-D(L)+E(L)/(G+SIGN(R,G))
            S=1.d0
            C=1.d0
            P=0.d0
            DO 14 I=M-1,L,-1
              F=S*E(I)
              B=C*E(I)
              IF(DABS(F).GE.DABS(G))THEN
                C=G/F
                R=dsqrt(C**2+1.d0)
                E(I+1)=F*R
                S=1.d0/R
                C=C*S
              ELSE
                S=F/G
                R=dsqrt(S**2+1.d0)
                E(I+1)=G*R
                C=1.d0/R  
                S=S*C
              ENDIF
              G=D(I+1)-P
              R=(D(I)-G)*S+2.d0*C*B
              P=S*R
              D(I+1)=G+P
              G=C*R-B
              DO 13 K=1,NV
                F=Z(K,I+1)
                Z(K,I+1)=S*Z(K,I)+C*F
                Z(K,I)=C*Z(K,I)-S*F
13            CONTINUE
14          CONTINUE
            D(L)=D(L)-P
            E(L)=G
            E(M)=0.d0
            GO TO 1
          ENDIF
15      CONTINUE
      ENDIF
      RETURN
      END

C***********************************************************************
           
      SUBROUTINE EIGSRT(D,V,W,N,NP)

C     REARRANGE ORDER AFTER INCREASING FREQUENCY D(I)

      implicit real*8 (a-h,o-z)
      dimension d(np),V(np,np),W(np,np)

      DO 10 I=1,N-1
        K=I
        P=D(I)
        DO 20 J=I+1,N
          IF(D(J).LE.P)THEN
            K=J
            P=D(J)
          ENDIF
20      CONTINUE
        IF(K.NE.I)THEN
          D(K)=D(I)
          D(I)=P
          DO 30 J=1,N
            P=V(J,I)
            V(J,I)=V(J,K)
            V(J,K)=P
            P=W(J,I)
            W(J,I)=W(J,K)
            W(J,K)=P
30        CONTINUE
        ENDIF
10    CONTINUE

      RETURN
      END

C***********************************************************************
