      PROGRAM VOF
C     +++ With Obstacles and With Velocity Perturbations
C     +++                    by
C     +++              R.W. Douglass
C     +++          7-'87, 2-'88, 5-'88, 6-'88
C     +++          7-2025
C     *** SOLA-VOF  VOLUME OF FLUID METHOD ***
C
      INCLUDE 'dint.h'
C
      real*8 Urms, Vrms, Ubar, Vbar
      CHARACTER(30) CNAME
      CHARACTER(20) fileName
      CHARACTER(12) icT
      INTEGER NX, NY

      common /new/ pltst1, pltst2, PLTDT1, PLTDT2

      COMMON /FL1/ NDIS,IDIS,IFL,IFR,JFB,JFT
C
      COMMON /FUND/FLENG,FTIME
C
      COMMON /OBS/ NOBS, IOMIN(NTOB), IOMAX(NTOB), JOMIN(NTOB),
     1 JOMAX(NTOB)
C
      COMMON /CONST/ EMF, EMF1, EM6, EM10, EP10, PI, RPD
C
      NAMELIST /XPUT/ DELT,VNU,ICYL,EPSI,GX,GY,UI,VI,VELMX,TWFIN,PRTDT
     1 ,PLTDT,OMG,ALPHA,IWL,IWR,IWT,IWB,IMOVY,AUTOT,FLHT,ISYMPLT,SIGMA
     2 ,ISURF10,CANGLE,CSQ,NMAT,RHOF,RHOFC,XPL,XPR,YPB,YPT,NPX,NPY
C
      NAMELIST /MSHSET/ NKX,XL,XC,NXL,NXR,DXMN,NKY,YL,YC,NYL,NYR,DYMN
C
      DATA EMF /1.0D-06/, EM6 /1.0D-06/, EM10 /1.0D-10/
      DATA EP10 /1.0D+10/
      PI = 4.D0*DATAN(1.D0)
      RPD = PI / 180.D0
C
C    READ THE INPUT DATA
C
      write(icT,'(I12)') JBAR2-2
      fileName = "int"//TRIM(adjustl(icT))//".in"

      OPEN(UNIT=1,FILE=fileName,STATUS='OLD')

      READ(1,120) CNAME
C
C      *** Select whether a plot-data file is to be saved.
C
      READ(1,"(/,15X,I5)") IPLOT
C
      READ(1,182) FLENG,FTIME
C
      ffleng = fleng

      READ( 1,110) VNU,VNUC,ICYL,EPSI,GX,GY,UI,VI,VELMX,IMOVY,OMG,ALPHA,
     1  IWL,IWR,IWT,IWB,CSQ,AUTOT,ISYMPLT,ISURF10,SIGMA,CANGLE,NMAT,RHOF
     2  ,RHOFC,FLHT,XPL,YPB,XPR,YPT,NPX,NPY,NOBS,NDIS,IDIS,VORIG,
     3  ISHVEL,SHVEL,DELT,TWFIN,PRTDT,pltst1,pltst2,PLTDT1,PLTDT2
C 
C     *** CONVERT INPUTS TO FUNDIMENTAL DIMENSIONS
C
      VNU = VNU * (FTIME/FLENG**2)
      VNUC = VNUC * (FTIME/FLENG**2)
      GX = GX * (FTIME**2/FLENG)
      GY = GY * (FTIME**2/FLENG)
      RHOF = RHOF * (FLENG**3)
      RHOFC = RHOFC * (FLENG**3)
      SIGMA = SIGMA * (FTIME**2)
      VORIG = VORIG * (FTIME/FLENG)
      write(6,"('VORIG:', 2X, F10.7)") VORIG
C
C     *** NOBS is the number of solid obstacles appearing within the
C     ***    computational domain.  Data must be added to the " ".IN file.
C     *** NDIS is used when stability computations are being done, so
C     ***    that if NDIS is not 0 superpose disturbances on the initial
C     ***    velocity distribution.  Data must be added to the " ".IN file.
C
      READ ( 1,140) NKX,NKY,NONUNIF
C      WRITE(6,"('NKX:' I5, 'NKY:' I5, 'NONUNIF:' I5)") NKX, NKY, NONUNIF
        DO I=1,NKX
          IF (I .EQ. 1) THEN
            READ ( 1,141) XC(I),NXL(I),NXR(I),DXMN(I)
          ELSE
            READ ( 1,142) XC(I),NXL(I),NXR(I),DXMN(I)
          END IF
        END DO 
        DO J=1,NKY
          IF (J .EQ. 1) THEN
            READ ( 1,141) YC(J),NYL(J),NYR(J),DYMN(J)
          ELSE
            READ ( 1,142) YC(J),NYL(J),NYR(J),DYMN(J)
          END IF
        END DO
        READ ( 1,150) (XL(I),I=1,NKX+1)
        READ ( 1,150) (YL(J),J=1,NKY+1)
C
C     *** CALCULATE VARIABLE MESH DATA
C
      CALL MESHSET
C
C     *** PRINT INITIAL INPUT DATA
C
      CALL PRT (1)
C
C     *** SET INITIAL CONDITIONS
C
      CALL SETUP
C
C      print '("IBAR: ",I5," IBAR2: ",I5," JBAR: ",I5, " JBAR2: ", I5)', 
C     1 IBAR, IBAR2, JBAR, JBAR2
C      print '("IMAX: ", I5, " IM1: ", I5," JMAX: ",I5," JM1: ",I5)',
C     1 IMAX, IM1, JMAX, JM1
C
C     *** SET INITIAL BOUNDARY CONDITIONS
C
      CALL BC
C
C     *** Set up the plot-data and area files.
C
      GO TO 20

c       this is the end of the things that are done only for new
c           runs.
C
C     *** START TIME CYCLE
C
C
   10 CONTINUE
      ITER=0
      FLG=1.D0
      FNOC=0.D0
C
C     *** EXPLICITLY APPROXIMATE NEW TIME-LEVEL VELOCITIES
C
      CALL TILDE
C
      IF (NMAT.EQ.2.AND.ISURF10.EQ.1) CALL TMS10
C
C     *** SET BOUNDARY CONDITIONS
C
      CALL BC
C
C     *** ITERATIVELY ADJUST CELL PRESSURE AND VELOCITY
C
      CALL PRESSIT
C
      IF (T.GT.EP10) GO TO 30
C
C     *** UPDATE FLUID CONFIGURATION
C
   20 CALL VFCONV
C
      IF (FLGC.GT.0.5D0) GO TO 90
C
C     *** SET BOUNDARY CONDITIONS
C
      CALL BC
C
C     *** MOVE MARKER PARTICLES
C
      CALL PARMOV
C
C     *** DETERMINE PRESSURE INTERPOLATION FACTOR AND NEIGHBOR
C     *** ALSO DETERMINE SURFACE TENSION PRESSURES AND
C     *** WALL ADHESION EFFECTS IN SURFACE CELLS
C
      CALL PETACAL
C
C     *** PRINT TIME AND CYCLE DATA
C
      CALL PRT (2)
   30 continue
C
      IF (iCYCLE.LE.0) GO TO 40
      IF (T+EM6.LT.TWPLT) GO TO 50
   40 CONTINUE                                        
      if(twplt+em6.gt.pltst2.and.pltdt.eq.pltdt1) pltdt = pltdt2
      twplt = twplt+pltdt
C
C     *** SAVE VELOCITY VECTOR, FREE SURFACE, AND MESH DATA FOR
C     *** FURTHER PROCESSING FOR PLOTTING ON DISSPLA.
C     *** OR USING PARAVIEW VTK 
C
C      WRITE(6,"('IPLOT:', I5)") IPLOT
      IF (IPLOT .eq. 1) THEN 
       CALL WRITEVTK()
      END IF
C
   50 CONTINUE
      IF (iCYCLE.LE.0) GO TO 60
      IF (T+EM6.LT.TWPRT) GO TO 70
      TWPRT=TWPRT+PRTDT
   60 CONTINUE
C
C     *** PRINT FIELD VARIABLE DATA TO OUTPUT FILE
C
           CALL PRT (3)
C
   70 CONTINUE
C
C     *** SET THE ADVANCE TIME ARRAYS INTO THE TIME-N ARRAYS
C
      DO I=1,IMAX
        DO J=1,JMAX
          UN(I,J)=U(I,J)
          VN(I,J)=V(I,J)
          U(I,J)=0.D0
          V(I,J)=0.D0
          PN(I,J)=P(I,J)
          FN(I,J)=F(I,J)
        END DO
      END DO

      NREGN=NREG

   90 CONTINUE
C
C     *** ADJUST DELT
C
      CALL DELTADJ
C
C     *** ADVANCE TIME
C
      T=T+DELT
      IF (T.GT.TWFIN) GO TO 100
C
C     *** ADVANCE CYCLE
C
      iCYCLE=iCYCLE+1
      IF (NFLGC.GE.25.OR.NOCON.GE.25) T=EP10

      GO TO 10
C
  100 CLOSE(UNIT=1)
C
      STOP
C
110   FORMAT(/2(/,15X,E16.5),/,15X,I2,/,15X,F15.6,/,15X,F10.4,/,15X,
     1  F10.4,/,15X,F10.4,/,15X,F10.4,/,15X,F10.4,/,15X,I2,/15X,F15.6,/,
     2  15X,F15.6,/,15X,I2,/,15X,I2,/,15X,I2,/,15X,I2,/,15X,F10.4,/,15X,
     3  F10.4,2(/,15X,I2),/,15X,F15.6,/,15X,F10.4,/,15X,I2,/,15X,F11.5,
     4  6(/,15X,F11.5),5(/,15X,I4),/,15X,F17.12,/,15X,I4,/,15X,F10.7,
     5  ///,7(/,16X,F15.6))
120   FORMAT(/,30A1)
140   FORMAT(/,3I5)
141   FORMAT(/,F5.2,2I5,1pE14.7)
142   FORMAT(F5.2,2I5,1pE14.7)
143   FORMAT(6(5X,1pE15.8))
150   FORMAT(//,16F5.2)
182   FORMAT(//,15X,F11.6,/,15X,F11.6)
      END

      SUBROUTINE BC
C
      INCLUDE 'dint.h'
C
      COMMON /CONST/ EMF, EMF1, EM6, EM10, EP10, PI, RPD
C
C     *** SET BOUNDARY CONDITIONS
C
      DO J=1,JMAX
       F(1,J)=F(2,J)
       F(IMAX,J)=F(IM1,J)
       P(1,J)=P(2,J)
       P(IMAX,J)=P(IM1,J)
C     Left Wall: IWL = 1 (rigid free-slip), 2 (rigid no-slip),
C                     3 (continuance), 4 (periodic in x), 
C                     5 (specified pressure)
       IF(IWL .EQ. 1) THEN
C      GO TO (10,20,30,40,30), IWL
        U(1,J)=0.D0
        V(1,J)=V(2,J)
      ELSE IF (IWL .EQ. 2) THEN
        U(1,J)=0.D0
        V(1,J)=-V(2,J)*DELX(1)/DELX(2)
      ELSE IF (IWL .EQ. 3 .OR. IWL .EQ. 5) THEN
        IF (ITER .LE. 0) THEN
         U(1,J)=U(2,J)*(X(2)*RX(1)*CYL+1.D0-CYL)
         V(1,J)=V(2,J)
        END IF
      ELSE IF (IWL .EQ. 4) THEN
        U(1,J)=U(IM2,J)
        V(1,J)=V(IM2,J)
        F(1,J)=F(IM2,J)
CC        BETA(1,J)=BETA(IM2,J)
      END IF
C   50 GO TO (60,70,80,90,80), IWR
      IF (IWR .EQ. 1) THEN
       U(IM1,J)=0.D0
       V(IMAX,J)=V(IM1,J)
      ELSE IF (IWR .EQ. 2) THEN
       U(IM1,J)=0.D0
       V(IMAX,J)=-V(IM1,J)*DELX(IMAX)/DELX(IM1)
      ELSE IF (IWR .EQ. 3 .OR. IWR .EQ. 5) THEN
       IF (ITER.LE.0)  THEN
        U(IM1,J)=U(IM2,J)*(X(IM2)*RX(IM1)*CYL+1.D0-CYL)
        V(IMAX,J)=V(IM1,J)
       END IF
      ELSE IF (IWR .EQ. 4)  THEN
       U(IM1,J)=U(2,J)
       V(IM1,J)=V(2,J)
       P(IM1,J)=P(2,J)
       PS(IM1,J)=PS(2,J)
       F(IM1,J)=F(2,J)
       V(IMAX,J)=V(3,J)
       F(IMAX,J)=F(3,J)
CC       BETA(IM1,J)=BETA(3,J)
CC       BETA(IMAX,J)=BETA(3,J)
      END IF
      END DO

      DO I=1,IMAX
       F(I,1)=F(I,2)
       F(I,JMAX)=F(I,JM1)
       P(I,1)=P(I,2)
       P(I,JMAX)=P(I,JM1)
C      GO TO (110,120,130,140,130), IWT
       IF (IWT .EQ. 1) THEN
        V(I,JM1)=0.D0
        U(I,JMAX)=U(I,JM1)
       ELSE IF (IWT .EQ. 2) THEN 
        V(I,JM1)=0.D0
        U(I,JMAX)=-U(I,JM1)*DELY(JMAX)/DELY(JM1)
       ELSE IF (IWT .EQ. 3 .OR. IWT .EQ. 5) THEN 
        IF (ITER.LE.0) THEN
         V(I,JM1)=V(I,JM2)
         U(I,JMAX)=U(I,JM1)
        END IF
       ELSE IF( IWT .EQ. 4) THEN
        V(I,JM1)=V(I,2)
        U(I,JM1)=U(I,2)
        P(I,JM1)=P(I,2)
        PS(I,JM1)=PS(I,2)
        F(I,JM1)=F(I,2)
        U(I,JMAX)=U(I,3)
        F(I,JMAX)=F(I,3)
       END IF
C  150 GO TO (160,170,180,190,180), IWB
       IF (IWB .EQ. 1) THEN
        V(I,1)=0.D0
        U(I,1)=U(I,2)
       ELSE IF (IWB .EQ. 2) THEN
        V(I,1)=0.D0
        U(I,1)=-U(I,2)*DELY(1)/DELY(2)
       ELSE IF(IWB .EQ. 3 .OR. IWB .EQ. 5) THEN
        IF (ITER.LE.0) THEN
         V(I,1)=V(I,2)
         U(I,1)=U(I,2)
        END IF
       ELSE IF (IWB .EQ. 4) THEN
        V(I,1)=V(I,JM2)
        U(I,1)=U(I,JM2)
        F(I,1)=F(I,JM2)
       END IF
      END DO

      IF (ISHVEL.EQ.1)THEN
       DO I=1,IMAX
        U(I,1) = SHVEL
        U(I,2) = SHVEL
        U(I,JMAX) = -SHVEL
        U(I,JM1) = -SHVEL
       END DO 
      END IF
C
C     *** FREE SURFACE AND SLOPED BOUNDARY CONDITIONS
C
      DO 450 I=2,IM1
       XRP=RDX(I)+RXI(I)/2.D0
       RXRP=1.D0/XRP
       XRM=RDX(I)-RXI(I)/2.D0
       IF (XRM.GT.0.D0) GO TO 210
       RXRM=0.D0
       GO TO 220
  210 CONTINUE
      RXRM=1.D0/XRM
  220 CONTINUE
      DO 450 J=2,JM1
      IF (BETA(I,J).GT.0.D0) GO TO 230
      BMR=0.D0
      BMT=0.D0
      BML=0.D0
      BMB=0.D0
      F(I,J)=0.D0
      P(I,J)=0.D0
      IF (BETA(I+1,J).GT.0.D0) BMR=1.D0
      IF (BETA(I,J+1).GT.0.D0) BMT=1.D0
      IF (BETA(I-1,J).GT.0.D0) BML=1.D0
      IF (BETA(I,J-1).GT.0.D0) BMB=1.D0
      BMTOT=BMR+BMT+BML+BMB
      IF (BMTOT.LE.0.D0) GO TO 450
      F(I,J)=(BMR*F(I+1,J)+BMT*F(I,J+1)+BML*F(I-1,J)+BMB*F(I,J-1))/BMTOT
      P(I,J)=(BMR*P(I+1,J)+BMT*P(I,J+1)+BML*P(I-1,J)+BMB*P(I,J-1))/BMTOT
      GO TO 450 
  230 CONTINUE
      IF (NMAT.EQ.2) GO TO 450
      IF (F(I,J).LT.EMF.OR.F(I,J).GT.EMF1) GO TO 450
      NFSB=0
      IF (F(I+1,J).LT.EMF) NFSB=NFSB+1
      IF (F(I,J+1).LT.EMF) NFSB=NFSB+2
      IF (F(I-1,J).LT.EMF) NFSB=NFSB+4
      IF (F(I,J-1).LT.EMF) NFSB=NFSB+8
      IF (NFSB.EQ.0) GO TO 450 
      IF (NFSB.GT.8) GO TO 240
      GO TO (250,260,270,280,290,300,310,320), NFSB
  240 NFSB1=NFSB-8
      GO TO (330,340,350,360,370,380,390), NFSB1
  250 U(I,J)=(U(I-1,J)-DELX(I)*RDY(J)*(V(I,J)-V(I,J-1)))*(1.D0-CYL)+CYL*
     1 (U(I-1,J)*XRM*RXRP-RDY(J)*RXRP*(V(I,J)-V(I,J-1)))
      GO TO 410
  260 V(I,J)=(V(I,J-1)-DELY(J)*RDX(I)*(U(I,J)-U(I-1,J)))*(1.D0-CYL)+CYL*
     1 (V(I,J-1)-DELY(J)*(XRP*U(I,J)-XRM*U(I-1,J)))
      GO TO 410
  270 U(I,J)=U(I-1,J)*(1.D0-CYL)+CYL*U(I-1,J)
      GO TO 260
  280 U(I-1,J)=(U(I,J)+DELX(I)*RDY(J)*(V(I,J)-V(I,J-1)))*(1.D0-CYL)+CYL*
     1 (U(I,J)*XRP*RXRM+RDY(J)*RXRM*(V(I,J)-V(I,J-1)))
      GO TO 410
  290 U(I-1,J)=U(I-1,J-1)
      GO TO 250
  300 U(I-1,J)=U(I,J)*(1.D0-CYL)+CYL*U(I,J)
      GO TO 260
  310 U(I-1,J)=U(I-1,J-1)
      U(I,J)=U(I,J-1)
      GO TO 260
  320 V(I,J-1)=(V(I,J)+DELY(J)*RDX(I)*(U(I,J)-U(I-1,J)))*(1.D0-CYL)+CYL*
     1 (V(I,J)+DELY(J)*(XRP*U(I,J)-XRM*U(I-1,J)))
      GO TO 410
  330 U(I,J)=U(I-1,J)*(1.D0-CYL)+CYL*U(I-1,J)
      GO TO 320
  340 V(I,J)=V(I-1,J)
      GO TO 320
  350 V(I,J)=V(I-1,J)
      V(I,J-1)=V(I-1,J-1)
      GO TO 250
  360 U(I-1,J)=U(I,J)*(1.D0-CYL)+CYL*U(I,J)
      GO TO 320
  370 U(I,J)=U(I,J+1)
      U(I-1,J)=U(I-1,J+1)
      GO TO 320
  380 V(I,J)=V(I+1,J)
      V(I,J-1)=V(I+1,J-1)
      GO TO 280
  390 U(I,J)=U(I-1,J)*(1.D0-CYL)+CYL*U(I-1,J)*XRM*RXRP
      V(I,J-1)=V(I,J)
      V(I,J+1)=V(I,J)
      GO TO 410
C
C     *** SET VELOCITIES IN EMPTY CELLS ADJACENT TO PARTIAL FLUID CELLS
C
  410 CONTINUE
      IF (FLG.GT.0.5D0.AND.ITER.GT.0) GO TO 450 
      IF (F(I+1,J).GT.EMF) GO TO 420
      IF (F(I+1,J+1).LT.EMF) V(I+1,J)=V(I,J)
      IF (F(I+1,J-1).LT.EMF) V(I+1,J-1)=V(I,J-1)
  420 IF (F(I,J+1).GT.EMF) GO TO 430
      IF (F(I+1,J+1).LT.EMF) U(I,J+1)=U(I,J)
      IF (F(I-1,J+1).LT.EMF) U(I-1,J+1)=U(I-1,J)
  430 IF (F(I-1,J).GT.EMF) GO TO 440
      IF (F(I-1,J+1).LT.EMF) V(I-1,J)=V(I,J)
      IF (F(I-1,J-1).LT.EMF) V(I-1,J-1)=V(I,J-1)
  440 IF (F(I,J-1).GT.EMF) CYCLE 
      IF (F(I+1,J-1).LT.EMF) U(I,J-1)=U(I,J)
      IF (F(I-1,J-1).LT.EMF) U(I-1,J-1)=U(I-1,J)
  450 CONTINUE
C
C     *** SPECIAL VELOCITY BOUNDARY CONDITIONS
C     *** Set up an inflow boundary for the jet along the top surface.
C
CRWD      IOM1 = IOMIN - 1
CRWD      DO 460 I=1,IOM1
CRWD      V(I,JM1) = VI
CRWD      V(I,JMAX) = VI
CRWD      U(I,JM1) = UI
CRWD      U(I,JMAX) = UI
CRWD      F(I,JM1) = 1.D0
CRWD      F(I,JMAX) = F(I,JM1)
CRWD  460 CONTINUE
      RETURN
C
      END
C
      SUBROUTINE CAVOVO
C
C
      INCLUDE 'dint.h'
C
      COMMON /CONST/ EMF, EMF1, EM6, EM10, EP10, PI, RPD
C
C
C     *** CALCULATE VOID VOLUMES
C
C     *** INITIALIZE VOID VOLUMES
C
      DO K=1,NVRM
       VOL(K)=0.D0
      END DO
C
C     *** COMPUTE VOID REGION VOLUMES
C
      DO 30 J=2,JM1
      DO 30 I=2,IM1
       INF=NF(I,J)
       IF (INF.EQ.0.OR.BETA(I,J).LT.0.D0) GO TO 30 
       VOLA=(1.D0-F(I,J))*DELX(I)*DELY(J)*(1.D0-CYL+CYL*2.D0*PI*XI(I))
       IF (INF.LE.5) THEN 
        INFR=NF(I+1,J)
        INFT=NF(I,J+1)
        INFL=NF(I-1,J)
        INFB=NF(I,J-1)
        INF=MAX0(INFR,INFT,INFL,INFB)
       END IF
       VOL(INF)=VOL(INF)+VOLA
   30 CONTINUE
      RETURN
      END
C
C     ******************************************************************
C     ******************************************************************
C
      SUBROUTINE CYLIND (R,H,ILOW,IHIGH)
C
C
      INCLUDE 'dint.h'
C
      REAL*4 RANDUM
      INTEGER SEED
C
      COMMON /FL1/ NDIS,IDIS,IFL,IFR,JFB,JFT
C
      COMMON /OBS/ NOBS, IOMIN(NTOB), IOMAX(NTOB), JOMIN(NTOB),
     1 JOMAX(NTOB)
C
      COMMON /CONST/ EMF, EMF1, EM6, EM10, EP10, PI, RPD
C
      real*8 AJN, EIG
      DIMENSION AJN(3), EIG(0:100,1:100)
C
C     *** Compute disturbances according to the cylindrical Rayleigh-
C     *** Taylor problem of Drazin and Reid (pg.31, Hydrodynamic Stability).
C
      ILOWP1 = ILOW + 1
C
C     Compute the eigenvalues of the disturbances.
C
      SEED = 1
      CALL EIGS(ILOW,ILOWP1,IHIGH,EIG)
      AHIGH = FLOAT(IHIGH - ILOW + 1)
      RANDUM = RAND(SEED)
      STHETA = DBLE(RANDUM)
      THETA = 2.D0*PI*STHETA
      IOBS1 = IOMIN(1) - 1
      DO 10 II=1,IOBS1
      XR = X(II)
      DO 20 JJ=1,JMAX
      IF (BETA(II,JJ) .LT. 0.D0) GO TO 20
      USUM = 0.D0
      VSUM = 0.D0
      Z = Y(JJ) - H
      R4R = 4.*R
      IF( DABS(Z) .GT. R4R) GO TO 20
      LP = 1
      IF (Z .LT. 0.D0) LP = 2
      Z =  DABS(Z)
        DO 30 IN=ILOW,IHIGH
        AN = FLOAT (IN)
        CON1 = DCOS (AN*THETA)
         DO 40 IM=ILOWP1,IHIGH
         RANDUM = RAND (0)
         AMP = DBLE(RANDUM)
         AMP = (AMP - .5D0)*2.D-3
         AR = R*AMP/AHIGH
         EV = EIG(IN,IM)/R
         RK = EV*XR
         IF (IN-1 .GE. 0) THEN
c         CALL  DBESJ (RK,IN-1,3,AJN,NZ)
         AJNP = (AJN(1) - AJN(3))/2.D0
         ELSE
c         CALL  DBESJ (RK,0,2,AJN,NZ)
         AJNP = -AJN(2)
         AJN(2) = AJN(1)
         END IF
         USUM = USUM + (-1.D0)**LP*AR*EV*AJNP*CON1*DEXP(-EV*Z)
         VSUM = VSUM + AR*EV*AJN(2)*CON1*DEXP(-EV*Z)
  40     CONTINUE
  30    CONTINUE
      IF (DABS(USUM).GE.1.D-10) THEN
        U(II,JJ) = U(II,JJ) + USUM
      END IF
      IF (DABS(VSUM).GE.1.D-10) THEN
        V(II,JJ) = V(II,JJ) + VSUM
      END IF
  20  CONTINUE
  10  CONTINUE
      RETURN
C
      END
C
C
C*********************************************************************
C*********************************************************************
C
C
      SUBROUTINE DATAPLT (NGRF)
C
C
      INCLUDE 'dint.h'

C      CHARACTER(6) PREFIX
      CHARACTER(20) fileName
C      CHARACTER(12) icT

C
      COMMON /OBS/ NOBS, IOMIN(NTOB), IOMAX(NTOB), JOMIN(NTOB),
     1 JOMAX(NTOB)
C
      COMMON /CONST/ EMF, EMF1, EM6, EM10, EP10, PI, RPD

      ngrf = 3
      open (unit=3,file=fileName,status='old',position='append')
c      open (unit=3,file='int300.out',status='old',position='append')
C
C      *** 
C      Write out the necessary data to complete the plots for
C      *** this run of SOLAVOF. Use an unformatted WRITE.
C
        IF (iCYCLE .EQ. 0) THEN
        WRITE (NGRF,100) IMAX, IM1, IMOVY, ISYMPLT, JMAX, JM1,
     1                   NMAT, NP, NOBS, ICYL, NONUNIF
        WRITE (NGRF,110) EM6, EP10, VELMX
        DO 10 I=1,IMAX
        WRITE (NGRF,120) DELX(I), X(I), XI(I), XP(I)
   10   CONTINUE
        DO 20 J=1,JMAX
        WRITE (NGRF,120) DELY(J), Y(J), YJ(J), YP(J)
   20   CONTINUE
        DO 25 I=1,NOBS
        WRITE (NGRF,105) IOMIN(I),IOMAX(I),JOMIN(I),JOMAX(I)
   25   CONTINUE
        END IF

        WRITE (NGRF,150) iCYCLE, T
        DO 5 I=1,NP
        WRITE (NGRF,101) XP(I), YP(I)
    5   CONTINUE

        IND = JMAX/5
        IDIF = JMAX - IND*5
        IF (IDIF .NE. 0) IND = IND + 1
        DO 50 K=1,IND
        JLO = (K-1)*5 + 1
        JHI = JLO + 4
        IF (JHI .GT. JMAX) JHI = JMAX
        DO 60 I=1,IMAX
        WRITE (NGRF,130) (F(I,J),J=JLO,JHI)
   60  CONTINUE
   50  CONTINUE
c        DO 70 K=1,IND
c        JLO = (K-1)*5 + 1
c        JHI = JLO + 4
c        IF (JHI .GT. JMAX) JHI = JMAX
c        DO 80 I=1,IMAX
c        WRITE (NGRF,130) (U(I,J),J=JLO,JHI)
c   80  CONTINUE
c   70  CONTINUE
c        DO 90 K=1,IND
c        JLO = (K-1)*5 + 1
c        JHI = JLO + 4
c        IF (JHI .GT. JMAX) JHI = JMAX
c        DO 95 I=1,IMAX
c        WRITE (NGRF,130) (V(I,J),J=JLO,JHI)
c   95  CONTINUE
c   90  CONTINUE

      close(unit = ngrf)
C
      RETURN
C
  100   FORMAT (11I5)
  101   FORMAT (2(E15.8,1X))
  105   FORMAT (4(2X,I5))
  110   FORMAT (3(E15.8,1X))
  120   FORMAT (4(E15.8,1X))
  130   FORMAT (5(E15.8,1X))
  150   FORMAT (I5,1X,E15.8)
      END
C
C
C   *****************************************************************
C   *****************************************************************
C
C
      SUBROUTINE DELTADJ
C
C
      INCLUDE 'dint.h'
C
      COMMON /CONST/ EMF, EMF1, EM6, EM10, EP10, PI, RPD
C
C
C     *** DELT (TIME STEP) ADJUSTMENT
C
      DELTN=DELT
      IF (FLGC.LT.0.5D0) GO TO 20
      T=T-DELT
      iCYCLE=iCYCLE-1
      DELT=DELT/2.D0
      DO 10 I=1,IMAX
      DO 10 J=1,JMAX
      P(I,J)=PN(I,J)
      F(I,J)=FN(I,J)
      U(I,J)=0.D0
      V(I,J)=0.D0
   10 CONTINUE
      NFLGC=NFLGC+1
   20 CONTINUE
      IF (AUTOT.LT.0.5D0.AND.FNOC.LT.0.5D0) GO TO 35
      DUMX=EM10
      DVMX=EM10
      IF (FNOC.GT.0.5D0) DELT=DELT/2.D0
      
C Find maximum velocity derivatives and problematic cells
      UMAX = 0.D0
      VMAX = 0.D0
      IMAX_U = 0
      JMAX_U = 0
      IMAX_V = 0
      JMAX_V = 0
      
      DO 30 I=2,IM1
      DO 30 J=2,JM1
      UDM= DABS(UN(I,J))/(XI(I+1)-XI(I))
      VDM= DABS(VN(I,J))/(YJ(J+1)-YJ(J))
      DUMX=DMAX1(DUMX,UDM)
      DVMX=DMAX1(DVMX,VDM)
      
C Track location of maximum velocities
      IF (UDM .GT. UMAX) THEN
        UMAX = UDM
        IMAX_U = I
        JMAX_U = J
      ENDIF
      IF (VDM .GT. VMAX) THEN
        VMAX = VDM
        IMAX_V = I
        JMAX_V = J
      ENDIF
   30 CONTINUE

C Calculate time step constraints BEFORE debug output
      DTMP=1.01D0
      IF (ITER.GT.25) DTMP=0.99D0
      DELTO=DELT*DTMP
      CON=0.25D0
      DELT_CFL_U = CON/DUMX
      DELT_CFL_V = CON/DVMX
      DELT_NEW = DMIN1(DELTO,DELT_CFL_U,DELT_CFL_V,DTVIS,DTSFT)
      IF (IMOVY.GT.0) DELT_NEW = DMIN1(DELT_NEW,PLTDT)

C NOW do the debug output with all variables properly defined
      IF (T .GT. 0.50D0) THEN
      WRITE(6,*) '=== DELTADJ DEBUG ==='
      WRITE(6,*) 'T=',T,' CYCLE=',iCYCLE
      WRITE(6,*) 'DUMX=',DUMX,' DVMX=',DVMX
      WRITE(6,*) 'DTVIS=',DTVIS,' DTSFT=',DTSFT
      WRITE(6,*) 'CON=',CON
      WRITE(6,*) 'DELT_CFL_U=',DELT_CFL_U,' DELT_CFL_V=',DELT_CFL_V
      WRITE(6,*) 'OLD DELT=',DELTN,' NEW DELT=',DELT_NEW
      
C Identify which constraint is most restrictive
      IF (DELT_NEW .EQ. DELT_CFL_U) THEN
        WRITE(6,*) 'TIME STEP LIMITED BY U-VELOCITY CFL'
      ELSEIF (DELT_NEW .EQ. DELT_CFL_V) THEN
        WRITE(6,*) 'TIME STEP LIMITED BY V-VELOCITY CFL'
      ELSEIF (DELT_NEW .EQ. DTVIS) THEN
        WRITE(6,*) 'TIME STEP LIMITED BY VISCOUS CONSTRAINT'
      ELSEIF (DELT_NEW .EQ. DTSFT) THEN
        WRITE(6,*) 'TIME STEP LIMITED BY SURFACE TENSION'
      ELSEIF (DELT_NEW .EQ. PLTDT) THEN
        WRITE(6,*) 'TIME STEP LIMITED BY PLOT OUTPUT'
      ELSE
        WRITE(6,*) 'TIME STEP LIMITED BY DELTO (PREVIOUS)'
      ENDIF
      
      WRITE(6,*) 'MAX U-VEL AT I=',IMAX_U,' J=',JMAX_U,' U=',
     1   UN(IMAX_U,JMAX_U)
      WRITE(6,*) 'MAX V-VEL AT I=',IMAX_V,' J=',JMAX_V,' V=',
     1   VN(IMAX_V,JMAX_V)
      WRITE(6,*) 'F AT MAX U CELL:',F(IMAX_U,JMAX_U)
      WRITE(6,*) 'F AT MAX V CELL:',F(IMAX_V,JMAX_V)
      
C Add a critical check for time step collapse
      IF (DELT_NEW .LT. 1.0D-8) THEN
        WRITE(6,*) '*** CRITICAL: TIME STEP TOO SMALL ***'
        WRITE(6,*) 'DELT would be:',DELT_NEW
        WRITE(6,*) 'This will likely cause crash!'
      ENDIF
      
      WRITE(6,*) '===================='
      ENDIF

C Actually set the new time step
C Set minimum allowable time step
      DELT_MIN = 1.0D-10  ! Adjust this value as needed
      
C Apply the calculated time step with safety check
      DELT = DELT_NEW
      
C Safety check to prevent catastrophic underflow
      IF (DELT .LT. DELT_MIN) THEN
        WRITE(6,*) '*** EMERGENCY: TIME STEP TOO SMALL ***'
        WRITE(6,*) 'Calculated DELT=',DELT
        WRITE(6,*) 'Setting to minimum:',DELT_MIN
        WRITE(6,*) 'T=',T,' CYCLE=',iCYCLE
        WRITE(6,*) 'DUMX=',DUMX,' DVMX=',DVMX
        WRITE(6,*) 'CFL_U=',CON/DUMX,' CFL_V=',CON/DVMX
        WRITE(6,*) 'DTVIS=',DTVIS,' DTSFT=',DTSFT
        
C Find and report the problematic cell
        UMAX_PROBLEM = 0.D0
        VMAX_PROBLEM = 0.D0
        DO I=2,IM1
        DO J=2,JM1
          UDM_CHK = DABS(UN(I,J))/(XI(I+1)-XI(I))
          VDM_CHK = DABS(VN(I,J))/(YJ(J+1)-YJ(J))
          IF (UDM_CHK .GT. UMAX_PROBLEM) THEN
            UMAX_PROBLEM = UDM_CHK
            I_PROBLEM = I
            J_PROBLEM = J
          ENDIF
        ENDDO
        ENDDO
        
        WRITE(6,*) 'PROBLEM CELL: I=',I_PROBLEM,' J=',J_PROBLEM
        WRITE(6,*) 'U=',UN(I_PROBLEM,J_PROBLEM)
        WRITE(6,*) 'F=',F(I_PROBLEM,J_PROBLEM)
        WRITE(6,*) 'P=',P(I_PROBLEM,J_PROBLEM)
        
C Option 1: Set to minimum and continue (may be unstable)
C        DELT = DELT_MIN
        
C Option 2: Stop gracefully (uncomment next line if preferred)
       STOP 'Time step too small - simulation terminated'
        
      ENDIF
     
   35 IF(DELT.EQ.DELTN .AND. NMAT.EQ.1) GO TO 50
      CTOS=DELT*RDTEXP
      COMG=DMIN1(CTOS**2,1.D0)
      OMG1=(OMG-1.D0)*COMG+1.D0
      DO 40 I=2,IM1
      DO 40 J=2,JM1
      IF (BETA(I,J).LT.0.D0) GO TO 40
      RHXR=(RHOFC+RHOD*F(I,J))*DELX(I+1)+(RHOFC+RHOD*F(I+1,J))*DELX(I)
      RHXL=(RHOFC+RHOD*F(I-1,J))*DELX(I)+(RHOFC+RHOD*F(I,J))*DELX(I-1)
      RHYT=(RHOFC+RHOD*F(I,J))*DELY(J+1)+(RHOFC+RHOD*F(I,J+1))*DELY(J)
      RHYB=(RHOFC+RHOD*F(I,J-1))*DELY(J)+(RHOFC+RHOD*F(I,J))*DELY(J-1)
      XX=DELT*RDX(I)*(2.D0/RHXL+2.D0/RHXR)+DELT*RDY(J)*
     ?   (2.D0/RHYT+2.D0/RHYB)
      RHOR=RHOF/(RHOFC+RHOD*F(I,J))
      BETA(I,J)=OMG1/(XX*COMG+RCSQ*RHOR/DELT)
   40 CONTINUE
   50 CONTINUE
      RETURN
      END

C
C
C   *****************************************************************
C   *****************************************************************
C
C
      SUBROUTINE DISTURB
C
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      COMMON /FL1/ NDIS,IDIS,IFL,IFR,JFB,JFT
C
C
C    *** Compute the initial disturbance velocity components.
C    *** The appropriate Rayleigh-Taylor disturbances are used.
C
C
C    *** R is the width of the disturbance region (x-direction)
C    *** H is the height above YMIN of the interface
C    *** ILOW is the lowest wave number to be used
C    *** IHIGH is the largest wave number to be used
C
      READ ( 1,110) R,H,ILOW,IHIGH
C
      IF (NDIS .EQ. 1) THEN
C
C     *** Use the disturbances for a planar problem.
C
      CALL PLANAR (R,ILOW,IHIGH)
C
      ELSE IF (NDIS .EQ. 2) THEN
C
C     *** Use the disturbances for a cylindrical problem.
C
      CALL CYLIND (R,H,ILOW,IHIGH)
C
      ELSE IF (NDIS .EQ. 3) THEN
C
C     *** Use a white noise initial velocity distribution.
C
      CALL NOISE
C
      ELSE IF (NDIS .EQ. 4) THEN
C
      CALL PLANA2(R,H,ILOW,IHIGH)
C
      ELSE IF (NDIS .EQ. 5) THEN
C
C     *** Use claude.ai (Anthropic.com) initial velocity distributions
C     ***   for finite x and y linear stability theory
C
      CALL RTINIT(H, R, IHIGH-ILOW)
C
      ELSE IF (NDIS .EQ. 6) THEN
C
C     *** Dalziel et al. (1999) interface height perturbation
C
      CALL DALZIEL_PERTURB(H, R, ILOW, IHIGH)
C
      END IF

      RETURN
C
  110 FORMAT(//,2F10.5,3I5)
      END
C
C
C     *****************************************************************
C     *****************************************************************
C
      SUBROUTINE EIGS (ILOW,ILOWP1,IHIGH,EIG)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
C     *** Compute the zeros (eigenvalues) of Jn'(km) = 0 .
C
      real*8 EIG, AJN
      DIMENSION EIG(0:100,1:100), AJN(3)
C
      PI = 4.D0* DTAN(1.D0)
      EPS = 1.E-6
      DO 10 I=0,100
      DO 10 J=1,100
      EIG(I,J) = 0.D0
  10  CONTINUE
      DO 20 N=ILOW,IHIGH
      AN = FLOAT(N)
      DO 30 M=ILOWP1,IHIGH
      IF (N .EQ. 0 .AND. M .EQ. 1) THEN
      EIG(N,M) = 0.D0
      GO TO 30
      END IF
      AM = FLOAT(M)
C
C     Estimate the zero.  Refer to Abramowitz and Stegun, Egn 9.5.13, pg 371,
C     1967.
C
      BPR = (AN/2.D0 + AM - .75D0)*PI
      EM = 4.D0*AN*AN
      CBPR = 8.D0*BPR
      CBPR3 = CBPR*CBPR*CBPR
      XP = BPR - (EM + 3.D0)/(8.D0*BPR) - 4.D0*(7.D0*EM*EM +
     1     82.D0*EM - 9.D0)/(3.D0*CBPR3)
C
C     Begin iteration procedure.
C
      IT = 0
  40  IT = IT + 1
      IF (IT .LE. 5000) THEN
      DO 41 LL=1,3
      AJN(LL) = 0.D0
   41 CONTINUE
      NZ = 0
c      CALL  DBESJ (XP,AN,3,AJN,NZ)
      TOP = XP*(AJN(1)*AN - XP*AJN(2))
      BOT = AN*(AN-1.D0)*AJN(1) - XP*(2.D0*AN+1.D0)*AJN(2) +
     1      XP*XP*AJN(3)
      DX = TOP/BOT
      XP = XP - DX
      ERR =  DABS(DX)
C      WRITE (2,110) N,M,IT,DX,XP
      IF (ERR .GT. EPS) GO TO 40
      EIG(N,M) = XP
      ELSE
ccccc      WRITE (2,100) DX,N,M
      STOP
      END IF
  30  CONTINUE
  20  CONTINUE
      RETURN
C
C  100 FORMAT (5X,'Root not found in 5000 iterations:',/,
C     1 10X,'DX = ',G14.7,' N,M = ',I5,',',I5)
C  110 FORMAT (1X,'N,M : ',I3,',',I3,' IT = ',I3,' DX = ',G14.7,
C     1 ' XP = ',G14.7)
C
      END
C
C
C*******************************************************************************
C*******************************************************************************
C
C
      SUBROUTINE LAVORE
C
      INCLUDE 'dint.h'
C
      COMMON /CONST/ EMF, EMF1, EM6, EM10, EP10, PI, RPD
C
C     *** LABEL VOID REGIONS - - VOID REGIONS ARE NF.EQ.6 AND ABOVE
C
      NNR=6
      NVR=6
      DO 30 J=2,JM1
      DO 30 I=2,IM1
      IF (NF(I,J).LT.6) GO TO 30
      INFB=NF(I,J-1)
      INFL=NF(I-1,J)
      IF (INFB.LT.6.AND.INFL.LT.6) GO TO 20
      IF (INFB.LT.6.OR.INFL.LT.6) GO TO 10
      NF(I,J)=MIN0(INFB,INFL)
      INRB=NR(INFB)
      INRL=NR(INFL)
      INRMN=MIN0(INRB,INRL)
      NR(INFB)=INRMN
      NR(INFL)=INRMN
      GO TO 30
   10 NF(I,J)=INFB
      IF (INFB.LT.6) NF(I,J)=INFL
      GO TO 30
   20 NF(I,J)=NVR
      NR(NVR)=NNR
      NVR=NVR+1
      NNR=NNR+1
   30 CONTINUE
C
C     *** REDEFINE REGION NUMBERS TO BE CONSECUTIVE
C
      NVR1=NVR-1
      NNR1=NNR-1
      KKN=7
      DO 50 KK=7,NNR1
      KFLG=0
      DO 40 K=7,NVR1
      IF (NR(K).NE.KK) GO TO 40
      NR(K)=KKN
      KFLG=1
   40 CONTINUE
      IF (KFLG.EQ.1) KKN=KKN+1
   50 CONTINUE
      NREG=KKN-6
C
C     *** REDEFINE VOID NUMBERS TO BE CONSECUTIVE IF NREG.GT.1
C
      DO 60 J=2,JM1
      DO 60 I=2,IM1
      INF=NF(I,J)
      IF (INF.LT.6) GO TO 60
      NF(I,J)=NR(INF)
   60 CONTINUE
      RETURN
      END
C
      SUBROUTINE MESHSET
C
      INCLUDE 'dint.h'
C
      COMMON /CONST/ EMF, EMF1, EM6, EM10, EP10, PI, RPD
      CHARACTER(30) fileName
      CHARACTER(11) PREFIX
      CHARACTER(12) cT1
      CHARACTER(12) cT2
      INTEGER nunit
C
C     *** MESH SETUP  (GENERATION)
C
      I=1
      J=1

      X(1)=XL(1)
      Y(1)=YL(1)

      DO 30 K=1,NKX
      DXML=(XC(K)-XL(K))/NXL(K)
      DXMR=(XL(K+1)-XC(K))/NXR(K)
      DXMN1=DXMN(K)
      NT=NXL(K)
      TN=FLOAT(NT)
      TN=DMAX1(TN,1.D0+EM6)
      DXMN(K)=DMIN1(DXMN1,DXML)
      CMC=(XC(K)-XL(K)-TN*DXMN(K))*TN/(TN-1.D0)
      IF (NT.EQ.1) CMC=0.D0
      BMC=XC(K)-XL(K)-CMC
      DO 10 L=1,NT
      I=I+1
      RLN=(FLOAT(L)-TN)/TN
   10 X(I)=XC(K)+BMC*RLN-CMC*RLN*RLN
      NT=NXR(K)
      TN=FLOAT(NT)
      TN=DMAX1(TN,1.D0+EM6)
      DXMN(K)=DMIN1(DXMN1,DXMR)
      CMC=(XL(K+1)-XC(K)-TN*DXMN(K))*TN/(TN-1.D0)
      IF (NT.EQ.1) CMC=0.D0
      BMC=XL(K+1)-XC(K)-CMC
      DO 20 L=1,NT
      I=I+1
      RLN=FLOAT(L)/TN
   20 X(I)=XC(K)+BMC*RLN+CMC*RLN*RLN
   30 CONTINUE

      IF (IWR.NE.4) GO TO 40
      I=I+1
      X(I)=X(I-1)+X(2)-X(1)
   40 CONTINUE

      DO 70 K=1,NKY
      DYML=(YC(K)-YL(K))/NYL(K)
      DYMR=(YL(K+1)-YC(K))/NYR(K)
      DYMN1=DYMN(K)
      NT=NYL(K)
      TN=FLOAT(NT)
      TN=DMAX1(TN,1.D0+EM6)
      DYMN(K)=DMIN1(DYMN1,DYML)
      CMC=(YC(K)-YL(K)-TN*DYMN(K))*TN/(TN-1.D0)
      IF (NT.EQ.1) CMC=0.D0
      BMC=YC(K)-YL(K)-CMC
      DO 50 L=1,NT
      J=J+1
      RLN=(FLOAT(L)-TN)/TN
   50 Y(J)=YC(K)+BMC*RLN-CMC*RLN*RLN
      NT=NYR(K)
      TN=FLOAT(NT)
      TN=DMAX1(TN,1.D0+EM6)
      DYMN(K)=DMIN1(DYMN1,DYMR)
      CMC=(YL(K+1)-YC(K)-TN*DYMN(K))*TN/(TN-1.D0)
      IF (NT.EQ.1) CMC=0.D0
      BMC=YL(K+1)-YC(K)-CMC
      DO 60 L=1,NT
      J=J+1
      RLN=FLOAT(L)/TN
   60 Y(J)=YC(K)+BMC*RLN+CMC*RLN*RLN
   70 CONTINUE
      IF (IWT.NE.4) GO TO 80
      J=J+1
      Y(J)=Y(J-1)+Y(2)-Y(1)
   80 CONTINUE

C      print '("I = NUMX: ",I3," J = NUMJ: ",I3)', I, J
      NUMX=I
      NUMY=J

      NUMXM1=NUMX-1
      NUMYM1=NUMY-1
      NUMXP1=NUMX+1
      NUMYP1=NUMY+1

      IF(IWR .EQ. 4) THEN
        IBAR = NUMX -2
        IMAX = IBAR + 3
      ELSE
        IBAR = NUMX - 1
        IMAX = IBAR + 2
      END IF

      IF(IWT .EQ. 4) THEN
        JBAR = NUMY - 2
        JMAX = JBAR + 3
      ELSE
        JBAR = NUMY - 1
        JMAX = JBAR + 2
      END IF

C      print '("IBAR: ",I3," JBAR: ",I3)', IBAR, JBAR
C      print '("IMAX: ",I3," JMAX: ",I3)', IMAX, JMAX
      IM1=IMAX-1
      JM1=JMAX-1
      IM2=IMAX-2
      JM2=JMAX-2
C      print '("IM1: ",I3," JM1: ",I3)', IM1, JM1
C      print '("IM2: ",I3," JM2: ",I3)', IM2, JM2
C
C     *** CALCULATE VALUES NEEDED FOR VARIABLE MESH
C
      DO 100 I=1,NUMX
      IF (X(I).EQ.0.D0) GO TO 90
      RX(I)=1.D0/X(I)
      GO TO 100
   90 RX(I)=0.D0
  100 CONTINUE
      DO 110 I=2,NUMX
      XI(I)=(X(I-1)+X(I))/2.D0
      DELX(I)=X(I)-X(I-1)
      RXI(I)=1.D0/XI(I)
  110 RDX(I)=1.D0/DELX(I)
      DELX(1)=DELX(2)
      XI(1)=XI(2)-DELX(2)
      RXI(1)=1.D0/XI(1)
      RDX(1)=1.D0/DELX(1)
      DELXA=DELX(NUMX)
      IF(IWR.EQ.4) DELXA=DELX(3)
      DELX(NUMXP1)=DELXA
      XI(NUMXP1)=XI(NUMX)+DELXA
      X(NUMXP1)=XI(NUMXP1)+DELX(NUMXP1)/2.D0
      RXI(NUMXP1)=1.D0/XI(NUMXP1)
      RDX(NUMXP1)=1.D0/DELX(NUMXP1)
      DO 120 I=2,NUMY
      YJ(I)=(Y(I-1)+Y(I))/2.D0
      RYJ(I)=1.D0/YJ(I)
      DELY(I)=Y(I)-Y(I-1)
      RDY(I)=1.D0/DELY(I)
  120 CONTINUE
      DELY(1)=DELY(2)
      RDY(1)=1.D0/DELY(1)
      YJ(1)=YJ(2)-DELY(2)
      RYJ(1)=1.D0/YJ(1)
      DELYA=DELY(NUMY)
      IF(IWT.EQ.4) DELYA=DELY(3)
      DELY(NUMYP1)=DELYA
      YJ(NUMYP1)=YJ(NUMY)+DELYA
      Y(NUMYP1)=YJ(NUMYP1)+DELY(NUMYP1)/2.D0
      RYJ(NUMYP1)=1.D0/YJ(NUMYP1)
      RDY(NUMYP1)=1.D0/DELY(NUMYP1)
C
C     Write out the mesh data.
C
      write(cT1,'(I12)') IBAR 
      write(cT2,'(I12)') JBAR 
      PREFIX = "Mesh"//TRIM(adjustl(cT1))//"x"//TRIM(adjustl(cT2))

      fileName = TRIM(PREFIX)//".vtk"
C      write(6,*) fileName

      OPEN(newunit=nunit, file=fileName)

      IF (IMOVY.EQ.0) GO TO 170
      WRITE ( nunit,190)
      DO I=1,NUMXP1
       IF(I .EQ. 1) WRITE ( nunit,191)
       WRITE (nunit,200) I,X(I),I,RX(I),I,DELX(I),I,RDX(I),I,XI(I),
     1   I,RXI(I)
      END DO 
      WRITE ( nunit,190)
      DO I=1,NUMYP1
        IF(I .EQ. 1) WRITE ( nunit,192)
        WRITE (nunit,210) I,Y(I),I,DELY(I),I,RDY(I),I,YJ(I),I,RYJ(I)
      END DO
  170 CONTINUE
      CLOSE(nunit)
C
C     *** TEST ARRAY SIZE
C
      IF (IMAX.LE.IBAR2.AND.JMAX.LE.JBAR2) GO TO 180
      WRITE (6,220)
C
      CALL EXIT
C
  180 CONTINUE
      RETURN
C
  190 FORMAT (///)
191   FORMAT(//,15X,/,15X,'----> X MESH <----',//)
  200 FORMAT (1X,2HX(,I2,2H)=,1PE12.5,2X,3HRX(,I2,2H)=,1PE12.5,
     1 2X,5HDELX(,I2,2H)=,1PE12.5,1X,4HRDX(,I2,2H)=,1PE12.5,2X,
     2 3HXI(,I2,2H)=,1PE12.5,2X,4HRXI(,I2,2H)=,1PE12.5)
192   FORMAT(///,15X,'---->Y MESH<----',//)
  210 FORMAT (1X,2HY(,I2,2H)=,1PE12.5,3X,5HDELY(,I2,2H)=,1PE12.5,
     1 3X,4HRDY(,I2,2H)=,1PE12.5,3X,3HYJ(,I2,2H)=,1PE12.5,3X,
     2 4HRYJ(,I2,2H)=,1PE12.5)
  220 FORMAT (41H  MESH SIZE GREATER THAN ARRAY DIMENSIONS)
      END
C
      SUBROUTINE NOISE
C
      INCLUDE 'dint.h'
C
      REAL*4 RANDUM
C
      COMMON /FL1/ NDIS,IDIS,IFL,IFR,JFB,JFT
C
C     *** Compute the initial disturbances using white noise from a uniformly
C     *** random distribution between 0. and 1.
C
      N = (NDIS-3)*2*IBAR2*JBAR2
      DO 5 I=1,N
      RANDUM = RAND(0)
 5    CONTINUE
      DO 10 II=1,IMAX
      DO 20 JJ=1,JMAX
      IF (BETA(II,JJ) .LT. 0.D0) GO TO 20
         RANDUM = RAND(0)
         UAMP = DBLE(RANDUM)
         UAMP = (UAMP - .5D0)*VORIG
         RANDUM = RAND(0)
         VAMP = DBLE(RANDUM)
         VAMP = (VAMP - .5D0)*VORIG
        IF (DABS(UAMP).GE.1.D-10) THEN
          U(II,JJ) = U(II,JJ) + UAMP 
        END IF
        IF (DABS(VAMP).GE.1.D-10) THEN
          V(II,JJ) = V(II,JJ) + VAMP
        END IF
  20   CONTINUE
  10  CONTINUE
C
        RETURN
        END
C
C
C   *****************************************************************
C   *****************************************************************
C
C
      SUBROUTINE PARMOV
C
      INCLUDE 'dint.h'
C     
      REAL*8 NORMX, NORMY
C
      COMMON /CONST/ EMF, EMF1, EM6, EM10, EP10, PI, RPD
C
C     *** MARKER PARTICLE MOVEMENT SECTION
C
      NPT=0
      NPN=0
      K=1
      KN=1
      PFLG=1.D0
      IPER=IM1
      IF(IWR.EQ.4) IPER=IM2
      JPER=JM1
      IF(IWT.EQ.4) JPER=JM2
   10 IF (NP.EQ.NPT) GO TO 150
C
C     *** CALCULATE U WEIGHTED VELOCITY OF PARTICLE
C
      I=IP(K)
      J=JP(K)
      IF (YP(K).GT.YJ(J)) GO TO 20
      HPX=X(I)-XP(K)
      HMX=DELX(I)-HPX
      HPY=YJ(J)-YP(K)
      NORMY=(DELY(J)+DELY(J-1))/2.D0
      HMY=NORMY-HPY
      UTOP=(U(I-1,J)*HPX+U(I,J)*HMX)*RDX(I)
      UBOT=(U(I-1,J-1)*HPX+U(I,J-1)*HMX)*RDX(I)
      UPART=(UTOP*HMY+UBOT*HPY)/NORMY
      GO TO 30
   20 HPX=X(I)-XP(K)
      HMX=DELX(I)-HPX
      HPY=YJ(J+1)-YP(K)
      NORMY=(DELY(J+1)+DELY(J))/2.D0
      HMY=NORMY-HPY
      UTOP=(U(I-1,J+1)*HPX+U(I,J+1)*HMX)*RDX(I)
      UBOT=(U(I-1,J)*HPX+U(I,J)*HMX)*RDX(I)
      UPART=(UTOP*HMY+UBOT*HPY)/NORMY
C
C     *** CALCULATE V WEIGHTED VELOCITY OF PARTICLE
C
   30 IF (XP(K).GT.XI(I)) GO TO 40
      NORMX=(DELX(I)+DELX(I-1))/2.D0
      RNORMX=1.D0/NORMX
      HPX=XI(I)-XP(K)
      HMX=NORMX-HPX
      HPY=Y(J)-YP(K)
      HMY=DELY(J)-HPY
      VTOP=(V(I-1,J)*HPX+V(I,J)*HMX)*RNORMX
      VBOT=(V(I-1,J-1)*HPX+V(I,J-1)*HMX)*RNORMX
      VPART=(VTOP*HMY+VBOT*HPY)*RDY(J)
      GO TO 50
   40 NORMX=(DELX(I)+DELX(I+1))/2.D0
      RNORMX=1.D0/NORMX
      HPX=XI(I+1)-XP(K)
      HMX=NORMX-HPX
      HPY=Y(J)-YP(K)
      HMY=DELY(J)-HPY
      VTOP=(V(I,J)*HPX+V(I+1,J)*HMX)*RNORMX
      VBOT=(V(I,J-1)*HPX+V(I+1,J-1)*HMX)*RNORMX
      VPART=(VTOP*HMY+VBOT*HPY)*RDY(J)
   50 XPART=XP(K)+UPART*DELT
      YPART=YP(K)+VPART*DELT
      IF (XPART.GT.X(I)) IP(KN)=IP(K)+1
      IF (XPART.LT.X(I-1)) IP(KN)=IP(K)-1
      IF (YPART.GT.Y(J)) JP(KN)=JP(K)+1
      IF (YPART.LT.Y(J-1)) JP(KN)=JP(K)-1
      XP(KN)=XPART
      YP(KN)=YPART
      IF(XP(KN).LT.X(1)) GO TO 90
      IF(YP(KN).LT.Y(1)) GO TO 100
      IF(XP(KN).GT.X(IPER)) GO TO 110
      IF(YP(KN).GT.Y(JPER)) GO TO 120
      GO TO 130
   90 IF(IWL.LE.2) GO TO 130
      IF(IWL.NE.4) GO TO 140
      XP(KN)=XP(KN)+X(IM2)-X(1)
      IP(KN)=IP(KN)+IM2-1
      GO TO 130
  100 IF(IWB.LE.2) GO TO 130
      IF(IWB.NE.4) GO TO 140
      YP(KN)=YP(KN)+Y(JM2)-Y(1)
      JP(KN)=JP(KN)+JM2-1
      GO TO 130
  110 IF(IWR.LE.2) GO TO 130
      IF(IWR.NE.4) GO TO 140
      XP(KN)=XP(KN)-X(IM2)+X(1)
      IP(KN)=IP(KN)-IM2+1
      GO TO 130
  120 IF(IWT.LE.2) GO TO 130
      IF(IWT.NE.4) GO TO 140
      YP(KN)=YP(KN)-Y(JM2)+Y(1)
      JP(KN)=JP(KN)-JM2+1
  130 KN=KN+1
      NPN=NPN+1
  140 K=K+1
      NPT=NPT+1
      PFLG=1.D0
      GO TO 10
  150 NP=NPN
      RETURN
      END
C
C
C   *****************************************************************
C   *****************************************************************
C
C
      SUBROUTINE PETACAL
C
      INCLUDE 'dint.h'
C
      COMMON /CONST/ EMF, EMF1, EM6, EM10, EP10, PI, RPD
C
C     *** DETERMINE THE PRESSURE INTERPOLATION FACTOR PETA
C     *** DETERMINE THE SURFACE TENSION PRESSURE AND
C     *** WALL ADHESION EFFECTS IN SURFACE CELLS
C
      DO 10 I=1,IMAX
      DO 10 J=1,JMAX
      NF(I,J)=0
      PS(I,J)=0.D0
   10 PETA(I,J)=1.0
      IPASS=0
      DO 150 I=2,IM1
      DO 150 J=2,JM1
      DTANTH(I,J)=EP10
      IF (BETA(I,J).LT.0.0) GO TO 150
      IF (F(I,J).LT.EMF) NF(I,J)=6
      IF (F(I,J).LT.EMF.OR.F(I,J).GT.EMF1) GO TO 150
      IF (F(I+1,J).LT.EMF) GO TO 20
      IF (F(I,J+1).LT.EMF) GO TO 20
      IF (F(I-1,J).LT.EMF) GO TO 20
      IF (F(I,J-1).LT.EMF) GO TO 20
      GO TO 150
   20 CONTINUE
C
C     *** CALCULATE THE PARTIAL DERIVATIVES OF F
C
      DXR=(DELX(I)+DELX(I+1))/2.D0
      DXL=(DELX(I)+DELX(I-1))/2.D0
      DYT=(DELY(J)+DELY(J+1))/2.D0
      DYB=(DELY(J)+DELY(J-1))/2.D0
      RXDEN=1.D0/(DXR*DXL*(DXR+DXL))
      RYDEN=1.D0/(DYT*DYB*(DYT+DYB))
      FL=F(I-1,J+1)
      IF (BETA(I-1,J+1).LT.0.D0.OR.(I.EQ.2.AND.IWL.LT.3)) FL=1.D0
      FC=F(I,J+1)
      IF (BETA(I,J+1).LT.0.D0) FC=1.D0
      FR=F(I+1,J+1)
      IF (BETA(I+1,J+1).LT.0.D0.OR.(I.EQ.IM1.AND.IWR.LT.3)) FR=1.D0
      AVFT=FL*DELX(I-1)+FC*DELX(I)+FR*DELX(I+1)
      FL=F(I-1,J-1)
      IF (BETA(I-1,J-1).LT.0.D0.OR.(I.EQ.2.AND.IWL.LT.3)) FL=1.D0
      FC=F(I,J-1)
      IF (BETA(I,J-1).LT.0.D0) FC=1.D0
      FR=F(I+1,J-1)
      IF (BETA(I+1,J-1).LT.0.D0.OR.(I.EQ.IM1.AND.IWR.LT.3)) FR=1.D0
      AVFB=FL*DELX(I-1)+FC*DELX(I)+FR*DELX(I+1)
      FL=F(I-1,J)
      IF (BETA(I-1,J).LT.0.D0.OR.(I.EQ.2.AND.IWL.LT.3)) FL=1.D0
      FR=F(I+1,J)
      IF (BETA(I+1,J).LT.0.D0.OR.(I.EQ.IM1.AND.IWR.LT.3)) FR=1.D0
      AVFCY=FL*DELX(I-1)+F(I,J)*DELX(I)+FR*DELX(I+1)
      FB=F(I,J-1)
      IF (BETA(I,J-1).LT.0.D0.OR.(J.EQ.2.AND.IWB.LT.3)) FB=1.D0
      FT=F(I,J+1)
      IF (BETA(I,J+1).LT.0.D0.OR.(J.EQ.JM1.AND.IWT.LT.3)) FT=1.D0
      AVFCX=FB*DELY(J-1)+F(I,J)*DELY(J)+FT*DELY(J+1)
      FB=F(I-1,J-1)
      IF (BETA(I-1,J-1).LT.0.D0.OR.(J.EQ.2.AND.IWB.LT.3)) FB=1.D0
      FC=F(I-1,J)
      IF (BETA(I-1,J).LT.0.D0) FC=1.D0
      FT=F(I-1,J+1)
      IF (BETA(I-1,J+1).LT.0.D0.OR.(J.EQ.JM1.AND.IWT.LT.3)) FT=1.D0
      AVFL=FB*DELY(J-1)+FC*DELY(J)+FT*DELY(J+1)
      FB=F(I+1,J-1)
      IF (BETA(I+1,J-1).LT.0.D0.OR.(J.EQ.2.AND.IWB.LT.3)) FB=1.D0
      FC=F(I+1,J)
      IF (BETA(I+1,J).LT.0.D0) FC=1.D0
      FT=F(I+1,J+1)
      IF (BETA(I+1,J+1).LT.0.D0.OR.(J.EQ.JM1.AND.IWT.LT.3)) FT=1.D0
      AVFR=FB*DELY(J-1)+FC*DELY(J)+FT*DELY(J+1)
C
C     *** BOUNDARY CONDITIONS FOR WALL ADHESION
C
      IF (ISURF10.EQ.0.D0.OR.CANGLE.EQ.0.D0) GO TO 60
      IF (BETA(I+1,J).GE.0.D0.AND.I.NE.IM1) GO TO 30
      AVFR=AVFCX+0.5D0*(DELX(I)+DELX(I+1))/DTANCA
      IF (F(I,J+1).LT.EMF.AND.F(I,J-1).GE.EMF) AVFT=AVFCY-0.5D0*(DELY(J)
     1 +DELY(J+1))*DTANCA
      IF (F(I,J-1).LT.EMF.AND.F(I,J+1).GE.EMF) AVFB=AVFCY-0.5D0*(DELY(J)
     1 +DELY(J-1))*DTANCA
   30 IF (BETA(I,J+1).GE.0.D0.AND.J.NE.JM1) GO TO 40
      AVFT=AVFCY+0.5D0*(DELY(J)+DELY(J+1))/DTANCA
      IF (F(I+1,J).LT.EMF.AND.F(I-1,J).GE.EMF) AVFR=AVFCX-0.5D0*(DELX(I)
     1 +DELX(I+1))*DTANCA
      IF (F(I-1,J).LT.EMF.AND.F(I+1,J).GE.EMF) AVFL=AVFCX-0.5D0*(DELX(I)
     1 +DELX(I-1))*DTANCA
   40 IF (BETA(I,J-1).GE.0.D0.AND.J.NE.2) GO TO 50
      AVFB=AVFCY+0.5D0*(DELY(J)+DELY(J-1))/DTANCA
      IF (F(I+1,J).LT.EMF.AND.F(I-1,J).GE.EMF) AVFR=AVFCX-0.5D0*(DELX(I)
     1 +DELX(I+1))*DTANCA
      IF (F(I-1,J).LT.EMF.AND.F(I+1,J).GE.EMF) AVFL=AVFCX-0.5D0*(DELX(I)
     1 +DELX(I-1))*DTANCA
   50 IF (BETA(I-1,J).GE.0.D0.AND.I.NE.2) GO TO 60
      IF (CYL.GT.0.5D0.AND.X(1).EQ.0.D0) GO TO 60
      AVFL=AVFCX+0.5D0*(DELX(I)+DELX(I-1))/DTANCA
      IF (F(I,J+1).LT.EMF.AND.F(I,J-1).GE.EMF) AVFT=AVFCY-0.5D0*(DELY(J)
     1 +DELY(J+1))*DTANCA
      IF (F(I,J-1).LT.EMF.AND.F(I,J+1).GE.EMF) AVFB=AVFCY-0.5D0*(DELY(J)
     1 +DELY(J-1))*DTANCA
   60 CONTINUE
      XTHM=3.D0*DMAX1(AVFT,AVFCY,AVFB)/(DELX(I-1)+DELX(I)+DELX(I+1))
      YTHM=3.D0*DMAX1(AVFL,AVFCX,AVFR)/(DELY(J-1)+DELY(J)+DELY(J+1))
      PFX=RXDEN*((AVFR-AVFCX)*DXL*DXL+(AVFCX-AVFL)*DXR*DXR)
      PFY=RYDEN*((AVFT-AVFCY)*DYB*DYB+(AVFCY-AVFB)*DYT*DYT)
      PF=PFX*PFX+PFY*PFY
      IF (PF.GT.EM10) GO TO 70
      NF(I,J)=5
      P(I,J)=(P(I+1,J)+P(I,J+1)+P(I-1,J)+P(I,J-1))/4.D0
      GO TO 150
   70 CONTINUE
C
C     *** DETERMINE THE PRESSURE INTERPOLATION CELL NF
C
      ABPFX= DABS(PFX)
      ABPFY= DABS(PFY)
      L=I
      M=J
      IF (ABPFY.GE.ABPFX) GO TO 80
      DXDYR=DELY(J)*RDX(I)
      PFMN=PFY
      NF(I,J)=2
      L=I+1
      DMX=DELX(I)
      DMIN=(DMX+DELX(I+1))/2.D0
      IF (PFX.GT.0.D0) GO TO 90
      NF(I,J)=1
      PFMN=-PFY
      L=I-1
      DMX=DELX(I)
      DMIN=(DMX+DELX(I-1))/2.D0
      GO TO 90
   80 CONTINUE
      DXDYR=DELX(I)*RDY(J)
      PFMN=-PFX
      NF(I,J)=4
      M=J+1
      DMX=DELY(J)
      DMIN=(DMX+DELY(J+1))/2.D0
      IF (PFY.GT.0.D0) GO TO 90
      NF(I,J)=3
      PFMN=PFX
      M=J-1
      DMX=DELY(J)
      DMIN=(DMX+DELY(J-1))/2.D0
   90 CONTINUE
      DTANTH(I,J)=PFMN
      ABDTAN= DABS(DTANTH(I,J))
C
C     *** DETERMINE THE CURVATURE AND SURFACE PRESSURE
C
      DFS=(0.5D0-F(I,J))*DMX
      IF (F(I,J).LT.0.5D0*ABDTAN*DXDYR) DFS=DMX*(1.D0+DXDYR*ABDTAN-
     1   DSQRT(8.D0*F(I,J)*DXDYR*ABDTAN))/2.D0
      IF (ISURF10.LT.1) GO TO 140
      NFC=NF(I,J)
      PXR=(AVFR-AVFCX)/DXR
      PXL=(AVFCX-AVFL)/DXL
      PYT=(AVFT-AVFCY)/DYT
      PYB=(AVFCY-AVFB)/DYB
      YDFS=-DFS
      IF(NFC.EQ.2 .OR. NFC.EQ.4) YDFS=DFS
      IF(NFC.GT.2) GO TO 100
      DXDN=DELY(J)
      XINB=YDFS+DTANTH(I,J)*DXDN/2.D0
      XINT=2.D0*YDFS-XINB
      GP1=PYT
      PX1=PXL
      IF(XINT.GT.0.D0) PX1=PXR
      IF(DABS(PX1).LT.DABS(GP1))GP1=DSIGN(1.D0,GP1)/( DABS(PX1)+EM10)
      GP2=PYB
      PX2=PXR
      IF(XINB.LT.0.D0) PX2=PXL
      IF( DABS(PX2).LT.DABS(GP2))GP2=DSIGN(1.D0,GP2)/(DABS(PX2)+EM10)
      GO TO 110
  100 DXDN=DELX(I)
      YINR=YDFS+DTANTH(I,J)*DXDN/2.D0
      YINL=2.D0*YDFS-YINR
      GP1=PXR
      PY1=PYT
      IF(YINR.LT.0.D0) PY1=PYB
      IF( DABS(PY1).LT.DABS(GP1))GP1=DSIGN(1.D0,GP1)/(DABS(PY1)+EM10)
      GP2=PXL
      PY2=PYB
      IF(YINL.GT.0.D0) PY2=PYT
      IF( DABS(PY2).LT.DABS(GP2))GP2=DSIGN(1.D0,GP2)/(DABS(PY2)+EM10)
  110 GP1D=1.D0+GP1*GP1
      GP2D=1.D0+GP2*GP2
      CURVXY=(GP2/DSQRT(GP2D)-GP1/DSQRT(GP1D))/DXDN
      CURVCYL=0.D0
      IF (CYL.LT.1.D0) GO TO 120
      XLITLR=XI(I)
      IF (NFC.EQ.1) XLITLR=X(I-1)+F(I,J)*DELX(I)
      IF (NFC.EQ.2) XLITLR=X(I)-F(I,J)*DELX(I)
      RLITLR=DMIN1(1.D0/XLITLR,RXI(2))
      TRIG= DABS(DSIN(DTAN(ABDTAN)))
      IF (NFC.LE.2) TRIG= DABS(DCOS(DTAN(ABDTAN)))
      CURVCYL=-CYL*TRIG* DSIGN(1.D0,PFX)*RLITLR
  120 CURV=CURVXY+CURVCYL
      PS(I,J)=SIGMA*CURV
      IF (XTHM.LT.1.D0.OR.YTHM.LT.1.D0) PS(I,J)=0.D0
  140 CONTINUE
C
C     *** CALCULATE PETA
C
      NFSB=0
      IF (F(I+1,J).LT.EMF.OR.I.EQ.IM1.OR.BETA(I+1,J).LT.0.D0) NFSB=
     1    NFSB+1
      IF (F(I,J+1).LT.EMF.OR.BETA(I,J+1).LT.0.D0) NFSB=NFSB+2
      IF (F(I-1,J).LT.EMF.OR.BETA(I-1,J).LT.0.D0) NFSB=NFSB+4
      IF (F(I,J-1).LT.EMF.OR.BETA(I,J-1).LT.0.D0) NFSB=NFSB+8
      IF (NFSB.EQ.15) PS(I,J)=0.D0
      IF (NMAT.EQ.2) GO TO 150
      PETA(I,J)=1.D0/(1.D0-DFS/DMIN)
      IF (L.EQ.1.OR.L.EQ.IMAX) PETA(I,J)=1.D0
      IF (M.EQ.1.OR.M.EQ.JMAX) PETA(I,J)=1.D0
      IF (BETA(L,M).LT.0.D0) PETA(I,J)=1.D0
  150 CONTINUE
C
      CALL LAVORE
C
      CALL CAVOVO
C
C     IF NECESSARY, DETERMINE PRESSURES PR FOR VOID REGIONS NF
C
      IF (NMAT.EQ.2) GO TO 300
C
C     *** SET PETA IN ADJACENT FULL CELL
C
      DO 290 J=1,JMAX
      DO 290 I=1,IMAX
      NFF=NF(I,J)
      IF (NFF.EQ.0.OR.BETA(I,J).LT.0.D0) GO TO 290
      IF (NFF.GT.5) GO TO 280
      L=I
      M=J
      GO TO (230,240,250,260,290), NFF
  230 L=I-1
      DMX=DELX(L)
      DMIN=(DMX+DELX(I))/2.D0
      GO TO 270
  240 L=I+1
      DMX=DELX(L)
      DMIN=(DMX+DELX(I))/2.D0
      GO TO 270
  250 M=J-1
      DMX=DELY(M)
      DMIN=(DMX+DELY(J))/2.D0
      GO TO 270
  260 M=J+1
      DMX=DELY(M)
      DMIN=(DMX+DELY(J))/2.D0
  270 CONTINUE
      IF (NF(L,M).GT.0) GO TO 290
      CTOS=DELT*RDTEXP
      COMG=DMIN1(CTOS**2,1.D0)
      BPD=1.D0/PETA(L,M)-BETA(L,M)*(1.D0-PETA(I,J))
     1*DELT/(DMIN*DMX)*(COMG/RHOF)
      PETA(L,M)= 1.D0/BPD
      GO TO 290
  280 CONTINUE
      P(I,J)=PR(NFF)
  290 CONTINUE
  300 CONTINUE
      RETURN
      END
C
C     ***************************************************************
C     ***************************************************************
C
      SUBROUTINE PLANA2 (R,H,ILOW,IHIGH)
C
C
      INCLUDE 'dint.h'
C
      REAL*4 RANDUM
C
      COMMON /CONST/ EMF, EMF1, EM6, EM10, EP10, PI, RPD
C
C     *** Compute the disturbances according to the classical Rayleigh-Taylor
C     *** stability problem outlined in Youngs, Physica, 12D (1984) 32-44
C     *** for a planer geometry.
C
C
      DO I = 1,IMAX
        DO J = 1,JMAX
          V(I,J) = 0
          U(I,J) = 0
        END DO
      END DO

      G = DSQRT(GX*GX + GY*GY)
      S = RHOF/RHOFC
      AHIGH = FLOAT(IHIGH-ILOW + 1)
      CON1 = PI/R
      CON2 = (S - 1.D0)/(S + 1.D0)
      CON3 = CON1*CON1*SIGMA/(G*(RHOF + RHOFC))
      DO II=1,IMAX
        XR = X(II)
        DO JJ=1,JMAX
          IF (BETA(II,JJ) .LT. 0.D0) THEN
            CYCLE
          END IF

          USUM = 0.D0
          VSUM = 0.D0
          Z = Y(JJ) - H
          IF (Z .EQ. 0.) THEN
            Z = 1.D-8
          END IF

          DO I=ILOW,IHIGH
            RANDUM = RAND(0)
            AMP = ((RANDUM - .5D0)*2.D-3)/AHIGH
            AR = R*AMP
            AI = FLOAT (I)
            AK = AI*PI/H
            USUM = USUM + AR*DCOS(AK*XR)*DEXP(-AK*Z)*Z/ABS(Z)
            VSUM = VSUM + AR*AK*DSIN(AK*XR)*DEXP(-AK*Z)
          END DO

          VSUM = - VSUM
          IF (DABS(USUM).GE.1.D-10) THEN
            U(II,JJ) = USUM
          END IF
          IF (DABS(VSUM).GE.1.D-10) THEN
            V(II,JJ) = - VSUM
          END IF
        END DO
      END DO
      RETURN
        END
C
C
C   *****************************************************************
C   *****************************************************************
      SUBROUTINE PLANAR (R,ILOW,IHIGH)
C
C
      INCLUDE 'dint.h'
C
      REAL*4 RANDUM
C
      COMMON /CONST/ EMF, EMF1, EM6, EM10, EP10, PI, RPD
C
C     *** Compute the disturbances according to the classical Rayleigh-Taylor
C     *** stability problem outlined in Youngs, Physica, 12D (1984) 32-44
C     *** for a planer geometry.
C
C
C      WRITE(6,*) 6HILOW: , ILOW, 7HIHIGH: , IHIGH, 3HR: , R

      PI = 3.141592654
      
      DO I = 1,IMAX
        DO J = 1,JMAX
          V(I,J) = 0.D0
          U(I,J) = 0.D0
        END DO
      END DO
C
C     *** Accumulate modes with random coefficients in [-1,+1].
C     *** VORIG scales the total sum so it controls max amplitude.
C
      DO I=ILOW,IHIGH
         RANDUM = 2.D0*RAND(0) - 1.D0
         DO II = 1,IMAX
           POS = 2.D0*DBLE(I)*PI*XI(II)/R
           VSUM = RANDUM*COS(POS)
           DO JJ = 1,JM1
             V(II,JJ) = V(II,JJ) + VSUM
           END DO
         END DO
      END DO
C
C     *** Scale by VORIG
C
      DO II = 1,IMAX
        DO JJ = 1,JM1
          V(II,JJ) = VORIG*V(II,JJ)
        END DO
      END DO
      RETURN
      END
C
C
C   *****************************************************************
C   *****************************************************************
C
C
      SUBROUTINE PRESSIT
C
      INCLUDE 'dint.h'
C
      COMMON /FL1/ NDIS,IDIS,IFL,IFR,JFB,JFT
C
      COMMON /CONST/ EMF, EMF1, EM6, EM10, EP10, PI, RPD
C
C     *** PRESSURE ITERATION
C
C     *** TEST FOR CONVERGENCE
C
   10 IF (FLG.EQ.0.) GO TO 140
      ITER=ITER+1
      ITMAX=1000
      IF (ITER.LT.ITMAX) GO TO 20
      FNOC=1.D0
      NOCON=NOCON+1
      GO TO 140
   20 FLG=0.D0
C
C     *** COMPUTE UPDATED CELL PRESSURE AND VELOCITIES
C
      DO 130 J=JPB,JPT
      DO 130 I=IPL,IPR
      IF (BETA(I,J).LT.0.D0) GO TO 130
      IF (NMAT.EQ.2) GO TO 80
      IF (F(I,J).LT.EMF) GO TO 130
      IF (NF(I,J).EQ.0) GO TO 80
C
C     *** CALCULATE PRESSURE FOR SURFACE CELLS
C
      NFF=NF(I,J)
      L=I
      M=J
      GO TO (30,40,50,60,130), NFF
   30 L=I-1
      GO TO 70
   40 L=I+1
      GO TO 70
   50 M=J-1
      GO TO 70
   60 M=J+1
   70 CONTINUE
      NFEL=NF(I-1,J)
      NFER=NF(I+1,J)
      NFEB=NF(I,J-1)
      NFET=NF(I,J+1)
      NFE=MAX0(NFEL,NFER,NFEB,NFET)
      PSURF=PS(I,J)+PR(NFE)
      PLM=P(L,M)
      IF (NF(L,M).NE.0.AND.BETA(I,J).GT.0.D0) PLM=PSURF
      DELP=(1.D0-PETA(I,J))*PLM+PETA(I,J)*PSURF-P(I,J)
      GO TO 90
   80 CONTINUE
      DIJ=RDX(I)*(U(I,J)-U(I-1,J))+RDY(J)*(V(I,J)-V(I,J-1))+CYL*RXI
     1 (I)*(U(I,J)+U(I-1,J))/2.D0
      RHOR=RHOF/(RHOFC+RHOD*F(I,J))
      DFUN=DIJ+RHOR*RCSQ*(P(I,J)-PN(I,J))/DELT
C
C     *** SET FLAG INDICATING CONVERGENCE
C
      IF ( DABS(DFUN).GE.EPSI) FLG=1.D0
      DELP=-BETA(I,J)*DFUN*PETA(I,J)
   90 CONTINUE
      P(I,J)=P(I,J)+DELP

C Add this at the end of the pressure update loop (after P(I,J)=P(I,J)+DELP)
      IF(T .GT. 8.28D0) THEN
      IF (ABS(P(I,J)) .GT. 1.0D10 .OR. P(I,J) .NE. P(I,J)) THEN
        WRITE(6,*) 'PRESSURE PROBLEM I=',I,' J=',J
        WRITE(6,*) 'P=',P(I,J),' DELP=',DELP,' DFUN=',DFUN
        WRITE(6,*) 'F=',F(I,J),' BETA=',BETA(I,J)
        STOP
      ENDIF
      ENDIF
      CTOS=DELT*RDTEXP
      COMG=DMIN1(CTOS**2,1.D0)
      DPTC=2.D0*DELT*DELP*COMG
      IF (BETA(I+1,J).LT.0.D0) GO TO 100
      RHOXR=(RHOFC+RHOD*F(I,J))*DELX(I+1)+(RHOFC+RHOD*F(I+1,J))*DELX(I)
      U(I,J)=U(I,J)+DPTC/RHOXR
  100 IF (BETA(I-1,J).LT.0.D0) GO TO 110
      RHOXL=(RHOFC+RHOD*F(I-1,J))*DELX(I)+(RHOFC+RHOD*F(I,J))*DELX(I-1)
      U(I-1,J)=U(I-1,J)-DPTC/RHOXL
  110 IF (BETA(I,J+1).LT.0.D0) GO TO 120
      RHOYT=(RHOFC+RHOD*F(I,J))*DELY(J+1)+(RHOFC+RHOD*F(I,J+1))*DELY(J)
      V(I,J)=V(I,J)+DPTC/RHOYT
  120 IF (BETA(I,J-1).LT.0.D0) GO TO 130
      RHOYB=(RHOFC+RHOD*F(I,J-1))*DELY(J)+(RHOFC+RHOD*F(I,J))*DELY(J-1)
      V(I,J-1)=V(I,J-1)-DPTC/RHOYB
  130 CONTINUE
      CALL BC
      GO TO 10
  140 CONTINUE
      RETURN
      END
C
C
C   *****************************************************************
C   *****************************************************************
C
C
      SUBROUTINE PRT (N)
C
      INCLUDE 'dint.h'
C
      COMMON /FL1/ NDIS,IDIS,IFL,IFR,JFB,JFT
C
      COMMON /FUND/ FLENG,FTIME
C
      COMMON /OBS/ NOBS, IOMIN(NTOB), IOMAX(NTOB), JOMIN(NTOB),
     1 JOMAX(NTOB)
C
      COMMON /CONST/ EMF, EMF1, EM6, EM10, EP10, PI, RPD
C
C     *** PRT--
C     *** PROVIDES FORMATTED WRITES TO A DATA FILE.
C
      IF (N .EQ. 1) THEN
C
C     *** PRT (1) WRITE OUT INITIAL DATA AND MESH DATA
C
      WRITE ( 2,240)
      WRITE ( 2,215) FLENG,FTIME
      WRITE ( 2,220) IBAR,JBAR,DELT,VNU,ICYL,EPSI,GX,GY,UI,VI,VELMX,
     1 TWFIN,PRTDT,PLTDT,OMG,ALPHA,IWL,IWR,IWT,IWB,IMOVY,AUTOT,FLHT,
     2 ISYMPLT,SIGMA,ISURF10,CANGLE,CSQ,NMAT,RHOF,RHOFC,NOBS,NDIS
C
C     *** WRITE THE VARIABLE MESH INPUT DATA
C
      WRITE ( 2,260) NKX
      DO 20 I=1,NKX
      WRITE ( 2,270) I,XL(I),XC(I),XL(I+1),NXL(I),NXR(I),DXMN(I)
   20 CONTINUE
      WRITE ( 2,280) NKY
      DO 30 I=1,NKY
      WRITE ( 2,275) I,YL(I),YC(I),YL(I+1),NYL(I),NYR(I),DYMN(I)
   30 CONTINUE
      RETURN
C
      ELSE IF (N .EQ. 2) THEN
C
C     *** PRT (2)  WRITE TIME STEP, CYCLE INFORMATION
C
      WRITE ( 2,210) ITER,T,DELT,iCYCLE,VCHGT
      WRITE ( 6,210) ITER,T,DELT,iCYCLE,VCHGT
      RETURN
C
      ELSE IF (N. EQ. 3) THEN
C
C     *** PRT (3)  WRITE FIELD VARIABLES TO DATA FILE
C
      IF (IMOVY.EQ.0) RETURN
C
      WRITE ( 2,240)
      WRITE ( 2,251)
      WRITE ( 2,251)
c      WRITE ( 2,250) CNAME
      WRITE ( 2,251)
      WRITE ( 2,251)
      WRITE ( 2,210) ITER,T,DELT,iCYCLE,VCHGT
      WRITE ( 2,290) NREG
      WRITE ( 2,300)
      KNR=NREG+5
      DO 100 K=6,KNR
      WRITE ( 2,310) K,VOL(K),PR(K)
  100 CONTINUE
      WRITE ( 2,190)
      DO I=1,IMAX
      DO J=1,JMAX
      IF (I .EQ. 1 .AND. J .GT. 1) THEN
      DIJ=RDX(I)*(U(I+1,J)-U(I,J))+RDY(J)*(V(I,J)-V(I,J-1))
     1  + CYL*RXI(I)*(U(I,J)+U(I+1,J))/2.D0
      ELSE IF (I .GT. 1 .AND. J .EQ. 1) THEN
      DIJ=RDX(I)*(U(I,J)-U(I-1,J))+RDY(J)*(V(I,J+1)-V(I,J))
     1  + CYL*RXI(I)*(U(I,J)+U(I-1,J))/2.D0
      ELSE IF (I .EQ. 1 .AND. J .EQ. 1) THEN
      DIJ=RDX(I)*(U(I,J)-U(I+1,J))+RDY(J)*(V(I,J+1)-V(I,J))
     1  + CYL*RXI(I)*(U(I,J)+U(I+1,J))/2.D0
      ELSE
      DIJ=RDX(I)*(U(I,J)-U(I-1,J))+RDY(J)*(V(I,J)-V(I,J-1))
     1  + CYL*RXI(I)*(U(I,J)+U(I-1,J))/2.D0
      END IF
      WRITE (2,200) I,J,U(I,J),V(I,J),P(I,J),F(I,J),PS(I,J),DIJ,NF(I,J)
     1  ,PETA(I,J)
      END DO 
      END DO

      RETURN
C
      END IF
C
  190 FORMAT (4X,1HI,5X,1HJ,9X,1HU,14X,1HV,15X,1HP,15X,1HF,12X,2HPS,13X,
     1 1HD,11X,2HNF,9X,4HPETA)
  200 FORMAT (2X,I3,3X,I3,6(3X,1PE12.5),3X,I3,3X,E12.5)
  210 FORMAT (6X,6HITER= ,I5,5X,6HTIME= ,1PE12.5,5X,6HDELT= ,1PE12.5,5X,
     1 7HiCYCLE= ,I7,3X,7HVCHGT= ,1PE12.5)
  215 FORMAT (5X,7HFLENG= ,E12.5,/,5X,7HFTIME= ,E12.5)
  220 FORMAT (//,6X,6HIBAR= ,I3,/,6X,6HJBAR= ,I3,/,6X,6HDELT= ,1PE12.5,
     1 /,8X,4HNU= ,E12.5,/,6X,6HICYL= ,I2,/,6X,6HEPSI= ,E12.5,/,8X,
     2 4HGX= ,E12.5,/,8X,4HGY= ,E12.5,/,8X,4HUI= ,E12.5,/,8X,4HVI= ,
     3 E12.5,/,5X,7HVELMX= ,E12.5,/,5X,7HTWFIN= ,E12.5,/,5X,7HPRTDT= ,
     4 E12.5,/,5X,7HPLTDT= ,E12.5,/,7X,5HOMG= ,E12.5,/,5X,7HALPHA= ,
     5 E12.5,/,8X,4HWL= ,I2,/,8X,4HWR= ,I2,/,8X,4HWT= ,I2,/,8X,4HWB= ,
     6 I2,/,5X,7HIMOVY= ,I4,/,5X,7HAUTOT= ,E12.5,/,6X,6HFLHT= ,E12.5,/,
     7 3X,9HISYMPLT= ,I2,/,5X,7HSIGMA= ,E12.5,/,3X,9HISURF10= ,I2,/,4X,
     8 8HCANGLE= ,E12.5,/,7X,5HCSQ= ,E12.5,/,6X,6HNMAT= ,I2,/,6X,
     9 6HRHOF= ,E12.5,/,6X,7HRHOFC= ,E12.5,/,6X,'NOBS= ',I5,/,6X,
     1 'NDIS= ',I5)
  240 FORMAT (1H1)
  251 FORMAT (18X,'******************************************')
  260 FORMAT (2X,5HNKX= ,I4)
  270 FORMAT(2X,8HMESH-X= ,I4,3X,4HXL= ,1PE12.5,3X,4HXC= ,E12.5,3X,
     1 4HXR= ,E12.5,3X,5HNXL= ,I4,3X,5HNXR= ,I4,3X,6HDXMN= ,E12.5)
  275 FORMAT(2X,8HMESH-Y= ,I4,3X,4HYL= ,1PE12.5,3X,4HYC= ,E12.5,3X,
     1 4HYR= ,E12.5,3X,5HNYL= ,I4,3X,5HNYR= ,I4,3X,6HDYMN= ,E12.5)
  280 FORMAT (2X,5HNKY= ,I4)
  290 FORMAT (2X,6HNREG= ,I4)
  300 FORMAT (15X,1HK,6X,6HVOL(K),9X,5HPR(K))
  310 FORMAT (13X,I3,2X,1PE12.5,3X,E12.5)
      END

C   END OF SUBROUTINE PRT

c=====================================================================
c     RT Instability initialization subroutine
c      Impose correct wall BC's
c      With help from claude (Anthropic AI)
c     
c     Initializes velocity field for Rayleigh-Taylor instability
c     for both inviscid and viscous cases, with or without surface 
c     tension.
c     
c     Based on linear stability analysis with appropriate eigenfunctions
c     for inviscid or viscous flow
c=====================================================================

       SUBROUTINE RTINIT(HB, L, NUMMODES)
      
C      IMPLICIT NONE
      
c     Include common block definitions
      INCLUDE 'dint.h'
      
c     Input parameters
      REAL*8 HB           ! Height of bottom domain (below interface)
      REAL*8 L           ! Horizontal domain length
      INTEGER NUMMODES   ! Number of modes to include
      
c     Boundary condition types:
c     1 => rigid free-slip
c     2 => rigid no-slip  
c     3 => continuative boundary (outflow)
c     4 => periodic boundary
c     5 => constant pressure boundary
      
      LOGICAL INVISCID   ! Flag for inviscid/viscous calculation
      
c     Local variables
      INTEGER I, J, N, NMAX
      INTEGER WL         ! Left boundary condition (1-5)
      INTEGER WR         ! Right boundary condition (1-5)  
      INTEGER WB         ! Bottom boundary condition (1-5)
      INTEGER WT         ! Top boundary condition (1-5)
      REAL*8 HT           ! Height of top domain (above interface)
      REAL*8 G            ! Gravitational acceleration
      REAL*8 KC, KMAX, KX, K2
      REAL*8 Z, AN, ALPHA, A0
      REAL*8 NUEFF, KNU
      REAL*8 TWOPI, PI
      REAL*8 DRHO, PHI1, PHI2, DPHI1, DPHI2
      REAL*8 KY, DX
      REAL*8 SINH_KB, SINH_KT, COSH_KB, COSH_KT
      REAL*8 X_COORD, PHASE_X_W, PHASE_X_U
      LOGICAL LATERAL_PERIODIC, VERTICAL_NOSLIP
      
c     Set BC values
      WL = IWL
      WR = IWR
      WB = IWB
      WT = IWT

c     Set constants
      PI = 3.14159265
      TWOPI = 2.0 * PI
      KY = 0.0           ! 2D case, no y-dependence
      ALPHA = 1.5        ! Power law for mode amplitudes
      A0 = 0.01          ! Base amplitude (1% of wavelength)

C     Calculate HT
      HT = Y(JBAR+1) - Y(1) - HB
      
c     Calculate density difference
      DRHO = RHOF - RHOFC

C     Set dx, assuming all delx are equal!!!
      DX = DELX(2)

c     Calculate gravitational acceleration
      G = SQRT(GX*GX + GY*GY)

c     Calculate the effective viscosity      
      NUEFF = (RHOFC*VNUC + RHOF*VNU)/(RHOFC + RHOF)  ! Effective viscosity

      INVISCID = .FALSE.
      IF(NUEFF .LE. 1.D-6) THEN
        INVISCID = .TRUE.
      END IF

c     Check boundary condition compatibility
      LATERAL_PERIODIC = (WL .EQ. 4 .AND. WR .EQ. 4)
      VERTICAL_NOSLIP = (WB .EQ. 2 .AND. WT .EQ. 2)
      
      WRITE(*,*) 'RTINIT: Boundary conditions:'
      WRITE(*,*) '  WL=', WL, ' WR=', WR, ' WB=', WB, ' WT=', WT
      WRITE(*,*) '  Lateral periodic: ', LATERAL_PERIODIC
      WRITE(*,*) '  Vertical no-slip: ', VERTICAL_NOSLIP
 
c     Calculate wavenumber parameters
      IF (INVISCID) THEN
c        Inviscid case
         IF (SIGMA .GT. 0.0) THEN
c           With surface tension
            KC = DSQRT(G*DABS(DRHO)/SIGMA)      ! Cutoff wavenumber
            KMAX = DSQRT(G*DABS(DRHO)/(3*SIGMA))  ! Most unstable wavenumber
         ELSE
c           Without surface tension - use grid-based limits
            KC = PI/DX        ! Grid-based cutoff (Nyquist limit)
            KMAX = PI/(10*DX) ! A reasonable scale for initialization
         ENDIF
      ELSE
c        Viscous case
         IF (SIGMA .GT. 0.0) THEN
c           With surface tension
            KC = DSQRT(G*DABS(DRHO)/SIGMA)  ! First approximation
            KC = DSQRT(G*DABS(DRHO)/(SIGMA + 
     1           2*NUEFF**2*(RHOFC+RHOF)**2/KC))
            KMAX = DSQRT(G*DABS(DRHO)/(3*SIGMA + 
     1           4*NUEFF**2*(RHOFC+RHOF)**2/SIGMA))
         ELSE
c           Without surface tension
            KNU = DSQRT(G*DABS(DRHO)/(2*NUEFF*(RHOFC+RHOF)))
            KC = KNU  ! Viscous cutoff wavenumber
            KMAX = KNU/DSQRT(3.D0)  ! Most unstable wavenumber
         ENDIF
      ENDIF
      
c     Determine number of modes to include and wavenumber structure
      IF (LATERAL_PERIODIC) THEN
c        Periodic lateral boundaries: k_n = 2*pi*n/L
         NMAX = INT(KC*L/TWOPI)  ! Maximum mode number based on cutoff
      ELSE
c        Non-periodic lateral: k_n = pi*n/L
         NMAX = INT(KC*L/PI)     ! Maximum mode number based on cutoff
      ENDIF
      NMAX = MIN0(NMAX, NUMMODES)  ! Use fewer if specified

C     Set number of cells in x and y directions     
      NX = IBAR+2
      NY = JBAR+2
 
c     Initialize velocity fields to zero
      DO J = 1, JBAR2
         DO I = 1, IBAR2
            U(I,J) = 0.0
            V(I,J) = 0.0
         ENDDO
      ENDDO
      
c     Add contributions from each mode
      DO N = 1, NMAX
c        Set wavenumber based on lateral boundary conditions
         IF (LATERAL_PERIODIC) THEN
c           Periodic: k_n = 2*pi*n/L
            KX = N*TWOPI/L
         ELSE
c           Non-periodic: k_n = pi*n/L  
            KX = N*PI/L
         ENDIF
         
         K2 = KX*KX + KY*KY
         
c        Set amplitude for this mode (decaying with mode number)
         AN = A0 * (1.0/N)**ALPHA
         
c        Pre-compute hyperbolic functions for no-slip cases
         IF (VERTICAL_NOSLIP) THEN
            SINH_KB = SINH(KX*HB)
            SINH_KT = SINH(KX*HT)
            COSH_KB = COSH(KX*HB)
            COSH_KT = COSH(KX*HT)
         ENDIF
         
c        Loop through the grid (skip ghost cells)
         DO J = 2, NY-1
            Z = YJ(J) - HB  ! Transformed coordinate (-HB to HT)
            
            DO I = 2, NX-1
               X_COORD = XI(I)
               
c              Set horizontal structure based on lateral boundary conditions
               IF (LATERAL_PERIODIC) THEN
c                 Periodic boundaries: exp(i*kx*x) -> cos(kx*x) for real part
                  PHASE_X_W = DCOS(KX*X_COORD)
                  PHASE_X_U = DCOS(KX*X_COORD)
               ELSEIF (WL .EQ. 1 .AND. WR .EQ. 1) THEN
c                 Free-slip both sides: cos(kx*x) for w, sin(kx*x) for u
                  PHASE_X_W = DCOS(KX*X_COORD)
                  PHASE_X_U = DSIN(KX*X_COORD)
               ELSEIF (WL .EQ. 3 .AND. WR .EQ. 3) THEN
c                 Outflow both sides: cos(kx*x) structure
                  PHASE_X_W = DCOS(KX*X_COORD)
                  PHASE_X_U = DCOS(KX*X_COORD)
               ELSEIF (WL .EQ. 5 .AND. WR .EQ. 5) THEN
c                 Constant pressure both sides: cos(kx*x) structure
                  PHASE_X_W = DCOS(KX*X_COORD)
                  PHASE_X_U = DCOS(KX*X_COORD)
               ELSE
c                 Mixed or no-slip lateral: default to cos(kx*x)
                  PHASE_X_W = DCOS(KX*X_COORD)
                  PHASE_X_U = DCOS(KX*X_COORD)
               ENDIF
               
c              Determine which fluid region we're in
               IF (Z .LT. 0.0) THEN
c                 Lower fluid (fluid 1)
                  IF (VERTICAL_NOSLIP) THEN
c                    No-slip top/bottom: use sinh/cosh structure
                     PHI1 = SINH(KX*(Z + HB))/SINH_KB
                     DPHI1 = KX*COSH(KX*(Z + HB))/SINH_KB
                  ELSEIF (WB .EQ. 1 .AND. WT .EQ. 1) THEN
c                    Free-slip top/bottom: sinh/cosh structure
c                    (equivalent to sinh(k(z+HB))/sinh(kHB))
                     PHI1 = (DEXP(-KX*DABS(Z)) -
     &                      DEXP(KX*DABS(Z)-2*KX*HB))
     &                      /(1.0-DEXP(-2*KX*HB))
                     DPHI1 = KX*(DEXP(-KX*DABS(Z)) +
     &                       DEXP(KX*DABS(Z)-2*KX*HB))
     &                       /(1.0-DEXP(-2*KX*HB))
                  ELSE
c                    Mixed vertical BCs: use sinh/cosh as approximation
                     PHI1 = SINH(KX*(Z + HB))/SINH(KX*HB)
                     DPHI1 = KX*COSH(KX*(Z + HB))/SINH(KX*HB)
                  ENDIF
                  
c                 Add contribution to velocity components
                  U(I,J) = U(I,J) + AN * (-KX/DSQRT(K2)*DPHI1) * 
     &                      PHASE_X_U
                  V(I,J) = V(I,J) + AN * PHI1 * PHASE_X_W
                  
               ELSE
c                 Upper fluid (fluid 2)
                  IF (VERTICAL_NOSLIP) THEN
c                    No-slip top/bottom: use sinh/cosh structure
                     PHI2 = SINH(KX*(HT - Z))/SINH_KT
                     DPHI2 = -KX*COSH(KX*(HT - Z))/SINH_KT
                  ELSEIF (WB .EQ. 1 .AND. WT .EQ. 1) THEN
c                    Free-slip top/bottom: sinh/cosh structure
c                    (equivalent to sinh(k(HT-z))/sinh(kHT))
                     PHI2 = (DEXP(-KX*Z) -
     &                      DEXP(KX*Z-2*KX*HT))
     &                      /(1.0-DEXP(-2*KX*HT))
                     DPHI2 = -KX*(DEXP(-KX*Z) +
     &                        DEXP(KX*Z-2*KX*HT))
     &                        /(1.0-DEXP(-2*KX*HT))
                  ELSE
c                    Mixed vertical BCs: use sinh/cosh as approximation
                     PHI2 = SINH(KX*(HT - Z))/SINH(KX*HT)
                     DPHI2 = -KX*COSH(KX*(HT - Z))/SINH(KX*HT)
                  ENDIF
                  
c                 Add contribution to velocity components
                  U(I,J) = U(I,J) + AN * (-KX/DSQRT(K2)*DPHI2) * 
     &                      PHASE_X_U
                  V(I,J) = V(I,J) + AN * PHI2 * PHASE_X_W
               ENDIF
            ENDDO
         ENDDO
      ENDDO

c     Cap velocities at VORIG to ensure CFL stability
      IF (VORIG .GT. 0.0) THEN
         DO J = 2, NY-1
            DO I = 2, NX-1
               IF (U(I,J) .GT. VORIG) U(I,J) = VORIG
               IF (U(I,J) .LT. -VORIG) U(I,J) = -VORIG
               IF (V(I,J) .GT. VORIG) V(I,J) = VORIG
               IF (V(I,J) .LT. -VORIG) V(I,J) = -VORIG
            ENDDO
         ENDDO
         WRITE(*,*) 'RTINIT: Velocities capped at VORIG = ', VORIG
      ENDIF

c     Handle lateral ghost cells based on boundary conditions
      DO J = 1, NY
c        Left boundary (WL)
         IF (WL .EQ. 1) THEN
c           Free-slip: u=0, dv/dx=0
            U(1,J) = 0.0
            V(1,J) = V(2,J)
         ELSEIF (WL .EQ. 2) THEN
c           No-slip: u=0, v=0
            U(1,J) = -U(2,J)
            V(1,J) = -V(2,J)
         ELSEIF (WL .EQ. 3) THEN
c           Outflow: du/dx=0, dv/dx=0
            U(1,J) = U(2,J)
            V(1,J) = V(2,J)
         ELSEIF (WL .EQ. 4) THEN
c           Periodic: values from right side
            U(1,J) = U(NX-1,J)
            V(1,J) = V(NX-1,J)
         ELSEIF (WL .EQ. 5) THEN
c           Constant pressure: dp/dx=0 -> approximate as outflow
            U(1,J) = U(2,J)
            V(1,J) = V(2,J)
         ENDIF
         
c        Right boundary (WR)
         IF (WR .EQ. 1) THEN
c           Free-slip: u=0, dv/dx=0
            U(NX,J) = 0.0
            V(NX,J) = V(NX-1,J)
         ELSEIF (WR .EQ. 2) THEN
c           No-slip: u=0, v=0
            U(NX,J) = -U(NX-1,J)
            V(NX,J) = -V(NX-1,J)
         ELSEIF (WR .EQ. 3) THEN
c           Outflow: du/dx=0, dv/dx=0
            U(NX,J) = U(NX-1,J)
            V(NX,J) = V(NX-1,J)
         ELSEIF (WR .EQ. 4) THEN
c           Periodic: values from left side
            U(NX,J) = U(2,J)
            V(NX,J) = V(2,J)
         ELSEIF (WR .EQ. 5) THEN
c           Constant pressure: dp/dx=0 -> approximate as outflow
            U(NX,J) = U(NX-1,J)
            V(NX,J) = V(NX-1,J)
         ENDIF
      ENDDO
      
c     Handle vertical ghost cells based on boundary conditions
      DO I = 1, NX
c        Bottom boundary (WB)
         IF (WB .EQ. 1) THEN
c           Free-slip: v=0, du/dy=0
            U(I,1) = U(I,2)
            V(I,1) = 0.0
         ELSEIF (WB .EQ. 2) THEN
c           No-slip: u=0, v=0
            U(I,1) = -U(I,2)
            V(I,1) = -V(I,2)
         ELSEIF (WB .EQ. 3) THEN
c           Outflow: du/dy=0, dv/dy=0
            U(I,1) = U(I,2)
            V(I,1) = V(I,2)
         ELSEIF (WB .EQ. 4) THEN
c           Periodic: values from top side
            U(I,1) = U(I,NY-1)
            V(I,1) = V(I,NY-1)
         ELSEIF (WB .EQ. 5) THEN
c           Constant pressure: dp/dy=0 -> approximate as outflow
            U(I,1) = U(I,2)
            V(I,1) = V(I,2)
         ENDIF
         
c        Top boundary (WT)
         IF (WT .EQ. 1) THEN
c           Free-slip: v=0, du/dy=0
            U(I,NY) = U(I,NY-1)
            V(I,NY) = 0.0
         ELSEIF (WT .EQ. 2) THEN
c           No-slip: u=0, v=0
            U(I,NY) = -U(I,NY-1)
            V(I,NY) = -V(I,NY-1)
         ELSEIF (WT .EQ. 3) THEN
c           Outflow: du/dy=0, dv/dy=0
            U(I,NY) = U(I,NY-1)
            V(I,NY) = V(I,NY-1)
         ELSEIF (WT .EQ. 4) THEN
c           Periodic: values from bottom side
            U(I,NY) = U(I,2)
            V(I,NY) = V(I,2)
         ELSEIF (WT .EQ. 5) THEN
c           Constant pressure: dp/dy=0 -> approximate as outflow
            U(I,NY) = U(I,NY-1)
            V(I,NY) = V(I,NY-1)
         ENDIF
      ENDDO
      
      WRITE(*,*) 'RTINIT: NMAX = ', NMAX, ' modes included'
      WRITE(*,*) 'RTINIT: KC = ', KC, ' KMAX = ', KMAX
      IF (LATERAL_PERIODIC) THEN
         WRITE(*,*) 'RTINIT: Using periodic lateral wavenumbers'
      ELSE
         WRITE(*,*) 'RTINIT: Using non-periodic lateral wavenumbers'
      ENDIF
      IF (VERTICAL_NOSLIP) THEN
         WRITE(*,*) 'RTINIT: Using no-slip vertical structure'
      ELSE
         WRITE(*,*) 'RTINIT: Using free-slip vertical structure'
      ENDIF
      
      RETURN
      END

C     Utility subroutine to check BC compatibility for RT instability
      SUBROUTINE CHECK_RT_BCS(WL, WR, WB, WT, COMPATIBLE)
      IMPLICIT NONE
      INTEGER WL, WR, WB, WT
      LOGICAL COMPATIBLE
      
      COMPATIBLE = .TRUE.
      
c     Check for incompatible combinations
      IF (WL .EQ. 4 .AND. WR .NE. 4) THEN
         WRITE(*,*) 'Warning: Left periodic but right not periodic'
         COMPATIBLE = .FALSE.
      ENDIF
      
      IF (WR .EQ. 4 .AND. WL .NE. 4) THEN
         WRITE(*,*) 'Warning: Right periodic but left not periodic'
         COMPATIBLE = .FALSE.
      ENDIF
      
      IF (WB .EQ. 4 .AND. WT .NE. 4) THEN
         WRITE(*,*) 'Warning: Bottom periodic but top not periodic'
         COMPATIBLE = .FALSE.
      ENDIF
      
      IF (WT .EQ. 4 .AND. WB .NE. 4) THEN
         WRITE(*,*) 'Warning: Top periodic but bottom not periodic'
         COMPATIBLE = .FALSE.
      ENDIF
      
      RETURN
      END
C
C     *****************************************************************
C
      SUBROUTINE DALZIEL_PERTURB(HB, L, ILOW, IHIGH)
C
C     Dalziel et al. (1999) interface perturbation.
C     Modifies F field directly with perturbed interface.
C     Wavelengths 4*dx to 8*dx, sigma = 0.08*dx.
C     F=1 (heavy fluid) ABOVE interface, F=0 below.
C
      INCLUDE 'dint.h'
C
      REAL*8 HB, L
      INTEGER ILOW, IHIGH
C
      INTEGER NMIN, NMAX, NMODES, I, J, N
      REAL*8 DX, SIGMA_P, PI, TWOPI
      REAL*8 ETA, KX, YTOP, YBOT, FRAC
      REAL*8 AMP(200), PHS(200)
      REAL*8 U1, U2
      INTEGER ISEED
C
      PI = 3.14159265358979D0
      TWOPI = 2.0D0 * PI
      DX = DELX(2)
      SIGMA_P = 0.08D0 * DX
C
C     Mode range: wavelengths from 4*dx to 8*dx
      NMIN = INT(L / (8.0D0 * DX)) + 1
      NMAX = INT(L / (4.0D0 * DX))
      NMODES = NMAX - NMIN + 1
      IF (NMODES .LT. 1) THEN
         WRITE(*,*) 'DALZIEL_PERTURB: ERROR - no modes in range'
         RETURN
      END IF
      IF (NMODES .GT. 200) THEN
         WRITE(*,*) 'DALZIEL_PERTURB: ERROR - too many modes'
         NMODES = 200
      END IF
C
C     Uniform amplitude, random phases (portable LCG)
C     All modes have same amplitude; only phases are random.
C     Amplitude chosen so RMS of eta = sigma_p.
      ISEED = 123
      DO N = 1, NMODES
C       Uniform random for phase
        ISEED = MOD(ISEED * 1103515245 + 12345, 2147483647)
        PHS(N) = TWOPI * DBLE(ISEED) / 2147483647.0D0
C       Uniform amplitude for all modes
        AMP(N) = SIGMA_P * DSQRT(2.0D0 / DBLE(NMODES))
      END DO
C
      WRITE(*,*) 'DALZIEL_PERTURB: Dalziel 1999 interface perturbation'
      WRITE(*,*) '  dx =', DX, '  sigma =', SIGMA_P
      WRITE(*,*) '  Modes:', NMIN, 'to', NMAX, '(', NMODES, 'modes)'
      WRITE(*,*) '  Interface height HB =', HB
C
C     Modify F field
      DO I = 2, IM1
        ETA = 0.0D0
        DO N = 1, NMODES
          KX = TWOPI * DBLE(NMIN + N - 1) / L
          ETA = ETA + AMP(N) * DCOS(KX * XI(I) + PHS(N))
        END DO
C
        DO J = 2, JM1
          YTOP = YJ(J) + 0.5D0 * DELY(J)
          YBOT = YJ(J) - 0.5D0 * DELY(J)
C
          IF (YTOP .LE. HB + ETA) THEN
C           Cell entirely below interface: light fluid
            F(I,J) = 0.0D0
          ELSE IF (YBOT .GE. HB + ETA) THEN
C           Cell entirely above interface: heavy fluid
            F(I,J) = 1.0D0
          ELSE
C           Cell straddles interface: fractional fill
            FRAC = (YTOP - (HB + ETA)) / DELY(J)
            F(I,J) = DMAX1(0.0D0, DMIN1(1.0D0, FRAC))
          END IF
        END DO
      END DO
C
      RETURN
      END
C
C   *****************************************************************
C
C
      SUBROUTINE SETUP
C
      INCLUDE 'dint.h'
C
      common /new/ pltst1, pltst2, PLTDT1, PLTDT2

      COMMON /FL1/ NDIS,IDIS,IFL,IFR,JFB,JFT
C
      COMMON /OBS/ NOBS, IOMIN(NTOB), IOMAX(NTOB), JOMIN(NTOB),
     1 JOMAX(NTOB)
C
      COMMON /CONST/ EMF, EMF1, EM6, EM10, EP10, PI, RPD
C
C     *** COMPUTE CONSTANT TERMS AND INITIALIZE NECESSARY VARIABLES
C
C     *** SET PARAMETER STATEMENT VALUE INTO CONSTANT
C
      NVRM=NVOR
C
      CYL=FLOAT(ICYL)
      EMF1=1.D0-EMF
      T=0.D0
      ITER=0
      iCYCLE=0
      TWPRT=PRTDT
                            
      pltdt = pltdt1
      TWPLT = pltst1

      VCHGT=0.D0
      NOCON=0
      NFLGC=0
      FNOC=0.D0
      RCSQ=1.D0/(RHOF*CSQ)
      IF (CSQ.LT.0.D0) RCSQ=0.D0
      IF (NMAT.EQ.1) RHOFC=RHOF
      RHOD=RHOF-RHOFC
      IF (CANGLE.EQ.90.D0) CANGLE=CANGLE-EM6
      CANGLE=CANGLE*RPD
      DTANCA=DTAN(CANGLE)
      IPL=2
      IF (IWL.EQ.5) IPL=3
      IPR=IM1
      IF (IWR.EQ.5) IPR=IM2
      JPB=2
      IF (IWB.EQ.5) JPB=3
      JPT=JM1
      IF (IWT.EQ.5) JPT=JM2
C
C     *** COMPUTE INITIAL VOLUME FRACTION FUNCTION F IN CELLS
C     *** FOR NMAT = 1; THAT IS FOR A FREE-SURFACE PROBLEM
C     DO 40 I=1,IMAX
C     DO 30 J=2,JMAX
C     F(I,J)=1.D0
C     IF (FLHT.GT.Y(J-1).AND.FLHT.LT.Y(J)) F(I,J)=RDY(J)*(FLHT-Y(J-1))
C     IF (Y(J-1).GE.FLHT) F(I,J)=0.D0
C  30 CONTINUE
C     F(I,1)=F(I,2)
C  40 CONTINUE
C
C     *** GENERATE SPECIAL F-FUNCTION (FLUID) CONFIGURATION
C
C     *** Set up a two region initial domain. 
C     *** 
C     *** The indices IFL,IFR,JFB,JFT define the volume number ranges
C     *** which contain fluid with F = 1.  The numbers should be biased by
C     *** by 1 to account for the left and bottom fictitious volumes.
C
      IF (NMAT .EQ. 2) THEN
       READ ( 1,"(//,4I5)") IFL,IFR,JFB,JFT
C       print'("IFL: ", I5," IFR: ", I5," JFB: ", I5, " JFT: ",I5)',
C     1 IFL, IFR, JFB, JFT
       DO 41 I=1,IMAX
       DO 41 J=1,JMAX
        F(I,J) = 0.D0
  41   CONTINUE

       DO 42 I=IFL,IFR
       DO 42 J=JFB,JFT
        F(I,J) = 1.D0
  42   CONTINUE

      END IF     

C
C     *** CALCULATE DTVIS AND DTSFT
C
      DS=1.D+10
      DTVIS=1.D+10
      DTSFT=1.D+10
      DO 50 I=2,IM1
      DO 50 J=2,JM1
      DXSQ=DELX(I)*DELX(I)
      DYSQ=DELY(J)*DELY(J)
      RDSQ=DXSQ*DYSQ/(DXSQ+DYSQ)
      RNUMAX = DMAX1(VNU,VNUC)
      RDSQ=RDSQ/(3.D0*RNUMAX+1.D-10)
      DTVIS=DMIN1(DTVIS,RDSQ)
      DS=DMIN1(DELX(I),DELY(J),DS)
   50 CONTINUE
      SIGX=SIGMA
      RHOMN=DMIN1(RHOF,RHOFC)
      IF(SIGX.EQ.0.D0) SIGX=EM10
      DTM=DSQRT(RHOMN*DS*DS*DS/(SIGX*4.D0*(1.D0+CYL)))
      DTSFT=DMIN1(DTSFT,DTM)
C
C     *** CALCULATE BETA(I,J) FOR MESH
C
      RDTEXP= 2.D0*DSQRT( DABS(CSQ))/DS
      IF(CSQ.LT.0.D0) RDTEXP= 1.D+10
      CTOS=DELT*RDTEXP
      COMG= DMIN1(CTOS*CTOS,1.D0)
      OMG1=(OMG-1.D0)*COMG+1.D0
      DO 55 I=2,IM1
      DO 55 J=2,JM1
      RHXR=(RHOFC+RHOD*F(I,J))*DELX(I+1)+(RHOFC+RHOD*F(I+1,J))*DELX(I)
      RHXL=(RHOFC+RHOD*F(I-1,J))*DELX(I)+(RHOFC+RHOD*F(I,J))*DELX(I-1)
      RHYT=(RHOFC+RHOD*F(I,J))*DELY(J+1)+(RHOFC+RHOD*F(I,J+1))*DELY(J)
      RHYB=(RHOFC+RHOD*F(I,J-1))*DELY(J)+(RHOFC+RHOD*F(I,J))*DELY(J-1)
      XX=DELT*RDX(I)*(2.D0/RHXL+2.D0/RHXR)+DELT*RDY(J)*
     1   (2.D0/RHYT+2.D0/RHYB)
      RHOR=RHOF/(RHOFC+RHOD*F(I,J))
      BETA(I,J)=OMG1/(XX*COMG+RCSQ*RHOR/DELT)
   55 CONTINUE
C
C     *** SET BETA(I,J)= -1.0  IN OBSTACLE CELLS.
C
      IF (NOBS .EQ. 0) GO TO 58
      DO 56 K=1,NOBS
        IF (K .EQ. 1) THEN
          READ ( 1,230) IOMIN(K),IOMAX(K),JOMIN(K),JOMAX(K)
        ELSE
          READ ( 1,240) IOMIN(K),IOMAX(K),JOMIN(K),JOMAX(K)
        END IF
       DO 57 J=JOMIN(K),JOMAX(K)
       DO 57 I=IOMIN(K),IOMAX(K)
        BETA(I,J) = -1.D0
  57   CONTINUE
  56  CONTINUE
C
C     *** PRINT BETA(I,J)
C
   58 IF (IMOVY.EQ.0) GO TO 70
      WRITE ( 2,210)
      DO 60 J=1,JM1
      DO 60 I=1,IM1
      WRITE ( 2,220) I,J,BETA(I,J)
   60 CONTINUE
   70 CONTINUE
C
C     *** CALCULATE HYDROSTATIC PRESSURE
C
      DO I=2,IM1

        P(I,JMAX)=0.D0

        DO J=2,JM1
          JJ=JM1-J+2
          RHOYA=(RHOFC+RHOD*F(I,JJ))*DELY(JJ)/2.D0+(RHOFC+RHOD*F(I,JJ+1)
     1    )*DELY(JJ+1)/2.D0
          IF (NMAT.EQ.1) RHOYA=(DMIN1(F(I,JJ+1),0.5D0)*DELY(JJ+1)+
     1        DMAX1(0.D0,F(I,JJ)-0.5D0)*DELY(JJ))*RHOF
          P(I,JJ)=P(I,JJ+1)-GY*RHOYA
        ENDDO

      ENDDO
C
C     *** PARTICLE SET UP
C
      NP=NPY*(1+NPX)
      IF (NP.EQ.0) GO TO 160
      DXP=(XPR-XPL)/FLOAT(NPX)
      DYP=(YPT-YPB)/FLOAT(NPY)
      K=0
      DO 100 JN=1,NPY,2
      DO 100 IN=1,NPX
      K=K+1
      XP(K)=XPL+(FLOAT(IN)-0.5D0)*DXP
      YP(K)=YPB+(FLOAT(JN)-1.D0)*DYP
      IF (YP(K).GT.YPT) YP(K)=YPT
  100 CONTINUE
      DO 110 JN=2,NPY,2
      K=K+1
      XP(K)=XPL
      YP(K)=YPB+(FLOAT(JN)-1.D0)*DYP
      IF (YP(K).GT.YPT) YP(K)=YPT
      DO 110 IN=1,NPX
      K=K+1
      XP(K)=XPL+FLOAT(IN)*DXP
      YP(K)=YPB+(FLOAT(JN)-1.D0)*DYP
      IF (YP(K).GT.YPT) YP(K)=YPT
  110 CONTINUE
      NP=K
      DO 150 K=1,NP
      DO 120 I=2,IM1
      IF (XP(K).GE.X(I-1).AND.XP(K).LE.X(I)) IP(K)=I
      IF (X(I-1).GT.XPR) GO TO 130
  120 CONTINUE
  130 DO 140 J=2,JM1
      IF (YP(K).GE.Y(J-1).AND.YP(K).LE.Y(J)) JP(K)=J
      IF (Y(J-1).GT.YPT) GO TO 150
  140 CONTINUE
  150 CONTINUE
  160 CONTINUE
C
C     *** SET INITIAL SURFACE PRESSURE
C
      DO 170 J=2,JM1
      DO 170 I=2,IM1
      PS(I,J)=0.D0
  170 CONTINUE
C
C     *** SET INITIAL VELOCITY FIELD INTO U AND V ARRAYS
C
      DO 180 I=2,IM1
      DO 180 J=2,JM1
CRWD      V(I,J)=VI
CRWD      U(I,J)=UI
CRWD      IF (F(I,J).GT.EMF.OR.NMAT.EQ.2) GO TO 180
C
          V(I,J) = 0.D0
          U(I,J) = 0.D0
  180 CONTINUE
      IF (NOBS .NE. 0) THEN
      IOM1 = IOMIN(1) - 1
      DO 181 I=2,IOM1
      DO 181 J=JM1,JMAX
         V(I,J) = VI
         U(I,J) = UI
  181 CONTINUE
      END IF
C
C     SETUP OF IMPOSED SHEAR FORCES
C
      IF (ISHVEL.EQ.1) THEN
      DO 185 I = 1,IMAX
      U(I,1) = SHVEL
      U(I,2) = SHVEL
      U(I,JM1) = -SHVEL
      U(I,JMAX) = -SHVEL
 185  CONTINUE
      ULINE = SHVEL
      DELTVL = 2.0*SHVEL/FLOAT(JMAX-3)
      DO 186 J = 3,JMAX-2
      ULINE = ULINE -DELTVL
      WRITE(6,*) J,ULINE
      DO 186 I=1,IMAX
      U(I,J) = ULINE
 186  CONTINUE
      END IF      
C
C     *** If NDIS .NE. 0, compute the disturbance velocity components for
C     *** the initial distribution.
C     *** NDIS = 1 --> use planer Rayleigh-Taylor disturbances
C     *** NDIS = 2 --> use cylindrical Rayleigh-Taylor disturbances
C
      IF (NDIS .NE. 0) THEN
        CALL DISTURB
        CALL BC
      END IF
C
C     *** SET INITIAL VOID REGION QUANTITIES
C
      DO K=1,NVRM
       NR(K)=0
       PR(K)=0.D0
       VOL(K)=0.D0
      END DO
      RETURN
C
  210 FORMAT (1H1)
  220 FORMAT (2X,5HBETA(,I2,1H,,I2,2H)=,1PE14.7)
  230 FORMAT (//,4I5)
  240 FORMAT (4I5)
      END
C
C
C   *****************************************************************
C   *****************************************************************
C
C
      SUBROUTINE TILDE
C
      INCLUDE 'dint.h'
C
      COMMON /CONST/ EMF, EMF1, EM6, EM10, EP10, PI, RPD
C
C     *** COMPUTE TEMPORARY U AND V DEXPLICITLY
C
      DO 20 J=2,JM1
      DO 20 I=2,IM1
      U(I,J)=0.D0
      RDELX=1.D0/(DELX(I)+DELX(I+1))
      RDELY=1.D0/(DELY(J)+DELY(J+1))
      IF (F(I,J)+F(I+1,J).LT.EMF.AND.NMAT.EQ.1) GO TO 10
      IF (BETA(I,J).LT.0.D0.OR.BETA(I+1,J).LT.0.D0) GO TO 10
      SGU= DSIGN(1.D0,UN(I,J))
      DUDR=(UN(I+1,J)-UN(I,J))*RDX(I+1)
      DUDL=(UN(I,J)-UN(I-1,J))*RDX(I)
      RDXA=DELX(I)+DELX(I+1)+ALPHA*SGU*(DELX(I+1)-DELX(I))
      RDXA=1.D0/RDXA
      FUX=RDXA*UN(I,J)*(DELX(I)*DUDR+DELX(I+1)*DUDL+ALPHA*SGU*(DELX(I+1)
     1 *DUDL-DELX(I)*DUDR))
      VBT=(DELX(I)*VN(I+1,J)+DELX(I+1)*VN(I,J))*RDELX
      VBB=(DELX(I)*VN(I+1,J-1)+DELX(I+1)*VN(I,J-1))*RDELX
      VAV=(VBT+VBB)/2.D0
      DYT=(DELY(J)+DELY(J+1))/2.D0
      DYB=(DELY(J-1)+DELY(J))/2.D0
      DUDT=(UN(I,J+1)-UN(I,J))/DYT
      DUDB=(UN(I,J)-UN(I,J-1))/DYB
      SGV= DSIGN(1.D0,VAV)
      DYA=DYT+DYB+ALPHA*SGV*(DYT-DYB)
      FUY=(VAV/DYA)*(DYB*DUDT+DYT*DUDB+ALPHA*SGV*(DYT*DUDB-DYB*DUDT))
      UBDYT=(DELY(J)*UN(I,J+1)+DELY(J+1)*UN(I,J))/(DELY(J)+DELY(J+1))
      UBDYB=(DELY(J-1)*UN(I,J)+DELY(J)*UN(I,J-1))/(DELY(J)+DELY(J-1))
      DUDXSQ=2.D0*(UN(I-1,J)*RDX(I)/(DELX(I)+DELX(I+1))+UN(I+1,J)*
     1 RDX(I+1)/(DELX(I)+DELX(I+1))-UN(I,J)*RDX(I)*RDX(I+1))
      DUDYT=(UN(I,J+1)*DELY(J)*RDY(J+1)-UN(I,J)*DELY(J+1)*RDY(J)-UBDYT*
     1 (DELY(J)*RDY(J+1)-DELY(J+1)*RDY(J)))/(0.5D0*(DELY(J)+DELY(J+1)))
      DUDYB=(UN(I,J)*DELY(J-1)*RDY(J)-UN(I,J-1)*DELY(J)*RDY(J-1)-UBDYB*
     1 (DELY(J-1)*RDY(J)-DELY(J)*RDY(J-1)))/(0.5D0*(DELY(J-1)+DELY(J)))
      DUDYSQ=(DUDYT-DUDYB)*RDY(J)
      DUDXL=(UN(I,J)-UN(I-1,J))*RDX(I)
      DUDXR=(UN(I+1,J)-UN(I,J))*RDX(I+1)
      RXDUDX=RX(I)*(DELX(I+1)*DUDXL+DELX(I)*DUDXR)/(DELX(I)+DELX(I+1))
      RXSQU=UN(I,J)*RX(I)*RX(I)
      RNUCAL = F(I,J)*VNU + (1.0-F(I,J))*VNUC
      VISX=RNUCAL*(DUDXSQ+DUDYSQ+CYL*RXDUDX-CYL*RXSQU)
      RHOX=(RHOFC+RHOD*F(I,J))*DELX(I+1)+(RHOFC+RHOD*F(I+1,J))*DELX(I)
      U(I,J)=UN(I,J)+DELT*((P(I,J)-P(I+1,J))*2.D0/RHOX+GX-FUX-FUY+VISX)
C Add this at the end of each velocity calculation (after U(I,J)= and V(I,J)= lines)
       IF (ABS(U(I,J)) .GT. 1.0D3 .OR. U(I,J) .NE. U(I,J)) THEN
        WRITE(6,*) '*** U VELOCITY EXPLOSION ***'
        WRITE(6,*) 'T=',T,' CYCLE=',iCYCLE,' ITER=',ITER
        WRITE(6,*) 'CELL I=',I,' J=',J
        WRITE(6,*) 'U(I,J)=',U(I,J)
        WRITE(6,*) 'UN(I,J)=',UN(I,J)
        WRITE(6,*) 'F(I,J)=',F(I,J),' F(I+1,J)=',F(I+1,J)
        WRITE(6,*) 'P(I,J)=',P(I,J),' P(I+1,J)=',P(I+1,J)
        WRITE(6,*) 'DELT=',DELT
        WRITE(6,*) 'PRESSURE GRADIENT=',(P(I,J)-P(I+1,J))*2.D0/RHOX
        WRITE(6,*) 'GRAVITY=',GX
        WRITE(6,*) 'ADVECTION FUX=',FUX,' FUY=',FUY
        WRITE(6,*) 'VISCOUS=',VISX
        STOP 'U velocity explosion detected'
      ENDIF

   10 CONTINUE
      V(I,J)=0.D0
      IF (F(I,J)+F(I,J+1).LT.EMF.AND.NMAT.EQ.1) GO TO 20
      IF (BETA(I,J).LT.0.D0.OR.BETA(I,J+1).LT.0.D0) GO TO 20
      UBR=(DELY(J+1)*UN(I,J)+DELY(J)*UN(I,J+1))*RDELY
      UBL=(DELY(J+1)*UN(I-1,J)+DELY(J)*UN(I-1,J+1))*RDELY
      UAV=(UBR+UBL)/2.D0
      DXR=(DELX(I)+DELX(I+1))/2.D0
      DXL=(DELX(I)+DELX(I-1))/2.D0
      SGU= DSIGN(1.D0,UAV)
      DXA=DXR+DXL+ALPHA*SGU*(DXR-DXL)
      DVDR=(VN(I+1,J)-VN(I,J))/DXR
      DVDL=(VN(I,J)-VN(I-1,J))/DXL
      FVX=(UAV/DXA)*(DXL*DVDR+DXR*DVDL+ALPHA*SGU*(DXR*DVDL-DXL*DVDR))
      SGV= DSIGN(1.D0,VN(I,J))
      DYA=DELY(J+1)+DELY(J)+ALPHA*SGV*(DELY(J+1)-DELY(J))
      DVDT=(VN(I,J+1)-VN(I,J))*RDY(J+1)
      DVDB=(VN(I,J)-VN(I,J-1))*RDY(J)
      FVY=(VN(I,J)/DYA)*(DELY(J)*DVDT+DELY(J+1)*DVDB+ALPHA*SGV*(DELY(J+1
     1 )*DVDB-DELY(J)*DVDT))
      VBDYR=(DELX(I+1)*VN(I,J)+DELX(I)*VN(I+1,J))/(DELX(I)+DELX(I+1))
      VBDYL=(DELX(I)*VN(I-1,J)+DELX(I-1)*VN(I,J))/(DELX(I)+DELX(I-1))
      DVDXR=(VN(I+1,J)*DELX(I)*RDX(I+1)-VN(I,J)*DELX(I+1)*RDX(I)-VBDYR*
     1 (DELX(I)*RDX(I+1)-DELX(I+1)*RDX(I)))/(0.5D0*(DELX(I+1)+DELX(I)))
      DVDXL=(VN(I,J)*DELX(I-1)*RDX(I)-VN(I-1,J)*DELX(I)*RDX(I-1)-VBDYL*
     1 (DELX(I-1)*RDX(I)-DELX(I)*RDX(I-1)))/(0.5D0*(DELX(I)+DELX(I-1)))
      DVDXSQ=(DVDXR-DVDXL)*RDX(I)
      DVDYSQ=2.D0*(VN(I,J-1)*RDY(J)/(DELY(J+1)+DELY(J))-
     1 VN(I,J)*RDY(J+1)*RDY(J)+VN(I,J+1)*RDY(J+1)/(DELY(J+1)+DELY(J)))
      DVDXRX=(VBDYR-VBDYL)*RDX(I)*RXI(I)
      RNUCAL = F(I,J)*VNU + (1.0-F(I,J))*VNUC
      VISY=RNUCAL*(DVDXSQ+DVDYSQ+CYL*DVDXRX)
      RHOY=(RHOFC+RHOD*F(I,J))*DELY(J+1)+(RHOFC+RHOD*F(I,J+1))*DELY(J)
      V(I,J)=VN(I,J)+DELT*((P(I,J)-P(I,J+1))*2.D0/RHOY+GY-FVX-FVY+VISY)
C Safety check for V velocity  
      IF (ABS(V(I,J)) .GT. 1.0D3 .OR. V(I,J) .NE. V(I,J)) THEN
        WRITE(6,*) '*** V VELOCITY EXPLOSION ***'
        WRITE(6,*) 'T=',T,' CYCLE=',iCYCLE,' ITER=',ITER
        WRITE(6,*) 'CELL I=',I,' J=',J
        WRITE(6,*) 'V(I,J)=',V(I,J)
        WRITE(6,*) 'VN(I,J)=',VN(I,J)
        WRITE(6,*) 'F(I,J)=',F(I,J),' F(I,J+1)=',F(I,J+1)
        WRITE(6,*) 'P(I,J)=',P(I,J),' P(I,J+1)=',P(I,J+1)
        WRITE(6,*) 'DELT=',DELT
        WRITE(6,*) 'PRESSURE GRADIENT=',(P(I,J)-P(I,J+1))*2.D0/RHOY
        WRITE(6,*) 'GRAVITY=',GY
        WRITE(6,*) 'ADVECTION FVX=',FVX,' FVY=',FVY
        WRITE(6,*) 'VISCOUS=',VISY
        STOP 'V velocity explosion detected'
      ENDIF
   20 CONTINUE
      RETURN
      END
C
C
C   *****************************************************************
C   *****************************************************************
C
C
      SUBROUTINE TMS10
C
      INCLUDE 'dint.h'
C
      COMMON /CONST/ EMF, EMF1, EM6, EM10, EP10, PI, RPD
C
C     *** TWO MATERIAL SURFACE TENSION
C
C     *** NOTE: THIS ROUTINE INTRODUCES SOME NUMERICAL NOISE
C     ***       AND MAY BE REPLACED IN THE FUTURE
C
      DO 30 I=2,IM1
      DO 30 J=2,JM1
      IF (NF(I,J).EQ.0.OR.NF(I,J).GE.5.OR.BETA(I,J).LT.0.D0) GO TO 30
      WHTL=0.5D0
      WHTR=0.5D0
      WHTT=0.5D0
      WHTB=0.5D0
      IF (NF(I,J).GT.2) GO TO 10
      WHTL=1.D0-F(I,J)
      IF (NF(I,J).EQ.2) WHTL=1.D0-WHTL
      WHTR=1.D0-WHTL
      STFX=PS(I,J)*DELY(J)
      IF (NF(I,J).EQ.1) STFX=-STFX
      STFY=STFX*DTANTH(I,J)
      GO TO 20
   10 WHTB=1.D0-F(I,J)
      IF (NF(I,J).EQ.4) WHTB=1.D0-WHTB
      WHTT=1.D0-WHTB
      STFY=PS(I,J)*DELX(I)
      IF (NF(I,J).EQ.3) STFY=-STFY
      STFX=-STFY*DTANTH(I,J)
   20 CONTINUE
      RHOXR=(RHOFC+RHOD*F(I,J))*DELX(I+1)+(RHOFC+RHOD*F(I+1,J))*DELX(I)
      U(I,J)=U(I,J)+2.D0*DELT*WHTR*STFX/(RHOXR*DELY(J))
      RHOXL=(RHOFC+RHOD*F(I-1,J))*DELX(I)+(RHOFC+RHOD*F(I,J))*DELX(I-1)
      U(I-1,J)=U(I-1,J)+2.D0*DELT*WHTL*STFX/(RHOXL*DELY(J))
      RHOYT=(RHOFC+RHOD*F(I,J))*DELY(J+1)+(RHOFC+RHOD*F(I,J+1))*DELY(J)
      V(I,J)=V(I,J)+2.D0*DELT*WHTT*STFY/(RHOYT*DELX(I))
      RHOYB=(RHOFC+RHOD*F(I,J-1))*DELY(J)+(RHOFC+RHOD*F(I,J))*DELY(J-1)
      V(I,J-1)=V(I,J-1)+2.D0*DELT*WHTB*STFY/(RHOYB*DELX(I))
   30 CONTINUE
      RETURN
      END
C
C
C
C   *****************************************************************
C   *****************************************************************
C
C
      SUBROUTINE VFCONV
C
      INCLUDE 'dint.h'
C
      COMMON /CONST/ EMF, EMF1, EM6, EM10, EP10, PI, RPD
C
C     *** CONVECT THE VOLUME OF FLUID FUNCTION F
C
      IF (iCYCLE.LT.1) GO TO 40
      FLGC=0.D0
      DO 30 J=1,JM1
      DO 30 I=1,IM1
      VX=U(I,J)*DELT
      VY=V(I,J)*DELT
      ABVX= DABS(VX)
      ABVY= DABS(VY)
      IF (ABVX.GT.0.5D0*DELX(I).OR.ABVY.GT.0.5D0*DELY(J)) FLGC=1.D0
      IA=I+1
      ID=I
      IDM=MAX0(I-1,1)
      RB=X(I)
      RA=XI(I+1)
      RD=XI(I)
      IF (VX.GE.0.D0) GO TO 10
      IA=I
      ID=I+1
      IDM=MIN0(I+2,IMAX)
      RA=XI(I)
      RD=XI(I+1)
   10 CONTINUE
      IAD=IA
      IF (NF(ID,J).EQ.3.OR.NF(ID,J).EQ.4) IAD=ID
      IF (FN(IA,J).LT.EMF.OR.FN(IDM,J).LT.EMF) IAD=IA
      FDM=DMAX1(FN(IDM,J),FN(ID,J))
      FX1=FN(IAD,J)* DABS(VX)+DMAX1((FDM-FN(IAD,J))* DABS(VX)-
     1    (FDM-FN(ID,J))*DELX(ID),0.D0)
      FX=DMIN1(FX1,FN(ID,J)*DELX(ID))
      F(ID,J)=F(ID,J)-FX*RDX(ID)*(( DABS(RB/RD))*CYL+(1.D0-CYL))
      F(IA,J)=F(IA,J)+FX*RDX(IA)*(( DABS(RB/RA))*CYL+(1.D0-CYL))
      JA=J+1
      JD=J
      JDM=MAX0(J-1,1)
      IF (VY.GE.0.D0) GO TO 20
      JA=J
      JD=J+1
      JDM=MIN0(J+2,JMAX)
   20 CONTINUE
      JAD=JA
      IF (NF(I,JD).EQ.1.OR.NF(I,JD).EQ.2) JAD=JD
      IF (FN(I,JA).LT.EMF.OR.FN(I,JDM).LT.EMF) JAD=JA
      FDM=DMAX1(FN(I,JDM),FN(I,JD))
      FY1=FN(I,JAD)* DABS(VY)+DMAX1((FDM-FN(I,JAD))* DABS(VY)-
     1    (FDM-FN(I,JD))*DELY(JD),0.D0)
      FY=DMIN1(FY1,FN(I,JD)*DELY(JD))
      F(I,JD)=F(I,JD)-FY*RDY(JD)
      F(I,JA)=F(I,JA)+FY*RDY(JA)
   30 CONTINUE
   40 CONTINUE
      DO 80 J=2,JM1
      DO 80 I=2,IM1
      IF (BETA(I,J).LT.0.D0) GO TO 80
      VCHG=0.D0
      IF (F(I,J).GT.EMF.AND.F(I,J).LT.EMF1) GO TO 60
      IF (F(I,J).GE.EMF1) GO TO 50
      VCHG=F(I,J)
      F(I,J)=0.D0
      GO TO 60
   50 CONTINUE
      VCHG=-(1.D0-F(I,J))
      F(I,J)=1.D0
   60 CONTINUE
      VCHGT=VCHGT+VCHG*DELX(I)*DELY(J)*(XI(I)*2.D0*PI*CYL+(1.D0-CYL))
      IF (F(I,J).LT.EMF1) GO TO 80
      IF (F(I+1,J).LT.EMF) GO TO 70
      IF (F(I-1,J).LT.EMF) GO TO 70
      IF (F(I,J+1).LT.EMF) GO TO 70
      IF (F(I,J-1).LT.EMF) GO TO 70
      GO TO 80
   70 F(I,J)=F(I,J)-1.1D0*EMF
      VCHG=1.1D0*EMF
      VCHGT=VCHGT+VCHG*DELX(I)*DELY(J)*(XI(I)*2.D0*PI*CYL+(1.D0-CYL))
   80 CONTINUE
      RETURN
      END

      SUBROUTINE WRITEVTK ()
C
      INCLUDE 'dint.h'
C
      CHARACTER(11) PREFIX
      CHARACTER(30) fileName
      CHARACTER(12) cT1
      CHARACTER(12) cT2
      CHARACTER(12) cT3
C
      INTEGER nunit, IND, NX, NY, NZ, JLO, JHI, IDIF
C
C     *** WRITEVTK--
C     *** PROVIDES FORMATTED WRITES TO A VTK FILE FOR PARAVIEW VISUALIZATION.
C
C     Calculate the number of real x,y,z coordinates in the mesh

      NX = IBAR + 1
      NY = JBAR + 1
      NZ = 1

      write(cT1,'(I12)') IBAR 
      write(cT2,'(I12)') JBAR 
      PREFIX = "RT"//TRIM(adjustl(cT1))//"x"//TRIM(adjustl(cT2))
c      write(6,*) PREFIX

      write(cT3,'(I12)') INT(T*1000.)
      fileName = TRIM(PREFIX)//"-"//TRIM(adjustl(cT3))//".vtk"

C      write(6,*) filename

      OPEN(newunit=nunit, file=fileName)

      write(nunit,'(A26)') "# vtk DataFile Version 2.0"
      write(nunit,'(A12,/,A5,/,A25)') "R-T Fractals","ASCII",
     ?"DATASET RECTILINEAR_GRID"
      write(nunit,'(A10,1X,I5,1X,I5,1X,I5)') "DIMENSIONS", NX, NY, NZ

      write(nunit,'(A13,1X,I5,1X,A5)') "X_COORDINATES", NX, "float"
      
      X1 = NX
      IND = FLOOR(X1/5.)
      IDIF = MOD(NX,5)
C      write(nunit,'("NX: ", I3, " IND: ", I3, " IDIF: ", I3)') NX, 
C     1 IND, IDIF
      IF (IDIF .NE. 0) IND = IND + 1
C      write(6,'("IMAX: ", I3, " NX: ", I3," IDIF: ", I3, " IND: ", I3)') 
C     1 IMAX, NX, IDIF, IDIF
      do K=1,IND
        JLO = (K-1)*5 + 1
        JHI = JLO + 4
        IF (JHI .GT. IMAX-2) JHI = NX
C        write(nunit, '("JLO: ", I3," JHIGH: ", I3)') JLO, JHIGH
        write(nunit,'(5(E15.8,1X))') (X(J),J=JLO,JHI)
      end do

      write(nunit,'(A13,1X,I5,1X,A5)') "Y_COORDINATES", NY, "float"
      Y1 = NY
      IND = FLOOR(Y1/5.)
      IDIF = MOD(NY,5)
      IF (IDIF .NE. 0) IND = IND + 1
      do K=1,IND
        JLO = (K-1)*5 + 1
        JHI = JLO + 4
        IF (JHI .GT. JMAX-2) JHI = NY
        write(nunit,'(5(E15.8,1X))') (Y(J),J=JLO,JHI)
      end do

      write(nunit,'(A13,1X,I5,1X,A5)') "Z_COORDINATES", NZ, "float"
      write(nunit, '(E15.8)') 0.

      write(nunit,'(A10,1X,I8)') "CELL_DATA", (NX-1)*(NY-1)*NZ

      X1 = IBAR
      IND = FLOOR(X1/5.)
      IDIF = MOD(IBAR,5)

      write(nunit,'(A7,1X,A1,1X,A5,1X,I2)') "SCALARS", "F", "float", 1
      write(nunit,'(A20)') "LOOKUP_TABLE default"
      IF (IDIF .NE. 0) IND = IND + 1
      DO J=2,JBAR+1
      DO K=1,IND
      ILO = (K-1)*5 + 2
      IHI = ILO + 4
      IF (IHI .GT. IBAR+1) IHI = IBAR+1
      WRITE (nunit,'(5(E15.8,1X))') (F(I,J),I=ILO,IHI)
      END DO
      END DO

      write(nunit,'(A7,1X,A1,1X,A5,1X,I2)') "SCALARS", "P", "float", 1
      write(nunit,'(A20)') "LOOKUP_TABLE default"
      DO J=2,JBAR+1
      DO K=1,IND
      ILO = (K-1)*5 + 2
      IHI = ILO + 4
      IF (IHI .GT. IBAR+1) IHI = IBAR+1
      WRITE (nunit,'(5(E15.8,1X))') (P(I,J),I=ILO,IHI)
      END DO
      END DO

      write(nunit,'(A7,1X,A1,1X,A5,1X,I2)') "SCALARS", "U", "float", 1
      write(nunit,'(A20)') "LOOKUP_TABLE default"
      DO J=2,JBAR+1
      DO K=1,IND
      ILO = (K-1)*5 + 2
      IHI = ILO + 4
      IF (IHI .GT. IBAR+1) IHI = IBAR+1
      WRITE (nunit,'(5(E15.8,1X))') (U(I,J),I=ILO,IHI)
      END DO
      END DO

      write(nunit,'(A7,1X,A1,1X,A5,1X,I2)') "SCALARS", "V", "float", 1
      write(nunit,'(A20)') "LOOKUP_TABLE default"
      DO J=2,JBAR+1
      DO K=1,IND
      ILO = (K-1)*5 + 2
      IHI = ILO + 4
      IF (IHI .GT. IBAR+1) IHI = IBAR+1
      WRITE (nunit,'(5(E15.8,1X))') (V(I,J),I=ILO,IHI)
      END DO
      END DO

      CLOSE(nunit)

      END
C   END OF SUBROUTINE WRITEVTK
