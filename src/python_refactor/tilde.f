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
