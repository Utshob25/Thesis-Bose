C     Last change:  MR   17 Apr 2007    7:24 pm
      PROGRAM PR12
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(MMX=1500)
      DIMENSION SG00(1),RS0(1),Q(MMX),G(MMX)
      DIMENSION R(MMX),VQ1(MMX),VQ2(MMX),F11(MMX),F12(MMX),
     &      F22(MMX),VI11(MMX),VI12(MMX),VI22(MMX),
     &          V11(MMX),V12(MMX),V22(MMX),W(MMX),VQ3(MMX)
      OPEN(1,FILE='PSEUDOPT.DAT',STATUS='OLD')
      OPEN(51,FILE='FOR051.DAT',STATUS='UNKNOWN')
      OPEN(52,FILE='FOR052.DAT',STATUS='UNKNOWN')
      OPEN(53,FILE='FOR053.DAT',STATUS='UNKNOWN')
      OPEN(54,FILE='FOR054.DAT',STATUS='UNKNOWN')
      OPEN(25,FILE='BRETIU.DAT',STATUS='UNKNOWN')
C     TIN IS THE TEMPERATURE TO INITIATE THE CALCULATION
C     T IS THE FINAL TEMPERATURE
      READ(1,*)TIN,T
C     RC IS THE CORE RADIUS. IT MUST BE IN A.U.
      READ(1,*)RC1,RC2
      READ(1,*)Z1,Z2
C     DENA IS THE NUMBER DENSITY (IN ANGSTROMS)
      READ(1,*)DENA, C1
C     MAXQ IS THE MAXIMUN Q NUMBER
C     DRIN IS THE DR IN ANSTRONGS. THE VALUE WILL BE MODIFIED
C     ACCORDING TO ISTEP
      READ(1,*)MAXQ,DRIN
C     MAXR IS THE MAXIMUN R NUMBER FOR PRINT OUT
      READ(1,*)MAXR
C     ISTEP IS THE STEPSIZE FOR PRINTING THE R. ALSO GIVES THE DR
C     USED IN THE CALCULATIONS.
      READ(1,*)ISTEP
      READ(1,*)DELTA,ERR
      READ(1,*)AA1,BB1
      READ(1,*)AA2,BB2
C  SOFTNESS PARAMETER OF THE BS POTENTIAL
      READ(1,*) AX
      CLOSE(1)
c      WRITE(25,455)T,DENA,RC1,Z1,AX
  455 FORMAT(1X,'T=',F9.2,2X,'DENA=',F7.5,2X,'RC1=',F6.3,2X,'Z1=',F5.2,
     1  2X,'AX=',F6.3/)
     
           
C     READ(12,*) AX
C     CLOSE(12)
      WRITE (52,*)' INTERATOMIC POTENTIALS USING BS PSEUDOPOT.'
        WRITE (52,*)' '
      WRITE(52,221) T
 221  FORMAT(1X,'FINAL TEMPERATURE=',E15.6)
      WRITE(52,231)RC1
 231   FORMAT(1X,'RC1=',2X,E15.6,'A.U.')
      WRITE(51,23)RC1
 23   FORMAT(1X,E15.6)
C
C     BOLTZMAN'S CONSTANT (IN  A.U./K)
C
      CKB=0.00000316674354
      CON=0.52917706
      CON3=CON*CON*CON
C     DENSITY IN A.U.
      DEN=DENA*CON3
      WRITE(51,15)DENA
 15   FORMAT(1X,E15.6)
      WRITE(52,151)DENA
 151  FORMAT(1X,'DENSITY=',2X,E15.6,'ANSGTROMS')

 16   FORMAT(1X,I5,2X,I5,2X,E15.6)
C
      WRITE(52,999)
  999 FORMAT('   ********** PR12 --- ANALYTIC VCORR --- ******')
      WRITE(52,888)
  888 FORMAT(//'              ****  GUNNARSON-LUNDQVIST ***')
      PAI=3.141592654
      PAI2=PAI*PAI
      PAI4=4.0*PAI
      C2=1.0-C1
      DENEL=DEN*(C1*Z1+C2*Z2)
      WRITE(52,5) C1,Z1,Z2
    5 FORMAT('  *** C1=',F7.3,5X,'Z1=',F5.1,5X,'Z2=',
     &       F5.1,'  ***'/)
C
C     FERMI WAVENUMBER
C
      QF=(3.0*PAI2*DENEL)**(1.0/3.0)
      RS=(3.0/(4.0*PAI*DENEL))**(1.0/3.0)
      QF2=2.0*QF
      VMQ=FLOAT(MAXQ)
      VQ0=-(2.0/3.0)*(QF*QF/2.0)
      WRITE(*,228)VQ0
  228  FORMAT(5X,F12.9) 
C
      PI2=2.0*PAI
      VISTEP=FLOAT(ISTEP)
      DRAU=DRIN/(VISTEP*CON)
C     DRAU=(DRIN/VISTEP)*CON
      VN=VMQ*QF*DRAU/PI2
      NDQ=INT(VN)
      VNN=FLOAT(NDQ)
C     DRAU=PI2*VNN/(VMQ*QF)
C     DR IN A.U.
      DR=DRAU
C     DQ=QF/VNN
C     DQ IS IN (A.U.)**(-1)
      DQ = PI2/(MAXQ*DR)
C
      N2QF=2*NDQ+1
C     DRAM IS DR IN ANSGTROMS
      DRAM=DR*CON
      WRITE(51,16)MAXR,ISTEP,DRAM
      WRITE(52,161)MAXR,ISTEP,DRAM
 161  FORMAT(1X,'MAXR=',I5,2X,'ISTEP=',I5,2X,'DR=',E15.6,'ANSGTROMS')


      WRITE(52,10) DEN,T,QF,RS
   10 FORMAT('  *** DEN=',E12.4,5X,'T=',F6.1,5X,'QF=',
     &       F10.5,5X,'RS=',F10.5)
      WRITE(52,12) DQ,N2QF,MAXQ
  12  FORMAT(/'  DQ= ',F10.5,5X,' N2QF= ',I5,5X,'MAXQ= ',I5)
C
C     VALUES OF SG(0)
C
      RS0(1)=RS
      CALL SGZERO(RS0,SG00)
      SG0=SG00(1)
      WRITE(52,20) SG0
   20 FORMAT(/'     **** SG(0)=',F10.5/)
C
C     G(Q) --- ICHIMARU - UTSUMI
C
      CALL COEF(RS,SG0,GAM0,B,C)
      WRITE(52,40) GAM0
   40 FORMAT(20X,'GAM0=',E15.6/)
C
C     Q -- Q IN A.U.
C     QV=Q/QF
C
      G(1)=0.0
      Q(1)=0.0
      DO 30 K=2,MAXQ
      Q(K)=FLOAT(K-1)*DQ
      QV=Q(K)/QF
      CALL IUG(GAM0,B,C,K,N2QF,QV,G0)
      G(K)=G0
   30 CONTINUE
C
C     CORRELATION POTENTIAL
C     GUNNARSON-LUNDQVIST
C      ANALYTIC FORM FOR VCORR --- GAUSSIAN  ----
C
C
C     CALCULATION OF INDIRECT PART
C
C     COMPUTATION OF THE FORM FACTOR
      VQ1(1)=PAI*PAI/QF
      VQ2(1)=VQ1(1)
C      SQ0=-(2.0/3.0)*(QF*QF/2.0)
C      WRITE(*,*)VQ0
      IS1=ISTEP+1
      DO 100 I=IS1,MAXR,ISTEP
      R(I)=FLOAT(I-1)*DR
      DO 200 J=1,MAXQ
      QR=Q(J)*R(I)
      SINQR=DSIN(QR)
      IF(J.EQ.1) GO TO 201
      GO TO 210
C
C     *** Q=0 ***
C
  201 F11(1)=-1.0
      F12(1)=-1.0
      F22(1)=-1.0
      GO TO 200
C
C     *** Q > 0 ***
C
  210 SQR=SINQR/QR
      Q2=Q(J)*Q(J)
      Q4=Q2*Q2
      X=Q(J)/QF2
C
C     *** PSEUDOPOTENTIAL ***
C
C     ASCROFT'S PSEUDOPOTENTIAL(AA1=AA2=0.0, BB1 & BB2 ANYTHING. OTHERWISE
C     YOU WILL GET YOUNG'S PSEUDOPOTENTIAL
C
C      V1=-DCOS(Q(J)*RC1)*ASHN(AA1,BB1,RC1,Q(J))
C      V2=-DCOS(Q(J)*RC2)*ASHN(AA2,BB2,RC2,Q(J))
C
C---------------->BRETONNET-SILBERT POTENTIAL
C     AX=RC1/RATIO
      RCA=RC1/AX
      
      B1=Z1/RC1*(1.0-2.0/RCA)*EXP(RCA)
      B2=2.0*Z1/RC1*(1.0/RCA-1.0)*EXP(RCA/2.0)
      AQ=AX*Q(J)
      AQS=AQ*AQ
      AQ2=2.0*AQ
      AQ2S=AQ2*AQ2
      QRC=Q(J)*RC1
      SNN=DSIN(QRC)
      CSS=DCOS(QRC)
      AQP=1.0+AQS
      AQM=1.0-AQS
      AQ2P=1.0+AQ2S
      AQ2M=1.0-AQ2S
      RCA2=RCA/2.0
      AJ1=2.0-EXP(-RCA)*((RCA*AQP+AQM)*SNN/AQ+(2.0+RCA*AQP)*CSS)
      AJ2=2.0-EXP(-RCA2)*((RCA2*AQ2P+AQ2M)*SNN/(2.0*AQ)+(2.0+RCA2*AQ2P)*
     1     CSS)
      AQPS=AQP*AQP
      AQ2PS=AQ2P*AQ2P
      W1=PAI4*DEN
      W2=B1*AJ1/AQPS+8.0*B2*AJ2/AQ2PS
      W3=W1*Z1/Q2
      W(J)=W1*AX**3*W2-W3*CSS

      V1=Q2*AX**3/Z1*(B1*AJ1/AQPS+8.0*B2*AJ2/AQ2PS)-CSS
      V2=V1
C	     --------------END OF BRETONNET'S POTENTIAL

      IF(J.EQ.N2QF) GO TO 215
      GO TO 220
C
C     *** Q=2*QF ***
C
  215 X0=-0.5/(PAI*QF)
      EPS1=X0/(X*X-(1.0-G(J))*X0)
      F11(J)=SQR*V1*V1*EPS1
      F12(J)=SQR*V1*V2*EPS1
      F22(J)=SQR*V2*V2*EPS1
C     CALCULATION OF THE FORM FACTORS
      EPS=1.0-(PAI4/Q2)*(-0.5*QF/PAI2)*(1.0-G(J))
      VQ1(J)=V1*(PAI4/Q2)/EPS/VQ1(1)
      VQ2(J)=V2*(PAI4/Q2)/EPS/VQ2(1)
      GO TO 200
C
C     *** Q >< 2*QF ***
C
  220 ABX=ABS((1.0+X)/(1.0-X))
      X0=-(1.0/(PAI*QF))*(0.5+((1.0-X*X)/(4.0*X))
     &       *DLOG(ABX))
      EPS1=X0/(X*X-(1.0-G(J))*X0)
      F11(J)=SQR*V1*V1*EPS1
      F12(J)=SQR*V1*V2*EPS1
      F22(J)=SQR*V2*V2*EPS1
C     CALCULATION OF THE FORM FACTORS
      EPS=1.0-(PAI4/Q2)*(QF*QF/PAI)*X0*(1.0-G(J))
C     FOR ATOMIC FORM FACTORS:
      EPS2=1.0-(PAI4/Q2)*(QF*QF/PAI)*X0
      VQ1(J)=V1*DEN*(PAI4/Q2)/EPS2/VQ1(1)
      VQ2(J)=V2*DEN*(PAI4/Q2)/EPS2/VQ2(1)
C      VQ3(J)=V1*DEN*(PAI4/Q2)*Z1/EPS2
      VQ3(J)=W(J)/EPS2
C      VQ4(J)=V2*DEN*(PAI4/Q2)*Z2/EPS2	
  200 CONTINUE
C
      CALL TRAP(MAXQ,F11,DQ,FSUM11)
      CALL TRAP(MAXQ,F22,DQ,FSUM22)
      CALL TRAP(MAXQ,F12,DQ,FSUM12)
C
C     *** VIND --- INDIRECT PART ***
C     *** VDIR --- DIRECT COULOMB ***
C
      VI11(I)=(2.0/PAI)*FSUM11
      VI12(I)=(2.0/PAI)*FSUM12
      VI22(I)=(2.0/PAI)*FSUM22
C           WILLS-HARRISON  POTENTIAL------------

C      ZD = 8.5
C      RD = 1.342
C      AP1 = 3.1415
C      AP2 = AP1*AP1
C      DBAND = ZD*225.0/AP2*RD**6/R(I)**8
C     >      - ZD*(1.0-ZD/10.0)*28.1/AP1*RD**3/R(I)**5

      V11(I)=Z1*Z1*(1.0/R(I)+VI11(I))
      V12(I)=Z1*Z2*(1.0/R(I)+VI12(I))
      V22(I)=Z2*Z2*(1.0/R(I)+VI22(I))
C
  100 CONTINUE
      WRITE(52,510)
C 510 FORMAT(/'   R   ',13X,'V11(R)',14X,'V12(R)',14X,'V22(R)',
C     &       14X,'V0RD(R)'//)
  510 FORMAT(/,4X,'  I  ',8X,'R (A.U.)',12X,'V11(R)/KBT',
     &       12X,'V22(R)',12X,'V12(R)'//)
      DO 500 J=IS1,MAXR,ISTEP
      VKB11=V11(J)/CKB
      VKB22=V22(J)/CKB
      VKB12=V12(J)/CKB
C      R(JN)=R(J)*CON
      WRITE(51,520)J,VKB11,VKB22,VKB12
  520 FORMAT(I5,1X,3(5X,E15.6))
      VKB11T=VKB11/T
      VKB22T=VKB22/T
      VKB12T=VKB12/T
C     R(JN)=R(J)*CON
      WRITE(52,521)J,VKB11T,VKB22T,VKB12T
  521 FORMAT(1X,I5,1X,4(5X,E15.6))

  500 CONTINUE
C
C     LISTING OF THE FORM FACTORS. WE USE IT FOR COMPUTING RESISTIVITY
C     SO WE ONLY NEED A FEW VALUES
      WRITE(54,609)RC1
  609 FORMAT(5X,'RC=',E20.10,'A.U.',/)
      WRITE(54,6090)
 6090 FORMAT(5X,'FORM FACTOR OF THE BS PSEUDOPOTENTIAL',/)
      DO 600 I = 1,200
      WRITE(53,610)I,Q(I)/CON,VQ3(I)
      WRITE(54,610)I,VQ1(I)
  610 FORMAT(I5,5X,F10.5,5X,F11.6,5X,F11.6)
  600 CONTINUE

c      WRITE(25,466)
c  466 FORMAT(5X,'R(A)',9X,'V(R)(AU)')
      DO 444 JN = IS1,MAXR,ISTEP
      R(JN)=R(JN)*CON
      WRITE(25,555) R(JN),V11(JN)
  555 FORMAT(2X,F7.3,2X,E14.6)
  444 CONTINUE

      STOP
      END
C
C
C     ***************************************************
C
      SUBROUTINE SGZERO(RS0,SG00)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION RS0(1),SG00(1),RS(7),SG0(7)
      DATA RS/1.0,2.0,3.0,4.0,5.0,6.0,10.0/
      DATA SG0/0.279,0.181,0.128,0.094,0.070,0.052,0.011/
C
      CALL S3N(RS,SG0,7,RS0,SG00,1)
      RETURN
      END
C
C
C     *****************************************************
C
      SUBROUTINE COEF(RS,SG0,GAM0,B,C)
      IMPLICIT REAL*8 (A-H,O-Z)
      DATA A/0.029/,B0/0.0621814/,B1/9.81379/,B2/2.82224/,
     &     B3/0.736411/
      PAI=3.141592654
      AL=(4.0/(9.0*PAI))**(1.0/3.0)
      X=DSQRT(RS)
      GAM0=0.25+(PAI*AL*B0/48.0)*(6.0*X*X*(1.0+B1*X)/
     &     (1.0+B1*X+B2*X*X+B3*X*X*X)+X**4*(2.0*B2+
     &     (B1*B2+3.0*B3)*X+2.0*B1*B3*X*X)/(1.0+B1*X+
     3     B2*X*X+B3*X*X*X)**2)
      B=(9.0/16.0)*GAM0-(3.0/64.0)*(1.0-SG0)-(16.0/15.0)*A
      C=-(3.0/4.0)*GAM0+(9.0/16.0)*(1.0-SG0)-(16.0/5.0)*A
      RETURN
      END
C
C
C     ******************************************************
C
      SUBROUTINE IUG(GAM0,B,C,K,N2QF,Q,G0)
C     CALCULATION OF LOCAL FIELD CORRECTION
C     ICHIMARU-UTSUMI
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DATA A/0.029/,B0/0.0621814/,B1/9.81379/,B2/2.82224/,
     &     B3/0.736411/
C
      Q2=Q*Q
      Q4=Q2*Q2
      IF(K.EQ.1) GO TO 100
      IF(K.EQ.N2QF) GO TO 200
      AQ=DABS((2.0+Q)/(2.0-Q))
      G0=A*Q4+B*Q2+C+(A*Q4+(B+(8.0/3.0)*A)*Q2-C)*
     1   ((4.0-Q2)/(4.0*Q))*DLOG(AQ)
      GO TO 300
  100 G0=0.0
      GO TO 300
  200 G0=A*Q4+B*Q2+C
  300 RETURN
      END
C
C
C     **************************************************
C
C
      SUBROUTINE S3N(X,Y,N,XX,YY,NN)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(10),Y(10),XX(10),YY(10)
      DIMENSION SM(10),H(10),ALM(10),AMU(10)
      DIMENSION C(10),P(10),Q(10),U(10)
      N1=N-1
      AMU(1)=0.0
      ALM(N)=0.0
      C(1)=0.0
      C(N)=0.0
C
      DO 110 I=2,N
  110 H(I)=X(I)-X(I-1)
      DO 120 I=2,N1
      ALM(I)=H(I+1)/(H(I)+H(I+1))
  120 AMU(I)=1.0-ALM(I)
      DO 130 I=2,N1
  130 C(I)=3.0*(ALM(I)*(Y(I)-Y(I-1))/H(I)+AMU(I)*(Y(I+1)-Y(I))/H(I+1))
      P(1)=2.0
      Q(1)=-AMU(1)/P(1)
      U(1)=C(1)/P(1)
      DO 140 K=2,N
      P(K)=ALM(K)*Q(K-1)+2.0
      Q(K)=-AMU(K)/P(K)
  140 U(K)=(C(K)-ALM(K)*U(K-1))/P(K)
      SM(N)=U(N)
      DO 150 K=1,N1
      K1=N1-K+1
  150 SM(K1)=Q(K1)*SM(K1+1)+U(K1)
      DO 160 I=1,NN
      XXI=XX(I)
      DO 170 K=2,N
      IF(XXI.GT.X(K)) GO TO 170
      J1=K
      GO TO 180
  170 CONTINUE
  180 J=J1-1
      SMJ=SM(J)
      SMJ1=SM(J1)
      YJ=Y(J)
      YJ1=Y(J1)
      HJ1=H(J1)
      XJ1=X(J1)-XXI
      XJ=XXI-X(J)
      HJ2=HJ1*HJ1
      HJ3=HJ2*HJ1
      YY(I)=SMJ*XJ1*XJ1*XJ/HJ2-SMJ1*XJ*XJ*XJ1/HJ2+YJ*XJ1*XJ1*
     &      (2.0*XJ+HJ1)/HJ3+YJ1*XJ*XJ*(2.0*XJ1+HJ1)/HJ3
  160 CONTINUE
      RETURN
      END
C
C     ***************************************************
C
C
      SUBROUTINE TRAP(M,FN,DQ,SUM)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION FN(M)

      TEND=FN(1)+FN(M)
      M1=M-1
      SUM=0.0

      DO 100 I=2,M1
        SUM=SUM+FN(I)
  100 CONTINUE

      TS=TEND*0.50+SUM
      SUM=TS*DQ

      RETURN
      END
C
C
      FUNCTION ASHN(A,B,R,Q)
      IMPLICIT REAL*8 (A-H,O-Z)
C
      Q2=Q*Q
      B2=B*B
      X=A*Q2/(Q2+B2)
      Y=1.0+TAN(Q*R)*B/Q
      Z=DEXP(-B*R)
      ASHN=1.0+X*Y*Z
      RETURN
      END
