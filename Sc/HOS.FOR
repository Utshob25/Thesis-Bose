C     Last change:  BIV   8 Jan 2011   11:30 am
       PROGRAM  HOS
c--------------------------------------------------------------------*   
c This program calculates the EAM effective pair potentials in a     *
c form (units) suitable for the VMHNC program.                       *
c Input: embedding energies (F) are in eV; A in ev; alp in (1/Angs); *
c Rc in Angs.; bta in (1/Angs.) ; gama=2*bta;                        *
c XBAR is the average equilibrium electron density at site i         *
c due to the neighbours (without atom i)                             * 
c--------------------------------------------------------------------*

      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER*30 BAR
      PARAMETER(MMX=3000)
      DIMENSION R(MMX),SAI(MMX), RR(MMX),SAIA(MMX)
      
      OPEN(1, FILE='HOS.DAT',  STATUS='OLD')
      OPEN(3, FILE='PSEUDOPT.DAT',STATUS='OLD')
      OPEN(51, FILE='FOR051.DAT',STATUS='unknown')
      OPEN(53,FILE='FOR053.DAT' ,STATUS = 'unknown')
      OPEN(66, FILE='FOR066.DAT',STATUS='unknown')
      OPEN(67, FILE='BRETEAM.DAT',STATUS='unknown')

      READ(1,*) X0,X1,X2,X3,X4
      READ(1,*) F0,F1,F2,F3,F4
      READ(1,*) XBAR
      READ(1,*) A,ALP
      READ(1,*) RC
      READ(1,*) BTA,GAMA
      READ(1,*) Z
      CLOSE(1)

      READ(3,*) TIN, T
c   -------------------     Temperature in Kelvin
      READ(3,*) RC1,RC2
      READ(3,*) Z1,Z2
c   --------------------   RC1,RC2,Z1,Z2 are not used in the calculations
      READ(3,*) DENA,C1
c-----------------------      DENA-- ion no den in Angstrom units 
      READ(3,*) MAXQ,DRIN
c   -------------------------- DRIN is the stepsize in r-space in Angs.
      READ(3,*) MAXR
      READ(3,*) ISTEP
      READ(3,*) DELTA,ERR
      READ(3,*) AA1,BB1
      READ(3,*) AA2,BB2
      CLOSE(3)

c -------------------------calculation of QF
      PAI=3.141592654
      PAI2=PAI*PAI
      PAI4=4.0*PAI
      C2=1.0-C1
      CON = 0.52917706
      CON3=CON*CON*CON
      DEN=DENA*CON3
      DENEL=DEN*(C1*Z1+C2*Z2)
      QF=(3.0*PAI2*DENEL)**(1.0/3.0)

C________________________________________


      vistep = float(istep)
      con = 0.52917706
      drau = drin/(vistep*con)
      dr = drau
      dram = dr*con
      WRITE(51,26) RC1
      WRITE(51,26) DENA
   26 FORMAT(1X,E15.6)
      WRITE(51,27) MAXR,ISTEP,DRAM
   27 FORMAT(1X,I5,2X,I5,2X,E15.6)

      PAI = 3.14159265

      CALL FFUN(X0,X1,X2,X3,X4,F0,F1,F2,F3,F4,XBAR,QX,F11,F12)
C      TYPE*,QX,F11,F12
      WRITE(*,*)QX,F11,F12   
c-----------program for the pair potential-------

C      IS1 = ISTEP+1

      DO 10 I = 1,MAXR,ISTEP
      R(I) = FLOAT(I-1)*DRIN
      D1 = R(I)
      CALL REPOT(D1,RC,A,ALP,FAI,FAI1,FAI2)

      CALL DENFUN(D1,RC,BTA,GAMA,ROA,RO1,RO2)

      ATSAI = 2.0*F11*ROA + F12*ROA*ROA

      SAI(I) = FAI + ATSAI

C  ------------------------- CONVERSION TO THE ATOMIC UNIT--------
      RR(I) = R(I)*CON
      SAIA(I) = SAI(I)*0.036768


C------------------------BOLTZMANN CONSTANT IN a.u.-------
      BK = 3.16674354E-6

      VKB11 = SAIA(I)/BK
      VKB22 = VKB11
      VKB12 = VKB11
      
      WRITE(51,55)I, VKB11,VKB22,VKB12
      WRITE(53,55)I, VKB11,VKB22,VKB12
      WRITE(66, 57)R(I),VKB11
   10 CONTINUE

      A=1.892
      B=2.398
      C=3.499

C      WRITE(66,'(A)')BAR
C      WRITE(66,666)A,B,C
C      WRITE(66,*)1
C      WRITE(66,777)Z1,QF,MAXR-1
C  666 FORMAT(1X,3F6.3)
C  777 FORMAT(1X,2F7.4,I5)

      DO 100 J = 1,MAXR
C     WRITE(66,56) R(J),SAI(J)
      WRITE(67,57) R(J),SAIA(J)
  100 CONTINUE

   55 FORMAT(I5,1X,3(5X,E15.6))
   56 FORMAT(1x,2E15.8)
C   57 FORMAT(1x,E15.6,2x,E15.6,2X,E15.6)
  57   FORMAT(1X,F7.3,2X,2(E14.6,1x))
C   56 FORMAT(F10.6,1X,E20.12)
      STOP
      END

C ------------PROGRAM TO CALCULATE FACTORIAL

      SUBROUTINE FAC(N,RE)
      IMPLICIT REAL*8 (A-H,O-Z)
       
      IFAC = 1

      DO 10 J = 1,N
      IFAC = IFAC*J
   10 CONTINUE
      RE = FLOAT(IFAC)
      RETURN
      END

c----SUBROUTINE FOR EAM FUNCTION F(rho) AND ITS 1ST AND 2ND DERIVATIVES--


      SUBROUTINE FFUN(X0,X1,X2,X3,X4,F0,F1,F2,F3,F4,XBAR,QX,F11,F12)

      IMPLICIT REAL*8 (A-H,O-Z)
        
      H0 = X1-X0
      H1 = X2-X1
      H2 = X3-X2
      H3 = X4-X3

      DY0 = (F1-F0)/H0
      DY1 = (F2-F1)/H1
      DY2 = (F3-F2)/H2
      DY3 = (F4-F3)/H3

      SIG0 = 0.0

      SNUM = 6.0*(DY2-DY1) - 12.0*(H1+H2)*(DY1-DY0)/H1
     >    + 3.0*H2*(H2*(DY1-DY0)-H1*(DY3-DY2))/(H1*(H2+H3))

      SDENU = H1- 4.0*(H1+H2)*(H0+H1)/H1
     >       + H2*H2*(H0+H1)/(H1*(H2+H3))

      SIG1 = SNUM/SDENU

      SIG2 = 6.0*(DY1-DY0)/H1 - 2.0*(H0+H1)/H1*SIG1

      SIG3 = H2*(H0+H1)/(H1*(H2+H3))*SIG1
     >     - 3.0/(H1*(H2+H3))*(H2*(DY1-DY0)-H1*(DY3-DY2))
      SIG4 = 0.0

C      WRITE(23,*) SIG0,SIG1,SIG2,SIG3,SIG4

      C10 = DY0 - H0*(SIG1+2.0*SIG0)/6.0
      C11 = DY1 - H1*(SIG2+2.0*SIG1)/6.0
      C12 = DY2 - H2*(SIG3+2.0*SIG2)/6.0
      C13 = DY3 - H3*(SIG4+2.0*SIG3)/6.0

C      WRITE(23,*) C10,C11,C12,C13

      C20 = SIG0/2.0
      C21 = SIG1/2.0
      C22 = SIG2/2.0
      C23 = SIG3/2.0

C       WRITE(23,*) C20,C21,C22,C23

      C30 = (SIG1-SIG0)/6.0/H0
      C31 = (SIG2-SIG1)/6.0/H1
      C32 = (SIG3-SIG2)/6.0/H2
      C33 = (SIG4-SIG3)/6.0/H3
      
C      WRITE(23,*) C30,C31,C32,C33

      x = XBAR

      IF( (X .LT. X1) .AND. (X .GT. X0) ) GO TO 5
      IF( (X .LT. X2) .AND. (X .GE. X1) ) GO TO 15
      IF( (X .LT. X3) .AND. (X .GE. X2) ) GO TO 25
      IF( (X .LT. X4) .AND. (X .GE. X3) ) GO TO 35

    5 QX = F0+C10*(X-X0)+C20*(X-X0)**2+C30*(X-X0)**3
      F11 = C10+2.0*C20*(X-X0)+3.0*C30*(X-X0)**2
      F12 = 2.0*C20+6.0*C30*(X-X0)

       GO TO 45

   15 QX = F1+C11*(X-X1)+C21*(X-X1)**2+C31*(X-X1)**3
      F11 = C11+2.0*C21*(X-X1)+3.0*C31*(X-X1)**2
      F12 = 2.0*C21+6.0*C31*(X-X1)

         GO TO 45

   25 QX = F2+C12*(X-X2)+C22*(X-X2)**2+C32*(X-X2)**3
      F11 = C12+2.0*C22*(X-X2)+3.0*C32*(X-X2)**2
      F12 = 2.0*C22+6.0*C32*(X-X2)

          GO TO 45

   35 QX = F3+C13*(X-X3)+C23*(X-X3)**2+C33*(X-X3)**3
      F11 = C13+2.0*C23*(X-X3)+3.0*C33*(X-X3)**2
      F12 = 2.0*C23+6.0*C33*(X-X3)

   45 XM = 0.0 
  
      RETURN
      END

C..................SUBROUTINE FOR REPULSIVE PART..............

      SUBROUTINE REPOT(R,RC,A,ALP,FAI,FAI1,FAI2)

      IMPLICIT REAL*8(A-H,O-Z)

      D = R
      B = RC
 
      IF(D .GE. B) GO TO 10
      X1 = EXP(-ALP*R)
      XC = EXP(-ALP*RC)
      FAIR = A*X1
      FAIRC = A*XC
      DFAIR = -A*ALP*XC
      FAI = FAIR-FAIRC+(RC-R)*DFAIR

      DF = -A*ALP*X1
       FAI1 = DF-DFAIR

       FAI2 = A*ALP**2*X1
       GO TO 20

   10 FAI = 0.0
      FAI1 = 0.0
      FAI2 = 0.0

   20 RETURN
      END
C-----------------------------------------------------

      SUBROUTINE PDENFUN(R,RC,BTA,RO,RO1,RO2)

      IMPLICIT REAL*8(A-H,O-Z)
      A = R
      B = RC
      IF(A .GT. B) GO TO 10
      X1 = EXP(-BTA*R)
      X2 = 512.0*EXP(-2.0*BTA*R)
      ROU = R**6*(X1+X2)

      XC1 = EXP(-BTA*RC)
      XC2 = 512.0*EXP(-2.0*BTA*RC)
      ROC = RC**6*(XC1+XC2)
      DROC = 6.0*RC**5*(XC1+XC2) - BTA*RC**6*(XC1+2.0*XC2)
      RO = ROU-ROC+(RC-R)*DROC

      DR1 = 6.0*R**5*(X1+X2) - BTA*R**6*(X1+2.0*X2)
      RO1 = DR1-DROC
      RO2 = 30.0*R**4*(X1+X2) - 12.0*BTA*R**5*(X1+2.0*X2)
     >       +BTA**2*R**6*(X1+4.0*X2)
       GO TO 20
   10 RO = 0.0
      RO1 = 0.0
      RO2 = 0.0
      
   20 RETURN
      END      

      SUBROUTINE DENFUN(R,RC,BTA,GAMA,ROS,ROS1,ROS2)
      
      IMPLICIT REAL*8(A-H,O-Z)

      A = R
      B = RC
      IF(A .GE. B) GO TO 10
      X1 = EXP(-BTA*R)
      X2 = 512.0*EXP(-GAMA*R)
      RO = R**6*(X1+X2)

      XC1 = EXP(-BTA*RC)
      XC2 = 512.0*EXP(-GAMA*RC)
      ROC = RC**6*(XC1+XC2)
      DRO = 6.0*RC**5*(XC1+XC2)
     >      -RC**6*(BTA*XC1+GAMA*XC2)
      XX = RC-R
      ROS = RO-ROC+XX*DRO

      YA = X1 + X2
      YB = BTA*X1 + GAMA*X2
      RO1 = 6.0*R**5*YA - R**6*YB 
      ROS1 = RO1 - DRO
      ROS2 = 30.0*R**4*YA - 12.0*R**5*YB
     >         + R**6*(BTA*BTA*X1+GAMA**2*X2)
      GO TO 20
   10 ROS = 0.0
      ROS1 = 0.0
      ROS2 = 0.0

C      WRITE(2,*) R,ROS

   20 RETURN
      END
C-----------------------------------------------
      SUBROUTINE  POTEN(R,RC,A,ALP,FAI,FAI1,FAI2)
    
      IMPLICIT  REAL*8(A-H,O-Z)

      IF(R .GT.RC) GO TO 10
      X1 = RC-R
      XA = X1**4
      XB = EXP(-ALP*R)
      FAI = A*XA*XB
      FAI1 = -4.0*A*X1**3*XB - ALP*FAI
      FAI2 = 12.0*A*X1**2*XB + 4.0*A*ALP*X1**3*XB - ALP*FAI1

      GO TO 20
   10 FAI = 0.0
      FAI1 = 0.0
      FAI2 = 0.0
   20 XXX = 0.0

      RETURN
      END
c--------------------------------------------

      SUBROUTINE POLPOT(R,RC,A,Z,FAI,FAI1,FAI2)

      IMPLICIT REAL*8(A-H,O-Z)

      X = R
      Y = RC
      IF(R .GT. RC) GO TO 10
      P = RC-R
      PZ = P**Z
      FAI = A*PZ/R
      Q1 = Z-1.0
      Q2 = Z-2.0

      FAI1 = -(A*Z*P**Q1/R+ FAI/R)
      FAI2 = A*Z*Q1*P**Q2/R + A*Z*P**Q1/R**2 - FAI1/R
     >       + FAI/R**2
      GO TO 20
   10 FAI = 0.0
      FAI1 = 0.0
      FAI2 = 0.0

   20 YY = 1.0
      RETURN
      END
