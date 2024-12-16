      PROGRAM EXCESS ENTROPY 
C     R IN ANGS
C     NUMBER DENSITY (DENA) IN /ANGS**3.0
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(N=175)
      DIMENSION R(N),R2(N),V(N),G(N),Q(N),SQ(N),H(N)
      DIMENSION F(N),V1(N)

      OPEN(1,FILE='STRUC1.DAT',STATUS='OLD')
      READ(1,*)(R(I),G(I),Q(I),SQ(I),H(I),V1(I),I=1,N)


      OPEN(3,FILE='PSEUDOPT.DAT',STATUS='OLD')
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

      OPEN(4,FILE='RES1.DAT',STATUS='UNKNOWN')
      
      H = R(4)-R(3)

      DO 16 I =1,N
      R2(I)=R(I)**2.0
      F(I)=G(I)*LOG(G(I))*R2(I)
   16 CONTINUE

      CALL SIMP(N,H,F,RINT)
   
      PI=3.1416
      CON=2.0*PI*DENA
      RES=CON*RINT

      PRINT*,RES
      WRITE(4,*)RES,DENA
 
      STOP
      END
C*******************************************************
      SUBROUTINE DIF(H,N,X,Y)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(N),Y(N)
      DO 10 I = 2,N-1
      Y(I) = (X(I+1)-X(I-1))/(2.0*H)
   10 CONTINUE
      Y(1) = (-X(3)+4.0*X(2)-3.0*X(1))/(2.0*H)
      Y(N) = (3.0*X(N)-4.0*X(N-1)+X(N-2))/(2.0*H)
      RETURN
      END
C*******************************************************
C SUBROUTINE FOR INTEGRATION USING SIMPSON'S ROLE
      SUBROUTINE SIMP(N,H,F,RINT)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION F(N)
      SUM1=0.0
      SUM2=0.0
      DO 77 I=2, N-1, 2
       	SUM1 = SUM1 + F(I)
   77 CONTINUE
      DO 99 J=3, N-2, 2
      	SUM2 = SUM2 + F(J)
   99 CONTINUE	
      RINT=(1.0/3.0)*H*(F(1) + 4.0*SUM1 + 2.0*SUM2 + F(N))
C      RHO=CON*TINT
      RETURN
      END

************************************************************
