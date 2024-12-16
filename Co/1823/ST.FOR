      PROGRAM SURFACE TENSION
CC    V(R) IN A.U, SURFACE TENSION (RES) IN N/M, R IN ANGS
CC    NUMBER DENSITY (DENA) IN /ANGS**3.0
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(N=199)
      DIMENSION R(N),R4(N),V(N),DV(N),G(N),Q(N),SQ(N),H(N)
      DIMENSION F(N),V1(N)
      OPEN(1,FILE='STRUC.DAT',STATUS='OLD')
      OPEN(2,FILE='STIN.DAT',STATUS='OLD')
      OPEN(3,FILE='BRET.DAT',STATUS='OLD')
    
      OPEN(4,FILE='INT.DAT',STATUS='UNKNOWN')
      READ(1,*)(R(I),G(I),Q(I),SQ(I),H(I),V1(I),I=1,N)
      READ(2,*)DENA
      READ(3,*)(R(I),V(I),I=1,N)
      H = R(4)-R(3)
      CALL DIF(H,N,V,DV)

      DO 16 I =1,N
      R4(I)=R(I)**4.0
      F(I)=DV(I)*G(I)*R4(I)
   16 CONTINUE

      DO 60 I = 1,N
C      WRITE(4,25) R(I),DV(I),F(I)
      WRITE(4,*) R(I),R4(I),G(I),F(I),V(I),DV(I)
   60 CONTINUE
C   25 FORMAT(1X,3(F7.3,2X))

      CALL SIMP(N,H,F,RINT)
C      PRINT*,RINT
C      T=933.0
C      DENA=0.0320
      PI=3.1416
      CON=435.70
      PIB8=PI/8.0
      DENA2=DENA**2.0
      CON1=PIB8*DENA2
      RES=CON*RINT*CON1
      PRINT*,RES
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
