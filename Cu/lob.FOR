 
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION RS(7),SG0(7)
      DATA RS/1.0,2.0,3.0,4.0,5.0,6.0,10.0/
      DATA SG0/0.279,0.181,0.128,0.094,0.070,0.052,0.011/
C
      
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
      PRINT*,SG0
      GAM0=0.25+(PAI*AL*B0/48.0)*(6.0*X*X*(1.0+B1*X)/
     &     (1.0+B1*X+B2*X*X+B3*X*X*X)+X**4*(2.0*B2+
     &     (B1*B2+3.0*B3)*X+2.0*B1*B3*X*X)/(1.0+B1*X+
     3     B2*X*X+B3*X*X*X)**2)
      B=(9.0/16.0)*GAM0-(3.0/64.0)*(1.0-SG0)-(16.0/15.0)*A
      C=-(3.0/4.0)*GAM0+(9.0/16.0)*(1.0-SG0)-(16.0/5.0)*A
      RETURN
      END
STOP
      END
