      PROGRAM DIFFUSION CONSTANT and Viscosity USING Hard sphere self consistently 

      IMPLICIT REAL*8 (A-H,O-Z)
      

      OPEN(1,FILE="SEINVHS.DAT",STATUS="OLD") 
      READ(1,*)T,ADEN,ETA,ATM
      PI=3.141592654
      BC=1.38
      CONE=0.01292
      COND=7.738
      SIG=((6.0/PI)*(ETA/ADEN))**(.333)
      
      D1=8.0*ADEN*SIG
      D2=5.0*(PI)
      
      D3=(PI*ATM*BC*T)**(.5)
      D4=ETA*(1.0-(ETA/2.0))
      D5=(1.0-ETA)**(3.0)
      VIS=(D1/D2)*(D3*D4)*(1.0/D5)*CONE
      
      C1=(1.0/16.0)
      C2=((PI*BC*T)/ATM)**(.5)
      C3=SIG*(1.0-(ETA/2.0))
      C4=ETA*((1.0-ETA))**(3.0)
      C5=0.33/ETA
      DIF=(COND*(C1*C2*C3*(1.0/C4)))+C5
      
      PRINT*,DIF,VIS

      STOP
      END
C
