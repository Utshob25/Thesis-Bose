c     PROGRAM SURFACE TENSION, VISCOSITY AND SELF DIFFUSION COEFFICIENT 
      IMPLICIT REAL*8(A-H,O-Z)
 
      OPEN(1,FILE="SEINVTY.DAT",STATUS="OLD")       
cc    ATM IS ATOMIC MASS IN amu
      READ(1,*)T,ADEN,ETA,ATM
c     *******************************************************
c     * HARD SPHERE SURFACE TENSION CALCULATION OF REISS'S  *
c     * EQUATION                                            *
c     *******************************************************

cc    SURFACE TENSION IN N/m

      PI=3.141592654
      BC=1.38
      DNDT1=DNDT*1E-6
      PI6=PI/6.0 
cc    ETA=(PI/6.0)*(ADEN)*((SIG)**3.0)
      SIG=((6.0/PI)*(ETA/ADEN))**(.333)

      SIGM=SIG*(1.126*(1-0.112*(T/TM)**0.5))

      DSIGDT=-(0.063*SIGM)/((T*TM)**0.5)

      ETASTR1=(2.0*ETA*(4.0-ETA))+(1.0-ETA)**4.0
      ETASTR2=(2.0-3.0*ETA+ETA**3.0)
      ETASTR4=4.0*((1.0+2.0*ETA)**2.0)*((1.0-ETA)**4.0)

      ETARUHUL=(ETASTR1*ETASTR2)/ETASTR4

      GAMMARUHUL=BC*T*ADEN*ETARUHUL*0.001*SIG

      ETASTR5=9.0*(ETA**2.0)*(1.0+ETA)
      ETASTR6=2.0*PI*((1.0-ETA)**3.0)
      

      ETAREISS=ETASTR5/ETASTR6
      GAMMAREISS=(BC*T*ETAREISS*0.001)/(SIG**2.0)
  
      ETASTR7=4.0*((1.0-ETA)**4.0)

      ETARATAN=ETASTR2/ETASTR7

      GAMMARATAN=BC*T*ADEN*ETARATAN*0.001*SIG

c     *****************************************************************
c     *              JELLIUM METAL SURFACE ENERGY                     *
c     *****************************************************************       
      AS=-0.57193
      A4=0.97866
      AX=-0.33878
      AC=-1.2754

c     EL IS THE ELECTRON DENSITY PROFILE PARAMETER  WITH UNIT OF a.u

      EL= 0.994
      CON=(3.0*ADEN*PI**2.0)**(1.0/3.0)
      
      TENJEL1= (3.0/10.0)*(CON**2.0)*ADEN*AS*EL
      TENJEL2= (0.693/72.0)*(ADEN/EL)
      TENJEL3= (1.0/(1620.0*PI**2.0))*CON*(A4/EL**3.0)
      TENJEL4= -(3.0/4.0*PI)*CON*ADEN*AX*EL
      TENJEL5= -(0.0155/3.0)*ADEN*AC*Le+(PI/2.0)*(ADEN*EL)**2.0*EL
      TENJEL=TENJEL1 + TENJEL2 +TENJEL3+ TENJEL4+ TENJEL5
              




c     CON2 IS CONSTANT WHICH IS EQUALS TO (16/15)*(1.66*1E-27/1E-23)**0.5 
c     THIS PROGRAM CALCULATES VISCOSITY (mPa.s)
      
      CON2=13.74304
      VISCORUHUL=GAMMARUHUL*CON2*((ATM/(BC*T))**0.5)
      VISCOREISS=GAMMAREISS*CON2*((ATM/(BC*T))**0.5)
      VISCORATAN=GAMMARATAN*CON2*((ATM/(BC*T))**0.5) 
c     SELF DIFFUSION COEFFICIENT IN ((m2/s)*1E-9)
c     CON3 IS CALCULETED FROM (15/(32*PI))*((1E-23)**1.5/(1.66E-27)**0.5)

      CON3=(15/(32*PI*1E-10))*((1E-23)**1.5/(1.66E-27)**0.5)

      DIFFRUHUL=(CON3*((BC*T)**1.5/ATM**0.5))/(SIG*GAMMARUHUL)
      DIFFREISS=(CON3*((BC*T)**1.5/ATM**0.5))/(SIG*GAMMAREISS)
      DIFFRATAN=(CON3*((BC*T)**1.5/ATM**0.5))/(SIG*GAMMARATAN) 


      PRINT*,GAMMARUHUL,GAMMAREISS,GAMMARATAN
      PRINT*,VISCORUHUL,VISCOREISS,VISCORATAN
      PRINT*,DIFFRUHUL,DIFFREISS,DIFFRATAN
      PRINT*,ETA,TENJEL,SIG
      STOP
      END
