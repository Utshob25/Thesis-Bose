
      IMPLICIT REAL*8(A-H,O-Z)
 
            
      READ(*,*)L,W 
      
      CALL area(L,W,A)
       

      STOP
      END

      SUBROUTINE area(L,W,A)
      IMPLICIT REAL*8(A-H,O-Z)
      h=1
      B=CALL E(L+h,B,R)
      A=L*W
      

      PRINT*,A
      
      END SUBROUTINE
    


      
      SUBROUTINE E(L+h,B,R)
      IMPLICIT REAL*8(A-H,O-Z)
      
            
      B=A*A
      PRINT*,R
      END SUBROUTINE


