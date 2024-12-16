
      IMPLICIT REAL*8(A-H,O-Z)
 
            
      READ(*,*)L,W 
      
      CALL area(L,W,A)
       
cc     CALL E(A,R)
      STOP
      END

      SUBROUTINE area(L,W,A)
      IMPLICIT REAL*8(A-H,O-Z)
      h=0.0001
      A=L*W
      

      PRINT*,A
       
      A=(L+h)*W
      CALL E(A,R)
      END SUBROUTINE
    


      
      SUBROUTINE E(A,R)
      IMPLICIT REAL*8(A-H,O-Z)
      
            
      R=A*A
      PRINT*,R
      END SUBROUTINE


