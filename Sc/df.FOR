
      IMPLICIT REAL*8(A-H,O-Z)
 
            
      READ(*,*)L,W 
      
      CALL area(L,W,A)
       

      STOP
      END

      SUBROUTINE area(L,W,A)
      IMPLICIT REAL*8(A-H,O-Z)
      h=1
      A=L*W
      

      PRINT*,A
       
      
      CALL E(L+h,W,A,R)
      END SUBROUTINE
    


      
      SUBROUTINE E(L,W,A,R)
      IMPLICIT REAL*8(A-H,O-Z)
      
            
      R=A*A
      PRINT*,R
      END SUBROUTINE


