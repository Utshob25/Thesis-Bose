program main
    implicit none
    real :: length, width, area

    length = 4.0
    width = 5.0

    call calculate_area(length, width, area)
    write(*,*) 'Area:', area

    call calculate_cube_of_area(area)
      write(*,*) 'Area cube:',cube
end program main

subroutine calculate_area(length, width, area)
    real, intent(in) :: length, width
    real, intent(out) :: area

    area = length * width
end subroutine calculate_area

subroutine calculate_cube_of_area(area,cube)
    real, intent(in) :: area
    real, intent(out) ::cube
    cube = area**3
end subroutine calculate_cube_of_area
