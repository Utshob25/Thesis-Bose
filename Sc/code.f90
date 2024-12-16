program circle
  implicit none
  real :: r, pi

  ! Assign a default value for pi
  pi = 3.14159

  ! Prompt the user for the radius
  print *, "Enter the radius of the circle:"
  read *, r

  ! Calculate and print the circumference
  print *, "The circumference of the circle is:", 2.0 * pi * r

end program circle
