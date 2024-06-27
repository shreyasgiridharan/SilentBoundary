# SilentBoundary
A FE implementation of absorbing/silent boundary; water-solid mixture implementation

The main subroutines are written in main_driver_program.f90. Use that as a starting point.

The file interface_files.f90 contains subroutines to call a plaxis dll. #Removed here for copyright reasons.

An example input file is provided in input.dat.

The original program accepts input file only with the name 'input.dat'. You can however change that to whatever you want in main.f90.

You will need GiD to pre- and post process the simulation.

If you encounted any errors, bugs or improvements, please contact shreyas.giridharan@gmail.com.
