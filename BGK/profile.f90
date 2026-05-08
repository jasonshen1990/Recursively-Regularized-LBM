program main
implicit none

integer i, j, k
integer, parameter :: nx = 199
integer istep
integer inttemp
integer, parameter :: nsteps = 2692, navg=800
real*8, dimension(nsteps, nx, 6) :: profile
real*8, dimension(nx, 6) :: profile_avg

integer :: ierr
character(len=256) :: line

open(11, file='data/Profiles.txt', status='old', action='read', iostat=ierr)
if (ierr /= 0) then
    write(*,*) "Error opening file: ", ierr
    write(*,*) "File: data/Profiles.txt"
    stop
end if

write(*,*) "Reading data..."
do istep = 1, nsteps
    read(11,*,iostat=ierr) inttemp
    if (ierr /= 0) then
        write(*,*) "Error reading step number at step", istep
        exit
    endif
    write(*,*) "Processing step:", inttemp
    
    do j = 1, nx
        read(11,*,iostat=ierr) profile(istep, j, :)
        if (ierr /= 0) then
            write(*,*) "Error reading data at step", istep, "line", j
            backspace(11)
            read(11,'(A)') line
            write(*,*) "Problematic line: ", trim(line)
            exit
        endif
    end do
    if (ierr /= 0) exit
end do
close(11)

if (ierr == 0) then
    do j = 1, 6
        do i = 1, nx
            profile_avg(i,j) = sum(profile(nsteps-navg+1:nsteps, i, j))
        end do
    end do
    
    profile_avg = profile_avg/real(navg)
    
    open(22, file='profile_forplt_.txt')
    do i = 1, nx
        write(22,460) profile_avg(i, :)
    end do
    close(22)
    
    write(*,*) "Successfully processed", nsteps, "steps"
else
    write(*,*) "Processing stopped due to errors"
endif

460 format(2x,6(1pe15.6))

end program