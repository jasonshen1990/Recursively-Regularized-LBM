program main
implicit none

integer i, j, k
integer, parameter :: nx = 64
integer istep
integer inttemp
integer, parameter :: nsteps = 1000, navg=400
real*8, dimension(nsteps, nx, 12) :: profile
real*8, dimension(nx, 12) :: profile_avg

integer :: ierr
character(len=256) :: line

! real*8 :: Re, visc, H, u0, ystar
! real*8, dimension(12) :: temp_row

! Re = 180.0d0
! visc = 0.000622d0
! H = nx / 2
! u0 = Re * visc / H
! ystar = visc / u0

open(12, file='data/Profiles.txt', status='old', action='read', iostat=ierr)
if (ierr /= 0) then
    write(*,*) "Error opening file: ", ierr
    write(*,*) "File: data/Profiles.txt"
    stop
end if

write(*,*) "Reading data..."
do istep = 1, nsteps
    read(12,*,iostat=ierr) inttemp
    if (ierr /= 0) then
        write(*,*) "Error reading step number at step", istep
        exit
    endif
    write(*,*) "Processing step:", inttemp
    
    do j = 1, nx
        read(12,*,iostat=ierr) profile(istep, j, :)
        if (ierr /= 0) then
            write(*,*) "Error reading data at step", istep, "line", j
            backspace(12)
            read(12,'(A)') line
            write(*,*) "Problematic line: ", trim(line)
            exit
        endif
    end do
    if (ierr /= 0) exit
end do
close(12)

if (ierr == 0) then
    do j = 1, 12
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

    ! open(22, file='profile_forplt_.txt')
    !     do i = 1, nx
    !         ! 将该行数据存入临时数组
    !         temp_row = profile_avg(i, :)

    !         ! 第一列 + 0.5
    !         temp_row(1) = temp_row(1) + 0.5d0
    !         ! 第二列 + 0.5/ystar
    !         temp_row(2) = temp_row(2) + 0.5d0 / ystar
            
    !         ! 写入文件
    !         write(22,460) temp_row
    !     end do
    ! close(22)
    
    write(*,*) "Successfully processed", nsteps, "steps"
else
    write(*,*) "Processing stopped due to errors"
endif

460 format(2x,12(1pe15.6))

end program