program main
!#Written by Daniel Buscombe, various times in 2011-2013
!#while at
!#School of Marine Science and Engineering, University of Plymouth, UK
!#and now:
!#Grand Canyon Monitoring and Research Center, U.G. Geological Survey, Flagstaff, AZ 

!#Please contact:
!#dbuscombe@usgs.gov

!#to report bugs and discuss the code, algorithm, collaborations

!#For the latest code version please visit:
!#https://github.com/dbuscombe-usgs

!#See also the project blog: 
!#http://dbuscombe-usgs.github.com/
  implicit none

  character ( len = 32 ) :: input1
  character ( len = 32 ) :: input2
  character ( len = 32 ) :: input3
  character ( len = 32 ) :: input4
  integer ( kind = 4 ) numgrains
  integer ( kind = 4 ) model
  integer ( kind = 4 ) batch
  integer ( kind = 4 ) it_max
  character ( len = 80 ) :: file_out_name

  real, dimension(8)     :: xp, zp, yp 

  xp=(/ -.5, .5, .5, -.5, -.5, .5, .5, -.5 /) 
  yp=(/ -.5, -.5, .5, .5, -.5, -.5, .5, .5 /) 
  zp=(/ -.5, -.5, -.5, -.5, .5, .5, .5, .5 /) 

  call timestamp ( )

  call get_command_argument(1, input1)
  read(input1,*) numgrains

  call get_command_argument(2, input2)
  read(input2,*) model

  call get_command_argument(3, input3)
  read(input3,*) batch

  call get_command_argument(4, input4)
  read(input4,*) it_max

  call get_command_argument(5, file_out_name)

  call csbinproc3d (numgrains, xp, zp, yp)
  !call gen_points (numgrains, model, batch, it_max, file_out_name)
  call timestamp ( )

  stop
end


subroutine csbinproc3d (numgrains, xp, zp, yp)

!*****************************************************************************
! generates points in 3D.

  ! passed variables
  real, dimension(8)     :: xp, zp, yp 

  ! internal variables
  real, dimension(numgrains)     :: x, y, z 
  real, dimension(1)     :: minx, maxx, miny, maxy, minz, maxz, cx, cy, cz
  real, dimension(1)     :: xt, yt, zt 
  integer ( kind = 4 ) i
  
  integer,parameter :: seed = 86456
          
  !find the maximum and the minimum for a 'box' around
  !the region. Will generate uniform on this, and throw
  !out those points that are not inside the region.
  minx = minval(xp)
  maxx = maxval(xp)
  miny = minval(yp)
  maxy = maxval(yp)
  minz = minval(zp)
  maxz = maxval(zp)
  cx = maxx-minx
  cy = maxy - miny
  cz = maxz - minz

  call srand(seed)

  i = 1
  do while (i <= numgrains)
       xt=rand()*cx + minx
       yt=rand()*cy + miny
       zt=rand()*cz + minz
       !print *, xt
       i=i+1
  end do

  !while i <= n
  !	xt = rand(1)*cx + minx;
  !	yt = rand(1)*cy + miny;
  !  zt = rand(1)*cz + minz;
  !  k = inhull([xt,yt,zt],[xp,yp,zp]);
  !	if k == 1
  !		% it is in the region
  !		x(i) = xt;
  !		y(i) = yt;
  !      z(i) = zt;
  !		i = i+1;
  !	end
  !end

  return
end







