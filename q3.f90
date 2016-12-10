program q3code

implicit none

character*100 :: buffer, fname, cellname   
character*10  :: PhaseName	   
integer :: nargc, iargc		    
integer, parameter :: maxneigh=100 
double precision :: rcut, rcut_sq, Sep_2 

!!!!!!!!!!!!!!DCD INPUT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
character*4              :: car4
character*80             :: car(10)
integer                  :: nstart, nsanc, nset, ntitle, count
integer                  :: charm, namin, i5(5), i9(9)
real*4                   :: DD
logical                  :: dcd_has_cell
real*8                   :: cell(6)
integer	                 :: nat, iostat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

integer			 			:: nof,cart, Phase, Nl
integer			 			:: ii, m, i, j, k
integer, allocatable     	:: list_pos(:), neigh(:)
double precision, parameter :: pi=4.0d0*datan(1.0d0), Na=6.02214129*10**23d0  
real*4, allocatable      	:: x(:), y(:), z(:)
real*8, allocatable	 		:: pos(:,:), q_3Tot(:), S_ij(:), S_i(:), S_j(:), q_3TotAv(:)
real*8			 			:: la, lz, y3, y2, y1, y0, y2_pos, NormM 
real*8			 			:: ThickNorm, dist_x, dist_y, dist_z
real*8			 			:: R, y1zxy, delta_thick, r3
real*8			 			:: P_0, P_1, P_2, P_3, P_m1, P_m2, P_m3, cos_theta, sin_theta, rho
real			 			:: densityIce
complex*16, allocatable  	:: q_m3(:), q_m2(:), q_m1(:), q_0(:), q_1(:), q_2(:), q_3(:) 
complex*16, allocatable  	:: q_m3_mag(:), q_m2_mag(:), q_m1_mag(:), q_0_mag(:), q_1_mag(:), q_2_mag(:), q_3_mag(:)  
complex*16, allocatable	 	:: q_m3Av(:), q_m2Av(:), q_m1Av(:), q_0Av(:), q_1Av(:), q_2Av(:), q_3Av(:)
complex, parameter	 		:: imag=(0,-1) !!!!!


nargc=iargc()

if (nargc.lt.4) then
   write(*,*) "Usage:"
   write(*,*) "code.x90 [number of frames] [some .dcd traj] [cell file] [cutoff]"
   stop
else
   call getarg(1,buffer)
   read(buffer,*) nof
   call getarg(2,buffer)
   read(buffer,*) fname 
   call getarg(3,buffer)
   read(buffer,*) cellname
   call getarg(4,buffer)
   read(buffer,*) rcut
   
endif

rcut_sq=rcut**2.0d0


!!!!!!!!!!!!!! READ IN DCD !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
open (unit=999,file=trim(fname),form="unformatted") 
  
read(999) car4, nset, nstart, nsanc, i5, namin, DD, i9, charm   ! nset is total num of frames 
read(999) ntitle, (car(i),i=1,ntitle)
read(999) nat							! number of (oxygen) atoms
 
write(*,*), 'nat', nat, nstart, nsanc, nset

if (nof.eq.-1) nof=nset
	if(i9(1)==1) then
		dcd_has_cell = .true.
	else 
		write(*,'(a57,i8)') "A dcd with pbc cell is required" 
		stop
	endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



cart=3

allocate(neigh(nat))
allocate(list_pos(maxneigh*nat))
allocate(x(nat))
allocate(y(nat))
allocate(z(nat))
allocate(pos(nat,cart))
allocate(S_ij(nat),S_i(nat),S_j(nat)) ! For Dot product method
allocate(q_m3(nat),q_m2(nat),q_m1(nat),q_0(nat),q_1(nat),q_2(nat),q_3(nat),q_3Tot(nat),q_3TotAv(nat))
allocate(q_m3Av(nat),q_m2Av(nat),q_m1Av(nat),q_0Av(nat),q_1Av(nat),q_2Av(nat),q_3Av(nat))
allocate(q_m3_mag(nat),q_m2_mag(nat),q_m1_mag(nat),q_0_mag(nat),q_1_mag(nat),q_2_mag(nat),q_3_mag(nat))


!!!!!!!!!!!!!! READ CELL DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
open(unit=200, file=trim(cellname),status='unknown')

read(200,*)cell(1),cell(4),cell(2),cell(5),cell(6),cell(3)
close(200)
write(*,*) cell(1),cell(4),cell(2),cell(5),cell(6),cell(3)

if ( cell(1) .gt. cell(2) .and. cell(1) .gt. cell(3)) then !xprism
	write(*,*), 'X-Prism'         	
	Phase=1 
	PhaseName='Xprism'
	la=cell(2)/10
else if (cell(2) .gt. cell(1) .and. cell(2) .gt. cell(3)) then !Y Basal
	write(*,*), 'Y-Basal'
	Phase=2 
	PhaseName='Ybasal'
	la=cell(1)/10 !Units: nm
else 
   	write(*,*), 'Error Not Basal Or Prism!. Lx,Ly,Lz:', cell(1), cell(2), cell(3)
   	stop        
end if


!!!!!!!!!!!!!! NORMALISE REAL SPHERICAL HARMONICS !!!!!!!!!!!!!!!!!!!
y3=dsqrt(35./pi)/8.
y2=dsqrt(105./(2.*pi))/4.
y1=dsqrt(21./pi)/8.
y0=dsqrt(7./pi)/4.

NormM=4.*pi/7. ! 4pi/2l+1 l = 3 


densityIce=0.9167						! g/cm3
densityIce=densityIce*1E-21 			! g/cm3---> g/nm3
densityIce=densityIce/18.01528     		! *mass h20 g/mol--->mol/nm3
densityIce=densityIce*Na		 		! ---> nm-3
ThickNorm=2*densityIce*la*cell(3)/10	! cell(3)--> nm



!!!!!!!!!!!!!! OUTPUT FILES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
open(unit=19, file='q3Porfile.dat', status='unknown')		 
open(unit=22, file='ThicknessTime.dat', status='unknown')	


!!!!!!!!!!!!!! READ IN DCD COORDINATES FOR EACH FRAME !!!!!!!!!!!!!
do ii=1,nof
	Nl=0
	q_m3Av(:)= cmplx(0.0 ,0.0)
	q_m2Av(:)= cmplx(0.0 ,0.0)
	q_m1Av(:)= cmplx(0.0 ,0.0)
	q_0Av(:)= cmplx(0.0 ,0.0)
	q_1Av(:)= cmplx(0.0 ,0.0)
	q_2Av(:)= cmplx(0.0 ,0.0)
	q_3Av(:)= cmplx(0.0 ,0.0)
	q_m3(:)= cmplx(0.0 ,0.0)
	q_m2(:)= cmplx(0.0 ,0.0)
	q_m1(:)= cmplx(0.0 ,0.0)
	q_0(:)= cmplx(0.0 ,0.0)
	q_1(:)= cmplx(0.0 ,0.0)
	q_2(:)= cmplx(0.0 ,0.0)
	q_3(:)= cmplx(0.0 ,0.0)
	y1zxy=0
	neigh(:)=0
	q_3Tot(:)=0
	S_ij(:)=0
	
	if(dcd_has_cell) then
		read(999,iostat=iostat) 
	    if(iostat/=0) exit
   	endif


	read(999,iostat=iostat)(x(m),m=1,nat) 
    if(iostat/=0) exit
    read(999,iostat=iostat)(y(m),m=1,nat)
    if(iostat/=0) exit
    read(999,iostat=iostat)(z(m),m=1,nat)
    if(iostat/=0) exit

	pos(:,1)=x
	pos(:,2)=y
	pos(:,3)=z


	!!!!!!!!!!!!!! Nearest neighbours !!!!!!!!!!!!!
	neigh(:)=0

	do i=1,nat-1
		do j=i+1,nat
			dist_x=abs(pos(i,1)-pos(j,1))
			dist_y=abs(pos(i,2)-pos(j,2))
			dist_z=abs(pos(i,3)-pos(j,3))

			!!!!!!!!!!!!!! pbc conditions !!!!!!!!!!!!!
			if (dist_x .gt. 0.5*cell(1) .and. Phase .eq. 2) then ! if Ybasal 
				dist_x=cell(1)-dist_x
			end if	
			if (dist_y .gt. 0.5*cell(2) .and. Phase .eq. 1) then ! if Xprism
				dist_y=cell(2)-dist_y
			end if
			if (dist_z .gt. 0.5*cell(3)) then
				dist_z=cell(3)-dist_z
			endif         	 
				
			Sep_2=dist_x**2 + dist_y**2 + dist_z**2
			R=dsqrt(Sep_2)	 
 			
			if (Sep_2 .lt. rcut_sq) then 	! if neighbours
				neigh(i)=neigh(i)+1			!neigh has num o atoms rows
				neigh(j)=neigh(j)+1    
                list_pos(neigh(i)+ maxneigh*(i-1))=j  ! j is a neighbour. maxneigh makes sure there are enough rows		 
                list_pos(neigh(j)+ maxneigh*(j-1))=i  ! atom j has neighbour i. It can have up to 100 neigbours
        	end if 
  		end do
	end do


	!!!!!!!!!!!!!! Calculate spherical harmonics for all neighbours !!!!!!!!!!!!!
	do i=1,nat
			
		do k=1,neigh(i)
    
			j=list_pos(k+maxneigh*(i-1))
		
			dist_x=pos(i,1)-pos(j,1) 
			dist_y=pos(i,2)-pos(j,2)
			dist_z=pos(i,3)-pos(j,3) 

			!!!!!!!!!!!!!! pbc conditions !!!!!!!!!!!!!
			if (abs(dist_x) .gt. 0.5*cell(1) .and. Phase .eq. 2) then ! Ybasal 
				if (dist_x .gt. 0) then 
					dist_x=cell(1)-dist_x
				else 
					dist_x=cell(1)+dist_x
				end if
			end if	
			if (abs(dist_y) .gt. 0.5*cell(2) .and. Phase .eq. 1) then !if Xprism
				if (dist_y .gt. 0) then 
					dist_y=cell(2)-dist_y
				else
					dist_y=cell(2)+dist_y
				end if							
			end if
			if (abs(dist_z) .gt. 0.5*cell(3)) then
				if (dist_z .gt. 0) then 
					dist_z=cell(3)-dist_z
				else 
					dist_z=cell(3)+dist_z
				end if
			endif  
					

			Sep_2=dist_x**2 + dist_y**2 + dist_z**2
			R=dsqrt(Sep_2)

			dist_x=dist_x/R	
			dist_y=dist_y/R	
			dist_z=dist_z/R	
					
			r3=1.
			write(*,*), i, dist_x					
				

			!!!!!!!!!!!!!! Calculate Spherical Harmonics !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			!!!!!!!!!!!!!! Use real and imaginary components !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			!!!!!!!!!!!!!! Project vectors connecting molecules onto spherical harmonics !!!!!!!!

			q_m3(i)= (y3*(dist_x-imag*dist_y)**3)/r3 + q_m3(i)
			q_m2(i)= y2*((dist_x-imag*dist_y)**2)*dist_z/r3 + q_m2(i)
			y1zxy= (4*dist_z**2 - dist_x**2 -dist_y**2)/r3 			
			q_m1(i)= y1*(dist_x-imag*dist_y)*y1zxy + q_m1(i)			
			q_0(i)= y0*dist_z*(2*dist_z**2 -3*dist_x**2 -3*dist_y**2)/r3 + q_0(i)
			q_1(i)= (-y1)*(dist_x+imag*dist_y)*y1zxy + q_1(i)
			q_2(i)= y2*((dist_x+imag*dist_y)**2)*dist_z/r3 + q_2(i)
			q_3(i)= ((-y3)*(dist_x +imag*dist_y)**3)/r3 + q_3(i)	


		end do ! neighbours of atom i
				
					
		q_m3(i)= q_m3(i)/neigh(i)
		q_m2(i)= q_m2(i)/neigh(i)
		q_m1(i)= q_m1(i)/neigh(i)
		q_0(i)= q_0(i)/neigh(i)
		q_1(i)= q_1(i)/neigh(i)
		q_2(i)= q_2(i)/neigh(i)
		q_3(i)=q_3(i)/neigh(i)

	end do ! atom i

	q_m3Av(:)=q_m3(:)
	q_m2Av(:)=q_m2(:)
	q_m1Av(:)=q_m1(:)
	q_0Av(:)=q_0(:)
	q_1Av(:)=q_1(:)
	q_2Av(:)=q_2(:)
	q_3Av(:)=q_3(:)


	!!!!!!!!!!!!!! Now have q3m values for m=-3,3 for every O atom !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	count=0
	do i=1,nat

		do k=1,neigh(i)
		
			j=list_pos(k+maxneigh*(i-1))
					
			q_m3Av(i)=q_m3Av(i)+q_m3(j)
			q_m2Av(i)=q_m2Av(i)+q_m2(j)
			q_m1Av(i)=q_m1Av(i)+q_m1(j)
			q_0Av(i)=q_0Av(i)+q_0(j)
			q_1Av(i)=q_1Av(i)+q_1(j)
			q_2Av(i)=q_2Av(i)+q_2(j)
			q_3Av(i)=q_3Av(i)+q_3(j)	
										
			S_i(i)=q_m3(i)*conjg(q_m3(j)) +q_m2(i)*conjg(q_m2(j))&
				& + q_m1(i)*conjg(q_m1(j)) +q_0(i)*conjg(q_0(j)) +q_1(i)*conjg(q_1(j)) &
				& + q_2(i)*conjg(q_2(j)) +q_3(i)*conjg(q_3(j))
					
			S_j(j)=q_m3(j)*conjg(q_m3(i)) +q_m2(j)*conjg(q_m2(i))&
				& + q_m1(j)*conjg(q_m1(i)) +q_0(j)*conjg(q_0(i)) +q_1(j)*conjg(q_1(i)) &
				& + q_2(j)*conjg(q_2(i)) +q_3(j)*conjg(q_3(i))

			S_ij(i)=S_i(i)/sqrt(S_i(i)*S_j(i)) + S_ij(i)

		end do
		
		q_m3Av(i)=q_m3Av(i)/(neigh(i)+1)
		q_m2Av(i)=q_m2Av(i)/(neigh(i)+1)
		q_m1Av(i)=q_m1Av(i)/(neigh(i)+1)
		q_0Av(i)=q_0Av(i)/(neigh(i)+1)
		q_1Av(i)=q_1Av(i)/(neigh(i)+1)
		q_2Av(i)=q_2Av(i)/(neigh(i)+1)
		q_3Av(i)=q_3Av(i)/(neigh(i)+1)

	
	end do

	write(*,*), 'count', count


	!!!!!!!!!!!!!! Calculate <q3m> and <q3> !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!! Calculate QLL Thickness from q3 threshold !!!!!!!!!!!!!!!!!!!!!!!
	do i=1,nat
			
			
		q_3Tot(i)= dsqrt((abs(q_m3(i))**2 + abs(q_m2(i))**2 + abs(q_m1(i))**2 + abs(q_0(i))**2 &
			& + abs(q_1(i))**2 + abs(q_2(i))**2 + abs(q_3(i))**2)*NormM)
			
		q_3TotAv(i)= dsqrt((abs(q_m3Av(i))**2 + abs(q_m2Av(i))**2 + abs(q_m1Av(i))**2 + abs(q_0Av(i))**2 &
			& + abs(q_1Av(i))**2 + abs(q_2Av(i))**2 + abs(q_3Av(i))**2)*NormM)			

			
		if (q_3Tot(i) .lt. 0.7) then 
			Nl=Nl+1
		end if

		write(19,*)ii,pos(i,1),pos(i,2),pos(i,3),q_3Tot(i) ,q_3TotAv(i)	!,q_m2(i),q_m1(i),q_0(i)!,q_1(i),q_2(i),q_3(i) !S_ij(i) !!ii is timestep. File with q3 and coords for each timestep
	end do

		
	delta_thick=(Nl/ThickNorm)*10 !QLL thickness at timestep ii in Angstrom 
	write(22,'(i8,1f10.4)')ii,delta_thick   ! timestep QLL thickness in Angstrom
		
end do !frames

end program q3code






