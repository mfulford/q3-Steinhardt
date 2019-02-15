PROGRAM q6q4nc

implicit none

character*100 :: buffer, command, fname, boxname, cellname
character*3, allocatable :: sym(:)
character*8 :: cdum
character*10  :: PhaseName, temp_string, cutoffname
integer :: cut_int
integer :: nargc, iargc, dump, cart, nat, nof, i, j, ns, idum, ii, k, m, m1, m2, m3, ncr, siz, iat, jat
integer :: count, ncr_tmp, ncrit, crit, l, al, be, lwork, info, time_ncrit, ntap, iostat,  Phase, n_liquid
integer, parameter :: mxvic=100
integer, allocatable    :: llist(:), nvic(:)
integer, allocatable :: graph_solid_connect(:,:), followgraph(:), cr_list(:)
double precision :: ts, xdf, ydf, zdf, rsqdf, xdfk,ydfk,zdfk, rdf, u, connect, dotik, tframe, dens, radg, esse, q4ave, q6ave
double precision :: dotpave
double precision :: q3_0n, maxX, minX, maxY,minY, maxZ, minZ
double precision, parameter :: thresh_bond=0.5d0, pi=4.0d0*datan(1.0d0), sigma=3.405d0, Na=6.02214129*10**23d0
double precision, allocatable :: pos(:,:), q6(:), q4(:),q3(:), q3AV(:), q4AV(:), q6AV(:), dotp(:), icell(:)
double precision, allocatable :: rad_sph(:),rad_cub(:),q3AV_tot(:)
double precision, allocatable :: coord(:,:), alpha(:), normq4(:), normq6(:),normq3(:),posSlice(:), q4AV_tot(:), q6AV_tot(:)
double precision, allocatable :: nbond(:), rcm(:,:), tin(:,:), mtemp(:,:), time_rcm(:,:)
double precision, allocatable :: eigen(:), work(:), asp(:,:), pol1(:), pol2(:), S_ij(:), S_i(:), S_j(:), S_ij_norm(:), Hist(:,:)
double precision, external :: wigner3j
double precision :: rcut, rcutsq, norm_aveq4, norm_aveq6, tmp, tmp2,ThickNorm, delta_thick, lz, boxslice, La, Lb, thickness
complex*16, allocatable :: q6vec(:,:), q4vec(:,:),q3vec(:,:)
complex*16, allocatable :: q3_m3i(:),q3_m2i(:),q3_m1i(:),q3_0i(:),q3_1i(:),q3_2i(:),q3_3i(:)
complex*16, allocatable :: q4_m4i(:), q4_m3i(:),q4_m2i(:),q4_m1i(:),q4_0i(:),q4_1i(:),q4_2i(:),q4_3i(:), q4_4i(:)
complex*16, allocatable :: q6_m6i(:), q6_m5i(:), q6_m4i(:), q6_m3i(:),q6_m2i(:),q6_m1i(:),q6_0i(:),q6_1i(:),q6_2i(:)
complex*16, allocatable :: q6_3i(:), q6_4i(:), q6_5i(:), q6_6i(:)
complex*16, allocatable :: q3_m3Av(:),q3_m2Av(:),q3_m1Av(:),q3_0Av(:),q3_1Av(:),q3_2Av(:),q3_3Av(:)
complex*16, allocatable :: q4_m4Av(:),q4_m3Av(:),q4_m2Av(:),q4_m1Av(:),q4_0Av(:),q4_1Av(:),q4_2Av(:),q4_3Av(:),q4_4Av(:)
complex*16, allocatable :: q6_m6Av(:), q6_m5Av(:), q6_m4Av(:),q6_m3Av(:),q6_m2Av(:),q6_m1Av(:),q6_0Av(:),q6_1Av(:)
complex*16, allocatable :: q6_2Av(:),q6_3Av(:),q6_4Av(:),q6_5Av(:),q6_6Av(:)

complex*16              :: ave_q6vec(-6:6),ave_q4vec(-4:4),ave_q3vec(-3:3)
real  :: NC, NC2, NC_sum, densityIce, p, Nl_sum, Nl2, Nl1, rmin, rmax
double precision, allocatable :: delta_thickRange(:), delta_sum(:),q3_neigh(:), coordAv(:), ErrorHist(:,:,:), Error(:,:)
! DCD stuff
character*4              :: car4
character*80             :: car(10)
integer                  :: nstart, nsanc, nset, ntitle
integer                  :: charm, namin, i5(5), i9(9)
real*4                   :: DD
real*4, allocatable      :: x(:), y(:), z(:)
logical                  :: dcd_has_cell
real*8                   :: cell(6)
integer                  :: count1, count2, counta, countb, countc, count3, count4, count5, count6,r
integer                  :: n_avhist, phase_id, temp_name, time_count, t, time_count1, time_count2, time_count3


nargc=iargc()

if (nargc.lt.3) then
   write(*,*) "Usage:"
   write(*,*) "analysis.x [some .dcd traj] [cellname] [cutoff]"
   stop
else
   call getarg(1,buffer)
   read(buffer,*) fname
   call getarg(2,buffer)
   read(buffer,*) cellname
   call getarg(3,buffer)
   read(buffer,*) rcut

endif

write(*,*), "What Phase?"
write(*,*), "Type 1 for ice"
write(*,*), "Type 2 for water"
read(*,*), phase_id
if (phase_id .eq. 1) then
  PhaseName="ice"
else if (phase_id .eq. 2) then
  PhaseName="water"
end if



write (cutoffname,'(I2.0)') int(rcut*10.0)


write(*,*) 'CUT-OFF Name, value:', cutoffname, rcut

rcutsq=rcut**2.0d0


!!! Reading DCD file
  open (unit=11,file=trim(fname),form="unformatted")
  ! Header
  read(11) car4, nset, nstart, nsanc, i5, namin, DD, i9, charm
  read(11) ntitle, (car(i),i=1,ntitle)
  read(11) nat

  nof = nset


  if(i9(1)==1) then
     dcd_has_cell = .true.
  else
     write(*,'(a57,i8)') "A dcd with pbc cell is required"
     stop
  endif

  ! Three dimensional system

write(*,*) nset
cart=3
lwork=cart*cart-1


allocate(work(lwork))
allocate(pos(nat,cart),sym(nat),nvic(nat),llist(mxvic*nat),tin(cart,cart),mtemp(cart,cart))
allocate(dotp(nat),icell(cart*cart),coord(cart+1,mxvic))
allocate(alpha(mxvic*nat),q4vec(-4:4,nat),q6vec(-6:6,nat),q3vec(-3:3,nat),nbond(nat),normq6(nat),normq4(nat),normq3(nat))
allocate(cr_list(nat),eigen(cart),rcm(nat,cart),pol1(cart),pol2(cart))
allocate(x(nat), y(nat), z(nat))
allocate(S_ij(nat),S_i(nat),S_j(nat),S_ij_norm(nat),posSlice(nat))
allocate(coordAv(nat))

allocate(q6(nat),q4(nat),q3(nat),q4AV(nat),q3AV(nat), q6AV(nat))
allocate(q3_m3i(nat),q3_m2i(nat),q3_m1i(nat),q3_0i(nat),q3_1i(nat),q3_2i(nat),q3_3i(nat))
allocate(q4_m4i(nat), q4_m3i(nat),q4_m2i(nat),q4_m1i(nat),q4_0i(nat),q4_1i(nat),q4_2i(nat),q4_3i(nat), q4_4i(nat))
allocate(q6_m6i(nat), q6_m5i(nat), q6_m4i(nat), q6_m3i(nat),q6_m2i(nat),q6_m1i(nat),q6_0i(nat),q6_1i(nat))
allocate(q6_2i(nat),q6_3i(nat), q6_4i(nat), q6_5i(nat), q6_6i(nat))

allocate(q3AV_tot(nat), q4AV_tot(nat), q6AV_tot(nat))

allocate(q3_m3Av(nat), q3_m2Av(nat), q3_m1Av(nat), q3_0Av(nat), q3_1Av(nat), q3_2Av(nat),q3_3Av(nat))
allocate(q4_m4Av(nat), q4_m3Av(nat), q4_m2Av(nat), q4_m1Av(nat), q4_0Av(nat), q4_1Av(nat), q4_2Av(nat),q4_3Av(nat), q4_4Av(nat))
allocate(q6_m6Av(nat), q6_m5Av(nat), q6_m4Av(nat), q6_m3Av(nat), q6_m2Av(nat), q6_m1Av(nat),q6_0Av(nat), q6_1Av(nat), &
         & q6_2Av(nat), q6_3Av(nat), q6_4Av(nat), q6_5Av(nat), q6_6Av(nat))




open(unit=341, file='q3or_'//trim(PhaseName)//'_cut'//trim(cutoffname)//'.dat',status='unknown')
open(unit=342,file='q3m_'//trim(PhaseName)//'_cut'//trim(cutoffname)//'.dat',status='unknown')
open(unit=343,file='q3tot_'//trim(PhaseName)//'_cut'//trim(cutoffname)//'.dat',status='unknown')


open(unit=441, file='q4or_'//trim(PhaseName)//'_cut'//trim(cutoffname)//'.dat',status='unknown')
open(unit=442,file='q4m_'//trim(PhaseName)//'_cut'//trim(cutoffname)//'.dat',status='unknown')
open(unit=443,file='q4tot_'//trim(PhaseName)//'_cut'//trim(cutoffname)//'.dat',status='unknown')

open(unit=641, file='q6or_'//trim(PhaseName)//'_cut'//trim(cutoffname)//'.dat',status='unknown')
open(unit=642,file='q6m_'//trim(PhaseName)//'_cut'//trim(cutoffname)//'.dat',status='unknown')
open(unit=643,file='q6tot_'//trim(PhaseName)//'_cut'//trim(cutoffname)//'.dat',status='unknown')


open (unit=200,file=trim(cellname),status='unknown')

read(200,*)cell(1), cell(2), cell(3)
close(200)

 cell(4)=0.0
 cell(5)=0.0
 cell(6)=0.0


if (cell(1) .gt. cell(2) .and. cell(1) .gt. cell(3)) then
La=cell(2)
Lb=cell(3)
else if (cell(2) .gt. cell(1) .and. cell(2) .gt. cell(3)) then
La=cell(1)
Lb=cell(3)
else if (cell(3) .gt. cell(1) .and. cell(3) .gt. cell(2)) then
La=cell(1)
Lb=cell(2)
end if

write(*,*) La, Lb




count1=0
count2=0
count3=0
count4=0
count5=0
count6=0


do ii=1,nset
    NC=0
    NC2=0
    Nl2=0
    Nl1=0
    time_count1=0
    time_count2=0
    time_count3=0
	! Initialize averages

	! read atoms of frame ii
    if(dcd_has_cell) then
      read(11,iostat=iostat)
    if(iostat/=0) exit
    endif

    read(11,iostat=iostat)(x(m),m=1,nat)
    if(iostat/=0) exit
    read(11,iostat=iostat)(y(m),m=1,nat)
    if(iostat/=0) exit
    read(11,iostat=iostat)(z(m),m=1,nat)
    if(iostat/=0) exit



   pos(:,1)=x
   pos(:,2)=y
   pos(:,3)=z

   ! Compute Density
   dens=dble(nat)/(cell(1)*cell(2)*cell(3))
   icell(:)=0.0d0
   icell(1)=cell(1)
   icell(5)=cell(2)
   icell(9)=cell(3)



      ! Neighbours
      nvic(:)=0
      llist(:)=0
      do i=1,nat
         do j=i+1,nat
            xdf=pos(i,1)-pos(j,1)
            ydf=pos(i,2)-pos(j,2)
            zdf=pos(i,3)-pos(j,3)
            !!call images(cart,0,1,1,icell,xdf,ydf,zdf)
		!!!pbc conditions

			xdf=xdf -cell(1)*nint(xdf/cell(1))
			ydf=ydf -cell(2)*nint(ydf/cell(2))
			zdf=zdf -cell(3)*nint(zdf/cell(3))



            rsqdf = xdf**2+ydf**2+zdf**2
            if (rsqdf.lt.rcutsq) then
                nvic(i) = nvic(i) + 1
                llist(nvic(i)+mxvic*(i-1)) = j
                nvic(j) = nvic(j) + 1
                llist(nvic(j)+mxvic*(j-1)) = i
            endif
         enddo
      enddo

      ! q6&q4
      do i=1,nat
         do k=1,nvic(i)
            j=llist(k+mxvic*(i-1))
            xdfk=pos(i,1)-pos(j,1)
            ydfk=pos(i,2)-pos(j,2)
            zdfk=pos(i,3)-pos(j,3)
            !!call images(cart,0,1,1,icell,xdfk,ydfk,zdfk)

			!!!pbc conditions
			xdfk=xdfk -cell(1)*nint(xdfk/cell(1))
			ydfk=ydfk -cell(2)*nint(ydfk/cell(2))
			zdfk=zdfk -cell(3)*nint(zdfk/cell(3))

	        rdf=dsqrt(xdfk**2.0d0+ydfk**2.0d0 +zdfk**2.0d0)
            coord(1,k)=xdfk/rdf
            coord(2,k)=ydfk/rdf
            coord(3,k)=zdfk/rdf
            coord(4,k)=rdf
            u=rdf/rcut
            alpha(k+mxvic*(i-1)) = (rdf-rcut)**2
         enddo

         call q3at(i,nvic(i),coord,mxvic,rcut,nat,q3vec(:,i),normq3(i),rcut,nbond(i), q3_m3i(i),q3_m2i(i),q3_m1i(i)&
	 	         &,q3_0i(i), q3_1i(i), q3_2i(i), q3_3i(i))

         call q4at(i,nvic(i),coord,mxvic,rcut,nat,q4vec(:,i),normq4(i),rcut,nbond(i),q4_m4i(i),q4_m3i(i),q4_m2i(i),q4_m1i(i)&
	 	         &,q4_0i(i),q4_1i(i),q4_2i(i),q4_3i(i),q4_4i(i))

         call q6at(i,nvic(i),coord,mxvic,rcut,nat,q6vec(:,i),normq6(i),rcut,nbond(i),q6_m6i(i), q6_m5i(i), q6_m4i(i), q6_m3i(i), &
         & q6_m2i(i), q6_m1i(i), q6_0i(i), q6_1i(i), q6_2i(i), q6_3i(i), q6_4i(i), q6_5i(i), q6_6i(i))


         q3(i) = dsqrt(4.0d0*pi/7.0d0)*normq3(i)
         q4(i) = dsqrt(4.0d0*pi/9.0d0)*normq4(i)
         q6(i) = dsqrt(4.0d0*pi/13.0d0)*normq6(i)

        enddo

    q3_m3Av(:)=q3_m3i(:)
	q3_m2Av(:)=q3_m2i(:)
	q3_m1Av(:)=q3_m1i(:)
	q3_0Av(:)=q3_0i(:)
	q3_1Av(:)=q3_1i(:)
	q3_2Av(:)=q3_2i(:)
	q3_3Av(:)=q3_3i(:)

	q4_m4Av(:)=q4_m4i(:)
	q4_m3Av(:)=q4_m3i(:)
	q4_m2Av(:)=q4_m2i(:)
	q4_m1Av(:)=q4_m1i(:)
	q4_0Av(:)=q4_0i(:)
	q4_1Av(:)=q4_1i(:)
	q4_2Av(:)=q4_2i(:)
	q4_3Av(:)=q4_3i(:)
	q4_4Av(:)=q4_4i(:)

    q6_m6Av(:)=q6_m6i(:)
    q6_m5Av(:)=q6_m5i(:)
    q6_m4Av(:)=q6_m4i(:)
	q6_m3Av(:)=q6_m3i(:)
	q6_m2Av(:)=q6_m2i(:)
	q6_m1Av(:)=q6_m1i(:)
	q6_0Av(:)=q6_0i(:)
	q6_1Av(:)=q6_1i(:)
	q6_2Av(:)=q6_2i(:)
	q6_3Av(:)=q6_3i(:)
	q6_4Av(:)=q6_4i(:)
    q6_5Av(:)=q6_5i(:)
    q6_6Av(:)=q6_6i(:)


	do i=1,nat

        q3AV_tot(i)=q3(i)
		q4AV_tot(i)=q4(i)
        q6AV_tot(i)=q6(i)
		coordAv(i)=nvic(i)

		do k=1,nvic(i)
			j=llist(k+mxvic*(i-1))


			coordAv(i)=coordAv(i)+nvic(j)

            q3_m3Av(i)=q3_m3Av(i)+q3_m3i(j)
			q3_m2Av(i)=q3_m2Av(i)+q3_m2i(j)
			q3_m1Av(i)=q3_m1Av(i)+q3_m1i(j)
			q3_0Av(i)=q3_0Av(i)+q3_0i(j)
			q3_1Av(i)=q3_1Av(i)+q3_1i(j)
			q3_2Av(i)=q3_2Av(i)+q3_2i(j)
			q3_3Av(i)=q3_3Av(i)+q3_3i(j)


			q4_m4Av(i)=q4_m4Av(i)+q4_m4i(j)
			q4_m3Av(i)=q4_m3Av(i)+q4_m3i(j)
			q4_m2Av(i)=q4_m2Av(i)+q4_m2i(j)
			q4_m1Av(i)=q4_m1Av(i)+q4_m1i(j)
			q4_0Av(i)=q4_0Av(i)+q4_0i(j)
			q4_1Av(i)=q4_1Av(i)+q4_1i(j)
			q4_2Av(i)=q4_2Av(i)+q4_2i(j)
			q4_3Av(i)=q4_3Av(i)+q4_3i(j)
			q4_4Av(i)=q4_4Av(i)+q4_4i(j)


            q6_m6Av(i)=q6_m6Av(i)+q6_m6i(j)
            q6_m5Av(i)=q6_m5Av(i)+q6_m5i(j)
            q6_m4Av(i)=q6_m4Av(i)+q6_m4i(j)
			q6_m3Av(i)=q6_m3Av(i)+q6_m3i(j)
			q6_m2Av(i)=q6_m2Av(i)+q6_m2i(j)
			q6_m1Av(i)=q6_m1Av(i)+q6_m1i(j)
			q6_0Av(i)=q6_0Av(i)+q6_0i(j)
			q6_1Av(i)=q6_1Av(i)+q6_1i(j)
			q6_2Av(i)=q6_2Av(i)+q6_2i(j)
			q6_3Av(i)=q6_3Av(i)+q6_3i(j)
			q6_4Av(i)=q6_4Av(i)+q6_4i(j)
            q6_5Av(i)=q6_5Av(i)+q6_5i(j)
            q6_6Av(i)=q6_6Av(i)+q6_6i(j)

            q3AV_tot(i)=q3AV_tot(i)+q3(j)
			q4AV_tot(i)=q4AV_tot(i)+q4(j)
            q6AV_tot(i)=q6AV_tot(i)+q6(j)
		end do

		coordAv(i)=coordAv(i)/(nvic(i)+1)

        q3_m3Av(i)=q3_m3Av(i)/(nvic(i)+1)
		q3_m2Av(i)=q3_m2Av(i)/(nvic(i)+1)
		q3_m1Av(i)=q3_m1Av(i)/(nvic(i)+1)
		q3_0Av(i)=q3_0Av(i)/(nvic(i)+1)
		q3_1Av(i)=q3_1Av(i)/(nvic(i)+1)
		q3_2Av(i)=q3_2Av(i)/(nvic(i)+1)
		q3_3Av(i)=q3_3Av(i)/(nvic(i)+1)

		q4_m4Av(i)=q4_m4Av(i)/(nvic(i)+1)
		q4_m3Av(i)=q4_m3Av(i)/(nvic(i)+1)
		q4_m2Av(i)=q4_m2Av(i)/(nvic(i)+1)
		q4_m1Av(i)=q4_m1Av(i)/(nvic(i)+1)
		q4_0Av(i)=q4_0Av(i)/(nvic(i)+1)
		q4_1Av(i)=q4_1Av(i)/(nvic(i)+1)
		q4_2Av(i)=q4_2Av(i)/(nvic(i)+1)
		q4_3Av(i)=q4_3Av(i)/(nvic(i)+1)
		q4_4Av(i)=q4_4Av(i)/(nvic(i)+1)



        q6_m6Av(i)=q6_m6Av(i)/(nvic(i)+1)
        q6_m5Av(i)=q6_m5Av(i)/(nvic(i)+1)
        q6_m4Av(i)=q6_m4Av(i)/(nvic(i)+1)
		q6_m3Av(i)=q6_m3Av(i)/(nvic(i)+1)
		q6_m2Av(i)=q6_m2Av(i)/(nvic(i)+1)
		q6_m1Av(i)=q6_m1Av(i)/(nvic(i)+1)
		q6_0Av(i)=q6_0Av(i)/(nvic(i)+1)
		q6_1Av(i)=q6_1Av(i)/(nvic(i)+1)
		q6_2Av(i)=q6_2Av(i)/(nvic(i)+1)
		q6_3Av(i)=q6_3Av(i)/(nvic(i)+1)
		q6_4Av(i)=q6_4Av(i)/(nvic(i)+1)
        q6_5Av(i)=q6_5Av(i)/(nvic(i)+1)
        q6_6Av(i)=q6_6Av(i)/(nvic(i)+1)

        q3AV(i)= dsqrt((abs(q3_m3Av(i))**2 + abs(q3_m2Av(i))**2 + abs(q3_m1Av(i))**2 + abs(q3_0Av(i))**2 &
		& + abs(q3_1Av(i))**2 + abs(q3_2Av(i))**2 + abs(q3_3Av(i))**2))*dsqrt(4.0d0*pi/7.0d0)

		q4AV(i)= dsqrt((abs(q4_m4Av(i))**2 + abs(q4_m3Av(i))**2 + abs(q4_m2Av(i))**2 + abs(q4_m1Av(i))**2 + &
                      & abs(q4_0Av(i))**2 + abs(q4_1Av(i))**2 + abs(q4_2Av(i))**2 + abs(q4_3Av(i))**2 + &
                      & abs(q4_4Av(i))**2))*dsqrt(4.0d0*pi/9.0d0)


        q6AV(i)= dsqrt((abs(q6_m6Av(i))**2 + abs(q6_m5Av(i))**2 + abs(q6_m4Av(i))**2 + abs(q6_m3Av(i))**2 + &
                      & abs(q6_m2Av(i))**2 + abs(q6_m1Av(i))**2 + abs(q6_0Av(i))**2  + abs(q6_1Av(i))**2 + &
                      & abs(q6_2Av(i))**2 + abs(q6_3Av(i))**2 + abs(q6_4Av(i))**2 + abs(q6_5Av(i))**2 + &
                      & abs(q6_6Av(i))**2))*dsqrt(4.0d0*pi/13.0d0)


        q3AV_tot(i)=q3AV_tot(i)/(nvic(i)+1)
        q4AV_tot(i)=q4AV_tot(i)/(nvic(i)+1)
        q6AV_tot(i)=q6AV_tot(i)/(nvic(i)+1)

        write(341,'(f10.4)') q3(i)
        write(342,'(f10.4)') q3Av(i)
        write(343,'(f10.4)') q3AV_tot(i)

		write(441,'(f10.4)') q4(i)
        write(442,'(f10.4)') q4Av(i)
        write(443,'(f10.4)') q4AV_tot(i)

        write(641,'(f10.4)') q6(i)
        write(642,'(f10.4)') q6Av(i)
        write(643,'(f10.4)') q6AV_tot(i)

end do !nat

end do !frames

end program q6q4nc



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine images (imcon,idnode,mxnode,natm,cell,xxx,yyy,zzz)

implicit real*8 (a-h,o-z)

dimension xxx(*),yyy(*),zzz(*)
dimension cell(9),rcell(9)

   call invert(cell,rcell,det)
   if(abs(det).lt.1.d-6) stop "zero determinant cell matrix"

   do i=idnode+1,natm,mxnode

      ssx=(rcell(1)*xxx(i)+rcell(4)*yyy(i)+rcell(7)*zzz(i))
      ssy=(rcell(2)*xxx(i)+rcell(5)*yyy(i)+rcell(8)*zzz(i))
      ssz=(rcell(3)*xxx(i)+rcell(6)*yyy(i)+rcell(9)*zzz(i))

      xss=ssx-nint(ssx)
      yss=ssy-nint(ssy)
      zss=ssz-nint(ssz)

      xxx(i)=(cell(1)*xss+cell(4)*yss+cell(7)*zss)
      yyy(i)=(cell(2)*xss+cell(5)*yss+cell(8)*zss)
      zzz(i)=(cell(3)*xss+cell(6)*yss+cell(9)*zss)

   end do

return
end subroutine images

!!!!!!!

subroutine invert(a,b,d)

real*8 a,b,d,r

dimension a(9),b(9)

b(1)=a(5)*a(9)-a(6)*a(8)
b(2)=a(3)*a(8)-a(2)*a(9)
b(3)=a(2)*a(6)-a(3)*a(5)
b(4)=a(6)*a(7)-a(4)*a(9)
b(5)=a(1)*a(9)-a(3)*a(7)
b(6)=a(3)*a(4)-a(1)*a(6)
b(7)=a(4)*a(8)-a(5)*a(7)
b(8)=a(2)*a(7)-a(1)*a(8)
b(9)=a(1)*a(5)-a(2)*a(4)

d=a(1)*b(1)+a(4)*b(2)+a(7)*b(3)
r=0.d0
if(abs(d).gt.0.d0)r=1.d0/d

b(1)=r*b(1)
b(2)=r*b(2)
b(3)=r*b(3)
b(4)=r*b(4)
b(5)=r*b(5)
b(6)=r*b(6)
b(7)=r*b(7)
b(8)=r*b(8)
b(9)=r*b(9)
return
end subroutine invert

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine q6at(iat,nvic,coord,mxv,r_c,nox,q6,normq6,rsoft,nbond, q_m6i, q_m5i, q_m4i,q_m3i,q_m2i,q_m1i,q_0i, &
                & q_1i,q_2i,q_3i,q_4i, q_5i, q_6i)
  implicit none

  integer, parameter                       :: dbl=kind(1.d0)
  real(dbl)                                :: normq6, nbond
  complex(dbl)                             :: q6(-6:6), q_m6i, q_m5i, q_m4i,q_m3i,q_m2i,q_m1i,q_0i,q_1i,q_2i,q_3i,q_4i, q_5i, q_6i
!*************************************************************************
  integer                                  :: m,i,j,k,ii,ix, mxv
  integer, intent(in)                      :: iat, nox, nvic
  real(dbl)                                :: rnorm(4), num
  real(dbl), dimension(4,*)                :: coord
  logical, save                            :: first=.true.
!*************************************************************************
  real(dbl)                                :: coeff_poly(0:6)
  real(dbl)                                :: u,rr
  real(dbl), save                          :: normaliz(0:6)
  real(dbl)                                :: pol_ass, r_c, frr
  real(dbl)                                :: alpha
  real(dbl)                                :: fact1(0:6),fact2(0:6)
  complex(dbl)                             :: com1, comm1(0:6)
  complex(dbl)                             :: yy(-6:6)
  real(dbl), external                      :: deriv_poly_6
  real(dbl)                                :: nnqk
  real(dbl)                                :: rsoft
  real(dbl), parameter :: sigma = 3.405 !LJ Sigma of Ar
  real(dbl), parameter :: pi = 3.141592653589793d0
!   data r_c /4.0/  !  LJ fluid
!   data r_c /2.6/  !  water

  data fact1 /720.0, 120.0, 24.0, 6.0, 2.0, 1.0, 1.0/
  data fact2 /720.,5040.,40320.,362880.,3628800.,39916800.,479001600./
  data coeff_poly /-0.3125,0.0,6.5625,0.0,-19.6875,0.0,14.4375/
!*************************************************************************
  if (first) then
    do m = 0, 6
      normaliz(m)=sqrt((13.0*fact1(m))/(4.0*pi*fact2(m)))*(-1.0d0)**m
    end do
  endif
  first = .false.

  nbond = 0.d0
  q6(:)  = cmplx(0.0 ,0.0)
  k          = iat
  do ii=1, nvic
     rnorm(:) = coord(:,ii)
     rr = rnorm(4)
     !
     u = rr/rsoft
     nbond = nbond + 1
     com1=cmplx(rnorm(1),rnorm(2))

     do m=0,6
        pol_ass=deriv_poly_6(coeff_poly,m,rnorm(3))
        yy(m)=normaliz(m)*pol_ass*com1**m
        yy(-m)=(-1.0d0)**m*conjg(yy(m))
     end do
     do m= -6,6
        q6(m)=q6(m)+yy(m)
     enddo



  enddo       !  loop on neighbors
  q6 = q6/nbond

  q_m6i=q6(-6)
  q_m5i=q6(-5)
  q_m4i=q6(-4)
  q_m3i=q6(-3)
  q_m2i=q6(-2)
  q_m1i=q6(-1)
  q_0i=q6(0)
  q_1i=q6(1)
  q_2i=q6(2)
  q_3i=q6(3)
  q_4i=q6(4)
  q_5i=q6(5)
  q_6i=q6(6)


  nnqk = q6(0)*q6(0)
  do m = 1,6
     nnqk = nnqk + 2.d0*(real(q6(m))**2 + dimag(q6(m))**2)
  enddo
  normq6 = dsqrt(nnqk)
  return
end subroutine q6at

!!!!!!

function deriv_poly_6(coeff_poly,order,x) result(res)
  integer, parameter                   :: dbl=kind(1.d0)
  real(dbl), dimension(0:6), intent(in):: coeff_poly
  integer, intent(in)                  :: order
  real(dbl)                            :: x,xi
  real(dbl)                            :: res
  integer                              :: i,j,fact

  res = 0.0_dbl
  xi  = 1.0_dbl
  do i=order,ubound(coeff_poly,1)
     ! fact= i!/(i-order)!
     fact=1
     do j=(i-order+1),i
        fact=fact*j
     end do
     res=res+coeff_poly(i)*dble(fact)*xi
     xi=xi*x
  end do

  return

end function deriv_poly_6

recursive function fact(n) result(res)
   implicit none
   integer     :: n
   integer(8)  :: res

   if ( n >= 1) then
        res = n * fact(n-1)
   else
        res = 1
   end if

end function fact

subroutine q4at(iat,nvic,coord,mxv,r_c,nox,q4,normq4,rsoft,nbond,q_m4i,q_m3i,q_m2i,q_m1i,q_0i,q_1i,q_2i,q_3i,q_4i)
  implicit none

  integer, parameter                       :: dbl=kind(1.d0)
  real(dbl)                                :: normq4, nbond
  complex(dbl)                             :: q4(-4:4), q_m3i,q_m2i,q_m1i,q_0i,q_1i,q_2i,q_3i, q_m4i, q_4i
!*************************************************************************
  integer                                  :: m,i,j,k,ii,ix, mxv
  integer, intent(in)                      :: iat, nox, nvic
  real(dbl)                                :: rnorm(4), num
  real(dbl), dimension(4,*)                :: coord
  logical, save                            :: first=.true.
!*************************************************************************
  real(dbl)                                :: coeff_poly(0:4)
  real(dbl)                                :: u,rr
  real(dbl), save                          :: normaliz(0:4)
  real(dbl)                                :: pol_ass, r_c
  real(dbl)                                :: alpha
  real(dbl)                                :: fact1(0:4),fact2(0:4)
  complex(dbl)                             :: com1, comm1(0:4)
  complex(dbl)                             :: yy(-4:4)
  real(dbl), external                      :: deriv_poly_4
  real(dbl)                                :: nnqk
  real(dbl)                                :: rsoft
  real(dbl), parameter :: sigma = 3.405 !LJ Sigma of Ar
  real(dbl), parameter :: pi = 3.141592653589793d0
  data fact1 /24.0, 6.0, 2.0, 1.0, 1.0/
  data fact2 /24.0, 120.0, 720.0, 5040.0, 40320.0/
  data coeff_poly /0.375, 0.0, -3.75, 0.0, 4.375/
!*************************************************************************
  if (first) then
    do m = 0, 4
      normaliz(m)=sqrt((9.0*fact1(m))/(4.0*pi*fact2(m)))*(-1.0d0)**m
    end do
  endif
  first = .false.

  nbond = 0.d0
  q4(:)  = cmplx(0.0 ,0.0)
  k          = iat
  do ii=1, nvic
     rnorm(:) = coord(:,ii)
     nbond = nbond + 1
     com1=cmplx(rnorm(1),rnorm(2))
     do m=0,4
        pol_ass=deriv_poly_4(coeff_poly,m,rnorm(3))
        yy(m)=normaliz(m)*pol_ass*com1**m
        yy(-m)=(-1.0d0)**m*conjg(yy(m))
     end do
     do m= -4,4
        q4(m)=q4(m)+yy(m)
     enddo
  enddo       !  loop on neighbors
  q4 = q4/nbond

  q_m4i=q4(-4)
  q_m3i=q4(-3)
  q_m2i=q4(-2)
  q_m1i=q4(-1)
  q_0i=q4(0)
  q_1i=q4(1)
  q_2i=q4(2)
  q_3i=q4(3)
  q_4i=q4(4)


  nnqk = q4(0)*q4(0)
  do m = 1,4
     nnqk = nnqk + 2.d0*(real(q4(m))**2 + dimag(q4(m))**2)
  enddo
  normq4 = dsqrt(nnqk)
  return
end subroutine q4at

!!!!!!

function deriv_poly_4(coeff_poly,order,x) result(res)
  integer, parameter                   :: dbl=kind(1.d0)
  real(dbl), dimension(0:4), intent(in):: coeff_poly
  integer, intent(in)                  :: order
  real(dbl)                            :: x,xi
  real(dbl)                            :: res
  integer                              :: i,j,fact

  res = 0.0_dbl
  xi  = 1.0_dbl
  do i=order,ubound(coeff_poly,1)
     ! fact= i!/(i-order)!
     fact=1
     do j=(i-order+1),i
        fact=fact*j
     end do
     res=res+coeff_poly(i)*dble(fact)*xi
     xi=xi*x
  end do

  return
end function deriv_poly_4


subroutine q3at(iat,nvic,coord,mxv,r_c,nox,q3,normq3,rsoft,nbond,q_m3i,q_m2i,q_m1i,q_0i,q_1i,q_2i,q_3i)
  implicit none

  integer, parameter                       :: dbl=kind(1.d0)
  real(dbl)                                :: normq3, nbond
  complex(dbl)                             :: q3(-3:3), q_m3i,q_m2i,q_m1i,q_0i,q_1i,q_2i,q_3i
!*************************************************************************
  integer                                  :: m,i,j,k,ii,ix, mxv
  integer, intent(in)                      :: iat, nox, nvic
  real(dbl)                                :: rnorm(4), num
  real(dbl), dimension(4,*)                :: coord
  logical, save                            :: first=.true.
!*************************************************************************
  real(dbl)                                :: coeff_poly(0:3)
  real(dbl)                                :: u,rr
  real(dbl), save                          :: normaliz(0:3)
  real(dbl)                                :: pol_ass, r_c
  real(dbl)                                :: alpha
  real(dbl)                                :: fact1(0:3),fact2(0:3)
  complex(dbl)                             :: com1, comm1(0:3)
  complex(dbl)                             :: yy(-3:3)
  real(dbl), external                      :: deriv_poly_3
  real(dbl)                                :: nnqk
  real(dbl)                                :: rsoft
  real(dbl), parameter :: sigma = 3.405 !LJ Sigma of Ar
  real(dbl), parameter :: pi = 3.141592653589793d0
  data fact1 /6.0, 2.0, 1.0, 1.0/
  data fact2 /6.0, 24.0, 120.0, 720.0/
  data coeff_poly /0.0, -1.5, 0.0, 2.5/
!*************************************************************************
  if (first) then
    do m = 0, 3
      normaliz(m)=sqrt((7.0*fact1(m))/(4.0*pi*fact2(m)))*(-1.0d0)**m
    end do
  endif
  first = .false.

  nbond = 0.d0
  q3(:)  = cmplx(0.0 ,0.0)
  k          = iat
  do ii=1, nvic
     rnorm(:) = coord(:,ii)
     nbond = nbond + 1
     com1=cmplx(rnorm(1),rnorm(2))
     do m=0,3
        pol_ass=deriv_poly_3(coeff_poly,m,rnorm(3))
        yy(m)=normaliz(m)*pol_ass*com1**m
        yy(-m)=(-1.0d0)**m*conjg(yy(m))
     end do
     do m= -3,3
        q3(m)=q3(m)+yy(m)
     enddo
  enddo       !  loop on neighbors
  q3 = q3/nbond

  q_m3i=q3(-3)
  q_m2i=q3(-2)
  q_m1i=q3(-1)
  q_0i=q3(0)
  q_1i=q3(1)
  q_2i=q3(2)
  q_3i =q3(3)

  nnqk = q3(0)*q3(0)
  do m = 1,3
     nnqk = nnqk + 2.d0*(real(q3(m))**2 + dimag(q3(m))**2)
  enddo
  normq3 = dsqrt(nnqk)
  return
end subroutine q3at

!!!!!!

function deriv_poly_3(coeff_poly,order,x) result(res)
  integer, parameter                   :: dbl=kind(1.d0)
  real(dbl), dimension(0:3), intent(in):: coeff_poly
  integer, intent(in)                  :: order
  real(dbl)                            :: x,xi
  real(dbl)                            :: res
  integer                              :: i,j,fact

  res = 0.0_dbl
  xi  = 1.0_dbl
  do i=order,ubound(coeff_poly,1)
     ! fact= i!/(i-order)!
     fact=1
     do j=(i-order+1),i
        fact=fact*j
     end do
     res=res+coeff_poly(i)*dble(fact)*xi
     xi=xi*x
  end do

  return
end function deriv_poly_3
