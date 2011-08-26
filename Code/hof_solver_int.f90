program hof_solver
implicit none
!INTERFACE
  ! FUNCTION Factorial(n)
    ! INTEGER :: Factorial
   !  INTEGER, INTENT(IN) :: n
  ! END FUNCTION Factorial
 ! 
!END INTERFACE
!This program will deal with interactions now
!This program goals are:
!1. To read in the basis states generated from octave code and 
!    make sure they are the correct basis states.
!2. To use the Bose-Hubbard Interaction Hamiltonian equations
!3. to build up the H matrix
!3. To solve the H matrix 

!4. After that is complete, I will put loops in to allow the
!program to generate a dataset for the different k space
! configurations.

!Inputs Needed:

! 1. Number of Particles we have in our system

! 2. The Value of Q.

! The first step is to read in the basis states

!Define Variables
integer, parameter :: i8 = 8
integer(kind=i8) :: bignum1, bignum2, bignum3, occupation, ROWS, END_K
character(len=100) :: BUFFER
double precision, dimension(:,:), allocatable :: yy
integer :: k,TOO, FROM, t, N, counter, ENDTEST,i, total, Q, e, m, create, anhil
integer :: TESTBREAK = 0, test, ii(1), mm(1), start, bb, tt, pp, P, jjj
integer, dimension(:), allocatable ::protobasis, transtest, oldlist
integer, dimension(:,:), allocatable :: xx
double precision :: tjingy, hh(4), jj(4), nrml, G
COMPLEX*16, dimension(:,:), allocatable :: H, hopr
COMPLEX*16, dimension(:), allocatable :: hsclr
COMPLEX :: J =(0,1)
double precision :: prefac, kx, ky, c, d, f, PI=4.0*atan(1.0)
integer :: v(4), LWORK
!DEFINE MATRIX SOLVER VARIABLES
integer :: DUMMY(1,1), ok
COMPLEX*16 , dimension(:), allocatable :: eig, work
double precision, dimension(:), allocatable :: eig_zheev, zheev_rwork
!Get ARGS


call getarg(1, BUFFER)
read(BUFFER,*) N

call getarg(2,BUFFER)
read(BUFFER,*) Q

call getarg(3, BUFFER)
read(BUFFER,*) G

call getarg(4, BUFFER)
read(BUFFER,*) END_K
!read in basis state information
!need to know how much memory to allocate for this
!it will be a HUGE array-each row will be a basis state
!each row will be Q cells long
!and there will be (n+k-1)choose(k-1) many rows.

ROWS = nchoosek_mult(N+Q-1,N)
if(ROWS .lt. 0) THEN
   write(*,*) "nchoosek overflow error"
   stop
endif
allocate(yy(ROWS,Q)) 
write(*,*) 'Brows;', ROWS
!OPEN INPUT FILE
!NOTE: you NEED to make sure that the input file generated
!is the same as the N, Q you are feeding in. To be safe, I will
!create a shell script to run the executable from that should
!take care of it, but it is something to watch out for

!Try and test for it?

open(UNIT = 3, FILE = "inff1.dat",action='READ', status='old') 
READ(3,*) ((yy(i,k), k =1,Q), i = 1,ROWS)
close(3)
!now, should have a large array that stores my basis states.
write(*,*) "Begin xx"
allocate(xx(ROWS,Q))
write(*,*) "END xx"
do i =1,ROWS
   do k= 1,Q
      xx(i,k) = int(yy(i,k))
   ENDDO
ENDDO
!This part may be optional, but I am converting the double array
!to an integer array just to be safe.
!I think I can take this down though, as I am now generating 
!my basis states using a fortran program, and I think
!it output them as ints. So I should just be able to read them
!in as ints now.

!DEBUGGING
!do i = 1,ROWS
!   irint *, (xx(i,k), k=1,Q)
!ENDDO
deallocate(yy)



!allocate memory for matrix building and solving
write(*,*) 'begin h'
allocate(H(ROWS, ROWS))
write(*,*) "END H"
allocate(WORK(1))
allocate(transtest(Q))
allocate(eig_zheev(ROWS))
allocate(zheev_rwork(3*ROWS-2))
!TEST LOOPS
!do P = 0, Q !disabled P for now, working with LARGE N, Q
P = 1.0
   do tt = 0, END_K-1
      kx = 2.0*PI/dble(Q)/dble(END_K-1)*dble(tt)
      do bb = 0,END_K-1
         ky = 2.0*PI/dble(END_K-1)*dble(bb)
         write(*,*) ky
         !Now build H matrix
         !set everything to zero (otherwise, when ZHEEV is called
         !there is junk in WORK, and I find that that messes things
         !up
         !  write(*,*) "BEGIN SET TO ZERO"
         H(:,:) = 0
         ! write(*,*) "H"
         WORK = 0.0
         !write(*,*) 'three'
         eig_zheev = 0.0
         zheev_rwork = 0.0
         !write(*,*) "SET TO ZERO"
         do i = 1, ROWS
            do k = 1, ROWS
               prefac = 1.
               
               transtest(:) = 0
               transtest = xx(i,:) - xx(k,:)
               !so, there are only two cases that we are intersted in
               !when i = k, transtest =[zeroes] and we are on a diagonal. 
               !We know that this transition will only come from the cos and interactions
               ! so we can include that right off the bat.
               IF(i .eq. k) THEN
                  !immediatly we know cos and interactions is all we care about (xx(i.m) stores the prefactor (# of P in Q)) 
                  do m = 1,Q
                     if(xx(i,m) .ne. 0)THEN
										 !interactions have been added
                                         !making quick changes to test current
                                         !cal. NEED to change back!
                        H(i,k) = H(i,k) + xx(i,m)*2.0*cos(kx+2.0*PI*dble(P)*dble(m)/dble(Q))+g*dble(N)*dble(N-1)
                     ENDIF
                  ENDDO
                  
                  !Now, I need to start testing for the other terms.
                  !See notes for a detailed explanation of algorithm
               ENDIF
                  !****************TRANSITION TEST***************************
                  create = 0
                  anhil = 0
                  TESTBREAK = 0
                  do m = 1,Q      
                     if( (transtest(m) .eq. 1) .and. (create .eq. 0) ) THEN
                        create = m
                     ELSEIF( (transtest(m) .eq. -1) .and. (anhil .eq. 0) ) THEN
                        anhil = m
                     ELSEIF(transtest(m) .eq. 0) THEN
                        TESTBREAK = 0
                     ELSE
                        TESTBREAK = 1
                        exit
                     ENDIF
                  ENDDO
                  
                  IF(TESTBREAK .eq. 1) THEN
                     cycle ! cycle next vale in  do k
                  ENDIF
                  
                  !**********************END*TRANSTION*TEST******************
                  !If the TESTBREAK == 1, then some part of transtest is .ne. 0 and
                  !the two states, xx(i,:) and xx(k,:) are seperated by more than
                  !one energy move, so the transition is not allowed. I put in TESTBREAK
                  !so that if this occurs we can skip the rest of the loop.
                  !NOTE: The rest of the loop is NOT robust enough to decect this case
                  !by itself. If this part is taken out, then it will allow any case
                  !where either create .ne. 0 or anhil .ne. 0 or some
                  !other part of transtest .eq. 0 and will result in an
                  !incorrect Hamiltonian matrix.
                  
                  
                  !Start testing for allowability of hamiltonian terms
                  IF( (create - anhil .eq. 1) .or. (create -(anhil-Q) .eq. 1) ) THEN
                     prefac = 1
                     do m = 1,Q     
                        if(transtest(m) .eq. 1) THEN
                           prefac=prefac*sqrt( dble(xx(i,m)) )
                           !     write(*,*) 'first prefact', prefac
                        ELSEIF(transtest(m) .eq. -1) THEN
                           prefac = prefac*sqrt( dble(xx(k,m)) )
                           !      write(*,*) 'second prefac', prefac
                        ENDIF
                     eNDDO
                     
                     H(i,k) = H(i,k)+prefac*exp(J*ky)
                     
                     !we know that exp(-j) is all we care about
                  ENDIF
                  
                  IF( (anhil - create .eq. 1) .or. (anhil - (create-Q) .eq. 1) ) THEN
                     !we know that exp(J) is all we care about
                     prefac = 1
                     do m = 1,Q   
                        if(transtest(m) .eq. 1) THEN
                           prefac=prefac*sqrt( dble(xx(i,m)) )
                        ELSEIF(transtest(m) .eq. -1) THEN
                           prefac = prefac*sqrt( dble(xx(k,m)) )
                        ENDIF
                     ENDDO
                     
                     H(i,k) = H(i,k)+ prefac*exp(-J*ky)
                     
                  ENDIF
            ENDDO !end K
         
         
         ENDDO !end i
         
         !******************END PREFACTOR CALCULATION*********************************************
         
         !The matrix should now be built
         !   write(*,*) dble(P)/dble(Q)
         !  do i = 1,ROWS
         !     write(*,*) H(i,:)
         !  ENDDO
         !need to test and run this code!
         !deallocate(yy)
         !H is now built.
         !Now, I can try to solve it
         
         !Query ZHEEV to find best workspace
         LWORK = -1
         call ZHEEV('V', 'U', ROWS, H, ROWS, eig_zheev, WORK, LWORK, zheev_rwork, ok)
         LWORK = WORK(1)
         write(*,*) 'WORK'
         deallocate(WORK)
         allocate(WORK(LWORK))
         
         !now solve with ZHEEV
         call ZHEEV('V', 'U', ROWS, H, ROWS, eig_zheev, WORK, LWORK,zheev_rwork,ok)
         
         
         
         write(*,*) "STATUS OF ZHEEV", ok
         !WRITE results to file
         !Open file for output
         open(unit = 4, FILE = 'out_solv.dat')
         
         do i = 1,ROWS
            write(4,*) dble(P)/dble(Q), '  ', real( eig_zheev(i))
         enddo
         
		 open(unit = 12, FILE='energy_spec.dat')

		 do i = 1, ROWS
			Write(12,*) dble(P)/dble(Q),'@',kx,'@',ky,'@',eig_zheev(i),'@',(H(jjj,i),'@', jjj=1,ROWS) 
         write(*,*) "END END"
         enddo 


         !Now, write two files that contain the imaginary and real parts, for
         !ease of processing.

         open(unit = 14, FILE='real_eig.dat')
         open(unit= 15, FILE = 'imag_eig.dat')

         do i =1,ROWS
			write(14,*) dble(P)/dble(Q),kx,ky,eig_zheev(i),(real(H(jjj,i)), jjj=1,ROWS) 
			write(15,*) dble(P)/dble(Q),kx,ky,eig_zheev(i),(AIMAG(H(jjj,i)), jjj=1,ROWS) 
         enddo

      ENDDO! end bb(ky)
   enddo ! end kx
!enddo! end p (disabled for now)
STOP

CONTAINS
!--------Factorial-----------------------------
!
!  Function to calculate factorials
!basis_gen.sh
!----------------------------------------------

  FUNCTION Factorial(n) RESULT(Fact)
    
    Implicit None
    integer, parameter :: I8 = 8
    Integer(kind=I8) :: fact, count
    Integer, Intent(IN) :: n
    fact = 1
    count = n
    do while( count .gt. 1)
       fact = fact*count
       write(*,*) fact
       count = count-1
    ENDDO
  END FUNCTION Factorial

FUNCTION nchoosek_mult(n_in,k_in)
  IMPLICIT NONE
  integer, intent(in) :: n_in, k_in
  integer, parameter :: I8=8
  double precision :: nchoosek_mult
  integer :: i,n,k 
  nchoosek_mult = 1
  n = n_in
  k = k_in
  if (k > n) THEN
     nchoosek_mult = 0
     return
  ENDIF
  
  if(k > n/2) THEN
     k = n-k
  ENDIF

  do i = 1, k
     nchoosek_mult = nchoosek_mult *(n-k+i)/i
     if (nchoosek_mult .gt. huge(0)) then
        nchoosek_mult = -99
        !ERROR RETURN
        !nchoosek_mult will be cast as int()
        !and int() can only go as high as huge(0)
        return
     endif
  ENDDO

  nchoosek_mult = nchoosek_mult+0.5
  return
ENDFUNCTION nchoosek_mult
 ENDPROGRAM hof_solver
