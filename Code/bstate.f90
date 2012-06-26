program bstate
implicit none
CHARACTER*100 :: BUFFER
INTEGER, parameter :: I8 = 8
INTEGER(kind = I8) :: Brows
INTEGER :: N,Q, start, add_col, i, k,p, TEST
INTEGER, allocatable, dimension(:,:) ::work, B, Bnew
INTEGER :: temp

!READ IN ARGS
call getarg(1,BUFFER)
read(BUFFER,*) N

call getarg(2,BUFFER)
read(BUFFER,*) Q


Brows = nchoosek_mult(N+Q-1,N)
!write(*,*) Brows
if(Brows .lt. 0) then
   write(*,*) "Nchoosek overflow error"
   STOP
endif
allocate(B(Brows, Q))
allocate(Bnew(Brows,Q))



B = 0
Bnew = 0
do k= 0, N-1
   start = 1
   B(:,:) = Bnew(:,:)
   do i = 1, Q
      temp = nchoosek_mult(k+i-1,k)
     ! write(*,*) 'temp', temp
      add_col = Q+1-i
      allocate(work(temp,Q))
      work(:,:) = B(1:temp,:)
      work(:,add_col) = work(:,add_col) + 1
!      write(*,*) work
 !     write(*,*) "Now B"
  !    write(*,*) "START:",start,"END",start+temp-1
      Bnew(start:start+temp-1,:) = work(:,:)
      !write(*,*) Bnew
     ! Write(*,*)"NEW"
     ! write(*,*) B
      !write(*,*) Bnew(1:temp,:)
      start = start+temp
      deallocate(work)
   enddo
enddo
!write(*,*) "END LOOP"
!write(*,*) B(:,:)
OPEN(UNIT=4, FILE = 'out.dat')
do i = 1, Brows
   write(4,*) Bnew(i,:)
   !write(4,*)
enddo
close(4)
deallocate(B)
deallocate(Bnew)



 



STOP

CONTAINS

  
  FUNCTION Factorial(n) RESULT(Fact)
    
    Implicit None
    integer, parameter :: I8 = 8
    Integer(kind=I8) :: fact, count
    Integer, Intent(IN) :: n

    fact = 1
    count = n
    do while( count .gt. 1)
       fact = fact*count

       count = count-1
    ENDDO
    if(count .eq. 0) THEN
       fact = 1
    endif

  END FUNCTION Factorial

       
  FUNCTION nchoosek(n,k) 
    IMPLICIT NONE
    integer, intent(in):: n, k
    integer, parameter :: I8 = 8
    INTEGER(kind=I8):: nchoosek
    nchoosek = Factorial(n)
    !write(*,*) 'n ,k', n, k
   ! write(*,*) nchoosek
    nchoosek=nchoosek/Factorial(k)
    !write(*,*) nchoosek, Factorial(k)
    nchoosek=nchoosek/Factorial(n-k)
    !write(*,*) nchoosek
    return
  END FUNCTION nchoosek

  
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
 !    write(*,*) nchoosek_mult
  ENDDO


  nchoosek_mult = nchoosek_mult+0.5
  return
ENDFUNCTION nchoosek_mult
end program bstate
