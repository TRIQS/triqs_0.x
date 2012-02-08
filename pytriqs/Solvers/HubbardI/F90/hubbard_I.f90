MODULE hubbard_I_data
  REAL(KIND=8) :: B
  TYPE occup
     LOGICAL :: ifdiag, run
     INTEGER :: ndeg 
     INTEGER :: n
     !COMPLEX(KIND=8), DIMENSION(:,:), ALLOCATABLE :: Hn
     !REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: En 
     !INTEGER, DIMENSION(:), ALLOCATABLE  :: st_n 
     !INTEGER, DIMENSION(:), ALLOCATABLE  :: narr
     COMPLEX(KIND=8), pointer :: Hn(:,:)
     REAL(KIND=8), pointer :: En(:) 
     INTEGER, pointer  :: st_n(:) 
     INTEGER, pointer  :: narr(:)
  ENDTYPE occup
END MODULE hubbard_I_data



SUBROUTINE gf_HI_fullU(GF,Tail,e0f,ur,umn,ujmn,zmsb,nlm,Iwmax,nmom,ns,atocc,atmag,temp,verbosity)

!
! Computes atomic GF with the full 4-index U 8.10.2007
! Finite temperature version
! (by L.V. Pourovskii)
!
!  Last change 23.10.2007
!
! Included tail calculation for the triqs package
! M. Aichhorn 10-2009
!
! GF /output/ - atomic Green's function
! e0f /input/ - atomic level position e0f_mm' C*_m C_m'
! ur  /input/ - 4-index U
! umn /input/ - 2-index reduced Coulomb interation U(m,m',m,m')
! ujmn /input/ - 2-index Coulomb interation U(m,m',m,m')-U(m,m',m',m)
! zmsb(Iwmax)/input/ - complex energy mesh
! nlm, ns /input/ - orbital and spin degeneracy
! atocc, atmag /output/ - occupancy and magnetic moment of the atom
! temp /input/ - temperature
! verbosity - 0: no text output, 1: basics, 2: all
     
  USE hubbard_I_data  
  IMPLICIT NONE

! Input/output variables      
  INTEGER, INTENT(in) :: nlm, ns, Iwmax, nmom
  COMPLEX(KIND=8), INTENT(in) :: e0f(nlm*ns,nlm*ns)
  COMPLEX(KIND=8), INTENT(in) :: zmsb(Iwmax)
  REAL(KIND=8), INTENT(in) :: umn(nlm,nlm), ujmn(nlm,nlm)
  REAL(KIND=8), INTENT(in) :: ur(nlm,nlm,nlm,nlm)
  REAL(KIND=8), INTENT(in) :: temp
  integer, intent(in) :: verbosity
  COMPLEX(KIND=8), INTENT(out) :: GF(nlm*ns,nlm*ns,Iwmax)
  Complex(kind=8), intent(out) :: Tail(nmom,nlm*ns,nlm*ns)
  REAL(KIND=8), INTENT(out) :: atocc, atmag
! Local
  TYPE(occup), DIMENSION(0:nlm*ns) :: N_occ
  integer, allocatable :: arr(:)
  REAL(KIND=8) :: U(nlm*ns,nlm*ns,nlm*ns,nlm*ns)
  COMPLEX(KIND=8) :: zener
  real(8) :: Z, nomin, denom, E_B, tresh, maxexp, norm, ge
  real(8) :: fsign, Eground, Zterm
  real(8), allocatable :: E_A(:), ummss(:,:), occ(:), docc(:), ener(:)
  integer, allocatable :: nground(:)
  INTEGER, PARAMETER :: numexp=650
  integer :: i, j, m, m1, is, is1, iom, ls
  integer :: k, l, ideg, ie, i1, k1, Nat, NN
  integer :: iloc 
  INTEGER :: nso, nstate
  REAL(KIND=8), EXTERNAL :: factor
  REAL(KIND=8), PARAMETER :: tol=1d-7
  LOGICAL :: zerotemp, Efirst, Left, Right
      
  WRITE(*,'()')
  IF(temp < 1d-7) THEN
     zerotemp=.TRUE.
     if (verbosity>0) WRITE(*,*)'H.-I calculations are run at T=0'
  ELSE
     zerotemp=.FALSE.
     if (verbosity>0) WRITE(*,*)'H.-I calculations are run at T= ',temp
  ENDIF

  nso    = 2 * nlm
  nstate = 2**nso

! Diagonalize H_at with the 2-index U to estimate the ground state
! occupancy
      
    
  ALLOCATE( occ(nso), arr(0:nso-1) )
  FORALL( i = 0:nso-1 ) arr(i) = i
    
  allocate( ummss(nso,nso), docc(nso), ener(nso) )
  DO i=1,nso
     ener(i)     = REAL(e0f(i,i))
  ENDDO
    
  ummss(1:nlm,1:nlm)         = ujmn
  ummss(nlm+1:nso,nlm+1:nso) = ujmn
  ummss(1:nlm,nlm+1:nso)     = umn
  ummss(nlm+1:nso,1:nlm)     = umn
    
  allocate( E_A(0:nstate-1),nground(0:nstate-1) )
!
!    Initialize energy state E_A, A={ n_i sigma } and calculate Z
!
  do i = 0,nstate - 1
     occ = merge( 1.d0, 0.d0, btest( i, arr(0:nso-1) ) )
     E_A(i) = dot_product( ener, occ ) 
     docc = matmul( occ, ummss ) 
     E_A(i) = E_A(i) + dot_product( docc, occ ) / 2.d0
  enddo
  ge=MINVAL(E_A)
  nground=0
  atocc=0d0
  atmag=0d0
  DO i=0,nstate-1
     IF(ABS(ge-E_A(i)) < 1d-9) THEN
        nground(i)=1
        occ = merge( 1.d0, 0.d0, btest( i, arr(0:nso-1) ) )
        atocc=atocc+SUM(occ)
        atmag=atmag+SUM(occ(1:nlm))-SUM(occ(nlm+1:nso))
     ENDIF
  ENDDO
  norm=SUM(nground)
  atocc=atocc/norm
  atmag=atmag/norm
  if (verbosity>0) write(*,'(/,a,f12.5)')'Atomic occupancy with 2-ind U :',atocc 
  if (verbosity>0) write(*,'(/,a,f13.7)')'Ground state energy with 2-ind U :',ge

!  Set up and diagonalize H_at matrices for N-1, N, N+1
 
  CALL vertex4ind(Ur,U,nlm,ns)

  Nat=NINT(atocc)
  iloc=0
  DO i=0,nso
     N_occ(i)%ifdiag=.FALSE.
     N_occ(i)%run=.FALSE.
     N_occ(i)%ndeg=0
  ENDDO
  IF(Nat==nso) THEN
     N_occ(nso-1:nso)%run=.TRUE.
  ELSEIF(Nat==0) THEN
     N_occ(0:1)%run=.TRUE.
  ELSE
     N_occ(Nat-1:Nat+1)%run=.TRUE.
  ENDIF

  DO WHILE(iloc==0)

     DO i=0,nso
        if (verbosity>0) write(*,'(a,I7)')'==> Starting N = ',i
        IF(.NOT.N_occ(i)%run.OR.N_occ(i)%ifdiag) CYCLE
        N_occ(i)%ifdiag=.TRUE.
        N_occ(i)%n=NINT(factor(nso)/factor(nso-i)/factor(i))
        NN=N_occ(i)%n
        ALLOCATE(N_occ(i)%Hn(NN,NN),N_occ(i)%En(NN))
        ALLOCATE(N_occ(i)%st_n(NN),N_occ(i)%narr(0:nstate-1))
        CALL diagH(N_occ(i)%Hn,N_occ(i)%En,e0f,U,N_occ(i)%st_n,arr,N_occ(i)%narr,nso,nstate,i,N_occ(i)%n,verbosity)
        if (verbosity>1) WRITE(*,'(a,I7,a,F14.7)')'The lowest energy for N= ',i,' is ',N_occ(i)%En(1)
        if (verbosity>1) write(*,'(a,I7,a,/)')'i = ',i,' done! <=='

     ENDDO

! If occupancy is constrained do not seach for the true GS

     Efirst=.TRUE.
     DO i=0,nso
        IF(N_occ(i)%ifdiag.AND.Efirst) THEN
           Eground=N_occ(i)%En(1)
           Nat=i
           Efirst=.FALSE.
        ELSEIF(N_occ(i)%ifdiag) then
           if(N_occ(i)%En(1)<Eground) THEN
              Eground=N_occ(i)%En(1)
              Nat=i
           endif
        ENDIF
     ENDDO
     
     IF((Nat.NE.0).and.(.NOT.N_occ(Nat-1)%ifdiag)) THEN
        N_occ(Nat-1)%run=.TRUE.
     ELSEIF((Nat.NE.nso) .and.(.NOT.N_occ(Nat+1)%ifdiag)) THEN
        N_occ(Nat+1)%run=.TRUE.
     ELSE
        iloc=1 
     ENDIF

  ENDDO

  atocc=0d0
  atmag=0d0
  IF(zerotemp) THEN
     ! T=0 case
     DO i=1,N_occ(Nat)%n-1
        IF(ABS(N_occ(Nat)%En(i)-N_occ(Nat)%En(i+1))>tol) EXIT
     ENDDO
     N_occ(Nat)%ndeg=i
     Z=i
     DO i=1,N_occ(Nat)%ndeg
        DO k=1,N_occ(Nat)%n
           occ = merge(1.d0,0.d0,btest(N_occ(Nat)%st_n(k),arr(0:nso-1)))
           atocc=atocc+SUM(occ)*N_occ(Nat)%Hn(k,i)*CONJG(N_occ(Nat)%Hn(k,i))
           atmag=atmag+(SUM(occ(1:nlm))-SUM(occ(nlm+1:nso)))*N_occ(Nat)%Hn(k,i)*CONJG(N_occ(Nat)%Hn(k,i))
        ENDDO
     ENDDO

  ELSE
     ! Finite T case
     Eground=N_occ(Nat)%En(1)
     Z=0d0
     DO k=0,nso
        IF(.NOT.N_occ(k)%ifdiag) CYCLE
        DO i=1,N_occ(k)%n
           Zterm=EXP((Eground-N_occ(k)%En(i))/temp)
           IF(Zterm < tol) EXIT
           Z=Z+Zterm
           N_occ(k)%ndeg=i
           DO l=1,N_occ(k)%n
              occ = merge(1.d0,0.d0,btest(N_occ(k)%st_n(l),arr(0:nso-1)))
              atocc=atocc+SUM(occ)*N_occ(k)%Hn(l,i)*CONJG(N_occ(k)%Hn(l,i))*Zterm
              atmag=atmag+(SUM(occ(1:nlm))-SUM(occ(nlm+1:nso)))*N_occ(k)%Hn(l,i)*CONJG(N_occ(k)%Hn(l,i))*Zterm
           ENDDO
        ENDDO
     ENDDO
     
  ENDIF

  atocc=atocc/Z
  atmag=atmag/Z
  IF(zerotemp) THEN
     if (verbosity>0) WRITE(*,'(a,i2,a,i5,a,f13.6)') &
          &'The ground state has occupancy ',Nat,&
          &', degeneracy ',N_occ(Nat)%ndeg,' and energy ',N_occ(Nat)%En(1)
  ELSE
     if (verbosity>0) WRITE(*,'(/,a,i2,a,f13.6)') &
          &'The ground state has occupancy ',Nat, &
          &' and energy ',Eground
     if (verbosity>0) WRITE(*,'(a,i5,a)') &
          &'Transitions from  ',N_occ(Nat)%ndeg, &
          &' atomic states are included in GF'
        
     if (verbosity>0) WRITE(*,'(a,f13.6)')'Z = ',Z
  ENDIF
  if (verbosity>0) write(*,'(/,a,f12.5)')'Atomic occupancy  :',atocc            
  if (verbosity>0) write(*,'(/,a,f12.5)')'Atomic mag. mom.  :',atmag

  OPEN(450,file='ATOMIC_LEVELS')
  DO i=1,N_occ(Nat)%n
     WRITE(450,*)i,N_occ(Nat)%En(i)
  ENDDO
  CLOSE(450)

! Compute the Green's function

  if (verbosity>0) WRITE(*,'(/,a)')'Start GF calculations'
  ! Diagonalize H for additional occupancies if needed for GF
  DO k=0,nso
     IF(N_occ(k)%ndeg > 0) THEN
        IF(k>0.AND..NOT.N_occ(k-1)%ifdiag) N_occ(k-1)%run=.TRUE.
        IF(k<nso.AND..NOT.N_occ(k+1)%ifdiag) N_occ(k+1)%run=.TRUE.
     ENDIF
  ENDDO

  DO i=0,nso
     IF(.NOT.N_occ(i)%ifdiag.AND.N_occ(i)%run) THEN
        N_occ(i)%run=.TRUE.
        N_occ(i)%n=NINT(factor(nso)/factor(nso-i)/factor(i))
        NN=N_occ(i)%n
        ALLOCATE(N_occ(i)%Hn(NN,NN),N_occ(i)%En(NN))
        ALLOCATE(N_occ(i)%st_n(NN),N_occ(i)%narr(0:nstate-1))
        CALL diagH(N_occ(i)%Hn,N_occ(i)%En,e0f,U,N_occ(i)%st_n,arr,N_occ(i)%narr,nso,nstate,i,N_occ(i)%n,verbosity)
        if (verbosity>1) WRITE(*,'(/,a,I7,a,F14.7)')'The lowest energy for N= ',i,' is ',N_occ(i)%En(1)
     ENDIF
  ENDDO

  GF=(0d0,0d0) 
  Tail=(0d0,0d0)
  DO i=0,nso
     LEFT=.FALSE.
     RIGHT=.FALSE.
     k=i; l=i
     IF(i>0.AND.N_occ(i)%ndeg>0) THEN
        k=i-1
        LEFT=.TRUE.
     ENDIF
! Fix 15.11.2011
!    IF(i<nso.AND.N_occ(i)%ndeg>0.AND.N_occ(i+1)%ndeg==0) THEN
     IF(i<nso.AND.N_occ(i)%ndeg>0) THEN
        l=i+1
        RIGHT=.TRUE.
     ENDIF
     
     IF(LEFT.OR.RIGHT) then
        CALL add_to_GF_N(GF,Tail,arr,nso,nmom,Nat,Iwmax,l-k+1,N_occ(k:l),zmsb,Z,Eground,temp,zerotemp,LEFT,RIGHT)
     endif
  ENDDO
  IF (verbosity>1) THEN
      WRITE(*,*)'Notmalization of atomic GF:'
      DO m=1,nso 
        WRITE(*,*)Tail(1,m,m)
      ENDDO
  ENDIF

  deallocate( occ, E_A, ummss, docc, ener, arr, nground )
  DO i=0,nso
     IF(N_occ(i)%ifdiag) DEALLOCATE(N_occ(i)%Hn,N_occ(i)%En,N_occ(i)%st_n,N_occ(i)%narr)
  ENDDO
           
  RETURN
END SUBROUTINE gf_HI_fullU

      
REAL(KIND=8) function factor(N)

! factorial of N

  IMPLICIT NONE
  INTEGER, INTENT(in) :: N
  INTEGER :: i

  factor=1d0
  IF(N==0) RETURN
  DO i=1,N
     factor=factor*i
  ENDDO
  RETURN
END function factor


SUBROUTINE diagH(H,E,e0f,U,st,arr,narr,nso,nstate,Nat,n,verbosity)

! Initilize and diagonalize Hat for occupancy Nat

  IMPLICIT NONE
  INTEGER, INTENT(in) :: nso, nstate, Nat, n, verbosity
  INTEGER, INTENT(in) :: arr(nso)
  COMPLEX(KIND=8), INTENT(in) :: e0f(nso,nso)
  REAL(KIND=8), INTENT(in) :: U(nso,nso,nso,nso)
  ! Output
  COMPLEX(KIND=8), INTENT(out) :: H(n,n)
  INTEGER, INTENT(out) :: st(n)
  REAL(KIND=8), INTENT(out) :: E(n)
  ! Locals
  INTEGER :: occ(nso), narr(0:nstate-1)
  INTEGER :: i, k, l, p, q, m, m1, m2, m3, m4, INFO
  INTEGER :: ks, ls, ps, qs
  REAL(KIND=8) :: fsign
  REAL(KIND=8), DIMENSION(3*n-2) :: rwork
  COMPLEX(KIND=8), DIMENSION(n*n) :: WORK

  k=0
  narr=-1
  DO i=0,nstate-1
     occ = merge( 1.d0, 0.d0, btest( i, arr(1:nso) ) )
     IF(SUM(occ)==Nat) THEN
        k=k+1
        st(k)=i
        narr(i)=k
     ENDIF
  ENDDO
  IF(k.NE.n) THEN
     WRITE(*,*)'Error in diagH: array size is different from number of states'
     STOP
  ENDIF
  H=(0d0,0d0)
  if (verbosity>0) WRITE(*,'(a,i7)')'Set up the Hamiltonian for N = ',Nat
  ! Set up one-particle term
  DO m=1,nso
     DO m1=1,nso
        IF(ABS(e0f(m,m1)) < 1d-9) CYCLE
        DO i=1,n
           occ = NINT(merge( 1.d0, 0.d0, btest( st(i), arr(1:nso) ) ))
           k=ibclr(st(i),m1-1)
           IF(k==st(i)) CYCLE
           occ(m1)=0
           ks=SUM(occ(1:m1-1))
           l=ibset(k,m-1)
           IF(l==k) CYCLE
           ls=SUM(occ(1:m-1))
           fsign=(-1d0)**(ks+ls)
           H(narr(l),i)=H(narr(l),i)+e0f(m,m1)*fsign
        ENDDO
     ENDDO
  ENDDO
  ! Add interaction term
  DO m=1,nso
     DO m1=1,nso
        DO m2=1,nso
           DO m3=1,nso
              IF(ABS(U(m,m1,m2,m3)) < 1d-9) CYCLE
              DO i=1,n
                 occ = NINT(merge( 1.d0, 0.d0, btest( st(i), arr(1:nso) ) ))
                 k=ibclr(st(i),m2-1)
                 IF(k==st(i)) CYCLE
                 occ(m2)=0
                 ks=SUM(occ(1:m2-1))
                 l=ibclr(k,m3-1)
                 IF(l==k) CYCLE
                 occ(m3)=0
                 ls=SUM(occ(1:m3-1))
                 p=ibset(l,m1-1)
                 IF(p==l) CYCLE
                 occ(m1)=1
                 ps=SUM(occ(1:m1-1))
                 q=ibset(p,m-1)
                 IF(q==p) CYCLE
                 qs=SUM(occ(1:m-1))
                 fsign=(-1d0)**(ks+ls+ps+qs)
                 H(narr(q),i)=H(narr(q),i)+0.5d0*U(m,m1,m2,m3)*fsign
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

! Diagonalize H

  if (verbosity>1) WRITE(*,'(a,i7)')'The Hamiltonian is set up'
  CALL ZHEEV('V','U',n,H,n,E,WORK,n*n,RWORK,INFO)
  IF(INFO.NE.0) THEN
     WRITE(*,*)'diagH : error in the matrix diagonalization'
     STOP
  ENDIF
  if (verbosity>1) WRITE(*,'(a)')'The Hamiltonian is diagonalized'

  RETURN
END SUBROUTINE diagH


subroutine sigma_atomic_fullU(GF,Sigma_mat,e0f,zmsb,nlm,Iwmax,ns)

! Computes Sigma from GF !!! Relativistic version

  implicit none

  ! Input variables
  INTEGER, INTENT(in) :: Iwmax, nlm, ns
  COMPLEX(KIND=8),INTENT(in) :: e0f(nlm*ns,nlm*ns)
  COMPLEX(KIND=8),INTENT(in) :: zmsb(Iwmax)
  COMPLEX(KIND=8),INTENT(in) :: GF(nlm*ns,nlm*ns,Iwmax)
  ! Output variables
  COMPLEX(KIND=8),INTENT(out) :: Sigma_mat(nlm*ns,nlm*ns,Iwmax)
  ! Locals
  integer :: iom, m, m1, is, nlms
  integer :: ind, ind1, INFO
  INTEGER, DIMENSION(nlm*ns*nlm*ns) :: IPIV
  COMPLEX(KIND=8), DIMENSION(nlm*ns) :: WORK
  COMPLEX(KIND=8), DIMENSION(nlm*ns,nlm*ns) :: g0, g

  nlms=ns*nlm
  DO iom = 1,Iwmax
     g0=-e0f
     DO m = 1,nlms
        g0(m,m)=g0(m,m)+zmsb(iom)
     ENDDO
     g=GF(:,:,iom)
!------- LAPACK inversion of GF
     CALL ZGETRF(nlms,nlms,g,nlms,IPIV,INFO)
     IF(info.ne.0) STOP'in SPTFLEX ZGETRF Info NE 0'
        
     CALL ZGETRI(nlms,g,nlms,IPIV,WORK,nlms,INFO)
     IF(info.ne.0) STOP'in SPTFLEX ZGETRI Info NE 0'

     DO ind=1,nlm*ns
        DO ind1=1,nlm*ns
           Sigma_mat(ind,ind1,iom)=g0(ind,ind1)-g(ind,ind1)
        ENDDO
     ENDDO
         
  ENDDO

end subroutine sigma_atomic_fullU
         



SUBROUTINE vertex4ind(U,Uc,nlm,ns)

  IMPLICIT NONE
  ! input/output    
  REAL(KIND=8), INTENT(in), DIMENSION(nlm,nlm,nlm,nlm) :: U
  INTEGER, INTENT(in):: nlm, ns
  REAL(KIND=8), INTENT(out), DIMENSION(nlm*ns,nlm*ns,nlm*ns,nlm*ns)  :: Uc
  ! Locals
  INTEGER :: is1, is2, is3, is4
  !-------------- Vertex 
  Uc=0d0
  DO is1=1,ns
     DO is2=1,ns
        Uc((is1-1)*nlm+1:is1*nlm,(is2-1)*nlm+1:is2*nlm,(is1-1)*nlm+1:is1*nlm,(is2-1)*nlm+1:is2*nlm)=U(:,:,:,:)
     ENDDO
  ENDDO
  RETURN
END SUBROUTINE vertex4ind


SUBROUTINE add_to_GF_N(GF,Tail,arr,nso,nmom,Nat,Iwmax,num,N_occ,zmsb,Z,Eground,temp,zerotemp,LEFT,RIGHT)

! Adds to the GF the contribution associated to transitions from the states with
! occupancy N

  USE hubbard_I_data  
  IMPLICIT NONE

  INTEGER, INTENT(in) :: nso, Nat, Iwmax, num,nmom
  LOGICAL, INTENT(in) :: zerotemp, LEFT, RIGHT      
  REAL(KIND=8), INTENT(in) :: Z, Eground, temp
  COMPLEX(KIND=8), INTENT(in) :: zmsb(Iwmax)
  INTEGER, INTENT(in) :: arr(0:nso-1)     
  TYPE(occup), INTENT(in) :: N_occ(num)

  COMPLEX(KIND=8), INTENT(inout) :: GF(nso,nso,Iwmax)
  complex(kind=8), intent(inout) :: Tail(nmom,nso,nso)

  COMPLEX(KIND=8) :: Gstore(nso,nso,Iwmax)
  complex(kind=8) :: Tailstore(nmom,nso,nso)
  COMPLEX(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: anmat,crmat  
  COMPLEX(KIND=8) :: zener
  INTEGER :: i, k, k1, num1, m, m1, l, ls, ie      
  INTEGER :: occ(nso)
  REAL(KIND=8) :: fsign, ecoff  
     
  Gstore=GF
  GF=(0d0,0d0)
  Tailstore=Tail
  Tail=(0d0,0d0)

  IF(LEFT) THEN
     ! Matrix elements |<N-1|d_m|N>|
     ALLOCATE(anmat(nso,N_occ(1)%n,N_occ(2)%ndeg))
     anmat=0d0
     DO i=1,N_occ(2)%ndeg
        DO k=1,N_occ(2)%n
           DO m=1,nso
              occ = merge(1.d0,0.d0,btest(N_occ(2)%st_n(k),arr(0:nso-1)))
              l=ibclr(N_occ(2)%st_n(k),m-1)
              IF(l==N_occ(2)%st_n(k)) CYCLE
              ls=SUM(occ(1:m-1))
              fsign=(-1d0)**ls
              DO k1=1,N_occ(1)%n
                 anmat(m,k1,i)=anmat(m,k1,i)+CONJG(N_occ(2)%Hn(k,i))*N_occ(1)%Hn(N_occ(1)%narr(l),k1)*fsign
              ENDDO
           ENDDO
        ENDDO
     ENDDO
     ! Compute contribution to GF
     DO i=1,N_occ(2)%ndeg              ! sum over contributing ground states
        DO k=1,N_occ(1)%n              ! sum over excited states 
           IF(zerotemp) THEN
              ecoff=1d0
           ELSEIF(k >N_occ(1)%ndeg) THEN
              ecoff=EXP((Eground-N_occ(2)%En(i))/temp)
           ELSE
              ecoff=EXP((Eground-N_occ(1)%En(k))/temp)+EXP((Eground-N_occ(2)%En(i))/temp)
           ENDIF
           DO ie=1,Iwmax
              zener=1d0/(zmsb(ie)-N_occ(2)%En(i)+N_occ(1)%En(k))*ecoff
              DO m=1,nso
                 DO m1=1,nso
                    GF(m,m1,ie)=GF(m,m1,ie)+anmat(m,k,i)*CONJG(anmat(m1,k,i))*zener
                 ENDDO
              ENDDO
           ENDDO

           ! Tails:
           do ie=1,nmom   ! number of moments 
              do m=1,nso
                 do m1=1,nso
                    Tail(ie,m,m1) = Tail(ie,m,m1) + ecoff * anmat(m,k,i)*CONJG(anmat(m1,k,i)) * &
                         &(N_occ(2)%En(i)-N_occ(1)%En(k))**(ie-1)
                 enddo
              enddo
           enddo

        ENDDO
     ENDDO
     DEALLOCATE(anmat)
  ENDIF

  IF(RIGHT) THEN
     ! Matrix elements |<N+1|d_m^+|N>|
     num1=num-1
     ALLOCATE(crmat(nso,N_occ(num)%n,N_occ(num1)%ndeg))
     crmat=0d0
     DO i=1,N_occ(num1)%ndeg
        DO k=1,N_occ(num1)%n
           DO m=1,nso
              occ = merge(1.d0,0.d0,btest(N_occ(num1)%st_n(k),arr(0:nso-1)))
              l=ibset(N_occ(num1)%st_n(k),m-1)
              IF(l==N_occ(num1)%st_n(k)) CYCLE
              ls=SUM(occ(1:m-1))
              fsign=(-1d0)**ls
              DO k1=1,N_occ(num)%n
                 crmat(m,k1,i)=crmat(m,k1,i)+CONJG(N_occ(num1)%Hn(k,i))*N_occ(num)%Hn(N_occ(num)%narr(l),k1)*fsign
              ENDDO
           ENDDO
        ENDDO
     ENDDO
     ! Compute contribution to GF
     DO i=1,N_occ(num1)%ndeg
     !   DO k=1,N_occ(num)%n
     ! Fix 15.11.2011: only states above ndeg are included
        DO k=N_occ(num)%ndeg+1,N_occ(num)%n
           IF(zerotemp) THEN
              ecoff=1d0
           ELSEIF(k >N_occ(num)%ndeg) THEN
              ecoff=EXP((Eground-N_occ(num1)%En(i))/temp)
           ELSE
              ecoff=EXP((Eground-N_occ(num)%En(k))/temp)+EXP((Eground-N_occ(num1)%En(i))/temp)
           ENDIF
           DO ie=1,Iwmax
              zener=1d0/(zmsb(ie)-N_occ(num)%En(k)+N_occ(num1)%En(i))*ecoff
              DO m=1,nso
                 DO m1=1,nso
                    GF(m,m1,ie)=GF(m,m1,ie)+CONJG(crmat(m,k,i))*crmat(m1,k,i)*zener
                 ENDDO
              ENDDO
           ENDDO

           ! Tails:
           do ie=1,nmom   ! number of moments
              do m=1,nso
                 do m1=1,nso
                    Tail(ie,m,m1) = Tail(ie,m,m1) + ecoff * CONJG(crmat(m,k,i))*crmat(m1,k,i) * &
                         &(N_occ(num)%En(k)-N_occ(num1)%En(i))**(ie-1)
                 enddo
              enddo
           enddo


        ENDDO
     ENDDO
     DEALLOCATE(crmat)
  ENDIF
  GF=Gstore+GF/Z
  Tail = Tailstore + Tail/Z
  RETURN
END SUBROUTINE add_to_GF_N
