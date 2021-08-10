module constants  
implicit none 
!=======================================================================
!                             parameters
   double precision , parameter :: pi = 3.1415926536    
   double precision , parameter :: epsil=1.d-5, hc=197.33d0 !MeV.fm
   double precision , parameter :: a0  = - 18.7 !fm 
   double precision , parameter :: r1  = - 0.8819 !fm^-1
   double precision , parameter :: a1  = - 62.951      !fm^3
   double precision , parameter :: gamma1=-r1/2., gamma0=1./a0  ! fm⁻1
   double precision , parameter :: mn = 939.56533/hc ,mp = 938.27203/hc          !fm^-1
   double precision , parameter :: kR  =   dsqrt(2.d0/a1/r1) ! fm^-1 

 
end module constants 


!=======================================================================
!                         Main Program                      
!=======================================================================
program three_body_bound_state
use constants  
! This code compute three-body bound states with halo EFT -> core-n-n
! To compute three-body bound energy  B3[core-n-n]

 integer      ::  ii,i,j,kk ! to do loops 
 integer , parameter ::  nn=60, niter=60
 double precision :: q(nn), wq(nn), B(niter), y(nn), dy(nn) 
 double precision :: taun ,Xnc ,Xcn, Xnn ,tauc , knc, knn, kcn
 double precision :: lambda,q0,c,d 
 double precision :: matriz(2*nn,2*nn) ,det(niter) , dd 
 double precision :: mc, A,   B3,qi,qj 


      open(unit=30,file='B3xLambda.dat', status='unknown')
    
!=======================================================================
!                               inputs
      print*,''
      print*,' Hazard a guess  initial  value for 3B bound energy  [MeV]'
      read*,B(1)
      B(1)=B(1)/hc
      B(2)=B(1)*1.01
       print*,' core mass number'
      read*,mc
      mc = mc*mn          !fm^-1
      A  =  (mc/mn)
!=======================================================================
 


do 2 kk=1,175                    ! varia lambda

lambda=650./hc   + 2*kk/hc                      ! cutoff  fm^-1 





call legauss(-1.d0,1.d0,nn,y,dy,1.d-15)

q0 = 0./hc
do i=1,nn          ! constroi o vetor q(nn) que  varia de 0 a lambda
  
  c=(q0+lambda)*0.5d0
  d=(lambda-q0)*0.5d0
!  
     q(i)=c+d*y(i)
     wq(i)=d*dy(i)
!     
enddo   
! 
!
det=0.d0
    do 1 ii=1,niter ! loop para calcular determinante ii


!  B = energia de ligacao  
         if (ii.gt.2) then 
    B(ii)=B(ii-1)-(B(ii-1)-B(ii-2))/( (det(ii-1))- (det(ii-2)))* (det(ii-1))    ! newton-rapshon para a energia
  

         endif
            B3=B(ii)    ! guardando o valor da energia 
            print*, 'B3=', b3
!
!========================= Matriz M ========================================
     
!
!
matriz=0.d0
       do 4 i=1,nn             !q_i         !varia a linha
!       
             qi=q(i)
!
         do 5 j=1,nn           !q_j          !varia a coluna
             qj=q(j)  
!
        
        kcn =   4.*pi*qj**2*Xcn(qi,qj, A, B3)*taun(qj,A,B3 )*wq(j)
        
        knc =   4.*pi*qj**2*Xnc(qi,qj, A, B3)*tauc(qj,A,B3)*wq(j)
        
        knn =   4.*pi*qj**2*Xnn(qi,qj, A, B3)*taun(qj,A, B3)*wq(j)

 
! 
!
!     Bloco (1,1)                 ! matriz identidade 
!
            matriz(i,j)       = 0.
   if (i==j)matriz(i,j)       = -1.
! 
!     Bloco(1,2)                           bloco CN
!           
            matriz(i,j+nn)    = 2.*kcn 
!
!     Bloco(2,1)                           bloco NC
!
           matriz(i+nn,j)    = knc
!
!     Bloco(2,2)                           bloco NN
!
           matriz(i+nn,j+nn) = knn                        
       if (i==j) matriz(i+nn,j+nn) = knn - 1.    
!
!        
5        continue           !laço de q i    linha
4      continue             !laço de q_j    coluna
!
!
  

          call duminv(matriz,2*nn,2*nn,dd)        !subrotina duminv retorna determinante
         
          det(ii)=dd 
   
!         
!   

       
      

        if(ii.eq.niter) then
        print*, 'nao convergiu para ', niter, 'interacoes', 'det=', dd, 'lambda=',lambda*hc, 'B3=',B3*hc
        stop
        endif

        print*,  ii,  lambda*hc,  'det=', dd, 'E =', B(ii)*hc 
        if(ii.gt.2.and.dabs((B(ii-1)-B(ii))/B(ii-1)).lt.epsil)then 
        print*, 'det=', dd, 'lambda=',lambda*hc, 'Convergiu para E =', B3*hc 
        exit
        

        
       endif
!
!
!=======================================================================         
!         
1    continue              ! determinante ii 
 
write(30,*)  lambda*hc , B3*hc 

2    continue!lambda laço kk
!
!
end program 
!=======================================================================
!=======================================================================
!                              functions

!-----------------------------------------------------------------------
      double precision function taun(q,A,EB)
      use constants         
        double precision ::  q , kn, A,EB 
         
          Kn = dsqrt(2.*A/(A+1.)*(mn*EB + (A+2.)/2./(A+1.)*q**2))
          taun = -1./(4.*pi**2*mn*gamma1)*(A+1.)/A /(Kn**2+kR**2)

        return  
        end
        
!-----------------------------------------------------------------------
      double precision function tauc(q ,A,EB)
      use constants        
        double precision ::  q  , kc, EB ,A
        
          Kc = dsqrt(mn*EB  +  (A+2.)/4./A*q**2)
          tauc =  1./(2.*pi**2*mn) / (gamma0-Kc) 
          
        return  
        end        

!-----------------------------------------------------------------------
      double precision function Xcn(q1,q2,A,EB)
      use constants 
        double precision :: q1,q2, Q0cn, Q1cn, EB, A ,zcn!,intpoli0,intpoli1
        
        zcn = -1/q1/q2*(mn*EB + (A+1.)/(2.*A)*(q1**2) + q2**2  )  
         Q0cn= 0.5*log((zcn+1.)/(zcn-1.))
         Q1cn= 0.5*zcn*log((zcn+1.)/(zcn-1.)) - 1. 
        !Q0cn=0.5*intpoli0(zcn)
        !Q1cn=0.5*intpoli1(zcn)

        Xcn = - dsqrt(2.d0)*mn*(A/(A+1.)/q1 * Q0cn  + 1./q2 *Q1cn)
                  
      return
      end
! ----------------------------------------------------------------------
      double precision function Xnc(q1,q2,A,EB)
      use constants 
        double precision :: q1,q2, Q0nc, Q1nc, EB, A,znc!,intpoli0,intpoli1
        znc = -1/q1/q2*(mn*EB + q1**2 + (A+1.)/(2.*A)*(q2**2))      
         Q0nc= 0.5*log((znc+1.)/(znc-1.))
         Q1nc= 0.5*znc*log((znc+1.)/(znc-1.)) - 1.
         
           
        !Q0nc= 0.5*intpoli0(znc)
        !Q1nc= 0.5*intpoli1(znc) 
        Xnc = - dsqrt(2.d0)*mn*(A/(A+1.)/q2 * Q0nc  + 1./q1 *Q1nc)  
      return
      end
! ----------------------------------------------------------------------
      double precision function Xnn(q1,q2,A,EB)
      use constants 
        double precision :: q1,q2, Q0nn, Q1nn, Q2nn,EB, A, znn!, intpoli0, intpoli1,intpoli2
      
        znn = -A/q1/q2*(mn*EB + (A+1.)/(2.*A)*(q1**2+q2**2))
         Q0nn= 0.5*log((znn+1.)/(znn-1.))
         Q1nn= 0.5*znn*log((znn+1.)/(znn-1.))-1.
         Q2nn= 0.5*(-0.5 + 3./2.*znn**2) *  log((znn+1)/(znn-1))-3./2.*znn
         
        !Q0nn= 0.5*intpoli0(znn)
        !Q1nn= 0.5*intpoli1(znn) 
        !Q2nn= 0.5*intpoli2(znn) 
         
        Xnn =  A*mn*( (A**2+2.*A+3.)/(A+1.)**2*Q0nn + 2./(A+1.)*(q1**2+q2**2)/q1/q2*Q1nn + Q2nn)
       
      return
      end      


! ----------------------------------------------------------------------
      double precision function intpoli0(z)   
       integer   :: n 
       double precision x(12),dx(12)   ,z 
      
        n=12
       call legauss(-1.d0,1.d0,n,x ,dx,1.d-15)  
     
       intpoli0=0.
       do i=1,n
       intpoli0=intpoli0 + 1.d0/(z-x(i))*dx(i)
       enddo
       
      return
      end
      

! ----------------------------------------------------------------------
      double precision function intpoli1(z) 
       integer   :: n   
       double precision x(12),dx(12)   ,z 
      
         n=12
         call legauss(-1.d0,1.d0,n,x ,dx,1.d-15) 
        !intpoli1=   z*(log(z+1.)-log(z-1.) )                              ! para |(z)| >1
      
        intpoli1=0.
        do i=1,n
         intpoli1=intpoli1 + z/(z-x(i))*dx(i)
        enddo 
      return
      end
      
      
! ----------------------------------------------------------------------
      double precision function intpoli2(z) 
      integer   :: n 
      double precision x(12),dx(12)   ,z , poli2
      
       n=12
      call legauss(-1.d0,1.d0,n,x ,dx,1.d-15) 
 
      intpoli2=0.
      do i=1,n
      poli2 = 0.5*(3.*x(i)**2 - 1.)
      intpoli2=intpoli2 + poli2/(z-x(i))*dx(i)
      enddo
      
      return
      end
      

!=======================================================================
!=======================================================================
!
! 
      subroutine legauss(xs,xl,n,x,dx,zz)
!
!  xs e xl - limites inferior e superior para os pontos de gauss
!  n - número de pontos de gauss
!  x - pontos de gauss
!  dx - peso
!  xx - fator de convergência
!
      implicit double precision (a-h,o-z)
      integer n
      dimension x(n),dx(n)
      if(n)10,10,20
10    write(5,600) n
      write(2,600) n
600  format(1h ,i10,' rejeitados ptos.leg-gauss')
      return
20   if(n-2) 30,40,40
30   x(1)=0.d0
      dx(1)=.5d0
      go to 140
40   i=1
      g=-1.d0
      ic=(n+1)/2
50   s=g
      t=1.d0
      u=1.d0
      v=0.d0
      do 60 k=2,n
      a=k
      fact1=(2.d0*a-1.d0)/a
      fact2=(a-1.d0)/a
      p=fact1*g*s-fact2*t
      dp=fact1*(s+g*u)-fact2*v
      t=s
      s=p
      v=u
60   u=dp
      sum=0.d0
      if(i-1)90,90,70
70   im1=i-1
      do 80 k=1,im1
80   sum=sum+1.d0/(g-x(k))
90   test=g
      g=g-p/(dp-p*sum)
      r=dabs(test-g)
      if(r.lt.zz)goto 100
      goto 50
100  r=n
      x(i)=g
      dx(i)=2.d0/r/t/dp
      if(ic-i)120,120,110
110  fim1=im1
g=g-(dp-p*sum)/((2.d0*g*dp-a*(a+1.d0)*p)/(1.d0-g*g)-2.d0*dp*sum-p*sum**2+fim1*p)
      i=i+1
      goto 50
120  k0=2*ic-n+2*(n/2)+1
      ic=ic+1
      do 130 i=ic,n
      k=k0-i
      x(i)=-x(k)
130  dx(i)=dx(k)
140  fact1=(xl-xs)/2.d0
      fact2=(xl+xs)/2.d0
      do 150 i=1,n
      dx(i)=dx(i)*fact1
150  x(i)=x(i)*fact1+fact2
      return
      end
      



!=======================================================================
!=======================================================================
!      *******SUBROTINA DUMINV ********* 
!      ****CALCULA O DETERMINANTE*******
!      ***E INVERTE A MATRIZ A(I,J)*****
!	 *********************************
! 
!
! 
      SUBROUTINE duminv(A,N,NN,D) 
!
!	A - MATRIZ ASER INVERTIDA 
!     N - TAMANHO DA MATRIZ 
!     NN - TAMANHO MAXIMO DA MATRIZ 
!     D - DETERMINANTE
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
!      COMPLEX*16 A,D,BIGA,HOLD 
      DIMENSION A(N*N),L(N),M(N)
! 
!     PROCURA DO MAIOR ELEMENTO 
! 
      D=1.0 
      NK=-NN
      DO 80 K=1,N 
      NK=NK+NN
      L(K)=K
      M(K)=K
      KK=NK+K 
      BIGA=A(KK)
      DO 20 J=K,N 
      IZ=NN*(J-1) 
      DO 20 I=K,N 
      IJ=IZ+I 
!10    IF(cDABS(BIGA)-cDABS(A(IJ)))15,20,20
10    IF(DABS(BIGA)-DABS(A(IJ)))15,20,20
15    BIGA=A(IJ)
      L(K)=I
      M(K)=J
20    CONTINUE
! 
!      TROCA LINHAS 
! 
      J=L(K)
      IF(J-K)35,35,25 
25    KI=K-NN 
      DO 30 I=1,N 
      KI=KI+NN
      HOLD=-A(KI) 
      JI=KI-K+J 
      A(KI)=A(JI) 
30    A(JI)=HOLD
! 
!       TROCA   COLUNAS 
! 
35    I=M(K)
      IF(I-K)45,45,38 
38    JP=NN*(I-1) 
      DO 40 J=1,N 
      JK=NK+J 
      JI=JP+J 
      HOLD=-A(JK) 
      A(JK)=A(JI) 
40    A(JI)=HOLD
! 
! DIVIDE A COLUNA POR (-PIVOT) - O VALOR DO PIVO ESTA CONTIDO EM BIGA 
!  
!45    IF(cDABS(BIGA))48,46,48
45    IF(DABS(BIGA))48,46,48
46    D=0.0 
      RETURN
48    DO 55 I=1,N 
      IF(I-K)50,55,50 
50    IK=NK+I 
      A(IK)=A(IK)/(-BIGA) 
55     CONTINUE 
! 
!       REDUZ A MATRIZ
! 
      DO 65 I=1,N 
      IK=NK+I 
      HOLD=A(IK)
      IJ=I-NN 
      DO 65 J=1,N 
      IJ=IJ+NN
      IF(I-K)60,65,60 
60    IF(J-K)62,65,62 
62    KJ=IJ-I+K 
      A(IJ)=HOLD*A(KJ)+A(IJ)
65    CONTINUE
! 
!      DIVIDE A LINHA PELO PIVO 
! 
      KJ=K-NN 
      DO 75 J=1,N 
      KJ=KJ+NN
      IF(J-K)70,75,70 
70    A(KJ)=A(KJ)/BIGA
75    CONTINUE
! 
!      PRODUTO DO PIVO
!!
      D=D*BIGA
! 
!      TROCA O PIVO PELO RECIPROCO
! 
      A(KK)=1./BIGA 
80     CONTINUE 
! 
!      TROCA FINAL DE LINHA E COLUNA
! 
      K=N 
100   K=K-1 
      IF(K)150,150,105
105   I=L(K)
      IF(I-K)120,120,108
108   JQ=NN*(K-1) 
      JR=NN*(I-1) 
      DO 110 J=1,N
      JK=JQ+J 
      HOLD=A(JK)
      JI=JR+J 
      A(JK)=-A(JI)
110   A(JI)=HOLD
120   J=M(K)
      IF(J-K)100,100,125
125   KI=K-NN 
      DO 130 I=1,N
      KI=KI+NN
      HOLD=A(KI)
      JI=KI-K+J 
      A(KI)=-A(JI)
130   A(JI)=HOLD
      GO TO 100 
150   RETURN
      END 


!*******************************************************************************
! REID93: Updated Reid potential, regularized with a dipole form factor
! of 8 pion masses, in momentum space on LSJ basis.
! Reference: V.G.J. Stoks et al., Phys. Rev. C 49, 2950 (1994).
!-------------------------------------------------------------------------------
! This is f90 version 1.1, 20 November 1999
! WWW: http://nn-online.org
! Email: info@nn-online.org
!-------------------------------------------------------------------------------
! IN:  real(wp)     :: qi        = center of mass momentum initial state in MeV
!      real(wp)     :: qf        = center of mass momentum final   state in MeV
!      character(2) :: type      = reaction; 'PP', 'NN', 'NP', or 'PN'
!      character(3) :: name      = name of the partial wave (see below),
!                                  maximum total angular momentum j = 9
! OUT: real(wp)     :: vpot(2,2) = potential matrix in MeV^{-2}
!-------------------------------------------------------------------------------
! Integer variable 'wp' defines the precision of 'real' variables. It is by
! default put to double precision. This code is made with double-precision
! of 8 bytes in mind.  We do not know about the behavior of the code for
! other precisions, but do not expect anything anomalous.
!-------------------------------------------------------------------------------
! Defining the K-matrix as   2i*mu*q*K = (1-S)(1+S)^-1
! (i.e. for singlet channel   tan(delta) = -2*mu*q*K )
! the partial-wave Lippmann-Schwinger equation reads
!     K(q',q) = V(q',q) + 2/pi int dk k^2 V(q',k) G(q,k) K(k,q)
! with
!     G(q,k) = P / (E(q) - k^2/2/mu)
!     V(q',k) = 1 / (4*pi) * vpot(qi=k,qf=q')
! S.Y.M. convention is used, so additional minus sign in off-diagonal tensor
! potential (shows up in mixing parameter).
!-------------------------------------------------------------------------------
! Potential decomposition in momentum space plane-wave basis:
! V(qf,qi) =    Vc
!             + Vs   (sigma_1.sigma_2)        (only in one-pion-exchange)
!             + Vt   [(sigma_1.k)(sigma_1.k)-k2/3(sigma_1.sigma_1)]
!             + Vls  (i/2)(sigma_1+sigma_1).n
!
! k = qf - qi ,   q = (qf+qi)/2 ,   n = qi x qf = q x k
!-------------------------------------------------------------------------------
! The OPE-potential distinguishes between neutral and charged pion masses,
! and has coupling constants f0pi**2 = fcpi**2 = 0.075. The delta-function
! (smeared out due to the form factor) only contributes to the S waves.
!-------------------------------------------------------------------------------
! The variable 'name' contains the name of the partial wave in spectral
! notation:  singlets:                 1S0  1P1  1D2  1F3  1G4 ...
!            triplets uncoupled:       3P0  3P1  3D2  3F3  3G4 ...
!            triplets coupled:              3C1  3C2  3C3  3C4 ...
! where 3C1 denotes   3S1 - EPS1 - 3D1 channel
!       3C2 denotes   3P2 - EPS2 - 3F2 channel
!       ...
! Only capitals are recognized!
!*******************************************************************************

      subroutine reid93(qi,qf,type,name,vpot)

      implicit none

!     integer, parameter        :: wp = selected_real_kind(15)
      integer, parameter        :: wp = kind(1.d0)       ! double precision

      double complex, intent(in):: qi,qf
      character(2), intent(in)  :: type
      character(3), intent(in)  :: name
      double complex, intent(out):: vpot(2,2)

      integer, save             :: nchan,l,spin,j,iso
      integer, save             :: jmm,jm,jp,jpp,jmax
      real(wp), save            :: tjmm,tjm,tj,tjp,tjpp,tjj
      character(3), save        :: name0 = '***'
      integer                   :: ll,im
      complex(wp)                  :: eln(0:12)
      complex(wp)                  :: vpis(-1:1),vpit(-2:2)
      complex(wp)                  :: vc(-1:1),vl(-2:2),vt(-2:2)
      complex(wp)                  :: fac,vten,vpisl0,l2,m2,x,y
      complex(wp)                  :: qiqf,qi2,qf2,qi2f2,s2psi,cpsi2,spsi2

      real(wp), parameter       :: f0pi = 0.075_wp
      real(wp), parameter       :: fcpi = 0.075_wp
      real(wp), parameter       :: pi = 3.1415926535897932384626433832795028_wp
      real(wp), parameter       :: m_pi0 = 134.9739_wp
      real(wp), parameter       :: m_pic = 139.5675_wp
      real(wp), parameter       :: m_pis = m_pic
      real(wp), parameter       :: m_pis2 = m_pis**2
      real(wp), parameter       :: m_pi = (m_pi0+2.0*m_pic)/3.0
      real(wp), parameter       :: m_pi2 = m_pi**2
      real(wp), parameter       :: zero = 0.0_wp
      real(wp), parameter       :: half = 0.5_wp
      real(wp), parameter       :: one = 1.0_wp
      real(wp), parameter       :: two = 2.0_wp
      real(wp), parameter       :: three = 3.0_wp
      real(wp), parameter       :: four = 4.0_wp
      real(wp), parameter       :: five = 5.0_wp

      complex(wp)                  :: ul(6,-2:12) = zero
      complex(wp)                  :: ulc(-2:12) = zero

      real(wp), save            :: a(5,5), b(5,5)
      logical, save             :: first = .true.

! pp parameters (isospin = 1)
      real(wp), parameter       :: aa(25) =                                 &
      (/  0.1756084e0_wp,-0.1414234e2_wp, 0.1518489e3_wp,-0.6868230e3_wp,   &
          0.1104157e4_wp,-0.4224976e2_wp, 0.2072246e3_wp,-0.3354364e3_wp,   &
         -0.1989250e1_wp,-0.6178469e2_wp, 0.2912845e2_wp, 0.1511690e3_wp,   &
          0.8151964e1_wp, 0.5832103e2_wp,-0.2074743e2_wp,-0.5840566e0_wp,   &
         -0.1029310e2_wp, 0.2263391e2_wp, 0.2316915e2_wp,-0.1959172e1_wp,   &
         -0.2608488e1_wp, 0.1090858e2_wp,-0.4374212e0_wp,-0.2148862e2_wp,   &
         -0.6584788e0_wp                                                    &
      /)

! np parameters (isospin = 0)
      real(wp), parameter       :: bb(25) =                                 &
      (/ -0.2234989e2_wp, 0.2551761e3_wp,-0.1063549e4_wp, 0.1609196e4_wp,   &
         -0.3505968e1_wp,-0.4248612e1_wp,-0.5352001e1_wp, 0.1827642e3_wp,   &
         -0.3927086e3_wp, 0.5812273e2_wp,-0.2904577e1_wp, 0.3802497e2_wp,   &
          0.3395927e0_wp, 0.8318097e0_wp, 0.1923895e1_wp, 0.0913746e0_wp,   &
         -0.1274773e2_wp, 0.1458600e3_wp,-0.6432461e3_wp, 0.1022217e4_wp,   &
         -0.0461640e0_wp, 0.7950192e1_wp,-0.1925573e1_wp, 0.5066234e2_wp,   &
          0.8359896e1_wp                                                    &
      /)

      if (first) then
          a = reshape(aa, (/5,5/), order=(/2,1/) )
          b = reshape(bb, (/5,5/), order=(/2,1/) )
          first = .false.
      endif

      if (name /= name0) then
          name0 = name
          nchan = 1
          if (name(2:2) == 'C') nchan = 2
          if (name(1:1) == '1') spin = 0
          if (name(1:1) == '3') spin = 1
          read(name,'(2x,i1)') j
          if (j > 9) stop 'REID93: j > 9 not allowed'
          l = j
          if (name == '3P0') l = 1
          if (nchan == 2) l = j-1
          iso = mod(spin+l+1,2)
          jmm  = j-2
          jm   = j-1
          jp   = j+1
          jpp  = j+2
          tjmm = two*j-three
          tjm  = two*j-one
          tj   = two*j+one
          tjp  = two*j+three
          tjpp = two*j+five
          tjj  = sqrt(j*(j+one))
          jmax = j+2
      endif

      qi2 = qi**2
      qf2 = qf**2
      qiqf = qi*qf
      qi2f2 = qi2+qf2
      s2psi = two*qiqf/qi2f2
      spsi2 = qf2/qi2f2
      cpsi2 = qi2/qi2f2

!     Neutral OPE
      m2 = m_pi0**2
      l2 = 64.0*m2
      x = half*(qi2f2+m2)/qiqf
      y = half*(qi2f2+l2)/qiqf
      call xdip(x,y,jmax,eln)
      ul(1,0:jmax) = eln(0:jmax)
      fac = two*pi*f0pi/(m_pis2*qiqf)
      do ll=-1,1
          vpis(ll) = fac*ul(1,j+ll)*m_pi0**2/three
      enddo
      vpisl0 = four*pi*f0pi/(m_pis2*three) * (y-x)**2/(one-y**2)
      do ll=-2,2
          vpit(ll) = -fac*ul(1,j+ll)
      enddo

!     Charged OPE
      select case(type)
      case('NP','PN','np','pn')
          m2 = m_pic**2
          l2 = 64.0*m2
          x = half*(qi2f2+m2)/qiqf
          y = half*(qi2f2+l2)/qiqf
          call xdip(x,y,jmax,eln)
          ulc(0:jmax) = eln(0:jmax)
          fac = (four*iso-two)*two*pi*fcpi/(m_pis2*qiqf)
          do ll=-1,1
              vpis(ll) = fac*ulc(j+ll)*m_pic**2/three - vpis(ll)
          enddo
          vpisl0 = four*pi*fcpi/(m_pis2*three) *                             &
              (four*iso-two)*(y-x)**2/(one-y**2) - vpisl0
          do ll=-2,2
              vpit(ll) = -fac*ulc(j+ll) - vpit(ll)
          enddo
      end select

!     Other Yukawa's with multiples of average pion mass
      l2 = 64.0*m_pi2
      y = half*(qi2f2+l2)/qiqf
      do im=2,6
          x = half*(qi2f2+im**2*m_pi2)/qiqf
          call xdip(x,y,jmax,eln)
          ul(im,0:jmax) = eln(0:jmax)
      enddo

      fac = two*pi/qiqf

!     Potential for each partial wave separately
      select case(name)
      case('1S0')
          select case(type)
          case('PP','NN','pp','nn')
              vpot(1,1) = fac*( a(1,1)*ul(2,j)+a(1,2)*ul(3,j)+a(1,3)*ul(4,j) &
                               +a(1,4)*ul(5,j)+a(1,5)*ul(6,j) )
          case('NP','PN','np','pn')
              vpot(1,1) = fac*( b(1,1)*ul(3,j)+b(1,2)*ul(4,j)+b(1,3)*ul(5,j) &
                               +b(1,4)*ul(6,j) )
          end select
      case('1D2')
          vpot(1,1) = fac*(a(2,1)*ul(4,j)+a(2,2)*ul(5,j)+a(2,3)*ul(6,j) )
      case('1G4')
          vpot(1,1) = fac*a(2,4)*ul(3,j)
      case('3P0')
          vten = ul(3,jp)-half*s2psi*(tjpp*ul(3,j)+tj*ul(3,jpp))/tjp
          vpot(1,1) = fac*( a(3,1)*ul(3,jp)+a(3,2)*ul(5,jp)                  &
                           +a(2,5)*(-vten/9.0/m_pi2)*qi2f2/three )
      case('3F3')
          vpot(1,1) = fac*a(4,5)*ul(3,j)
      case('3P1')
          vten = ul(3,j)-half*s2psi*(tjp*ul(3,jm)+tjm*ul(3,jp))/tj
          vpot(1,1) = fac*( a(3,3)*ul(3,j)+a(3,4)*ul(5,j)                    &
                           +a(3,5)*(-vten/9.0/m_pi2)*qi2f2/three )
      case('3P3')
          vpot(1,1) = fac*a(4,5)*ul(3,j)
      case('1P1','1H5','1K7','1M9')
          vpot(1,1) = fac*( b(2,1)*ul(3,j)+b(2,2)*ul(4,j)+b(2,3)*ul(5,j)     &
                           +b(2,4)*ul(6,j) )
      case('1F3')
          vpot(1,1) = fac*(b(1,5)*ul(3,j)+b(2,5)*ul(5,j) )
      case('3D2')
          vten = ul(3,j)-half*s2psi*(tjp*ul(3,jm)+tjm*ul(3,jp))/tj
          vpot(1,1) = fac*( b(3,1)*ul(3,j)+b(3,2)*ul(5,j)                    &
                           +b(3,3)*(-vten/9.0/m_pi2)*qi2f2/three )
      case('3G4')
          vpot(1,1) = fac*b(3,4)*ul(3,j)
      case('1I6','1L8')     ! 1S0 PP parameters used
          vpot(1,1) = fac*( a(1,1)*ul(2,j)+a(1,2)*ul(3,j)+a(1,3)*ul(4,j)     &
                           +a(1,4)*ul(5,j)+a(1,5)*ul(6,j) )
      case default       ! spin=1, uncoupled j>5 and coupled
          select case(iso)
          case(1)
              do ll=-1,1
                  vc(ll) = fac*( a(4,1)*ul(3,j+ll)+a(4,2)*ul(4,j+ll)     &
                                +a(4,3)*ul(5,j+ll)+a(4,4)*ul(6,j+ll) )
              enddo
              do ll=-2,2
  vt(ll) = -fac*( a(5,1)*ul(4,j+ll)/16.0/m_pi2+a(5,2)*ul(6,j+ll)/36.0/m_pi2 )
              enddo
          case(0)
              do ll=-1,1
                  vc(ll) = fac*( b(4,1)*ul(2,j+ll)+b(4,2)*ul(3,j+ll)     &
                                +b(4,3)*ul(4,j+ll)+b(4,4)*ul(5,j+ll)     &
                                +b(4,5)*ul(6,j+ll) )
              enddo
              do ll=-2,2
                  vt(ll) = -fac*( b(3,5)*ul(4,j+ll)/16.0/m_pi2           &
                                 +b(5,5)*ul(6,j+ll)/36.0/m_pi2 )
              enddo
          end select

          select case(nchan)
          case(1)
              vl(:) = zero
              vpot(1,1) = vc(0) - qiqf*(vl(-1)-vl(1))/tj +               &
                          two/three*qi2f2*                               &
                          (vt(0)-half*s2psi*(tjp*vt(-1)+tjm*vt(1))/tj)
              vpot(1,2) = zero
              vpot(2,1) = zero
              vpot(2,2) = zero
          case(2)
              select case(name)
              case('3C1')
                  do ll=-2,2
                      vl(ll) = fac*( b(5,1)*ul(3,j+ll)/ 9.0/m_pi2        &
                                    +b(5,2)*ul(5,j+ll)/25.0/m_pi2 )
                  enddo
              case('3C2')
                  do ll=-2,2
                      vl(ll) = fac*( a(5,3)*ul(3,j+ll)/ 9.0/m_pi2        &
                                    +a(5,4)*ul(5,j+ll)/25.0/m_pi2 )
                  enddo
              case('3C3')
                  do ll=-2,2
                      vl(ll) = fac*( b(5,3)*ul(3,j+ll)/ 9.0/m_pi2        &
                                    +b(5,4)*ul(5,j+ll)/25.0/m_pi2 )
                  enddo
              case('3C4')
                  do ll=-2,2
                      vl(ll) = fac*a(5,5)*ul(3,j+ll)/ 9.0/m_pi2
                  enddo
              case default
                  vl(:) = zero
              end select

              vpot(1,1) = vc(-1) + qiqf*jm/tjm*(vl(-2)-vl(0)) +          &
                          two/three*qi2f2*jm/tj*                         &
                          (-vt(-1)+half*s2psi*(tjmm*vt(0)+tj*vt(-2))/tjm)
              vpot(1,2) = -two*qi2f2*tjj/tj*                             &
                          (-s2psi*vt(0)+cpsi2*vt(-1)+spsi2*vt(1))
              vpot(2,1) = -two*qi2f2*tjj/tj*                             &
                          (-s2psi*vt(0)+spsi2*vt(-1)+cpsi2*vt(1))
              vpot(2,2) = vc(1) - qiqf*jpp/tjp*(vl(0)-vl(2)) +           &
                          two/three*qi2f2*jpp/tj*                        &
                          (-vt(1)+half*s2psi*(tjpp*vt(0)+tj*vt(2))/tjp)
          end select
      end select

!     Add OPE potential
      select case(nchan)
      case(1)
          if (spin == 0) then
              vpot(1,1) = vpot(1,1) - three*vpis(0)
              if (l == 0) vpot(1,1) = vpot(1,1) - three*vpisl0
          elseif (l == j) then
              vpot(1,1) = vpot(1,1) + vpis(0) + two/three*qi2f2*             &
                          (vpit(0)-half*s2psi*(tjp*vpit(-1)+tjm*vpit(1))/tj)
          elseif (name == '3P0') then
              vpot(1,1) = vpot(1,1) + vpis(1) + two/three*qi2f2*jpp/tj*      &
                          (-vpit(1)+half*s2psi*(tjpp*vpit(0)+tj*vpit(2))/tjp)
          endif
      case(2)
          vpot(1,1) = vpot(1,1) + vpis(-1) + two/three*qi2f2*jm/tj*          &
                      (-vpit(-1)+half*s2psi*(tjmm*vpit(0)+tj*vpit(-2))/tjm)
          if (l == 0) vpot(1,1) = vpot(1,1) + vpisl0
          vpot(1,2) = vpot(1,2) - two*qi2f2*tjj/tj*                          &
                      (-s2psi*vpit(0)+cpsi2*vpit(-1)+spsi2*vpit(1))
          vpot(2,1) = vpot(2,1) - two*qi2f2*tjj/tj*                          &
                      (-s2psi*vpit(0)+spsi2*vpit(-1)+cpsi2*vpit(1))
          vpot(2,2) = vpot(2,2) + vpis(1) + two/three*qi2f2*jpp/tj*          &
                      (-vpit(1)+half*s2psi*(tjpp*vpit(0)+tj*vpit(2))/tjp)
      end select

      return

      end subroutine reid93

!     contains

!*******************************************************************************

      subroutine xdip(x,y,jmax,eln)

      implicit none

!     integer, parameter    :: wp = kind(1.d0)
      integer, parameter        :: wp = kind(1.d0)       ! double precision
      complex(wp), intent(in)  :: x,y
      logical, save         :: first = .true.
      integer, intent(in)   :: jmax
      real(wp)              :: p(12)
      complex(wp), intent(out) :: eln(0:12)
      integer, parameter    :: nogp = 64
      real(wp), save        :: zz(nogp),wz(nogp)
      complex(wp)              :: qx0,qy0,form
      integer               :: iz,l

      qx0 = 0.5*log((x+1.0)/(x-1.0))
      qy0 = 0.5*log((y+1.0)/(y-1.0))
      eln(0) = qx0-qy0+(y-x)/(1.0-y**2)
      if (jmax == 0) return
      if (first) then
          call gauss(-1.0_wp,1.0_wp,zz,wz,nogp)
          first = .false.
      endif

      eln(1:) = 0.0_wp
      do iz=1,nogp
          form = (y-x)/(y-zz(iz))
          form = form**2/(x-zz(iz))
!         Calculate Legendre polynomials first kind and integrate
          p(1) = zz(iz)
          p(2) = 1.5*zz(iz)*zz(iz)-0.5
          eln(1) = eln(1)+0.5*wz(iz)*form*p(1)
          eln(2) = eln(2)+0.5*wz(iz)*form*p(2)
          do l=3,jmax
              p(l) = ( (2*l-1)*zz(iz)*p(l-1) - (l-1)*p(l-2) ) / l
              eln(l) = eln(l)+0.5*wz(iz)*form*p(l)
          enddo
      enddo

      end subroutine xdip

!*******************************************************************************
! Calculate abscissas and weigths for Gauss-Legendre integration
! (n-points Gauss-Legendre quadrature).
!
! IN:    x1,x2     begin and end of integration range
!        n         number of integration points
! OUT:   x(1:n)    abscissas
!        w(1:n)    weights
!
! Literature: * Numerical Recipes in Fortran, paragraph 4.5 (1st edition)
!             * Abramowitz and Stegun, 25.4.29 and 25.4.30
!             * Mathews and Walker, paragraph 13.2
!             * J.J. de Swart, Lecture notes MTF (THEF-NIJM 84.04) par. F4.8
!
! Only for integrals with finite integration range. For improper integrals:
! perform a mapping first to an integral with a finite integration range.

      subroutine gauss(x1,x2,x,w,n)

      implicit none
!     integer, parameter      :: wp = kind(1.d0)
      integer, parameter        :: wp = kind(1.d0)       ! double precision
      real(wp), intent(in)    :: x1,x2
      integer, intent(in)     :: n
      real(wp), intent(out)   :: x(n),w(n)
      integer                 :: i,j,m
      real(wp)                :: xm,xl,z,z1,p1,p2,p3,pp
      real(wp), parameter     :: eps = 3e-14_wp
      real(wp), parameter     :: pi = 3.1415926535897932384626433832795028_wp
      real(wp), parameter     :: r0 = 0.0_wp
      real(wp), parameter     :: r025 = 0.25_wp
      real(wp), parameter     :: r05 = 0.5_wp
      real(wp), parameter     :: r1 = 1.0_wp
      real(wp), parameter     :: r2 = 2.0_wp

      m  = (n+1)/2
      xm = r05*(x2+x1)
      xl = r05*(x2-x1)
      do i=1,m
          z = cos(pi*(i-r025)/(n+r05))
          do
              p1 = r1
              p2 = r0
              do j=1,n
                  p3 = p2
                  p2 = p1
                  p1 = ((2*j-1)*z*p2-(j-1)*p3)/j
              enddo
              pp = n*(z*p1-p2)/(z*z-r1)
              z1 = z
              z  = z1-p1/pp
              if (abs(z-z1) <= eps) exit
          enddo
          x(i) = xm-xl*z
          x(n+1-i) = xm+xl*z
          w(i) = r2*xl/((r1-z*z)*pp*pp)
          w(n+1-i) = w(i)
      enddo

      end subroutine gauss

!*******************************************************************************

!     end subroutine reid93
