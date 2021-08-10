 
module constants  
implicit none 
!=======================================================================
!                             parameters
   double precision , parameter :: pi = 3.1415926536    
   double precision , parameter :: epsil=1.d-5, hc=197.33d0 !MeV.fm
   double precision , parameter :: a0  = - 18.7 !fm 
   double precision , parameter :: r1  = - 0.8819 !fm-1           
   double precision , parameter :: a1  = - 62.951 !fm3                   
   double precision , parameter :: gamma1 = -r1/2., gamma0 = 1./a0  !  
   double precision , parameter :: mn = 939.56533/hc ,mp = 938.27203/hc          !fm^-1
   double precision , parameter :: kR  =   dsqrt(2.d0/a1/r1) !                   fm^-1 
 
end module constants 


!=======================================================================
!                         Main Program                      
!=======================================================================
program three_body_bound_state
use constants   
! This code compute three-body bound states with halo EFT -> core-n-n       PRC 90, 044004 (2014) 

 integer      ::  ii,i,j,k ! to do loops 
 integer , parameter ::  nn=100, niter=60 , ntheta=16
 double precision :: q(nn), wq(nn),  H(niter),y(nn), dy(nn),x(ntheta), dx(ntheta)  
 double precision :: taun ,Xnc ,Xcn, Xnn ,tauc , knc, knn, kcn,FindDet
 double precision :: lambda,q0,c,d ,soma
 double precision :: matriz(nn,nn),mataux(nn,nn),matriz2(nn-1 ,nn-1 ),matriz2aux(nn-1,1) 
 double precision :: det(niter)  ,passo,Fcintp,Fnintp, Intncn, G0n,G0c,Fcnc,Fnnn,Fn1,Fc1,argcn,Fncn
 double precision :: mc,A,B3,qi,qj, psicmod(nn,nn), psinmod(nn,nn),F(nn),qq(nn),pp(nn)
 double precision :: chute, psic(nn,nn), psin(nn,nn),Fc(nn), Fn(nn),B(nn-1,1) ,H0,H0c,H0n
 double precision :: cmais,cmenos,nmais,nmenos,argnn,argnc, Fcmais,Fcmenos ,FnMenos,FnMais,somais,somenos
 double precision ::  pc,qc,pn,qn,  mnn,mnnc,mnc,mncn  
 character(len =4) :: lambc
    !  open(unit=55,file='H0xLambda.dat', status='unknown')      
      open(unit=51,file='input.in', status='old')
      open(unit=37,file='matriz.dat', status='unknown')
      open(unit=44,file='psic.dat', status='unknown')
      open(unit=45,file='psin.dat', status='unknown')

!=======================================================================
!                               inputs
 

      read(51,*) mc
      read(51,*) B3
      read(51,*) lambda
      read(51,*) chute ! de H0
close(51)
      mc =   mc*mn                         !fm^-1
      A  =  (mc/mn) 
      write(lambc,'(I4.4)') int(lambda) 
       
      lambda=lambda/hc
      B3 = B3/hc
      H(1)=chute
      print*, 'B3=', B3*hc, 'valor inicial para H0=', H(1) , ' lambda =', lambda*hc

 
      open(unit=33,file='Faddeev-components-Lamb='//lambc//'.dat', status='unknown') 
!=======================================================================
!                   Calculo da funcao H_0(Lambda)
        

       H(2)=H(1)*1.01


! do 2 kk=1,995                   ! varia lambda para calcalar a função H0  

 ! lambda=50./hc     + 10*kk/hc                      ! cutoff  fm^-1 

 

call gausslegend(-1.d0,1.d0,y,dy,nn)

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



         if (ii.gt.2) then  
    H(ii)=H(ii-1)-(H(ii-1)-H(ii-2))/( (det(ii-1))- (det(ii-2)))* (det(ii-1))    ! newton-rapshon para H0

         endif  
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
         
        knc =   2.*4.*pi*qj**2*Intncn(qi,qj, A, B3,q,wq,nn)*taun(qj,A, B3)*wq(j)          
 
        knn =   4.*pi*qj**2*Xnn(qi,qj, A, B3, H(ii),lambda)*taun(qj,A, B3)*wq(j)  
! 
!
!      !
            matriz(i,j)       = knn + knc
   if (i==j)matriz(i,j)       = matriz(i,j) -1.
!  
!        
5        continue           !laço de q i    linha
4      continue             !laço de q_j    coluna
!
! 
           
          mataux = matriz   
          det(ii)= FindDet(mataux, nn)
           
!         
!    

        if(ii.eq.niter) then
        print*, 'nao convergiu para ', niter, 'interacoes'
        stop
        endif
 
        print*,  ii,  lambda*hc,  'det=', det(ii), 'lambda =', lambda*hc, 'H0=',H(ii)  
        if(ii.gt.2.and.dabs((H(ii-1)-H(ii))/H(ii-1)).lt.epsil)then
        print*, 'det=', det(ii), 'Convergiu para H=',H(ii) , ' lambda =', lambda*hc  
        write(55,*) lambda*hc, H(ii)
        H0=H(ii)
        exit        
        endif
!
! 
!=======================================================================         
!         
1    continue              ! determinante ii 
  !2    continue       !lambda laço kk 

!======================================================================= 

!                        Calculo das funcoes espectadoras

!======================================================================= 
  
          
          write(37,*) nn
           do i=1,nn
           write(37,*) i, q(i)
           enddo
           do i=1,nn  
             do j=1,nn  
               write(37,*)  matriz(i ,j)        ! para resolver no metodo gauss-jordan  gj.f90
             enddo 
         enddo 
           
           
          
                
             
          do i=2,nn 
           B(i-1,1)=  -matriz(i,1)                   !matrix x F = -B
          enddo  
          
         do i=2, nn  
             do j=2, nn  
               matriz2(i-1,j-1) = matriz(i,j)   
             enddo 
         enddo 
          
          call GaussInv(matriz2, nn-1) ! inverte a matriz2  
          !print*, FindDet(matriz2, 2*nn)
             
         do i=2,nn 
           F(i)=0.
           do j=2,nn 
             F(i) = F(i) + matriz2(i-1,j-1)*B(j-1,1)
           enddo
         enddo
         
         !F(1)=1.
        do i=1,nn-1
           !write(10,*)q(i)*hc,F(i)  
        enddo
        
        ! multiplicando diretamente B x M^-1
        
          matriz2aux = matmul(matriz2,B)           ! F = - B x M^-1
          
          
          !------------------------ teste inverso se F=M^-1B, entao MF-B = 0      
          !  do i=1,2*nn  
          !   do j=1,2*nn  
          !    matriz2(i,j) = matriz(i,j)   
          !    if(j.eq.1)  matriz2(i ,j)=1.d-9 
          !   enddo 
          ! enddo
         
          ! do i=1,2*nn  
          !    print*, matmul(matriz2(i,:),matriz2aux)- first(i,1)
          ! enddo
          
           
         do  i=1, nn-1 
            Fn(i) = matriz2aux(i,1) 
         enddo   
!----------------------------------------------- Fc
            
           do i=1,nn-1 
             soma=0.
             do j=1,nn  
               soma = soma + 8.*pi*q(j)**2*wq(j) * Xcn(q(i),q(j),A,B3)*taun(q(j),A,B3)*Fn(j) 
             enddo
             Fc(i) =soma        ! Fc em q^-1 -> 1/fm 
             write(11,*)    q(i)*hc , Fc(i) 
             write(10,*)    q(i)*hc , Fn(i) 
            write(33,*)     q(i)*hc , Fn(i) *hc  /Fc(1), Fc(i) /Fc(1)
           enddo 
           
          print*,Fc(1) 
           
   
!=======================================================================

!                        Calculo das wave functions   (passo fixo)

!=======================================================================  
              Fn=Fn!*(hc)/Fc(1)  
              Fc=Fc!/Fc(1)
             call gausslegend(-1.d0,1.d0,x,dx,ntheta)
             mnn  = 1./2.    *mn
             mnnc = 2.*A/(A+2.) *mn
             mnc  = A/(A+1.) *mn
             mncn = (A+1)/(A+2.)  *mn
  

               
               passo = q(nn)/nn
               qc=1.d-5
               qn=qc 
           do i=1, nn                      !q
               pc=1.d-5
               pn=pc
            do j=1, nn                    !p
             qq(i)=qc
             pp(j)=pc
           !----------------------------- psic     
             somenos =0.
             somais =0.
            do k=1,ntheta 
              cmais=dsqrt( pc**2 + qc**2/4.  +   pc*qc *x(k))
              cmenos=dsqrt(pc**2 + qc**2/4.  -   pc*qc *x(k))
             call interpolation(nn,q,Fn,1,1, cmais,ixx,FcMais)
             call interpolation(nn,q,Fn,1,1, cmenos,ixx,FcMenos) 
             somenos = somenos + FcMenos*dx(k)
             somais  = somais + FcMais*dx(k)
            enddo  
              H0c=  pc**2/2./mnn  +  qc**2 /2./mnnc  

             call interpolation(nn,q,Fc,1,1,qc,ixx,Fcintp) 
             
               psic(i,j )=(Fcintp + somenos + somais) / ( B3 + H0c )
               
       !----------------------------- psin   
                     H0n=  pn**2/2./mnc  +  qn**2 /2./mncn  
                     
             somenos =0.
             somais =0.
               do k=1,ntheta 
               nmais=dsqrt( pn**2 + qn**2 /(A+1.)**2  +  2.*pn*qn /(A+1.)*y(k))
               nmenos=dsqrt(pn**2 + qn**2 *A**2/(A+1.)**2 - 2.*pn*qn*A/(A+1.)*y(k))
  
             call interpolation(nn,q,Fc,1,1, nmenos,ixx,FnMenos) 
             call interpolation(nn,q,Fn,1,1, nmais,ixx,FnMais)
             somenos = somenos + FnMenos*dx(k)
             somais  = somais + FnMais*dx(k)
             
             
               enddo
               
               
             call interpolation(nn,q,Fn,1,1,qn,ixx,Fnintp) 
               psin(i,j )= (Fnintp + somenos + somais ) / (  B3 + H0n )        !psi(q,p,)
             
            
               
            
              
             pc=pc+passo
             pn=pn+passo
            enddo  
             qc=qc+passo
             qn=qn+passo
            enddo 
            
            
            do i=1,40
             do j=1,40     
                write(44,*) qq(i)*hc, pp(j)*hc,   (psic(i,j))            !psi(q,p)
                write(45,*) qq(i)*hc, pp(j)*hc,   (psin(i,j)) 
            enddo
            enddo !
!
close(33) 
end program 

!=======================================================================
!=======================================================================
!                              functions
!-----------------------------------------------------------------------

      double precision function Intncn(q1,q2, A, B3,q,wq,n)
      use constants
       integer n
       double precision :: soma, tauc ,Xnc ,Xcn,q1,q2,q(n),wq(n),A,B3
       
       soma=0.
       do i=1,n  
           soma = soma + q(i)**2*wq(i)*Xnc(q1,q(i),A,B3)*tauc(q(i),A,B3)*Xcn(q(i),q2,A,B3)
       enddo
       Intncn = 4.*pi*soma

      return
      end
!-----------------------------------------------------------------------
      double precision function taun(q,A,EB)
      use constants         
        double precision ::  q , kn, A,EB 
         
          Kn   = dsqrt(2.*A/(A+1.)*(mn*EB + (A+2.)/2./(A+1.)*q**2))
          taun = -1./(4.*pi**2*gamma1 *mn)*(A+1.)/A /(Kn**2  +  kR**2)

        return  
        end
        
!-----------------------------------------------------------------------
      double precision function tauc(q ,A,EB)
      use constants        
        double precision ::  q  , kc, EB ,A
        
          Kc   = dsqrt(mn*EB  +  (A+2.)/4./A*q**2)
          tauc =  1./(2.*pi**2*mn) / (gamma0-Kc) 
          
        return  
        end        

!-----------------------------------------------------------------------
      double precision function Xcn(q1,q2,A,EB)
      use constants 
        double precision :: q1,q2, Q0cn, Q1cn, EB, A ,zcn 
        
          zcn = -1/q1/q2*(mn*EB + (A+1.)/(2.*A)*(q1**2) + q2**2  )  
         Q0cn = 0.5*log((zcn+1.)/(zcn-1.))
         Q1cn = 0.5*zcn*log((zcn+1.)/(zcn-1.)) - 1.  

        Xcn   = - dsqrt(2.d0)*mn *(A/(A+1.)/q1 * Q0cn  + 1./q2 *Q1cn)        
      return
      end
! ----------------------------------------------------------------------
      double precision function Xnc(q1,q2,A,EB)
      use constants 
        double precision :: q1,q2, Q0nc, Q1nc, EB, A,znc 
        znc  = -1/q1/q2*(mn*EB + q1**2 + (A+1.)/(2.*A)*(q2**2))      
        Q0nc = 0.5*log((znc+1.)/(znc-1.))
        Q1nc = 0.5*znc*log((znc+1.)/(znc-1.)) - 1. 
        Xnc  = - dsqrt(2.d0)*mn *(A/(A+1.)/q2 * Q0nc  + 1./q1 *Q1nc)  
      return
      end
! ----------------------------------------------------------------------
      double precision function Xnn(q1,q2,A,EB,H,lambda)
      use constants 
        double precision :: q1,q2, Q0nn, Q1nn, Q2nn,EB, A, znn 
       double precision :: H, lambda
         znn = -A/q1/q2*(mn*EB + (A+1.)/(2.*A)*(q1**2+q2**2))
         Q0nn= 0.5*log((znn+1.)/(znn-1.))
         Q1nn= 0.5*znn*log((znn+1.)/(znn-1.))-1.
         Q2nn= 0.5*(-0.5 + 3./2.*znn**2) *  log((znn+1)/(znn-1))-3./2.*znn 
         
        Xnn =  A *mn* ( (A**2+2.*A+3.)/(A+1.)**2*Q0nn + 2./(A+1.)*(q1**2+q2**2)/q1/q2*Q1nn + Q2nn)
        Xnn = Xnn -  q1*q2*H*mn  /lambda**2  
       
      return
      end      
  
!=======================================================================
!                              Subroutines
 
! --------------------------------------------------------------------
SUBROUTINE GaussInv(a,n)       ! Invert matrix by Gauss method
! --------------------------------------------------------------------
IMPLICIT NONE
INTEGER :: n
REAL(8) :: a(n,n)
! - - - Local Variables - - -
REAL(8) :: b(n,n), c, d, temp(n)
INTEGER :: i, j, k, m, imax(1), ipvt(n)
! - - - - - - - - - - - - - -
b = a
ipvt = (/ (i, i = 1, n) /)
DO k = 1,n
   imax = MAXLOC(ABS(b(k:n,k)))
   m = k-1+imax(1)
   IF (m /= k) THEN
      ipvt( (/m,k/) ) = ipvt( (/k,m/) )
      b((/m,k/),:) = b((/k,m/),:)
   END IF
   d = 1/b(k,k)
   temp = b(:,k)
   DO j = 1, n
      c = b(k,j)*d
      b(:,j) = b(:,j)-temp*c
      b(k,j) = c
   END DO
   b(:,k) = temp*(-d)
   b(k,k) = d
END DO
a(:,ipvt) = b
END SUBROUTINE GaussInv



!Function to find the determinant of a square matrix
!Author : Louisda16th a.k.a Ashwith J. Rego
!Description: The subroutine is based on two key points:
!1] A determinant is unaltered when row operations are performed: Hence, using this principle,
!row operations (column operations would work as well) are used
!to convert the matrix into upper traingular form
!2]The determinant of a triangular matrix is obtained by finding the product of the diagonal elements
!
  double precision function FindDet(matrix, n)
    IMPLICIT NONE
    INTEGER , INTENT(IN) :: n
     double precision , DIMENSION(n,n) :: matrix
     double precision  :: m, temp
    INTEGER  :: i, j, k, l
    LOGICAL :: DetExists = .TRUE.
    l = 1
    !Convert to upper triangular form
    DO k = 1, n-1
        IF (matrix(k,k) == 0) THEN
            DetExists = .FALSE.
            DO i = k+1, n
                IF (matrix(i,k) /= 0) THEN
                    DO j = 1, n
                        temp = matrix(i,j)
                        matrix(i,j)= matrix(k,j)
                        matrix(k,j) = temp
                    END DO
                    DetExists = .TRUE.
                    l=-l
                    EXIT
                ENDIF
            END DO
            IF (DetExists .EQV. .FALSE.) THEN
                FindDet = 0
                return
            END IF
        ENDIF
        DO j = k+1, n
            m = matrix(j,k)/matrix(k,k)
            DO i = k+1, n
                matrix(j,i) = matrix(j,i) - m*matrix(k,i)
            END DO
        END DO
    END DO
   
    !Calculate determinant by finding product of diagonal elements
    FindDet =   l 
    DO i = 1, n
        FindDet = FindDet * matrix(i,i)
    END DO
   
END FUNCTION FindDet



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

      subroutine gausslegend(x1,x2,x,w,n)

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

      end subroutine gausslegend

!*******************************************************************************


!******************  SUBROTINA de INTERPOLAÇÃO  ************************       
!23456789a123456789b123456789c123456789d123456789e123456789f123456789g12
!
! Copyright D.J. Jeffery, 2006jan01.
!
! interpolation.f does interpolations and extrapolations.
!
! Program interpolationt.f seems to confirm adequately all three modes
! in both spacing options and for both interpolation and extrapolation.
!
!   References:
!
! \bibitem[Arfken(1970)]{arfken1970} Arfken, G. 1970,
!        Mathematical Methods for Physicists
!        (New York:  Academic Press)
!        % The old, but still very good Arfken.
!        % See p. xiii-xiv of my own notes.
!
!     \bibitem[Metcalf et al.(2004)]{metcalf} Metcalf, M., Reid, J.,
!          \&~Cohen, M. 2004,
!          Fortran 95/2003 Explained
!          (Oxford:  Oxford University Press) (MRC)
!          % Pretty good reference book.
!
!23456789a123456789b123456789c123456789d123456789e123456789f123456789g12
!
      subroutine interpolation(nn,x,y,imode,idx,xx,ixx,yy)
      !use numerical
      implicit none                        
!
      real (kind=8) :: a
      real (kind=8) :: b
      integer :: i
      integer, intent(in) :: idx    ! 0 for equally spaced points;  1 for unequally spaced points.
      integer, intent(in) :: imode  ! 1 for linear, 2 for semi-log, 3 for log-log.
      integer, intent(out) :: ixx   ! Index of array space that contains xx.
      integer, intent(in) :: nn     ! Number of points in the array of data.
      real (kind=8), intent(in) :: x(nn) 
      real (kind=8), intent(in) :: xx
      real (kind=8), intent(in) :: y(nn)
      real (kind=8), intent(out) :: yy       ! The output interpolated value. 
!
      if(idx .eq. 0) then
          i=max( ceiling( (xx-x(1))/(x(2)-x(1)) ) + 1, 2)  ! +1 since there is no zero element. 
!                                                          ! See MRC-162 for ceiling.
!                                                          ! max is for extrapolation case. 
        else
          do i=2,nn                  ! Starts from 2 to handle extrapolation case.
           if(x(i) .ge. xx) exit     ! See MRC-61 for do-loop variable final value.
          end do
      end if
      ixx=min(i,nn)  ! min is for extrapolation case for too-large xx.
!
      if(imode .eq. 1) then
          a=(y(ixx)-y(ixx-1))/(x(ixx)-x(ixx-1))    ! linear
          b=y(ixx)-a*x(ixx)
          yy=a*xx+b
        else if(imode .eq. 2) then
          a=log(y(ixx)/y(ixx-1))/(x(ixx)-x(ixx-1)) ! semi-logarithmic 
          b=log(y(ixx))-a*x(ixx)
          yy=exp( a*xx+b )
        else
          a=log(y(ixx)/y(ixx-1))/log(x(ixx)/x(ixx-1)) ! logarithmic 
          b=log(y(ixx))-a*log(x(ixx))
          yy=exp( a*log(xx)+b )
      end if
!
      end subroutine interpolation 
 
