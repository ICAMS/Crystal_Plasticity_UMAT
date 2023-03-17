cccccccccccccccccccccccccccccccccccccccccccc
      subroutine caleulang(Mx,vv,ising)     
cccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer irdummy,ising,ising_1
      real(8) det
      real(8) vv(3)
      real(8) Mx(3,3),rm(3,3),um(3,3)

cccccccccccccccccccccccccccccccccccccccccccc
      ising_1=0
!      call PDECOMPOSITION(Mx,Um,Rm,ising_1)
      call polar_decomp(Mx,Um,Rm,ising_1)
      if(ising_1/=0)then
         write(6,*) 'pdecompsition of Mx is failure'
         call pm(Mx,3,3)
         ising=1
         return
      endif
      call icams_Q2Eang(transpose(Rm),vv(1),vv(2),vv(3))
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine icams_misori(v1,v2,ang)
ccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      real(8) v1(3),v2(3),ang,x1,x2
      real(8) QM1(3,3),QM2(3,3),dQM(3,3)
      real(8) pi2
      pi2 = datan(1.0d0)*2.0d0

      call icams_Eang2Q(v1(1),v1(2),v1(3),QM1)
      call icams_Eang2Q(v2(1),v2(2),v2(3),QM2)
      dQM=matmul(QM2,transpose(QM1))
      x1=dQM(1,1)+dQM(2,2)+dQM(3,3)
      x2=(x1-1.d0)/2
      if(dabs(x2)>1.d0)x2=1.d0*sign(1.d0,x2)
      ang=dabs(pi2-dasin(x2))
      return
      end

C****************************************************************
      subroutine pdecomposition(Mx,UMx,RMx,ising)
C****************************************************************
      implicit none
      integer ising
      real(8) Mx(3,3),ce(3,3),RMx(3,3),UMx(3,3),IUMx(3,3)
      real(8) eb1(3,3),eb2(3,3),eb3(3,3)
      real(8) ev1,ev2,ev3,det
      ising=0
      ce=matmul(transpose(Mx),Mx)
      call spectral(ce,ev1,ev2,ev3,eb1,eb2,eb3,ising)
      if(ev1<=0 .or. ev2<=0 .or. ev3<=0 .or. ising/=0)then
         write(6,*) 'eigen value of ce <0'
         print*, ev1,ev2,ev3
         ising=1
         return
      endif
      UMx=dsqrt(ev1)*eb1+dsqrt(ev2)*eb2+dsqrt(ev3)*eb3
      IUMx=1/dsqrt(ev1)*eb1+1/dsqrt(ev2)*eb2+1/dsqrt(ev3)*eb3
      RMx=matmul(Mx,IUMx)
      return 
      end

      subroutine icams_determ(a,det)
c********************************************************************
c     This routine calculates the determinant of a
c********************************************************************
      implicit none
      real(8) a(3,3),v1,v2,v3,det
      v1= a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2))
      v2= a(1,2)*(a(2,1)*a(3,3)-a(2,3)*a(3,1))
      v3= a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1))
      det= v1-v2+v3
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine icams_Eang2Q(p1,p,p2,QM)
ccccccccccccccccccccccccccccccccccccccccccccccccccc
c     rotate from X[100],Y[010],Z[001] to v1,v2,v3
ccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      real(8) QM(3,3)
      real(8) p1,p,p2
      real(8) c1,c,c2,s1,s,s2
      c1=dcos(p1)
      s1=dsin(p1)
      c =dcos(p )
      s =dsin(p )
      s2=dsin(p2)
      c2=dcos(p2)
      QM(1,1)=+c1*c2-s1*s2*c
      QM(1,2)=+s1*c2+c1*s2*c
      QM(1,3)=+s2*s
      QM(2,1)=-c1*s2-s1*c2*c
      QM(2,2)=-s1*s2+c1*c2*c
      QM(2,3)=+c2*s
      QM(3,1)=+s1*s
      QM(3,2)=-c1*s
      QM(3,3)=+c
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine icams_Q2Eang(QM,phi1,PHI,phi2)
cccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      real(8) QM(3,3)
      real(8) phi1,PHI,phi2
      real(8) sqhkl,squvw,sqhk,val
      real(8) Tol
      Tol=1.d-15


c---------------------------------------------
c
c             v1   v2    v3
c
c           | u   v2_1    h  |
c     QM =  | v   v2_2    k  |  with v2 = v3 x v1   
c           | w   v2_3    l  | 
c
c             100   010   001
c          X|  u   v2_1    h  |
c     QM = Y|  v   v2_2    k  |    
c          Z|  w   v2_3    l  | 
c
c---------------------------------------------


c---------------------------------------------
c
c    if the roation tensor is defined as following
c    QM_ij:= (dx/dX)_ij = dx_i/dX_j
c
c               X     Y     X
c          v1|  u   v2_1    h  |
c     QM'= v2|  v   v2_2    k  |    
c          v3|  w   v2_3    l  | 
c    the ratation must be take 
c
c    QM=transpose(QM')
c
c---------------------------------------------
      squvw=dsqrt(QM(1,1)**2+QM(2,1)**2+QM(3,1)**2)
      sqhkl=dsqrt(QM(1,3)**2+QM(2,3)**2+QM(3,3)**2)
      sqhk =dsqrt(QM(1,3)**2+QM(2,3)**2           )

      val=QM(3,3)/sqhkl
      if(dabs(val)>1.d0)val=1.d0*sign(1.d0,val)
      PHI=dacos(val)

      if(PHI < TOL) then
         phi2=0.0
         val=QM(1,1)/squvw
         if(QM(2,1) <= 0.d0) then
            if(dabs(val)>1.d0)val=1.d0*sign(1.d0,val)
            phi1=dacos(val)
         else
            if(dabs(val)>1.d0)val=1.d0*sign(1.d0,val)
c           phi1=2*pi-dacos(val)
            phi1=-dacos(val)
         endif
      else
c        val=QM(2,3)/sqhk
         val=QM(2,3)/dsin(PHI)
         if(QM(1,3) >= 0.d0) then
            if(dabs(val)>1.d0)val=1.d0*sign(1.d0,val)
            phi2=dacos(val)
         else
            if(dabs(val)>1.d0)val=1.d0*sign(1.d0,val)
c           phi2=2*pi-dacos(val)
            phi2=-dacos(val)
         endif
         val=-QM(3,2)/dsin(PHI)
         if(QM(3,1) >= 0.d0) then
            if(dabs(val)>1.d0)val=1.d0*sign(1.d0,val)
            phi1=dacos(val)
         else
            if(dabs(val)>1.d0)val=1.d0*sign(1.d0,val)
c           phi1=2*pi-dacos(val)
            phi1=-dacos(val)
         endif
      endif

      return
      end

cccccccccccccccccccccccccccccccccccccccccccc
      subroutine icams_angax2QM(ang,u,v,w,QM)
cccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      real(8) QM(3,3)
      real(8) s,c,u2,v2,w2,ang,u,v,w,x1

      x1=dsqrt(u**2+v**2+w**2)
      u2=u/x1
      v2=v/x1
      w2=w/x1
      s=dsin(ang)
      c=dcos(ang)

      QM(1,1)=(1-u2**2)*c+u2**2
      QM(2,2)=(1-v2**2)*c+v2**2
      QM(3,3)=(1-w2**2)*c+w2**2

      QM(1,2)=u2*v2*(1-c)+w2*s
      QM(2,1)=u2*v2*(1-c)-w2*s

      QM(1,3)=u2*w2*(1-c)-v2*s
      QM(3,1)=u2*w2*(1-c)+v2*s

      QM(2,3)=v2*w2*(1-c)+u2*s
      QM(3,2)=v2*w2*(1-c)-u2*s


      return
      end

c
c
C********************************************************************** 
      SUBROUTINE HI(M,HI1M,HI2M,HI3M)
C**** HAUPTINVARIANTEN HI1M, HI2M, HI3M DER 3X3 MATRIX M
      IMPLICIT NONE
      real(8) M(3,3),HI1M,HI2M,HI3M 
      HI1M = M(1,1)+M(2,2)+M(3,3)
      HI2M =(M(1,1)+M(2,2)+M(3,3))**2/2.d0
     &      -M(1,1)**2/2.d0
     &      -M(2,2)**2/2.d0
     &      -M(3,3)**2/2.d0
     &      -M(1,2)*M(2,1)
     &      -M(1,3)*M(3,1)
     &      -M(2,3)*M(3,2)
      HI3M =+M(1,1)*M(2,2)*M(3,3)
     &      +M(2,1)*M(1,3)*M(3,2)
     &      +M(3,1)*M(1,2)*M(2,3)
     &      -M(1,1)*M(2,3)*M(3,2)
     &      -M(2,1)*M(1,2)*M(3,3)
     &      -M(3,1)*M(1,3)*M(2,2)
      RETURN  
      END


C**********************************************************************
      SUBROUTINE NORM(M,ZM,SM,NO)
C**********************************************************************
      IMPLICIT NONE
      INTEGER I,J,ZM,SM
      DOUBLE PRECISION M(ZM,SM),NO  
      NO=0.d0
      DO I=1,ZM
      DO J=1,SM 
         NO=NO+M(I,J)**2.d0 
      END DO 
      END DO
      NO=dsqrt(NO)
      RETURN 
      END  


C**********************************************************************
      subroutine wm(M,ni,nj)
      implicit none
      integer ni,nj,i,j
      real(8) M(ni,nj)
      write(6,*) 
      do i=1,ni
         write(6,100) (m(i,j),j=1,nj)          
      enddo
      write(6,*) 
100   format(48e20.8)       
      return
      end


C**********************************************************************
      subroutine pm(M,ni,nj)
      implicit none
      integer ni,nj,i,j
      real(8) M(ni,nj)
      write(*,*) 
      do i=1,ni
         write(*,100) (m(i,j),j=1,nj)          
      enddo
      write(*,*) 
c100   format(48e20.8)        
100   format(48e12.4)       
c100   format(48f6.3)       
      return
      end

C**********************************************************************
      subroutine pv(M,ni)
      implicit none
      integer ni,i
      real(8) M(ni)
      write(*,100) (M(i),i=1,ni)          
c      write(6,*) 
c100   format(48e20.8)       
100   format(48e12.4)       
      return
      end


C******************************************** 
      subroutine gaussj(Amat,n,Bmat,ising)
C******************************************** 
      implicit none
      integer, intent(in) :: n
      integer, intent(out) :: ising
      real(8), intent(in) :: Amat(n,n)
      real(8), intent(out) ::  Bmat(n,n)
      integer   i,icol,irow,j,k,l,ll,i1,i2
      integer   indxc(n),indxr(n),ipiv(n)
      real(8)   vx(n)
      real(8)   big,dum,pivinv,x1,x2,x3
c-------------------------------------------
c      write(6,*) 'coming to jd'
c      call pm(A,3,3)
c      call flush(6)

      ising=0
      Bmat=Amat
      ipiv=0
c-----loop for the pivot procedures, from 1 to n
      do i=1,n
         big=0.d0
         do j=1,n
         do k=1,n
            if(ipiv(j).ne.1 .and. ipiv(k)==0)then
               if(dabs(Bmat(j,k)).ge.big)then
                  big=dabs(Bmat(j,k))
                  irow=j
                  icol=k
               endif
            endif
            if(ipiv(j).ne.1 .and. ipiv(k).gt.1)then
c              print*,'sigular matrix in gauss_jordan'
c              write(6,*)'sigular matrix in gauss_jordan'
c              call flush(6)
               ising=1
               return
            endif
         enddo
         enddo
c--------check whether the pivot element is zero or not
         if(big.le.1.0d-10)then
c	      print*,'sigular matrix in gauss_jordan'
c	      print*,'indices is:', irow,icol
c           write(6,*)'sigular matrix in gauss_jordan'
c           write(6,*)'indices is:', irow,icol
c           call flush(6)
            ising=2
            return
         endif
c-------------------------------------------------------
c        if one component is selected as pivot element
c        the second indice is important, so it is marked
c        from 0 to 1 in ipiv(:) array
c        after convert, it is the row number
c-------------------------------------------------------
         ipiv(icol)=ipiv(icol)+1
c--------record the row and collum number for ith pivot element
         indxr(i)=irow
         indxc(i)=icol
c--------change pivot element to diagonal position
         if(irow.ne.icol)then
            vx=Bmat(irow,:)
            Bmat(irow,:)=Bmat(icol,:)
            Bmat(icol,:)=vx
         endif
c--------eliminate the elements besides Bmat(icol,icol)
         pivinv=1.d0/Bmat(icol,icol)
         Bmat(icol,icol)=1.d0
         Bmat(icol,:)=Bmat(icol,:)*pivinv 
         do i2=1,n
            if(i2.ne.icol)then
               dum=Bmat(i2,icol)
               Bmat(i2,icol)=0.d0
               Bmat(i2,:)=Bmat(i2,:)-Bmat(icol,:)*dum
            endif
         enddo
      enddo

c-----after maximum pivot strategy elimination
c-----rearrage the left matrix
      do l=n,1,-1
         if(indxr(l).ne.indxc(l))then
            vx=Bmat(:,indxr(l))
            Bmat(:,indxr(l))=Bmat(:,indxc(l))
            Bmat(:,indxc(l))=vx
         endif
      enddo

      RETURN
      END


C***************************************************************C
C      SUBROUTINE SPOLAR (F,U,R)				        C
C							  	        C
C *performs the polar decomposition of the	               C
C   tensor F=RU using a Cayley-Hamilton theorem.		        C
C *Then, U=R^T F is found					        C
C***************************************************************C
      SUBROUTINE polar_decomp(F,U,R,ising)
      implicit none
      integer i,j,k,iflag1,iflag2,ising
      real(8) F(3,3),U(3,3),R(3,3),RT(3,3),C(3,3),CS(3,3),UINV(3,3)
      real(8) x1,x2,x3,f1,f2,c1,c2,c3,u1,u2,u3,b,b1,b2,b3,b4
C
      C=matmul(transpose(F),F)
      CS=matmul(C,C)

c-----1st, 2st, 3st invariant for tensor C
      C1=C(1,1)+C(2,2)+C(3,3)
      C2=(C1**2.d0-(CS(1,1)+CS(2,2)+CS(3,3)))/2
      C3=+C(1,1)*(C(2,2)*C(3,3)-C(2,3)*C(3,2))
     1	 -C(1,2)*(C(2,1)*C(3,3)-C(2,3)*C(3,1))
     2	 +C(1,3)*(C(2,1)*C(3,2)-C(2,2)*C(3,1))

c-----3st invariant for tensor U
      U3=dsqrt(C3)

      X1= 2.0**5.0 /27.0 * (2.0*C1**3.0-9.0*C1*C2+27.0*C3)
      X2= 2.**10.0 /27.0 * (4.0*C2**3.0 - C1**2.0*C2**2.0 +
     1	  4.0*C1**3.0*C3 - 18.0 * C1*C2*C3 + 27.0 * C3**2.0)

      IF(X2<0)X2=0
      IFLAG1=0
      IF(F1<0)IFLAG1=1
      F2=X1-dsqrt(X2)
      IFLAG2=0
      IF(F2<0)IFLAG2=1

      IF(IFLAG1==1) F1=-F1
      IF(IFLAG2==1) F2=-F2

      X3= -2.0/3.0*C1 + F1**(1.0/3.0) + F2**(1.0/3.0)
      IF(IFLAG1==1) X3= -2.0/3.0*C1 - F1**(1.0/3.0) + F2**(1.0/3.0)
      IF(IFLAG2==1) X3= -2.0/3.0*C1 + F1**(1.0/3.0) - F2**(1.0/3.0)

c-----1st, 2st invariant for tensor U
      B=-2.0*C1
      if(X3==B)then
         U1=dsqrt(C1+2.0*dsqrt(C2))
      else

         x1=dsqrt(2.0*C1+X3)
         if(x1==0)then
            ising=1
            return
         endif      
         U1= 0.5 * ( x1 + dsqrt(2.0*C1 - X3 + 16.0*dsqrt(C3)/x1) )
      endif
      U2=dsqrt(C2+2.0*U3*U1)

      B1= U3**2.0 * (U3+U1*C1) + U1**2.0 * (U1*C3+U3*C2)
      if(B1==0)then
         ising=1
         return
      endif      

      B2= U1*(U1*U2-U3) / B1
      B3=-(U1*U2-U3) * (U3+U1*C1) / B1
      B4= (U2*U3*(U3+U1*C1) + U1**2.0 * (U2*C2+C3))/B1

      UINV=B2*CS + B3*C
      Uinv(1,1)=Uinv(1,1)+B4
      Uinv(2,2)=Uinv(2,2)+B4
      Uinv(3,3)=Uinv(3,3)+B4

      R=matmul(F,UINV)
      U=matmul(transpose(R),F)

      RETURN
      END


C**********************************************************************
      SUBROUTINE spectral(M,EW1,EW2,EW3,EB1,EB2,EB3,ising)
C**********************************************************************
      implicit none
      integer ising
      real(8) M(3,3),EB1(3,3),EB2(3,3),EB3(3,3),EW1,EW2,EW3,
     &        HI1M,HI2M,HI3M,TOL,R,S,T,P,Q,RHO,PHI,Y1,Y2,Y3,D1,D2,D3,
     &        E(3,3),M1(3,3),M2(3,3),M3(3,3)
      real(8) x1,x2,x3,PI
      ising=0
      TOL=1.d-15
      PI=dacos(-1.d0)
      CALL HI(M,HI1M,HI2M,HI3M)
      R=-HI1M
      S= HI2M
      T=-HI3M
      P=S-R**2.d0/3.d0
      Q=2.d0/27.d0*R**3.d0-R*S/3.d0+T
      RHO=dsqrt(-3.d0*P**3.d0)/9.d0
      if(dabs(RHO)<=Tol)then
         ising=1
         return         
      endif
      !======================================
      !  modification for overcome overflow
      !======================================
      if(dabs(RHO)<=Tol) RHO=sign(1.d0,RHO)*TOL
      x1=-Q/RHO/2.d0
      if(dabs(x1)>1.d0) x1=1.d0*sign(1.d0,x1)
      PHI=dacos(x1)
      Y1=2*RHO**(1.d0/3.d0)*dcos(PHI/3.d0)
      Y2=2*RHO**(1.d0/3.d0)*dcos(PHI/3.d0+2.d0/3.d0*PI)
      Y3=2*RHO**(1.d0/3.d0)*dcos(PHI/3.d0+4.d0/3.d0*PI)
      EW1=Y1-R/3.d0
      EW2=Y2-R/3.d0
      EW3=Y3-R/3.d0
      E(1,:)=[1,0,0]
      E(2,:)=[0,1,0]
      E(3,:)=[0,0,1]
      EB1=0.d0
      EB2=0.d0
      EB3=0.d0
      IF(dabs(ew1-ew2)<tol .and. dabs(ew1-ew3)<tol)THEN
         EB1=E
      elseif(dabs(ew2-ew3)<TOL)THEN
         EB1=MATMUL(M-EW2*E,M-EW3*E)/((ew1-ew2)*(ew1-ew3))
         EB2=E-EB1
      elseif(dabs(ew1-ew3)<TOL)THEN
         EB2=MATMUL(M-EW1*E,M-EW3*E)/((ew2-ew1)*(ew2-ew3))
         EB1=E-EB2
      elseif(dabs(ew1-ew2)<TOL)THEN
         EB3=MATMUL(M-EW1*E,M-EW2*E)/((ew3-ew1)*(ew3-ew2))
         EB1=E-EB3
      else
         EB1=MATMUL(M-EW2*E,M-EW3*E)/((ew1-ew2)*(ew1-ew3))
         EB2=MATMUL(M-EW1*E,M-EW3*E)/((ew2-ew1)*(ew2-ew3))
         EB3=MATMUL(M-EW1*E,M-EW2*E)/((ew3-ew1)*(ew3-ew2))
      endif

      RETURN
      END


c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c +                                                               +
c +   Fast Fourier Transformation Code                            +
c +                                                               +
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine fourn(Datx,nn,ndim,isg)
      implicit none
      integer isg,ndim,nn(ndim)
      integer i1,i2,i2rev,i3,i3rev,ibit,idim,ifp1
      integer k2,n,nprev,nrem,ntot,ifp2,ip1,ip2,ip3,k1
      real(8) Datx(*)
      real(8) tempi,tempr
      real(8) theta,wi,wpi,wpr,wr,wtemp
      ntot=1
      do 11 idim=1,ndim
        ntot=ntot*nn(idim)
11    continue
      nprev=1
      do 18 idim=1,ndim
        n=nn(idim)
        nrem=ntot/(n*nprev)
        ip1=2*nprev
        ip2=ip1*n
        ip3=ip2*nrem
        i2rev=1
        do 14 i2=1,ip2,ip1
          if(i2.lt.i2rev)then

            do 13 i1=i2,i2+ip1-2,2
              do 12 i3=i1,ip3,ip2
                i3rev=i2rev+i3-i2
                tempr=Datx(i3)
                tempi=Datx(i3+1)
                Datx(i3)=Datx(i3rev)
                Datx(i3+1)=Datx(i3rev+1)
                Datx(i3rev)=tempr
                Datx(i3rev+1)=tempi
12            continue
13          continue
          endif
          ibit=ip2/2
1         if ((ibit.ge.ip1).and.(i2rev.gt.ibit)) then
            i2rev=i2rev-ibit
            ibit=ibit/2
          goto 1
          endif

          i2rev=i2rev+ibit
14      continue
        ifp1=ip1
2       if(ifp1.lt.ip2)then
          ifp2=2*ifp1
          theta=isg*6.28318530717959d0/(ifp2/ip1)
          wpr=-2.d0*dsin(0.5d0*theta)**2
          wpi=dsin(theta)
          wr=1.d0
          wi=0.d0
          do 17 i3=1,ifp1,ip1
            do 16 i1=i3,i3+ip1-2,2
              do 15 i2=i1,ip3,ifp2
                k1=i2
                k2=k1+ifp1
                tempr=sngl(wr)*Datx(k2)-sngl(wi)*Datx(k2+1)
                tempi=sngl(wr)*Datx(k2+1)+sngl(wi)*Datx(k2)
                Datx(k2)=Datx(k1)-tempr

                Datx(k2+1)=Datx(k1+1)-tempi
                Datx(k1)=Datx(k1)+tempr
                Datx(k1+1)=Datx(k1+1)+tempi
15            continue
16          continue
            wtemp=wr
            wr=wr*wpr-wi*wpi+wr
            wi=wi*wpr+wtemp*wpi+wi
17        continue
          ifp1=ifp2
        goto 2
        endif
        nprev=n*nprev
18    continue
      return
      end

c**********************************************
      subroutine icams_conv33to6(Am,ibx1,ibx2,Av)
c**********************************************
      implicit none
      integer i
      integer ibx1(9),ibx2(9)
      real(8) Am(3,3),Av(6)
      do i=1,6
         Av(i)=Am(ibx1(i),ibx2(i))
      enddo
      return
      end


c***********************************************
      subroutine icams_conv33to9(Am,ibx1,ibx2,Av)
c***********************************************
      implicit none
      integer i
      integer ibx1(9),ibx2(9)
      real(8) Am(3,3),Av(9)
      do i=1,9
         Av(i)=Am(ibx1(i),ibx2(i))
      enddo
      return
      end


c**************************************************
      subroutine icams_conv6to33(Av,ibx1,ibx2,Am)
c**************************************************
      implicit none
      integer i
      integer ibx1(9),ibx2(9)
      real(8) Am(3,3),Av(6)
      do i=1,6
         Am(ibx1(i),ibx2(i))=Av(i)
      enddo
      do i=7,9
         Am(ibx1(i),ibx2(i))=Av(i-3)
      enddo
      return
      end


c====================================================================================================
c     Computes all eigenvalues and eigenvectors of a real symmetric matrix a, which is of size n
c     by n, stored in a physical np by np array. On output, elements of a above the diagonal are
c     destroyed. d returns the eigenvalues of a in its first n elements. v is a matrix with the same
c     logical and physical dimensions as a, whose columns contain, on output, the normalized
c     eigenvectors of a. nrot returns the number of Jacobi rotations that were required.
c====================================================================================================
      subroutine jacobi(a,n,np,d,v,nrot,ising)
      implicit none
      integer n,np,nrot,nmax
      real(8) a(np,np),d(np),v(np,np)
      parameter (nmax=500)
      integer i,ip,iq,j,ising
      real(8) c,g,h,s,sm,t,tau,theta,tresh,b(nmax),z(nmax)
      ising=0
      do ip=1,n  !initialize to the identity matrix.
         do iq=1,n
            v(ip,iq)=0.
         enddo
         v(ip,ip)=1.
      enddo 
      do ip=1,n
         b(ip)=a(ip,ip) !initialize b and d to the diagonal of a.
         d(ip)=b(ip)
         z(ip)=0.  !this vector will accumulate terms of the form tapq
      enddo        !as in equation (11.1.14).
      nrot=0
      do i=1,50
         sm=0.
         do ip=1,n-1  !sum off-diagonal elements.
         do iq=ip+1,n
            sm=sm+dabs(a(ip,iq))
         enddo 
         enddo 
         if(sm==0.)return   
         !===================!
         ! sucessful return  !
         !===================!
         if(i.lt.4)then       
            tresh=0.2*sm/n**2  !...on the first three sweeps.
         else
            tresh=0.           !...thereafter.
         endif
         do ip=1,n-1
         do iq=ip+1,n
            g=100.*dabs(a(ip,iq))
            !after four sweeps, skip the rotation if the off-diagonal element is small.
            if((i.gt.4).and.(dabs(d(ip))+g==dabs(d(ip)))
     *         .and.(dabs(d(iq))+g==dabs(d(iq))))then
               a(ip,iq)=0.
            else if(dabs(a(ip,iq)).gt.tresh)then
               h=d(iq)-d(ip)
               if(dabs(h)+g==dabs(h))then
                  t=a(ip,iq)/h          !t = 1/(2*theta)
               else
                  theta=0.5*h/a(ip,iq)  !equation (11.1.10).
                  t=1./(dabs(theta)+dsqrt(1.+theta**2))
                  if(theta.lt.0.)t=-t
               endif
               c=1./dsqrt(1+t**2)
               s=t*c
               tau=s/(1.+c)
               h=t*a(ip,iq)
               z(ip)=z(ip)-h
               z(iq)=z(iq)+h
               d(ip)=d(ip)-h
               d(iq)=d(iq)+h
               a(ip,iq)=0.
               do j=1,ip-1  !case of rotations 1 = j < p.
                  g=a(j,ip)
                  h=a(j,iq)
                  a(j,ip)=g-s*(h+g*tau)
                  a(j,iq)=h+s*(g-h*tau)
               enddo 
               do j=ip+1,iq-1 !case of rotations p < j < q.
                  g=a(ip,j)
                  h=a(j,iq)
                  a(ip,j)=g-s*(h+g*tau)
                  a(j,iq)=h+s*(g-h*tau)
               enddo 
               do j=iq+1,n !case of rotations q < j = n.
                  g=a(ip,j)
                  h=a(iq,j)
                  a(ip,j)=g-s*(h+g*tau)
                  a(iq,j)=h+s*(g-h*tau)
               enddo
               do j=1,n
                  g=v(j,ip)
                  h=v(j,iq)
                  v(j,ip)=g-s*(h+g*tau)
                  v(j,iq)=h+s*(g-h*tau)
               enddo 
               nrot=nrot+1
            endif
         enddo 
         enddo 
         do ip=1,n
            b(ip)=b(ip)+z(ip)
            d(ip)=b(ip) !update d with the sum of tapq ,
            z(ip)=0.    !and reinitialize z.
         enddo
      enddo 

      write(*,*) 'too many iterations in jacobi, give up'
      ising=1

      return
      end

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+                                                           +
c+     calculate stress due to dislocation line              +
c+                                                           +
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine stress_dislocation(p1,p2,t,b,mu,nu,p,sm)
      implicit none
      integer i,j,k,l,m,n,ip
      real(8) p1(3),p2(3),p(3),t(3),b(3),mu,nu,sm(3,3),sv(6)
      real(8) r(2,3),rL(2,3),rP(2,3),rY(2,3)
      real(8) c_r(2),c_rL(2),c_rP(2),c_rY(2)
      real(8) stf(2,3,3,3)
      real(8) x1,x2,x3,mx1(3,3),mx2(3,3),vx(3)
      real(8) eij(3,3),eijk(3,3,3)
      real(8) c1,c2,c3,c4,c5
      real(8) z1,z2,z3,z4,z5,z6,z7

      eij(1,:)=[1,0,0]
      eij(2,:)=[0,1,0]
      eij(3,:)=[0,0,1]
      eijk=0
      eijk(1,2,3)= 1
      eijk(1,3,2)=-1
      eijk(2,1,3)=-1
      eijk(2,3,1)= 1
      eijk(3,1,2)= 1
      eijk(3,2,1)=-1

      r(1,:)=p-p1
      c_r(1)=dsqrt(sum(r(1,:)**2))
      c_rL(1)=dot_product(r(1,:),t)
      rP(1,:)=r(1,:)-c_rL(1)*t;          
      rY(1,:)=r(1,:)+c_r(1)*t
      c_rY(1)=dsqrt(sum(rY(1,:)**2))

      r(2,:)=p-p2
      c_r(2)=dsqrt(sum(r(2,:)**2))
      c_rL(2)=dot_product(r(2,:),t)
      rP(2,:)=r(2,:)-c_rL(2)*t;          
      rY(2,:)=r(2,:)+c_r(2)*t
      c_rY(2)=dsqrt(sum(rY(2,:)**2))

c      print*,c_r
c      print*,c_rL
c      print*,c_rY
c      read*

      stf=0
      do ip=1,2
         c1=mu/(3.14*c_rY(ip)**2)
         c2=1/(1-nu)
         vx=0
         do i=1,3
         do j=1,3
         do k=1,3
            vx(k)=vx(k)+eijk(k,i,j)*rY(ip,i)*t(j)
         enddo
         enddo
         enddo
         vx=vx/(2*(1-nu))
         c4=2/c_rY(ip)**2
         c5=c_rL(ip)/c_r(ip)
         do i=1,3
         do j=1,3
         do k=1,3
            z1=( dot_product(eijk(i,k,:),rY(ip,:))*t(j)
     &          +dot_product(eijk(j,k,:),rY(ip,:))*t(i))/2
            z2=( dot_product(eijk(i,k,:),t       )*rY(ip,j)
     &          +dot_product(eijk(j,k,:),t       )*rY(ip,i))/2
            z3=eij(i,j)+t(i)*t(j)
            z4=rP(ip,i)*rY(ip,j)+rP(ip,j)*rY(ip,i)
            z5=rY(ip,i)*rY(ip,j)
            stf(ip,i,j,k)=c1*(z1-c2*z2-vx(k)*(z3+c4*(z4+c5*z5)))         
         enddo
         enddo
         enddo
      enddo
      do i=1,3
      do j=1,3
         sm(i,j)=dot_product(stf(2,i,j,:)-stf(1,i,j,:),b)
      enddo
      enddo

      return
      end




c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c +                                                               +
c +   Crystal Plasticity Fast Fourier Transformation Code         +
c +                                                               +
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine smooth(nx,ny,nz,L_size,rve_X)
      implicit none
      integer i,j,k,k1,ising
      integer nx,ny,nz,ix,iy,iz,mx,my,mz
      integer, parameter :: ax=0
      integer, parameter :: ay=0
      integer, parameter :: az=0
      real(8) fct_smooth,L_size
      real(8) WSfft((nx+ax*2)*(ny+ay*2)*(nz+az*2)*2)
      real(8) rve_V((nx+ax*2),(ny+ay*2),(nz+az*2))
      real(8) rve_R((nx+ax*2),(ny+ay*2),(nz+az*2))
      real(8) rve_I((nx+ax*2),(ny+ay*2),(nz+az*2))
      real(8) rve_X(nx,ny,nz)
      real(8) x1,x2,x3,Fqc(3)

      mx=nx+ax*2
      my=ny+ay*2
      mz=nz+az*2

!      fct_smooth = 4*3.14**2*10.d0
!      fct_smooth = 4*3.14**2*20.d0
!      fct_smooth = 4*3.14**2*40.d0
!      fct_smooth = 4*3.14**2*100.d0
!      fct_smooth = 4*3.14**2*50.d0
!      fct_smooth = 4*3.14**2*30.d0

      fct_smooth = 4*3.14**2*60.d0
!      fct_smooth = 4*3.14**2*40.d0
!      fct_smooth = 4*3.14**2*40.d0*1.d-6

!      fct_smooth = 4*3.14**2*4.d19*L_size**2

c-----transfer stress to frequency space
      rve_V=0
      rve_V( ax+1:ax+nx, ay+1:ay+ny, az+1:az+nz ) = rve_X 

      WSfft=0               
      k1=0
      do iz=1,mz
      do iy=1,my
      do ix=1,mx
         k1=k1+1; WSfft(k1)=rve_V(ix,iy,iz)
         k1=k1+1; WSfft(k1)=0               
      enddo
      enddo
      enddo
      call fourn(WSfft,[mx,my,mz],3,1)
      WSfft=WSfft/(mx*my*mz)
      k1=0
      do iz=1,mz
      do iy=1,my
      do ix=1,mx
         k1=k1+1; rve_R(ix,iy,iz)=WSfft(k1) !--real part
         k1=k1+1; rve_I(ix,iy,iz)=WSfft(k1) !--image part
      enddo
      enddo
      enddo

c-----calculate equilibrium value in frequency space
      do ix=1,mx
      do iy=1,my
      do iz=1,mz
         if(ix<=mx/2+1)then
            Fqc(1) = (ix-1.)/mx
         else
            Fqc(1) = (ix-(mx+1.))/mx
         endif
         if(iy<=my/2+1)then
            Fqc(2) = (iy-1.)/my
         else
            Fqc(2) = (iy-(my+1.))/my
         endif
         if(iz<=mz/2+1)then
            Fqc(3) = (iz-1.)/mz
         else
            Fqc(3) = (iz-(mz+1.))/mz
         endif
         x1= 1 + fct_smooth*sum(Fqc**2)
         rve_R(ix,iy,iz) = rve_R(ix,iy,iz)/x1
         rve_I(ix,iy,iz) = rve_I(ix,iy,iz)/x1
      enddo
      enddo
      enddo

c-----transfer new value from frequency space to physical space
      k1=0
      do iz=1,mz
      do iy=1,my
      do ix=1,mx
         k1=k1+1; WSfft(k1)=rve_R(ix,iy,iz)
         k1=k1+1; WSfft(k1)=rve_I(ix,iy,iz)
      enddo
      enddo
      enddo
      call fourn(WSfft,[mx,my,mz], 3, -1)
      k1=0
      do iz=1,mz
      do iy=1,my
      do ix=1,mx
         k1=k1+1; rve_V(ix,iy,iz)=WSfft(k1) 
         k1=k1+1
      enddo
      enddo
      enddo

      rve_X = rve_V( ax+1:ax+nx, ay+1:ay+ny, az+1:az+nz ) 
c
      return
      end

c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c +                                                               +
c +   Output date to ovito                                        +
c +                                                               +
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine output(nx,ny,nz,sigR)
         implicit none
         integer ix,iy,iz,nx,ny,nz
         real(8) sigR(6,nx,ny,nz),v3(3)

         open(888,file='../stress_in.chkpt')

         write(888,'("#F A 1 1 1 3 1 1")')
         write(888,'("#C Np Ng Nl x y z sig11 sig22 sig33 
     &   sig12 sig13 sig23")')
         write(888,'("#X 1 0 0")')
         write(888,'("#Y 0 1 0")')
         write(888,'("#Z 0 0 1")')
         write(888,'("#E")')

         do ix=1,nx
         do iy=1,ny
         do iz=1,nz
            v3=[ix,iy,iz]*2.d0  
            write(888,'(3I10,100e12.3)') 
     &         1,1,1,v3,sigR(1,ix,iy,iz),sigR(2,ix,iy,iz),
     &         sigR(3,ix,iy,iz),sigR(4,ix,iy,iz),sigR(5,ix,iy,iz),
     &         sigR(6,ix,iy,iz)
         enddo
         enddo
         enddo
         
         close(888)

         return
         end
