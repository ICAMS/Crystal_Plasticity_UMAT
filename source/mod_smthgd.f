c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c +                                                               +
c +   Calculate the first and second order of strain gradient     +
c +                                                               +
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine  cal_Fp_grad(Nele,Ngpt,        ! FEM mesh size
     &                        Nftx,Nfty,Nftz,   ! FFT grid size
     &                        CD_ref,           ! reference Smooth fct 
     &                        detx,dety,detz,   ! FFT grid increment
     &                        fft_D0,           ! Smooth factor 
     &                        fft_xyz,          ! FFT grid coordinates
     &                        fft_Fqc,          ! FFT grid freuencies
     &                        fft_Inf,          ! FFT-->FEM connection 
     &                        fem_Inf,          ! FEM-->FFT connection
     &                        fem_xyz,          ! FEM gp coordinates
     &                        fem_Fp0,          ! original Fp
     &                        fem_Fp1,          ! smoothed Fp
     &                        fem_dFpdx,        ! 1st grad of smthed Fp
     &                        fem_ddFpddx)      ! 2st grad of smthed Fp

      implicit none
      integer Nele,Ngpt
      integer fem_Inf    (Nele,Ngpt,4)
      real(8) fem_xyz    (Nele,Ngpt,3)
      real(8) fem_Fp0    (Nele,Ngpt,3,3)
      real(8) fem_Fp1    (Nele,Ngpt,3,3)
      real(8) fem_dFpdx  (Nele,Ngpt,3,3,3)
      real(8) fem_ddFpddx(Nele,Ngpt,3,3,3,3)
      integer Nftx,Nfty,Nftz
      real(8) fft_D0     (Nftx,Nfty,Nftz)
      real(8) fft_Fqc    (Nftx,Nfty,Nftz,3)
      real(8) fft_Ieg    (Nftx,Nfty,Nftz,2)
      integer fft_Inf    (Nftx,Nfty,Nftz,3)
      real(8) fft_xyz    (Nftx,Nfty,Nftz,3)
      real(8) fft_V0     (Nftx,Nfty,Nftz)
      real(8) fft_V1     (Nftx,Nfty,Nftz)
      real(8) fft_dVdx   (Nftx,Nfty,Nftz,3)
      real(8) fft_ddVddx (Nftx,Nfty,Nftz,3,3)
      real(8) Lx(6),Ly(6),Lz(6),CDiff(3)
c
      real(8) V1,dVdx(3),ddVddx(3,3),Fqc(3)
      real(8) detx,dety,detz,x_min,x_tmp,y_tmp
      real(8) rsd,Crsd,CDI,CDM_fre,CDM_rgd,CD_ref
      real(8) x1,x2,x3,x_ie,x_ig
      integer i,j,k,ising,i1,i2,i3,i4,j1,j2,j3,j4
      integer ix,iy,iz,Iloop,Nloop,Nx,ip1,ip2,iex,igx,idx
c
c----------------------------------------------------
c-----Calcuate 1st and 2nd gradient of smoothed field
c----------------------------------------------------
      fft_V0=0
      fft_V1=0
      fem_Fp1=0
      fem_dFpdx=0
      fem_ddFpddx=0
      do ip1=1,3
      do ip2=1,3

        !--------------------------------------- 
        ! Mapping Fp_ij from FEM gps to FFT grids
        !--------------------------------------- 
         do ix=1,Nftx
         do iy=1,Nfty
         do iz=1,Nftz
            iex=fft_Inf(ix,iy,iz,1)
            igx=fft_Inf(ix,iy,iz,2)
            if(iex/=0.and.igx/=0) then
c               fft_V0(ix,iy,iz)=fem_Fp0(iex,igx,ip1,ip2)

               if(fft_Inf(ix,iy,iz,3)==1)then              !=====> internal
                  fft_V0(ix,iy,iz)=fem_Fp0(iex,igx,ip1,ip2)
               elseif(fft_Inf(ix,iy,iz,3)==2)then          !=====> rigid
                  if(ip1/=ip2) fft_V0(ix,iy,iz)=0.d0
                  if(ip1==ip2) fft_V0(ix,iy,iz)=1.d0
               else                                        !=====> free
                  fft_V0(ix,iy,iz)=fem_Fp0(iex,igx,ip1,ip2)*1.0
               endif
            endif


c       print*,ix,iy,iz,fft_Inf(ix,iy,iz,3),fft_V0(ix,iy,iz)

         enddo         
         enddo         
         enddo         

         !-------------------------------------------------------
         ! Calculate smoothed field Fp'_ij of component Fp_ij 
         !-------------------------------------------------------
         call cal_V_smth(CD_ref,detx,dety,detz,
     &                   Nftx,Nfty,Nftz,
     &                   fft_Fqc,fft_D0,
     &                   fft_V0,fft_V1)

         !-------------------------------------------------------
         ! Calculate 1st and 2nd gradiend of Fp'_ij
         !-------------------------------------------------------
         call cal_V_grad(fft_V0,fft_V1,
     &                   Nftx,Nfty,Nftz,
     &                   detx,dety,detz,
     &                   fft_dVdx,fft_ddVddx)


         !-------------------------------------------------------
         ! Mapping values from FFT grids to FEM gps
         !-------------------------------------------------------
         do i1=1,Nele
         do i2=1,Ngpt
         if(fem_Inf(i1,i2,1)/=0)then
            ix=fem_Inf(i1,i2,2)
            iy=fem_Inf(i1,i2,3)
            iz=fem_Inf(i1,i2,4)
            fem_Fp1(i1,i2,ip1,ip2)=fft_V1(ix,iy,iz)
            fem_dFpdx(i1,i2,ip1,ip2,:)=fft_dVdx(ix,iy,iz,:)
            fem_ddFpddx(i1,i2,ip1,ip2,:,:)=fft_ddVddx(ix,iy,iz,:,:)
         endif
         enddo
         enddo

c         !-------------------------------------------------------
c         ! Remove imaginary matrix
c         !-------------------------------------------------------
c         do ix=1,Nftx
c         do iy=1,Nfty
c         do iz=1,Nftz
c            x1=fft_xyz(ix,iy,iz,1)
c            x2=fft_xyz(ix,iy,iz,2)
c            x3=fft_xyz(ix,iy,iz,3)
c            if((x1<Lx(2).or.x1>Lx(3)).or.
c     &         (x2<Ly(2).or.x2>Ly(3)).or.
c     &         (x3<Lz(2).or.x3>Lz(3))) then
c               fft_V0    (ix,iy,iz)    =0
c               fft_V1    (ix,iy,iz)    =0
c               fft_dVdx  (ix,iy,iz,:)  =0
c               fft_ddVddx(ix,iy,iz,:,:)=0
c            endif
c         enddo
c         enddo
c         enddo
c         !-------------------------------------------------------
c         ! Make plot file for ovito
c         !-------------------------------------------------------
c         if(ip1==1.and.ip2==1)then
c            call output_ovito(fft_V0,fft_V1,fft_dVdx,fft_ddVddx,
c     &                        Nftx,Nfty,Nftz)
c         endif

      enddo
      enddo
c
      return
      end


c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c +                                                               +
c +   Calculate smoothed field by diffusion method                +
c +                                                               +
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine cal_V_smth(CD_ref,detx,dety,detz,Nftx,Nfty,Nftz,
     &                      fft_Fqc,fft_D0,fft_V0,fft_V1)
      implicit none
      integer Nftx,Nfty,Nftz,Nele,Ngpt
      real(8) WSfft  (Nftx*Nfty*Nftz*2)
      real(8) fft_Fqc(Nftx,Nfty,Nftz,3)
      real(8) fft_R0 (Nftx,Nfty,Nftz)
      real(8) fft_I0 (Nftx,Nfty,Nftz)
      real(8) fft_R1 (Nftx,Nfty,Nftz)
      real(8) fft_I1 (Nftx,Nfty,Nftz)
      real(8) fft_Rs (Nftx,Nfty,Nftz)
      real(8) fft_Is (Nftx,Nfty,Nftz)
      real(8) fft_D0 (Nftx,Nfty,Nftz)
      real(8) fft_V0 (Nftx,Nfty,Nftz)
      real(8) fft_V1 (Nftx,Nfty,Nftz)
      real(8) Lx(6),Ly(6),Lz(6)
      real(8) rsd,Crsd,CD_ref,Fqc(3)
      real(8) detx,dety,detz
c
      integer i,j,k,k1,ising,i1,i2,i3,i4,j1,j2,j3,j4
      integer ix,iy,iz,Iloop,Nloop,ip
      real(8) x1,x2,x3
c
      Nloop=500
      fft_R0=0
      fft_I0=0

c-----transfer V00 to frequency space
      WSfft=0               
      k1=0
      do iz=1,Nftz
      do iy=1,Nfty
      do ix=1,Nftx
         k1=k1+1; WSfft(k1)=fft_V0(ix,iy,iz)
         k1=k1+1; WSfft(k1)=0               
      enddo
      enddo
      enddo
      call fourn(WSfft,[Nftx,Nfty,Nftz],3,1)
      WSfft=WSfft/(Nftx*Nfty*Nftz)
      k1=0
      do iz=1,Nftz
      do iy=1,Nfty
      do ix=1,Nftx
         k1=k1+1; fft_R0(ix,iy,iz)=WSfft(k1) !--real part
         k1=k1+1; fft_I0(ix,iy,iz)=WSfft(k1) !--image part
      enddo
      enddo
      enddo

      fft_R1=fft_R0
      fft_I1=fft_I0
      Crsd=1.d-8
      rsd=Crsd*2
      Iloop=0

      do while(Iloop<Nloop) 
         Iloop=Iloop+1

c--------check convergence of (R1,I1) in hetergenious material
         rsd=0
         do ix=1,Nftx
         do iy=1,Nfty
         do iz=1,Nftz
            Fqc = fft_fqc(ix,iy,iz,:)
            x1 = 1 + fft_D0(ix,iy,iz)*sum(Fqc**2)
            rsd = rsd + dabs(fft_R1(ix,iy,iz)*x1 - fft_R0(ix,iy,iz))
     &                + dabs(fft_I1(ix,iy,iz)*x1 - fft_I0(ix,iy,iz))
         enddo
         enddo
         enddo
         rsd = rsd/(Nftx*Nfty*Nftz)

c         print '(2I5,10e15.5)',Iloop,Nloop,rsd,Crsd
c         read*

         if(rsd < Crsd) goto 102

c--------calculate sources
         rsd=0
         do ix=1,Nftx
         do iy=1,Nfty
         do iz=1,Nftz
            Fqc = fft_fqc(ix,iy,iz,:)
            x1 = (fft_D0(ix,iy,iz)-CD_ref)*sum(Fqc**2)
            fft_Rs(ix,iy,iz) = fft_R0(ix,iy,iz) - x1*fft_R1(ix,iy,iz)
            fft_Is(ix,iy,iz) = fft_I0(ix,iy,iz) - x1*fft_I1(ix,iy,iz)
         enddo
         enddo
         enddo

c--------calculate (R1,I1) for homegenious reference material
         do ix=1,Nftx
         do iy=1,Nfty
         do iz=1,Nftz
            Fqc = fft_fqc(ix,iy,iz,:)
            x1= 1 + CD_ref*sum(Fqc**2)
            fft_R1(ix,iy,iz) = fft_Rs(ix,iy,iz)/x1
            fft_I1(ix,iy,iz) = fft_Is(ix,iy,iz)/x1
         enddo
         enddo
         enddo

      enddo

102   continue

c-----Transfer (R1,I1) from frequency space to physical space
      k1=0
      do iz=1,Nftz
      do iy=1,Nfty
      do ix=1,Nftx
         k1=k1+1; WSfft(k1)=fft_R1(ix,iy,iz)
         k1=k1+1; WSfft(k1)=fft_I1(ix,iy,iz)
      enddo
      enddo
      enddo
      call fourn(WSfft,[Nftx,Nfty,Nftz], 3, -1)
      k1=0
      do iz=1,Nftz
      do iy=1,Nfty
      do ix=1,Nftx
         k1=k1+1; fft_V1(ix,iy,iz)=WSfft(k1) 
         k1=k1+1
      enddo
      enddo
      enddo
c
      return
      end




c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c +                                                               +
c +   Calculate first and second gradient                         +
c +                                                               +
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine cal_V_grad(fft_V0,fft_V1,
     &                      Nftx,Nfty,Nftz,
     &                      detx,dety,detz,
     &                      fft_dVdx,
     &                      fft_ddVddx)
      implicit none
      integer Nftx,Nfty,Nftz
      real(8) fft_V0    (Nftx,Nfty,Nftz)
      real(8) fft_V1    (Nftx,Nfty,Nftz)
      real(8) fft_dVdx  (Nftx,Nfty,Nftz,3)
      real(8) fft_ddVddx(Nftx,Nfty,Nftz,3,3)
      real(8) detx,dety,detz
c
      integer i,j,k
      integer ix,iy,iz,ix1,ix2,iy1,iy2,iz1,iz2
      real(8) x1,x2,x3

c--------------------
      fft_dVdx=0
c--------------------
      do iz=1,Nftz
      do iy=1,Nfty
      do ix=1,Nftx
         if(ix==1)then
            ix1=Nftx
            ix2=ix+1
         elseif(ix==Nftx)then
            ix1=Nftx-1
            ix2=1
         else
            ix1=ix-1
            ix2=ix+1
         endif
         if(iy==1)then
            iy1=Nfty
            iy2=iy+1
         elseif(iy==Nfty)then
            iy1=Nfty-1
            iy2=1
         else
            iy1=iy-1
            iy2=iy+1
         endif
         if(iz==1)then
            iz1=Nftz
            iz2=iz+1
         elseif(iz==Nftz)then
            iz1=Nftz-1
            iz2=1
         else
            iz1=iz-1
            iz2=iz+1
         endif

         if(Nftx>2)then
            x1=fft_V1(ix1, iy, iz  )
            x2=fft_V1(ix2, iy, iz  )
!            x1=fft_V0(ix1, iy, iz  )
!            x2=fft_V0(ix2, iy, iz  )
            fft_dVdx (ix,  iy, iz,1)=(x2-x1)/2/detx
         endif
         if(Nfty>2)then
            x1=fft_V1(ix, iy1, iz  )
            x2=fft_V1(ix, iy2, iz  )
!            x1=fft_V0(ix, iy1, iz  )
!            x2=fft_V0(ix, iy2, iz  )
            fft_dVdx (ix, iy,  iz,2)=(x2-x1)/2/dety
         endif
         if(Nftz>2)then
            x1=fft_V1(ix, iy, iz1  )
            x2=fft_V1(ix, iy, iz2  )
!            x1=fft_V0(ix, iy, iz1  )
!            x2=fft_V0(ix, iy, iz2  )
            fft_dVdx (ix, iy, iz, 3)=(x2-x1)/2/detz
         endif

      enddo
      enddo
      enddo

c-----------------------
      fft_ddVddx=0
c-----------------------
      do iz=1,Nftz
      do iy=1,Nfty
      do ix=1,Nftx
         if(ix==1)then
            ix1=Nftx
            ix2=ix+1
         elseif(ix==Nftx)then
            ix1=Nftx-1
            ix2=1
         else
            ix1=ix-1
            ix2=ix+1
         endif
         if(iy==1)then
            iy1=Nfty
            iy2=iy+1
         elseif(iy==Nfty)then
            iy1=Nfty-1
            iy2=1
         else
            iy1=iy-1
            iy2=iy+1
         endif
         if(iz==1)then
            iz1=Nftz
            iz2=iz+1
         elseif(iz==Nftz)then
            iz1=Nftz-1
            iz2=1
         else
            iz1=iz-1
            iz2=iz+1
         endif

         do i=1,3
            if(Nftx>2)then
               x1=fft_dVdx(ix1, iy, iz, i  )
               x2=fft_dVdx(ix2, iy, iz, i  )
               fft_ddVddx (ix,  iy, iz, i,1)=(x2-x1)/2/detx
            endif
            if(Nfty>2)then
               x1=fft_dVdx(ix, iy1, iz, i  )
               x2=fft_dVdx(ix, iy2, iz, i  )
               fft_ddVddx (ix,  iy, iz, i,2)=(x2-x1)/2/dety
            endif
            if(Nftz>2)then
               x1=fft_dVdx(ix, iy, iz1, i  )
               x2=fft_dVdx(ix, iy, iz2, i  )
               fft_ddVddx (ix, iy, iz,  i,3)=(x2-x1)/2/detz
            endif
         enddo

      enddo
      enddo
      enddo
c
      return
      end


c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c +                                                               +
c +   Output date to ovito                                        +
c +                                                               +
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine output_ovito(fft_V0,fft_V1,fft_dVdx,fft_ddVddx,
     &                        Nftx,Nfty,Nftz)
         implicit none
         integer Nftx,Nfty,Nftz
         real(8) fft_V0(Nftx,Nfty,Nftz)
         real(8) fft_V1(Nftx,Nfty,Nftz)
         real(8) fft_dVdx(Nftx,Nfty,Nftz,3)
         real(8) fft_ddVddx(Nftx,Nfty,Nftz,3,3)
c
         integer i,j,k
         integer ix,iy,iz
         real(8) x1,x2,x3,v3(3)
c
         open(923,file='smth.chkpt')
c
         write(923,'("#F A 1 1 1 3 1 1")')
         write(923,'("#C Np Ng Nl x y z v0 v1 g1 g2 g3")')
         write(923,'("#X 1 0 0")')
         write(923,'("#Y 0 1 0")')
         write(923,'("#Z 0 0 1")')
         write(923,'("#E")')
         do ix=1,Nftx
         do iy=1,Nfty
         do iz=1,Nftz
            v3=[ix,iy,iz]*2.d0  
            write(923,'(3I10,100e12.3)') 
     &         1,1,1,v3,
     &         fft_V0(ix,iy,iz),
     &         fft_V1(ix,iy,iz),
     &         fft_dVdx(ix,iy,iz,1:3)
         enddo
         enddo
         enddo
         close(923)
         return
      endsubroutine

c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c +                                                               +
c +   Fast Fourier Transformation Code                            +
c +                                                               +
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine fourn_tmp(Datx,nn,ndim,isg)
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

