c================================================================
c
c    weak coupling code
c
c================================================================
      module mod_mesh_grid_size
         implicit none
         integer, parameter :: Tnel=27000
         integer, parameter :: Tngp=8
         integer, parameter :: Tnfx=64
         integer, parameter :: Tnfy=64
         integer, parameter :: Tnfz=64
         integer, parameter :: Iwkcoup_trip=0  !-->with(1), without(0) 
         integer, parameter :: Iwkcoup_grad=0  !-->with(1), without(0) 
         integer, parameter :: Iwkcoup_int=0   !-->with(1), without(0)
         integer, parameter :: Iwkcoup_sup=0   !-->with(1), without(0)
         integer, parameter :: Iwkcoup_bk=3    !-->with(1)AFKH, with(2)CHKH, with(3)OWKH, without(0)
      endmodule mod_mesh_grid_size


c================================================================
c
c    weak coupling module: Grad
c    output for plastic deformation calculation: Rho_gnd, IVB_gnd, pk2i_gnd
c
c================================================================
      module mod_wkcoup_grad
         use mod_mesh_grid_size
         implicit none
         integer, parameter :: I_gnd_iso = 1       !-->with(1),without(0)
         integer, parameter :: I_gnd_kin = 1       !-->with(1),without(0)  
         real(8), parameter :: C_taui    = 0.01    !-->Taylor hardening parameter
         real(8), parameter :: L_size    = 1.0d-9  !-->unit m, length scale
         real(8), parameter :: C_unit    = 1.0d-6  !-->unit m, mesh
         real(8), parameter :: B_lattice = 2.5d-10 !-->unit m, Burgers vector 
         real(8), parameter :: CD_smooth = 3.d-15  
c--------strain gradient
         real(8) fem_Fp0        (Tnel,Tngp,3,3)
         real(8) fem_dFp        (Tnel,Tngp,3,3)
         real(8) fem_dFpp       (Tnel,Tngp,3,3)
         real(8) fem_dFpp_1gd   (Tnel,Tngp,3,3,3)
         real(8) fem_dFpp_2gd   (Tnel,Tngp,3,3,3,3)
c
         integer Icall_fem_fft_map
         integer Icall_ieig_grad(Tnel,Tngp)
         real(8) fem_xyz        (Tnel,Tngp,3)
         integer fem_Inf        (Tnel,Tngp,4)
         real(8) BOX_x(4),BOX_y(4),BOX_z(4)
         real(8) CD_parts(3),CD_ref
         integer fft_Inf        (Tnfx,Tnfy,Tnfz,3)
         real(8) fft_Fqc        (Tnfx,Tnfy,Tnfz,3)
         real(8) fft_xyz        (Tnfx,Tnfy,Tnfz,3)
         real(8) fft_CD         (Tnfx,Tnfy,Tnfz  )
         real(8) detx,dety,detz
      contains
c================================================================
c        grad: state variable ini
c================================================================
         subroutine sub_grad_ini
         implicit none
         integer i,j,k,l,m,n,ix1,ix2,ix3
         real(8) x1,x2,x3
         character(len=100) Fname,Tname
         Icall_ieig_grad=0
         Icall_fem_fft_map=0
         fem_Fp0=0
         fem_dFp=0
         fem_dFpp=0
         fem_dFpp_1gd=0
         fem_dFpp_2gd=0
         fem_Inf=0
         fem_xyz=0
         fft_Inf=0
         fft_xyz=0
         fft_Fqc=0
         fft_CD =0
         return
         endsubroutine

c================================================================
c        grad: strain gradient evolution
c================================================================
         subroutine sub_grad_evolution
         implicit none
         integer i,j,k,l,m,n,ix1,ix2,ix3,ii,np,ising,ising_1
         real(8) mx33_1(3,3),mx33_2(3,3),mx33_3(3,3)
         real(8) x1,x2,x3,dydx,ddyddx
         real(8) x_min,x_max,x_add
         real(8) y_min,y_max,y_add
         real(8) z_min,z_max,z_add
         real(8) Lenx,Leny,Lenz,x_tmp
         integer ix,iy,iz,i1,i2,j1,j2,iex,igx,idx

         CD_parts(1)=CD_smooth        ! inside  domain
         CD_parts(2)=CD_parts(1)   ! outside domain, rigid
         CD_parts(3)=CD_parts(1)   ! outside domain, free
         CD_ref = max( CD_parts(1),CD_parts(2),CD_parts(3) )
         
         x_min =+1.d50  
         x_max =-1.d50
         y_min =+1.d50  
         y_max =-1.d50
         z_min =+1.d50  
         z_max =-1.d50
         do i=1,Tnel
         do j=1,Tngp
         if(x_min>fem_xyz(i,j,1)) x_min=fem_xyz(i,j,1)
         if(y_min>fem_xyz(i,j,2)) y_min=fem_xyz(i,j,2)
         if(z_min>fem_xyz(i,j,3)) z_min=fem_xyz(i,j,3)
         if(x_max<fem_xyz(i,j,1)) x_max=fem_xyz(i,j,1)
         if(y_max<fem_xyz(i,j,2)) y_max=fem_xyz(i,j,2)
         if(z_max<fem_xyz(i,j,3)) z_max=fem_xyz(i,j,3)
         enddo
         enddo
         if(x_max==x_min) x_max=x_min+1
         if(y_max==y_min) y_max=y_min+1
         if(z_max==z_min) z_max=z_min+1
         x_add=(x_max-x_min)*1.d-1
         y_add=(y_max-y_min)*1.d-1
         z_add=(z_max-z_min)*1.d-1
         BOX_x=[x_min-x_add,x_min,x_max,x_max+x_add]
         BOX_y=[y_min-y_add,y_min,y_max,y_max+y_add]
         BOX_z=[z_min-z_add,z_min,z_max,z_max+z_add]

         if(Icall_fem_fft_map==0)then
            Icall_fem_fft_map=1
            Lenx = BOX_x(4)-BOX_x(1)
            Leny = BOX_y(4)-BOX_y(1)
            Lenz = BOX_z(4)-BOX_z(1)
            detx = Lenx/Tnfx
            dety = Leny/Tnfy
            detz = Lenz/Tnfz

            !-----------------------------------------
            !Frequency and coordinite for FFT grids
            !-----------------------------------------
            fft_Fqc=0
            do ix=1,Tnfx
            do iy=1,Tnfy
            do iz=1,Tnfz
               if(ix<=Tnfx/2+1)then
                  x1 = (ix-1.       )/Lenx
               else
                  x1 = (ix-(Tnfx+1.))/Lenx
               endif
               if(iy<=Tnfy/2+1)then
                  x2 = (iy-1.       )/Leny
               else
                  x2 = (iy-(Tnfy+1.))/Leny
               endif
               if(iz<=Tnfz/2+1)then
                  x3 = (iz-1.       )/Lenz
               else
                  x3 = (iz-(Tnfz+1.))/Lenz
               endif
               fft_Fqc(ix,iy,iz,:) = [x1,x2,x3]
               fft_xyz(ix,iy,iz,1) = BOX_x(1) + detx*(ix-1)
               fft_xyz(ix,iy,iz,2) = BOX_y(1) + dety*(iy-1)
               fft_xyz(ix,iy,iz,3) = BOX_z(1) + detz*(iz-1)
            enddo
            enddo
            enddo

            !----------------------------------------------------
            ! Connection between FFT and FEM and BCs
            !----------------------------------------------------
            fft_Inf=0
            do ix=1,Tnfx
            do iy=1,Tnfy
            do iz=1,Tnfz
               x1=fft_xyz(ix,iy,iz,1)
               x2=fft_xyz(ix,iy,iz,2)
               x3=fft_xyz(ix,iy,iz,3)
               x_min=1.d50; iex=0; igx=0
               do i1=1,Tnel
               do i2=1,Tngp
               if(fem_Inf(i1,i2,1)/=0)then
                  x_tmp=dsqrt(sum((fem_xyz(i1,i2,:)-[x1,x2,x3])**2))
                  if(x_min>x_tmp)then
                     x_min=x_tmp; iex=i1; igx=i2
                  endif                
               endif
               enddo
               enddo
               idx=1
               if( x1<BOX_x(2) ) idx=2
               if( x1>BOX_x(3) ) idx=3
               if( x2<BOX_y(2) ) idx=4
               if( x2>BOX_y(3) ) idx=5
               if( x3<BOX_z(2) ) idx=6
               if( x3>BOX_z(3) ) idx=7
               fft_Inf(ix,iy,iz,1:3)=[iex,igx,idx]
               fft_CD(ix,iy,iz)=CD_ref 
            enddo         
            enddo         
            enddo

            do i1=1,Tnel
            do i2=1,Tngp
            if(fem_Inf(i1,i2,1)/=0)then
               x_min=1.d50; iex=0; igx=0
               do ix=1,Tnfx
               do iy=1,Tnfy
               do iz=1,Tnfz
                  x1=fft_xyz(ix,iy,iz,1)
                  x2=fft_xyz(ix,iy,iz,2)
                  x3=fft_xyz(ix,iy,iz,3)
                  if(fft_Inf(ix,iy,iz,3)==1)then
                     x_tmp=dsqrt(sum((fem_xyz(i1,i2,:)-[x1,x2,x3])**2))
                     if(x_min>x_tmp)then
                        x_min=x_tmp; fem_Inf(i1,i2,2:4)=[ix,iy,iz]
                     endif                
                  endif
               enddo
               enddo
               enddo
            endif
            enddo
            enddo

         endif


         fem_dFpp=0
         fem_dFpp_1gd=0
         fem_dFpp_2gd=0
         call cal_Fp_grad(Tnel,Tngp,          ! FEM mesh size
     &                    Tnfx,Tnfy,Tnfz,     ! FFT grid size
     &                    CD_ref,             ! reference Smooth fct 
     &                    detx,dety,detz,     ! FFT grid increment
     &                    fft_CD,             ! Smooth factor 
     &                    fft_xyz,            ! FFT grid coordinates
     &                    fft_Fqc,            ! FFT grid freuencies
     &                    fft_Inf,            ! FFT-->FEM connection 
     &                    fem_Inf,            ! FEM-->FFT connection
     &                    fem_xyz,            ! FEM gp coordinates
     &                    fem_dFp,            ! dFp
     &                    fem_dFpp,           ! dFp_smth
     &                    fem_dFpp_1gd,       ! d(dFp_smth)_dx
     &                    fem_dFpp_2gd)       ! dd(dFp_smth)_ddx

         return
         endsubroutine

c================================================================
c        grad: additional hardening due to grad
c================================================================
         subroutine sub_grad_effect(iex,igx,IB1,IB2,
     &              Nslp_mx,Nslp,STFei26,
     &              vd_slp,vl_slp,vn_slp,
     &              Rho_gnd, IVB_gnd, pk2i_gnd)
            implicit none
            integer Nslp_mx,Nslp,iex,igx
            real(8) STFei26(6,6)
            real(8) vd_slp(Nslp_mx,3)
            real(8) vl_slp(Nslp_mx,3)
            real(8) vn_slp(Nslp_mx,3)
            integer IB1(9),IB2(9)
c
            real(8) Rho_gnd(Nslp_mx)
            real(8) IVB_gnd(Nslp_mx)
            real(8) pk2i_gnd(6)            
c
            integer i,j,k,m,n,is
            integer i1,i2,i3,i4,i5,i6
            integer j1,j2,j3,j4,j5,j6

            real(8) Fp0(3,3),dFpp(3,3),Fp(3,3)
            real(8) dFpp_1gd(3,3,3),dFpp_2gd(3,3,3,3)

            real(8) sm_R_gnd(3,3),sm_I_gnd(3,3)
            real(8) bvct(3),tvct(3),gvct(3),pvct1(3),pvct2(3),pvct(3)
            real(8) Xmu,Xnu,smij(3,3)
            real(8) x1,x2,x3
            real(8) IDT33(3,3),IDT333(3,3,3)
c
            IDT33=0
            do i=1,3
               IDT33(i,i)=1                           
            enddo         
            IDT333=0
            IDT333(1,2,3)=+1
            IDT333(1,3,2)=-1
            IDT333(2,1,3)=-1
            IDT333(2,3,1)=+1
            IDT333(3,1,2)=+1
            IDT333(3,2,1)=-1
c
            Rho_gnd=0
            IVB_gnd=0
            pk2i_gnd=0
c
            Fp0      = fem_Fp0       (iex,igx,:,:    ) 
            dFpp     = fem_dFpp      (iex,igx,:,:    ) 
            dFpp_1gd = fem_dFpp_1gd  (iex,igx,:,:,:  ) 
            dFpp_2gd = fem_dFpp_2gd  (iex,igx,:,:,:,:) 
            Fp=matmul(dFpp,Fp0)

c-----------elastic constants
            Xmu=STFei26(4,4)
            Xnu=0.3

c-----------GND density due to strain gradient (1/mm^2)
            Rho_gnd=0
            do is=1,Nslp
               x1=0
               x2=0
               do i1=1,3
               do i2=1,3
               do i3=1,3
               do i4=1,3
                  x1=x1-dFpp_1gd(i4,i2,i3)*IDT333(i1,i2,i3)
     &                 *vd_slp(is,i4)*vl_slp(is,i1)
                  x2=x2-dFpp_1gd(i4,i2,i3)*IDT333(i1,i2,i3)
     &                 *vd_slp(is,i4)*vd_slp(is,i1)
               enddo
               enddo
               enddo
               enddo
               Rho_gnd(is)=(dabs(x1)+dabs(x2))/B_lattice
            enddo

c-----------passing stress due to strain gradient (MPa)
            IVB_gnd=0
            do is=1,Nslp
               x1=dabs(sum(Rho_gnd(1:Nslp)))
               x2=Xmu*B_lattice*dsqrt(x1)*C_taui 
               IVB_gnd(is)=x2*I_gnd_iso
            enddo

c-----------Internal stress due to strain gradient (1/mm^2)
            pvct=0
            sm_R_gnd=0            
            do i=1,3
            do j=1,3
               if(i/=j)then
                  tvct = IDT33(i,:)               
                  gvct = IDT33(j,:)               
                  pvct1 = gvct*L_size - tvct*L_size/2
                  pvct2 = gvct*L_size + tvct*L_size/2
                  bvct=0
                  do k=1,3
                     x1=0
                     do i1=1,3
                     do i2=1,3
                     do i3=1,3
                     do i4=1,3
                        x1=x1-IDT333(i1,i2,i3)*gvct(i4)
     &                  *tvct(i1)*dFpp_2gd(k,i2,i3,i4)
                     enddo
                     enddo
                     enddo
                     enddo
                     bvct(k)=x1*L_size*L_size*L_size
                  enddo
                  smij=0
                  call stress_dislocation
     &            (pvct1,pvct2,tvct,bvct,Xmu,Xnu,pvct,smij)
                  sm_R_gnd = sm_R_gnd + smij            
               endif
            enddo
            enddo
            sm_I_gnd=matmul(Fp,matmul(sm_R_gnd,transpose(Fp)))
c
            do i=1,6
               x1=sm_I_gnd(IB1(i),IB2(i)) *10
               pk2i_gnd(i) = x1*I_gnd_kin
            enddo

            return
         endsubroutine
c
      endmodule mod_wkcoup_grad


c================================================================
c
c    weak coupling module: Trip
c    output for plastic deformation calculation: Ftrp, IFtrp, IVB_trp
c
c================================================================
      module mod_wkcoup_trip
         use mod_mesh_grid_size
         implicit none
         integer, parameter :: NStrp=12
         integer, parameter :: NTtrp=12
         real(8), parameter :: B_bcc=2.8665d-10
         real(8), parameter :: B_fcc=3.6800d-10
         real(8) fem_cs0    (Tnel,Tngp,6)
         real(8) fem_Fe0    (Tnel,Tngp,3,3)
         real(8) fem_cs     (Tnel,Tngp,6)
         real(8) fem_Fe     (Tnel,Tngp,3,3)
         real(8) fem_gm0    (Tnel,Tngp,NStrp)
         real(8) fem_gm     (Tnel,Tngp,NStrp)
         real(8) fem_Vtrp0  (Tnel,Tngp,NTtrp)
         real(8) fem_Vtrp0_e(Tnel,Tngp,NTtrp)
         real(8) fem_Vtrp0_s(Tnel,Tngp,NTtrp)
         real(8) fem_Vtrp   (Tnel,Tngp,NTtrp)
         real(8) fem_Vtrp_e (Tnel,Tngp,NTtrp)
         real(8) fem_Vtrp_s (Tnel,Tngp,NTtrp)
         real(8) Pm_fcc_twin(NStrp,NStrp)
         real(8) Bm_trp     (NTtrp,3,3)
         real(8) Rm_trp     (NTtrp,3,3)
         real(8) Cab        (NTtrp,2)
         real(8) Msmd       (NTtrp,3,3)
         real(8) V2smd      (NStrp,6)
      contains
c================================================================
c        trip: constants and state variable ini
c================================================================
         subroutine sub_trip_ini(IB1,IB2)
            implicit none
            integer IB1(9),IB2(9)
c            
            integer i,j,k,l,m,i1,i2,i3,i4,ising,j1,j2,j3,j4
            integer k1,k2,k3,k4,ic1,ic2,iex,igx,ip,np,is,js
            real(8) M2X3_1(3,3),x1,x2,x3,vx3_1(3),vx3_2(3),vx3_3(3)
            real(8) mx33_1(3,3),mx33_2(3,3),mx33_3(3,3),mx33_4(3,3)
            real(8) xm33_1(3,3),xm33_2(3,3),vx9(9),vx6(6)
            real(8) IDT33(3,3)
c
            IDT33=0
            IDT33(1,1)=1.d0
            IDT33(2,2)=1.d0
            IDT33(3,3)=1.d0
            
            x1 = B_bcc/B_fcc
            x2 = 10.26d0/180*dacos(-1.d0)

            Bm_trp=0
            Bm_trp(1,1,1)=x1
            Bm_trp(1,2,2)=x1*dsqrt(2.d0)
            Bm_trp(1,3,3)=x1*dsqrt(2.d0)
            Bm_trp(2,1,1)=x1*dsqrt(2.d0)
            Bm_trp(2,2,2)=x1
            Bm_trp(2,3,3)=x1*dsqrt(2.d0)
            Bm_trp(3,1,1)=x1*dsqrt(2.d0)
            Bm_trp(3,2,2)=x1*dsqrt(2.d0)
            Bm_trp(3,3,3)=x1

            Rm_trp=0
            Rm_trp(1,1,:)=[      1.d0,       0.d0,         0.d0]
            Rm_trp(1,2,:)=[      0.d0,  dcos(+x2),   -dsin(+x2)]
            Rm_trp(1,3,:)=[      0.d0,  dsin(+x2),    dcos(+x2)]
            Rm_trp(2,1,:)=[ dcos(+x2),       0.d0,   -dsin(+x2)]
            Rm_trp(2,2,:)=[      0.d0,       1.d0,         0.d0]
            Rm_trp(2,3,:)=[ dsin(+x2),       0.d0,    dcos(+x2)]
            Rm_trp(3,1,:)=[ dcos(+x2),  -dsin(+x2),        0.d0]
            Rm_trp(3,2,:)=[ dsin(+x2),   dcos(+x2),        0.d0]
            Rm_trp(3,3,:)=[      0.d0,        0.d0,        1.d0]
            Rm_trp(4,1,:)=[      1.d0,       0.d0,         0.d0]
            Rm_trp(4,2,:)=[      0.d0,  dcos(-x2),   -dsin(-x2)]
            Rm_trp(4,3,:)=[      0.d0,  dsin(-x2),    dcos(-x2)]
            Rm_trp(5,1,:)=[ dcos(-x2),       0.d0,   -dsin(-x2)]
            Rm_trp(5,2,:)=[      0.d0,       1.d0,         0.d0]
            Rm_trp(5,3,:)=[ dsin(-x2),       0.d0,    dcos(-x2)]
            Rm_trp(6,1,:)=[ dcos(-x2),  -dsin(-x2),        0.d0]
            Rm_trp(6,2,:)=[ dsin(-x2),   dcos(-x2),        0.d0]
            Rm_trp(6,3,:)=[      0.d0,        0.d0,        1.d0]

            Msmd( 1,:,:) = matmul(Rm_trp(2,:,:),Bm_trp(1,:,:))-IDT33
            Msmd( 2,:,:) = matmul(Rm_trp(5,:,:),Bm_trp(1,:,:))-IDT33
            Msmd( 3,:,:) = matmul(Rm_trp(3,:,:),Bm_trp(1,:,:))-IDT33
            Msmd( 4,:,:) = matmul(Rm_trp(6,:,:),Bm_trp(1,:,:))-IDT33
            Msmd( 5,:,:) = matmul(Rm_trp(1,:,:),Bm_trp(2,:,:))-IDT33
            Msmd( 6,:,:) = matmul(Rm_trp(4,:,:),Bm_trp(2,:,:))-IDT33
            Msmd( 7,:,:) = matmul(Rm_trp(3,:,:),Bm_trp(2,:,:))-IDT33
            Msmd( 8,:,:) = matmul(Rm_trp(6,:,:),Bm_trp(2,:,:))-IDT33
            Msmd( 9,:,:) = matmul(Rm_trp(1,:,:),Bm_trp(3,:,:))-IDT33
            Msmd(10,:,:) = matmul(Rm_trp(4,:,:),Bm_trp(3,:,:))-IDT33
            Msmd(11,:,:) = matmul(Rm_trp(2,:,:),Bm_trp(3,:,:))-IDT33
            Msmd(12,:,:) = matmul(Rm_trp(5,:,:),Bm_trp(3,:,:))-IDT33

            do is=1,NTtrp      
               do i=1,9
                  vx9(i)=Msmd(is,IB1(i),IB2(i))
               enddo
               V2smd(is,1:3)=vx9(1:3)
               V2smd(is,4:6)=vx9(4:6)+vx9(7:9)
            enddo                

            Pm_fcc_twin( 1,:) = [ 0,-1,-1,  0, 0, 0,  0, 0, 0,  0, 0, 0]
            Pm_fcc_twin( 2,:) = [+1, 0,-1,  0, 0, 0,  0, 0, 0,  0, 0, 0]
            Pm_fcc_twin( 3,:) = [+1,+1, 0,  0, 0, 0,  0, 0, 0,  0, 0, 0]
            Pm_fcc_twin( 4,:) = [ 0, 0, 0,  0,-1,-1,  0, 0, 0,  0, 0, 0]
            Pm_fcc_twin( 5,:) = [ 0, 0, 0, +1, 0,+1,  0, 0, 0,  0, 0, 0]
            Pm_fcc_twin( 6,:) = [ 0, 0, 0, +1,-1, 0,  0, 0, 0,  0, 0, 0]
            Pm_fcc_twin( 7,:) = [ 0, 0, 0,  0, 0, 0,  0,-1,-1,  0, 0, 0]
            Pm_fcc_twin( 8,:) = [ 0, 0, 0,  0, 0, 0, +1, 0,+1,  0, 0, 0]
            Pm_fcc_twin( 9,:) = [ 0, 0, 0,  0, 0, 0, +1,-1, 0,  0, 0, 0]
            Pm_fcc_twin(10,:) = [ 0, 0, 0,  0, 0, 0,  0, 0, 0,  0,+1,+1]
            Pm_fcc_twin(11,:) = [ 0, 0, 0,  0, 0, 0,  0, 0, 0, -1, 0,+1]
            Pm_fcc_twin(12,:) = [ 0, 0, 0,  0, 0, 0,  0, 0, 0, -1,-1, 0]

            Cab(1:12,1)=[+1, -4, +1,-4,-2, +8,-5, -2, -9,+3, +6,+3 ]
            Cab(1:12,2)=[+7,-10,-10,+7,-5,+11,+8,+11,-12,+6,-12,-9 ]

            do i=1,Tnel
            do j=1,Tngp
               fem_Vtrp0(i,j,:)=0
               fem_Vtrp(i,j,:)=0
               fem_Vtrp0_e(i,j,:)=0
               fem_Vtrp0_s(i,j,:)=0
               fem_Vtrp_e(i,j,:)=0
               fem_Vtrp_s(i,j,:)=0
               fem_cs0(i,j,:)=0
               fem_cs(i,j,:)=0
               fem_Fe0(i,j,:,:)=IDT33
               fem_Fe(i,j,:,:)=IDT33
               fem_gm0(i,j,:)=0
               fem_gm(i,j,:)=0
            enddo
            enddo
            return
         endsubroutine

c================================================================
c        grad: state variable evolution
c================================================================
         subroutine sub_trip_evolution
            implicit none
            integer i,j,k
            real(8) x1,x2,x3
            fem_Fe0=fem_Fe
            fem_cs0=fem_cs
            fem_gm0=fem_gm
            fem_Vtrp0=fem_Vtrp
            fem_Vtrp0_e=fem_Vtrp_e
            fem_Vtrp0_s=fem_Vtrp_s
            return
         endsubroutine 

c================================================================
c        grad: additional eigendeformation and hardening due to trip
c================================================================
         subroutine sub_trip_effect(iex,igx,IB1,IB2,
     &              Nslp_mx,Nslp,STFei26,dt1,
     &              Ftrp, IFtrp, IVB_trp)
            implicit none
            integer Nslp_mx,Nslp,iex,igx
            real(8) STFei26(6,6),dt1
            integer IB1(9),IB2(9)
c
            real(8) Ftrp(3,3)
            real(8) IFtrp(3,3)
            real(8) IVB_trp(Nslp_mx)
c            
            real(8) Vtrp(Nslp_mx)
            real(8) Vtrp_e(Nslp_mx)
            real(8) Vtrp_s(Nslp_mx)
            real(8) Vtrp0(Nslp_mx),TVtrp0,TVtrp
            real(8) Vtrp0_e(Nslp_mx)
            real(8) Vtrp0_s(Nslp_mx)
            integer i,j,k,l,m,i1,i2,i3,i4,ising,j1,j2,j3,j4
            integer k1,k2,k3,k4,ic1,ic2,ip,np,is,js
            real(8) x1,x2,x3,y1,y2
            real(8) M2X3_1(3,3),vx3_1(3),vx3_2(3),vx3_3(3),vx4_1(4)
            real(8) mx33_1(3,3),mx33_2(3,3),mx33_3(3,3),mx33_4(3,3)
            real(8) xm33_1(3,3),xm33_2(3,3)
            real(8) XMcs(3,3),XVcs(6)
            real(8) XMfe(3,3),IXMfe(3,3)
            real(8) XVgam(NStrp),XVtwin(NStrp)
            real(8) Rvf_shbd(NTtrp),Rvf_strs(NTtrp)
            real(8) Rvf_bken(NTtrp),Rvf_psbd(NTtrp)
            real(8) Umech(NTtrp)
            real(8), parameter :: T_akt      = 300.d0
            real(8), parameter :: par1       = 1.d-5
            real(8), parameter :: par2       = 1.d-5
            real(8), parameter :: par3       = 1.d-0
            real(8), parameter :: gam_min    = 1.d-10
            real(8), parameter :: vf_crt     = 1.d-5
            real(8), parameter :: dG_fcc_bcc = 1.d0
            real(8), parameter :: wp_c       = 1.2d0
            real(8), parameter :: wp_mn      = 1.5d0
            real(8), parameter :: wp_si      = 1.5d0
            real(8), parameter :: wp_al      = 0.3d-20
            real(8), parameter :: gam_crt    = 0.0224  !==> max strain for TRIP
            real(8), parameter :: c_nucl  = 1       !==> strain: with(1) nucleation or not(0)
            real(8), parameter :: c_grth  = 0       !==> stress: with(1) growth or not(0)
            real(8), parameter :: Vf0_trp = 1.d-3   !==> inital martensite volme fraction
            real(8), parameter :: c_pben  = 0       !==> PB energy, MPa
            real(8), parameter :: c_pbmb0 = 1.d-1   !==> PB mobility
            real(8), parameter :: c_hard  = 2.d3    !==> hardening coefficient of trip
            real(8), parameter :: B       = 5.0d0   !==> fitting parameter speed
            real(8), parameter :: R       = 8.314d0 !==> gas constant   
            real(8) c_bken                             !==> Bulk energy difference per volume, MPa
            real(8) c_pbmb	                           !==> transformation speed
            real(8) Ms_trp, T0_trp, G_crit, G_chem
            real(8) Xmu,Xnu
            real(8) IDT33(3,3)
c
            IDT33=0
            IDT33(1,1)=1.d0
            IDT33(2,2)=1.d0
            IDT33(3,3)=1.d0
c
            Vtrp=0
            Vtrp_e=0
            Vtrp_s=0
            Ftrp=IDT33
            IFtrp=IDT33
            IVB_trp=0
c
            Xmu=STFei26(4,4)
            Xnu=0.3
c            
            XVgam  = fem_gm0(iex,igx,1:NStrp)
            XVcs   = fem_cs0(iex,igx,:)
            XMfe   = fem_Fe0(iex,igx,:,:)
            XVtwin = matmul(Pm_fcc_twin,XVgam)*dsqrt(3.d0/4)
            do i=1,9
               j=i; if(j>6) j=j-3
               XMcs(IB1(i),IB2(i))=XVcs(j)
            enddo

            Vtrp0  (1:NTtrp) = fem_Vtrp0  (iex,igx,1:NTtrp)  !==> Trip volume fraction
            Vtrp0_e(1:NTtrp) = fem_Vtrp0_e(iex,igx,1:NTtrp)  !==> strain induced trip
            Vtrp0_s(1:NTtrp) = fem_Vtrp0_s(iex,igx,1:NTtrp)  !==> stress induced trip
            TVtrp0  = sum(Vtrp0(1:NTtrp))
            if(TVtrp0>1.d0) TVtrp0=1.d0
               
c-----------------------------------------------------------------------
            Ms_trp = 273.15d0+
     &      (539.d0-423.d0*wp_C-30.4d0*wp_Mn-7.5d0*wp_Si+30d0*wp_Al)
            T0_trp = 655.d0-(wp_C-1.1d0)*325   !494+((1.6-wp_C))*335
            G_crit = 334.d0
            G_chem = G_crit-((T_akt-MS_trp)*G_crit/(T0_trp-Ms_trp))
            c_bken =-G_chem
            c_pbmb = c_pbmb0*dexp(B*7.0*G_chem/(R*T_akt))
c-----------------------------------------------------------------------
         

c-----------Shear band intersection induced martensite nucleation
            do is=1,NTtrp
               i=int(dabs(Cab(is,1)))
               j=int(dabs(Cab(is,2)))
               x1=dsign(1.d0, Cab(is,1))
               x2=dsign(1.d0, Cab(is,2))
               y1=dsign(1.d0, XVtwin(i))
               y2=dsign(1.d0, XVtwin(j))
               if(x1==y1.and.x2==y2)then
                  x3=min(dabs(XVtwin(i)),dabs(XVtwin(j)))/gam_crt
               else
                  x3=0
               endif               
!               Vtrp_e(is)=Vtrp0_e(is)+c_nucl*(1-TVtrp0)*x3
               Vtrp_e(is)=c_nucl*x3*0.1
               if(Vtrp_e(is)<0) Vtrp_e(is)=0
               if(Vtrp_e(is)>1) Vtrp_e(is)=1
            enddo         


c-----------Stress induced martensite growth
            !==>mechanical energy driving part
            Rvf_strs=0
            ising=0
            call gaussj(XMfe,3,IXMfe,ising)
            if(ising==0)then
               do is=1,NTtrp
                  xm33_1=matmul(matmul(XMfe,Msmd(is,:,:)),IXMfe)
                  xm33_2=(xm33_1+transpose(xm33_1))/2
                  x1=0
                  do i=1,3
                  do j=1,3
                     x1=x1 + XMcs(i,j)*xm33_2(j,i)
                  enddo         
                  enddo         
                  if(x1<0) x1=0
                  Rvf_strs(is) = x1
               enddo
            endif
            !==>Bulk energy difference driving part
            Rvf_bken=0
            do is=1,NTtrp
               Rvf_bken(is) = -c_bken * 12*(Vtrp0(is)**2-Vtrp0(is)**3)
            enddo
            !==>Phase boundary effect (double well potential) driving part
            Rvf_psbd=0
            do is=1,NTtrp
               Rvf_psbd(is) = -c_pben * 
     &             (2*Vtrp0(is)-6*Vtrp0(is)**2+4*Vtrp0(is)**3)
            enddo
            !==>finial stress helped growth
            do is=1,NTtrp
               x3 = c_pbmb*(Vtrp0_e(is)+Vf0_trp)
     &            *(Rvf_strs(is)+Rvf_bken(is)+Rvf_psbd(is)) 
               Vtrp_s(is)=Vtrp0_s(is)+c_grth*(1-TVtrp0)*x3*dt1
               if(Vtrp_s(is)<0) Vtrp_s(is)=0
               if(Vtrp_s(is)>1) Vtrp_s(is)=1
            enddo

c-----------Total phase transformation
            Vtrp = Vtrp_e + Vtrp_s
            TVtrp=sum(Vtrp(1:NTtrp))
            if(TVtrp>1)then
               Vtrp_e = Vtrp_e / TVtrp
               Vtrp_s = Vtrp_s / TVtrp
               Vtrp   = Vtrp   / TVtrp
            endif

c-----------Eigen deformation and its inverse matrix of Trip
            Ftrp=IDT33
            do is=1,NTtrp
               Ftrp = Ftrp + Vtrp_s(is)*Msmd(is,:,:)
            enddo
            mx33_1=0
            do is=1,NTtrp
               mx33_1 = mx33_1 + Vtrp_e(is)*Msmd(is,:,:)
            enddo
            x1=(mx33_1(1,1)+mx33_1(2,2)+mx33_1(3,3))/3
            Ftrp=Ftrp+x1*IDT33
            ising=0
            call gaussj(Ftrp,3,IFtrp,ising)
            if(ising/=0) IFtrp=IDT33

c-----------Hardening due to Trip lammelles for slip systems

            IVB_trp(1:Nslp) = B_fcc*Xmu*(dexp(8.8d0*TVtrp0)-1)*c_hard 

c-----------Save data for material point(iex,igx)
            fem_Vtrp_e(iex,igx,:) = Vtrp_e
            fem_Vtrp_s(iex,igx,:) = Vtrp_s
            fem_Vtrp  (iex,igx,:) = Vtrp
c
            return
         endsubroutine

c
      endmodule mod_wkcoup_trip

c================================================================
c
c    weak coupling module: Interal stress (misfit stress)
c    output for plastic deformation calculation:pk2i_int(x,y,z,p) 
c
c================================================================
      module mod_wkcoup_int
         use mod_mesh_grid_size
         implicit none 
         integer ios
         character(len=200) :: Path=''                                     
         real(8) epsR0(24),sigR0(24),matrixA(24,24)
         real(8) epsR(24),sigR(24),epsRP(24),sigRP(24)
         real(8) eps_x0(6),eps_y0(6),eps_z0(6),eps_p0(6)
         real(8) eps_x(6),eps_y(6),eps_z(6),eps_p(6)                               
         real(8) fem_intp(Tnel,Tngp,6),fem_intx(Tnel,Tngp,6)
         real(8) fem_inty(Tnel,Tngp,6),fem_intz(Tnel,Tngp,6)
         real(8) fem_fx(Tnel,Tngp),fem_fy(Tnel,Tngp)
         real(8) fem_fz(Tnel,Tngp),fem_fpp(Tnel,Tngp)
         real(8) fem_Fpx(Tnel,Tngp,3,3),fem_Fpy(Tnel,Tngp,3,3)
         real(8) fem_Fpz(Tnel,Tngp,3,3),fem_Fpx0(Tnel,Tngp,3,3) 
         real(8) fem_Fpy0(Tnel,Tngp,3,3),fem_Fpz0(Tnel,Tngp,3,3)
         real(8) fem_Fppp(Tnel,Tngp,3,3),fem_Fppp0(Tnel,Tngp,3,3) 
         real(8) fem_epsx0(Tnel,Tngp,6),fem_epsy0(Tnel,Tngp,6)
         real(8) fem_epsz0(Tnel,Tngp,6),fem_epsp0(Tnel,Tngp,6)
         real(8) fem_epsx(Tnel,Tngp,6),fem_epsy(Tnel,Tngp,6)
         real(8) fem_epsz(Tnel,Tngp,6),fem_epsp(Tnel,Tngp,6)
      contains
c================================================================
c        int: constants and state variable ini
c================================================================
      subroutine sub_int_ini
         implicit none
         integer i,j,i1,i2
         real(8) mx33(3,3)
                   
         eps_p0 = [-0.0015, -0.0015, -0.0015, 0., 0., 0.]  
         eps_x0 = [0., 0., 0., 0., 0., 0.]
         eps_y0 = [0., 0., 0., 0., 0., 0.]
         eps_z0 = [0., 0., 0., 0., 0., 0.]
         
         epsR0(1:6)=eps_p0(:)
         epsR0(7:12)=eps_x0(:)
         epsR0(13:18)=eps_y0(:)
         epsR0(19:24)=eps_z0(:)

         call getoutdir(Path,ios)
         open(1,file=trim(Path)//'/Matrix.dat',status='old')

         do i=1,24
            read (1,'(100e12.3)') (matrixA(i,j),j=1,24)
         enddo         

         sigR0(:)=matmul(matrixA(:,:),epsR0(:))

         close(1) 

         do i1=1,Tnel
         do i2=1,Tngp
         fem_intp(i1,i2,:)=sigR0(1:6)
         fem_intx(i1,i2,:)=sigR0(7:12)
         fem_inty(i1,i2,:)=sigR0(13:18)
         fem_intz(i1,i2,:)=sigR0(19:24)

         fem_fx(i1,i2)= 0.126
         fem_fy(i1,i2)= 0.126
         fem_fz(i1,i2)= 0.126
         fem_fpp(i1,i2)=0.62

         mx33(1,:)=[1.,0.,0.]
         mx33(2,:)=[0.,1.,0.]
         mx33(3,:)=[0.,0.,1.]
         fem_Fpx0(i1,i2,:,:)=mx33
         fem_Fpy0(i1,i2,:,:)=mx33
         fem_Fpz0(i1,i2,:,:)=mx33
         fem_Fppp0(i1,i2,:,:)=mx33

         fem_epsx0(i1,i2,:)=0
         fem_epsy0(i1,i2,:)=0
         fem_epsz0(i1,i2,:)=0
         fem_epsp0(i1,i2,:)=0
         enddo
         enddo 

         return
         endsubroutine

c================================================================
c        int: state variable evolution
c================================================================
         subroutine sub_int_evolution
            implicit none
            
            fem_epsx0 = fem_epsx
            fem_epsy0 = fem_epsy
            fem_epsz0 = fem_epsz
            fem_epsp0 = fem_epsp


            fem_Fpx0  = fem_Fpx
            fem_Fpy0  = fem_Fpy
            fem_Fpz0  = fem_Fpz
            fem_Fppp0 = fem_Fppp

            return
         endsubroutine 

c================================================================
c        int: internal stress due to misfit eigenstrain
c================================================================
          subroutine sub_int_effect(iex,igx,IB1,IB2,fx,fy,fz,fpp,
     &                          pk2i_intx, pk2i_inty, pk2i_intz,
     &                          pk2i_intp)

         implicit none
         integer IB1(9),IB2(9),iex,igx       
         real(8) pk2i_intp(6),pk2i_intx(6)
         real(8) pk2i_inty(6),pk2i_intz(6)
         real(8) fx,fy,fz,fpp

            eps_p(:) = eps_p0(:) + fem_epsp0(iex,igx,:)        
            eps_x(:) = eps_x0(:) + fem_epsx0(iex,igx,:)
            eps_y(:) = eps_y0(:) + fem_epsy0(iex,igx,:)
            eps_z(:) = eps_z0(:) + fem_epsz0(iex,igx,:)

            epsR(1:6)=eps_p(:)
            epsR(7:12)=eps_x(:)
            epsR(13:18)=eps_y(:)
            epsR(19:24)=eps_z(:) 

            sigR(:)=matmul(matrixA(:,:),epsR(:)) 

         pk2i_intp(:)=sigR(1:6)
         pk2i_intx(:)=sigR(7:12)
         pk2i_inty(:)=sigR(13:18)
         pk2i_intz(:)=sigR(19:24)

         fx=fem_fx(iex,igx)
         fy=fem_fy(iex,igx)
         fz=fem_fz(iex,igx)
         fpp=fem_fpp(iex,igx)

         return
         endsubroutine

      endmodule mod_wkcoup_int

c================================================================
c
c    weak coupling module: super stress
c
c================================================================
      module mod_wkcoup_sup
         use mod_mesh_grid_size
         implicit none                                  
         real(8) fem_gamm0(Tnel,Tngp,48),fem_gamp0(Tnel,Tngp,48)
         real(8) fem_gamm(Tnel,Tngp,48),fem_gamp(Tnel,Tngp,48)
         real(8) fem_cl(Tnel,Tngp,48),fem_mm(Tnel,Tngp,48)
         real(8) fem_kw(Tnel,Tngp,48),fem_gamx(Tnel,Tngp)
         real(8) fem_gamy(Tnel,Tngp),fem_gamz(Tnel,Tngp)
         real(8) fem_chttc(Tnel,Tngp,48)
      contains
c================================================================
c        sup: constants and state variable ini
c================================================================
      subroutine sub_sup_ini
         implicit none
         integer i1,i2

         do i1=1,Tnel
         do i2=1,Tngp

         fem_gamm0(i1,i2,:)=0.d0
         fem_gamp0(i1,i2,:)=0.d0

         fem_gamm(i1,i2,:)=0.d0
         fem_gamp(i1,i2,:)=0.d0

         fem_cl(i1,i2,:)=0.d0
         fem_mm(i1,i2,:)=0.d0
         fem_kw(i1,i2,:)=0.d0

         enddo
         enddo 

         return
         endsubroutine

c================================================================
c        sup: state variable evolution
c================================================================
         subroutine sub_sup_evolution
            implicit none
            integer i1,i2,is
            real(8), parameter :: s1   = 65.d0
            real(8), parameter :: s2   = 8.d0
            real(8), parameter :: a1   = 400.d0
            real(8), parameter :: a2   = 10.d0
            real(8), parameter :: a3   = 0  !800.d0
            real(8), parameter :: a4   = 0  !80.d0
            real(8), parameter :: p1   = 0.5d0
            real(8), parameter :: p2   = 0  !0.4d0

            fem_gamx=0.d0
            fem_gamy=0.d0
            fem_gamz=0.d0
            fem_chttc=0.d0
            fem_cl=0.d0
            fem_mm=0.d0
            fem_kw=0.d0

            do i1=1,Tnel
            do i2=1,Tngp

            do is=1,12
              fem_gamx(i1,i2)=fem_gamx(i1,i2)+fem_gamm(i1,i2,is)
            enddo
            do is=13,24
              fem_gamy(i1,i2)=fem_gamy(i1,i2)+fem_gamm(i1,i2,is)
            enddo
            do is=25,36
              fem_gamz(i1,i2)=fem_gamz(i1,i2)+fem_gamm(i1,i2,is)
            enddo
         
            fem_chttc(i1,i2,37)=fem_gamm(i1,i2,2)+fem_gamm(i1,i2,3)
     &                   +fem_gamm(i1,i2,2+12)+fem_gamm(i1,i2,3+12)
     &                   +fem_gamm(i1,i2,2+24)+fem_gamm(i1,i2,3+24)
            fem_chttc(i1,i2,38)=fem_gamm(i1,i2,1)+fem_gamm(i1,i2,3)
     &                   +fem_gamm(i1,i2,1+12)+fem_gamm(i1,i2,3+12)
     &                   +fem_gamm(i1,i2,1+24)+fem_gamm(i1,i2,3+24)
            fem_chttc(i1,i2,39)=fem_gamm(i1,i2,2)+fem_gamm(i1,i2,1)
     &                   +fem_gamm(i1,i2,2+12)+fem_gamm(i1,i2,1+12)
     &                   +fem_gamm(i1,i2,2+24)+fem_gamm(i1,i2,1+24)
            fem_chttc(i1,i2,40)=fem_gamm(i1,i2,5)+fem_gamm(i1,i2,6)
     &                   +fem_gamm(i1,i2,5+12)+fem_gamm(i1,i2,6+12)
     &                   +fem_gamm(i1,i2,5+24)+fem_gamm(i1,i2,6+24)
            fem_chttc(i1,i2,41)=fem_gamm(i1,i2,4)+fem_gamm(i1,i2,6)
     &                   +fem_gamm(i1,i2,4+12)+fem_gamm(i1,i2,6+12)
     &                   +fem_gamm(i1,i2,4+24)+fem_gamm(i1,i2,6+24)
            fem_chttc(i1,i2,42)=fem_gamm(i1,i2,4)+fem_gamm(i1,i2,5)
     &                   +fem_gamm(i1,i2,4+12)+fem_gamm(i1,i2,5+12)
     &                   +fem_gamm(i1,i2,4+24)+fem_gamm(i1,i2,5+24)
            fem_chttc(i1,i2,43)=fem_gamm(i1,i2,8)+fem_gamm(i1,i2,9)
     &                   +fem_gamm(i1,i2,8+12)+fem_gamm(i1,i2,9+12)
     &                   +fem_gamm(i1,i2,8+24)+fem_gamm(i1,i2,9+24)
            fem_chttc(i1,i2,44)=fem_gamm(i1,i2,7)+fem_gamm(i1,i2,9)
     &                   +fem_gamm(i1,i2,7+12)+fem_gamm(i1,i2,9+12)
     &                   +fem_gamm(i1,i2,7+24)+fem_gamm(i1,i2,9+24)
            fem_chttc(i1,i2,45)=fem_gamm(i1,i2,7)+fem_gamm(i1,i2,8)
     &                   +fem_gamm(i1,i2,7+12)+fem_gamm(i1,i2,8+12)
     &                   +fem_gamm(i1,i2,7+24)+fem_gamm(i1,i2,8+24)
            fem_chttc(i1,i2,46)=fem_gamm(i1,i2,11)+fem_gamm(i1,i2,12)
     &                   +fem_gamm(i1,i2,11+12)+fem_gamm(i1,i2,12+12)
     &                   +fem_gamm(i1,i2,11+24)+fem_gamm(i1,i2,12+24)
            fem_chttc(i1,i2,47)=fem_gamm(i1,i2,12)+fem_gamm(i1,i2,10)
     &                   +fem_gamm(i1,i2,12+12)+fem_gamm(i1,i2,10+12)
     &                   +fem_gamm(i1,i2,12+24)+fem_gamm(i1,i2,10+24)
            fem_chttc(i1,i2,48)=fem_gamm(i1,i2,10)+fem_gamm(i1,i2,11)
     &                   +fem_gamm(i1,i2,10+12)+fem_gamm(i1,i2,11+12)
     &                   +fem_gamm(i1,i2,10+24)+fem_gamm(i1,i2,11+24)
            enddo
            enddo

            do i1=1,Tnel
            do i2=1,Tngp
 
            do is=1,12
                 fem_cl(i1,i2,is)=s1*dtanh(s2*fem_gamx(i1,i2))
            enddo
            do is=13,24
                 fem_cl(i1,i2,is)=s1*dtanh(s2*fem_gamy(i1,i2))
            enddo
            do is=25,36
                 fem_cl(i1,i2,is)=s1*dtanh(s2*fem_gamz(i1,i2))
            enddo
            do is=37,48
                 fem_cl(i1,i2,is)=0.d0
            enddo   

            do is=1,36
               fem_mm(i1,i2,is)=0.d0
               fem_kw(i1,i2,is)=0.d0
            enddo 

            do is=37,48           
            fem_mm(i1,i2,is)=a1*fem_chttc(i1,i2,is)**p1
     &    *dexp(-a2*(fem_gamx(i1,i2)+fem_gamy(i1,i2)+fem_gamz(i1,i2)))
            fem_kw(i1,i2,is)=a3*fem_gamp(i1,i2,is)**p2
     &           *dexp(-a4*fem_gamp(i1,i2,is))
            enddo
             
            enddo
            enddo
            
            
            fem_gamm0=fem_gamm
            fem_gamp0=fem_gamp

            return
         endsubroutine 

c================================================================
c        sup: IVB_cl, IVB_m, IVB_kw
c================================================================
          subroutine sub_sup_effect(iex,igx,IB1,IB2,
     &                              IVB_cl,IVB_m,IVB_kw)

         implicit none
         integer IB1(9),IB2(9),iex,igx,is       
         real(8) IVB_cl(48),IVB_m(48),IVB_kw(48)

         IVB_cl(:)=fem_cl(iex,igx,:)
         IVB_m(:)=fem_mm(iex,igx,:)
         IVB_kw(:)=fem_kw(iex,igx,:)

         return
         endsubroutine

      endmodule mod_wkcoup_sup

c================================================================
c
c    weak coupling module: Back stress 
c
c================================================================
      module mod_wkcoup_bk
         use mod_mesh_grid_size
         use mod_Austenite
         use GlobalValue
         implicit none 

         real(8) fem_bk(Tnel,Tngp,N_slp)
         real(8) fem_bk_ch(Tnel,Tngp,N_slp,3)
         real(8) fem_dbkdt(Tnel,Tngp,N_slp)
         real(8) fem_dbkdt_ch(Tnel,Tngp,N_slp,3)
         real(8) fem_dgmdt(Tnel,Tngp,N_slp)
                            
      contains			
c================================================================	

c================================================================
c        int: constants and state variable ini
c================================================================
      subroutine sub_bk_ini
         implicit none
          
         fem_bk=0         
         fem_dgmdt=0
         fem_dbkdt=0
         fem_bk_ch=0
         fem_dbkdt_ch=0

         return
         endsubroutine

c================================================================
c        int: state variable evolution
c================================================================
         subroutine sub_bk_evolution(dtime)
         implicit none
       
         integer i1,i2,is
         real(8) dtime
c          print*,'sub_bk_evo', Adir, Adyn, A1
c AFKH
          if (Iwkcoup_bk.eq.1) then
          do i1=1,Tnel
          do i2=1,Tngp          
          do is=1,N_slp
           fem_dbkdt(i1,i2,is) = Adir*fem_dgmdt(i1,i2,is)-    
     &     Adyn*fem_bk(i1,i2,is)*dabs(fem_dgmdt(i1,i2,is))

           fem_bk(i1,i2,is)=fem_bk(i1,i2,is)+fem_dbkdt(i1,i2,is)
     &     *dtime
          enddo
          enddo
          enddo
          endif

c CHKH
          if (Iwkcoup_bk.eq.2) then
          do i1=1,Tnel
          do i2=1,Tngp          
          do is=1,N_slp
          fem_dbkdt_ch(i1,i2,is,1) = A1*fem_dgmdt(i1,i2,is)-    
     &     B1*fem_bk_ch(i1,i2,is,1)*dabs(fem_dgmdt(i1,i2,is))

          fem_dbkdt_ch(i1,i2,is,2) = A2*fem_dgmdt(i1,i2,is)-    
     &     B2*fem_bk_ch(i1,i2,is,2)*dabs(fem_dgmdt(i1,i2,is))

          fem_dbkdt_ch(i1,i2,is,3) = A3*fem_dgmdt(i1,i2,is)-    
     &     B3*fem_bk_ch(i1,i2,is,3)*dabs(fem_dgmdt(i1,i2,is))

          fem_bk_ch(i1,i2,is,1:3)=fem_bk_ch(i1,i2,is,1:3)+
     &     fem_dbkdt_ch(i1,i2,is,1:3)*dtime

          fem_bk(i1,i2,is)=fem_bk_ch(i1,i2,is,1)+fem_bk_ch(i1,i2,is,2)
     &     +fem_bk_ch(i1,i2,is,3)
          enddo
          enddo
          enddo
          endif

c OWKH
          if (Iwkcoup_bk.eq.3) then
          do i1=1,Tnel
          do i2=1,Tngp          
          do is=1,N_slp
           fem_dbkdt(i1,i2,is) = Adir*fem_dgmdt(i1,i2,is)-    
     &     Adyn*(dabs(fem_bk(i1,i2,is))/(Adir/Adyn))**M_OW
     &     *fem_bk(i1,i2,is)*dabs(fem_dgmdt(i1,i2,is))

           fem_bk(i1,i2,is)=fem_bk(i1,i2,is)+fem_dbkdt(i1,i2,is)
     &     *dtime
          enddo
          enddo
          enddo
          endif

          return
         endsubroutine 

c================================================================
c        int: back stress due to Bauschinger effect
c================================================================
         subroutine sub_bk_effect(iex,igx,IB1,IB2,IVB_bk)
         implicit none
         integer IB1(9),IB2(9),iex,igx       
         real(8) IVB_bk(N_slp)

         IVB_bk(:)=fem_bk(iex,igx,:)

         return
         endsubroutine

      endmodule mod_wkcoup_bk
c================================================================
c
c    weak coupling module
c
c================================================================
      module mod_wkcoup
         use mod_wkcoup_grad
         use mod_wkcoup_trip
         use mod_wkcoup_int
         use mod_wkcoup_sup
         use mod_wkcoup_bk
         implicit none
      contains
c================================================================
c        state variable ini
c================================================================
         subroutine mod_wkcp_ini(IB1,IB2)
            implicit none
            integer IB1(9),IB2(9)
            call sub_grad_ini
            call sub_trip_ini(IB1,IB2)
            call sub_int_ini
            call sub_sup_ini
            call sub_bk_ini
            return
         endsubroutine
c================================================================
c        state variable evolution
c================================================================
         subroutine wkcoup_evolution(dtime)
            implicit none
            real(8) dtime
            if(Iwkcoup_grad/=0) call sub_grad_evolution
            if(Iwkcoup_trip/=0) call sub_trip_evolution
            if(Iwkcoup_int/=0)  call sub_int_evolution
            if(Iwkcoup_sup/=0)  call sub_sup_evolution
            if(Iwkcoup_bk/=0)   call sub_bk_evolution(dtime)
            return
         endsubroutine

c================================================================
c        wkcoup influence
c================================================================
         subroutine wkcoup_effect(iex,igx,Ialloy,IB1,IB2,
     &              Nslp_mx,Nslp,STFei26,dt1,
     &              vd_slp,vl_slp,vn_slp,
     &              Rho_gnd, IVB_gnd, pk2i_gnd,
     &              Ftrp, IFtrp, IVB_trp,fx,fy,fz,fpp,
     &              pk2i_intx, pk2i_inty, pk2i_intz,
     &              pk2i_intp,IVB_cl,IVB_m,IVB_kw,IVB_bk)

            implicit none
            integer Nslp_mx,Nslp,iex,igx,Ialloy
            real(8) STFei26(6,6),dt1
            real(8) vd_slp(Nslp_mx,3)
            real(8) vl_slp(Nslp_mx,3)
            real(8) vn_slp(Nslp_mx,3)
            integer IB1(9),IB2(9)
c
            real(8) Rho_gnd(Nslp_mx)
            real(8) IVB_gnd(Nslp_mx)
            real(8) pk2i_gnd(6)
c
            real(8) Ftrp(3,3)
            real(8) IFtrp(3,3)
            real(8) IVB_trp(Nslp_mx)

            real(8) fx,fy,fz,fpp
            real(8) pk2i_intp(6)
            real(8) pk2i_intx(6),pk2i_inty(6),pk2i_intz(6)
            real(8) IVB_cl(48),IVB_m(48),IVB_kw(48)
            real(8) IVB_bk(Nslp_mx)

            if(Iwkcoup_grad/=0) call sub_grad_effect(iex,igx,IB1,IB2,
     &                          Nslp_mx,Nslp,STFei26,
     &                          vd_slp,vl_slp,vn_slp,
     &                          Rho_gnd, IVB_gnd, pk2i_gnd)
     
            if(Iwkcoup_trip/=0.and.Ialloy==4) 
     &                          call sub_trip_effect(iex,igx,
     &                          IB1,IB2,Nslp_mx,Nslp,STFei26,dt1,
     &                          Ftrp, IFtrp, IVB_trp)

            if(Iwkcoup_int/=0.and.Ialloy==5) 
     &                          call sub_int_effect(iex,igx,
     &                          IB1,IB2,fx,fy,fz,fpp,
     &                          pk2i_intx, pk2i_inty, pk2i_intz,
     &                          pk2i_intp)
            if(Iwkcoup_sup/=0.and.Ialloy==5) 
     &                          call sub_sup_effect(iex,igx,
     &                          IB1,IB2,IVB_cl,IVB_m,IVB_kw)
            if(Iwkcoup_bk/=0) 
     &                          call sub_bk_effect(iex,igx,
     &                          IB1,IB2,IVB_bk)
            return
         endsubroutine

      endmodule mod_wkcoup






