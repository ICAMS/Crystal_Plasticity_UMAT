c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c +                                                                         +
c +   Material library: #1==>Aluminum                                       +
c +                                                                         +
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      module mod_Aluminum
         use GlobalValue
         implicit none
         integer,parameter:: N_slip = 12
         real(8),parameter:: c11    = 247.d3
         real(8),parameter:: c12    = 147.d3
         real(8),parameter:: c44    = 125.d3 ! c44=2*(c11-c12) iso
         real(8),parameter:: shrt0  = 1.d-6
         real(8),parameter:: pwfl   = 20.d0
         real(8),parameter:: pwhd   = 2.25d0
         real(8),parameter:: crss0  = 20.d0
         real(8),parameter:: crsss  = 1500.d0
         real(8),parameter:: hdrt0  = 60.d0
         real(8),parameter:: c_cpl   = 1.d0
         real(8),parameter:: c_oth   = 1.4d0
         real(8) Mstiff(6,6)
         real(8) HMij(N_slip,N_slip)
         real(8) Dvct(N_slip,3)
         real(8) Lvct(N_slip,3)
         real(8) Nvct(N_slip,3)
         real(8) Msmd(N_slip,3,3)
         real(8) SMsmd(N_slip,3,3)
         real(8) AMsmd(N_slip,3,3)
         real(8) V1smd(N_slip,6)
         real(8) V2smd(N_slip,6)
      contains
c        +----------------------------+
c        +   flow and hardening lows  +
c        +----------------------------+
         subroutine sub_constant_ini_Aluminum(IB1,IB2)
            implicit none
            integer i,j,is,js,IB1(9),IB2(9)
            real(8) x1,x2,x3

            Dvct( 1,:)=[ 0,  1, -1] ; Nvct( 1,:)=[ -1,  1,  1]
            Dvct( 2,:)=[ 1,  0,  1] ; Nvct( 2,:)=[ -1,  1,  1]
            Dvct( 3,:)=[ 1,  1,  0] ; Nvct( 3,:)=[ -1,  1,  1]
            Dvct( 4,:)=[ 0,  1, -1] ; Nvct( 4,:)=[  1,  1,  1]
            Dvct( 5,:)=[ 1,  0, -1] ; Nvct( 5,:)=[  1,  1,  1]
            Dvct( 6,:)=[ 1, -1,  0] ; Nvct( 6,:)=[  1,  1,  1]
            Dvct( 7,:)=[ 0,  1,  1] ; Nvct( 7,:)=[  1,  1, -1]
            Dvct( 8,:)=[ 1,  0,  1] ; Nvct( 8,:)=[  1,  1, -1]
            Dvct( 9,:)=[ 1, -1,  0] ; Nvct( 9,:)=[  1,  1, -1]
            Dvct(10,:)=[ 0,  1,  1] ; Nvct(10,:)=[  1, -1,  1]
            Dvct(11,:)=[ 1,  0, -1] ; Nvct(11,:)=[  1, -1,  1]
            Dvct(12,:)=[ 1,  1,  0] ; Nvct(12,:)=[  1, -1,  1]
c
            Mstiff(:,:)=0
            Mstiff(1,1)=c11
            Mstiff(2,2)=c11
            Mstiff(3,3)=c11
            Mstiff(4,4)=c44*2
            Mstiff(5,5)=c44*2
            Mstiff(6,6)=c44*2
            Mstiff(2,3)=c12
            Mstiff(3,2)=c12
            Mstiff(1,3)=c12
            Mstiff(3,1)=c12
            Mstiff(1,2)=c12
            Mstiff(2,1)=c12
c
            do is=1,N_slip
            Lvct(is,1)=Nvct(is,2)*Dvct(is,3)-Nvct(is,3)*Dvct(is,2)
            Lvct(is,2)=Nvct(is,3)*Dvct(is,1)-Nvct(is,1)*Dvct(is,3)
            Lvct(is,3)=Nvct(is,1)*Dvct(is,2)-Nvct(is,2)*Dvct(is,1)
            x1=dsqrt(sum(Nvct(is,:)**2))
            x2=dsqrt(sum(Dvct(is,:)**2))
            x3=dsqrt(sum(Lvct(is,:)**2))
            Nvct(is,:)=Nvct(is,:)/x1
            Dvct(is,:)=Dvct(is,:)/x2
            Lvct(is,:)=Lvct(is,:)/x3
               do i=1,3
               do j=1,3
                  Msmd(is,i,j)=Dvct(is,i)*Nvct(is,j)
               enddo
               enddo
               SMsmd(is,:,:)=( Msmd(is,:,:)
     &                       +transpose(Msmd(is,:,:)) )/2
               AMsmd(is,:,:)=( Msmd(is,:,:)
     &                       -transpose(Msmd(is,:,:)) )/2
               do i=1,6
                  V1smd(is,i)=
     &            SMsmd(is,ib1(i),ib2(i))
               enddo
               V2smd(is,1:3)=V1smd(is,1:3)
               V2smd(is,4:6)=V1smd(is,4:6)*2
            enddo
c
            do is=1,N_slip
            do js=1,N_slip
               x1=sum(dabs(Nvct(is,:)-Nvct(js,:)))
               if(x1<1.d-10)then
                  HMij(is,js)=c_cpl
               else
                  HMij(is,js)=c_oth
               endif
            enddo
            enddo
c
            return
         endsubroutine

c        c----------------------------c
c        c   get material parameters  c
c        c----------------------------c
         subroutine sub_get_parameter_Aluminum(
     &                   Nslp_mx,
     &                   Nslp,         
     &                   STFei26,
     &                   smdMi,smdSMi,smdAMi,
     &                   smdVi1,
     &                   smdVi2,
     &                   vd_slp,vl_slp,vn_slp,
     &                   refv_pk2i,refv_IVB,
     &                   IVB_ini)
            implicit none
            integer Nslp_mx,Nslp
            real(8) STFei26(6,6)
            real(8) smdMi(Nslp_mx,3,3) 
            real(8) smdSMi(Nslp_mx,3,3) 
            real(8) smdAMi(Nslp_mx,3,3) 
            real(8) smdVi1(Nslp_mx,6) 
            real(8) smdVi2(Nslp_mx,6) 
            real(8) vd_slp(Nslp_mx,3) 
            real(8) vl_slp(Nslp_mx,3) 
            real(8) vn_slp(Nslp_mx,3) 
            real(8) refv_pk2i,refv_IVB
            real(8) IVB_ini(Nslp_mx)
            Nslp               = N_slip
            STFei26            = Mstiff
            smdMi (1:Nslp,:,:) = Msmd
            smdSMi(1:Nslp,:,:) = SMsmd
            smdAMi(1:Nslp,:,:) = AMsmd
            smdVi1(1:Nslp,:  ) = V1smd
            smdVi2(1:Nslp,:  ) = V2smd
            vd_slp(1:Nslp,:  ) = Dvct
            vl_slp(1:Nslp,:  ) = Lvct
            vn_slp(1:Nslp,:  ) = Nvct
            refv_pk2i          = c44*1.d-6
            refv_IVB           = c44*1.d-6
            IVB_ini(1:Nslp)    = crss0 
            return
         endsubroutine

c        c----------------------------c
c        c   flow and hardening lows  c
c        c----------------------------c
         subroutine sub_flow_harden_Aluminum(
     &                 Nslp_mx,Iexp_loc,
     &                 IVB_wcp,IVB_bk,
     &                 tau,IVB,
     &                 dgmdt,
     &                 ddgmdt_dtau,
     &                 ddgmdt_dIVB,
     &                 dIVBdt,
     &                 ddIVBdt_ddgmdt,
     &                 ddIVBdt_dIVB,
     &                 ising)
            implicit none
            integer i,j,is,js,ising
            integer Nslp_mx,Iexp_loc
            real(8) tau(Nslp_mx)
            real(8) IVB(Nslp_mx)
            real(8) IVB_wcp(Nslp_mx)
            real(8) IVB_bk(Nslp_mx)
            real(8) IVB_eff(Nslp_mx)
            real(8) dgmdt(Nslp_mx)
            real(8) ddgmdt_dtau(Nslp_mx)
            real(8) ddgmdt_dIVB(Nslp_mx)
            real(8) dIVBdt(Nslp_mx)
            real(8) ddIVBdt_ddgmdt(Nslp_mx,Nslp_mx)
            real(8) ddIVBdt_dIVB(Nslp_mx,Nslp_mx)
            real(8) x1,x2,x3
c
            ising=0
            dgmdt=0
            ddgmdt_dtau=0
            ddgmdt_dIVB=0
            dIVBdt=0
            ddIVBdt_ddgmdt=0
            ddIVBdt_dIVB=0
c-----------resolved shear stress and resistence
            do is=1,N_slip
               x1=crss0*1.d-10
               x2=crsss
               if(IVB(is)<x1 .or. IVB(is)>x2)then
                  ising=112
                  return
               endif
            enddo
            
c-----------shear rate, derivative of shear rate w.r.t. pk2i,IVB
            do is=1,N_slip
               IVB_eff(is)=IVB(is)+IVB_wcp(is)
               dgmdt(is)=shrt0*(dabs(tau(is)-IVB_bk(is))/IVB_eff(is))
     &                          **pwfl*dsign(1.d0,(tau(is)-IVB_bk(is)))
               if(Iexp_loc/=1)then
                  ddgmdt_dtau(is)=pwfl/IVB_eff(is)*shrt0
     &            *(dabs(tau(is)-IVB_bk(is))/IVB_eff(is))**(pwfl-1)
                  ddgmdt_dIVB(is)=-pwfl*dgmdt(is)/IVB_eff(is)
               endif
            enddo
c--------evolution rate, derivative of evolution rate w.r.t. pk2i,IVB
            do is=1,N_slip
            do js=1,N_slip
               x1=1-IVB(js)/crsss
               dIVBdt(is)=dIVBdt(is) + HMij(is,js)
     &         *hdrt0*dabs(dgmdt(js))*x1**pwhd 
               if(Iexp_loc/=1)then
                  ddIVBdt_ddgmdt(is,js)=HMij(is,js)*hdrt0
     &            *dsign(1.d0,tau(js))*x1**pwhd 
                  ddIVBdt_dIVB(is,js)=HMij(is,js)*hdrt0
     &            *ddgmdt_dIVB(js)*dsign(1.d0,tau(js))*x1**pwhd 
     &            -HMij(is,js)*hdrt0*dabs(dgmdt(js))
     &            *x1**(pwhd-1)*pwhd/crsss
               endif
            enddo
            enddo

            return
         endsubroutine
      endmodule mod_Aluminum
      
      
      
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c +                                                                         +
c +   Material library: #2==>Copper                                         +
c +                                                                         +
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      module mod_Copper
         use GlobalValue
         implicit none
         integer,parameter:: N_slip = 12
         real(8),parameter:: c11    = 247.d3
         real(8),parameter:: c12    = 147.d3
         real(8),parameter:: c44    = 125.d3 ! c44=2*(c11-c12) iso
         real(8),parameter:: shrt0  = 1.d-3
         real(8),parameter:: pwfl   = 20.d0
         real(8),parameter:: pwhd   = 2.25d0
         real(8),parameter:: crss0  = 20.d0
         real(8),parameter:: crsss  = 117.d0
         real(8),parameter:: hdrt0  = 180.d0
         real(8),parameter:: c_cpl   = 1.d0
         real(8),parameter:: c_oth   = 1.4d0
         real(8) Mstiff(6,6)
         real(8) HMij(N_slip,N_slip)
         real(8) Dvct(N_slip,3)
         real(8) Lvct(N_slip,3)
         real(8) Nvct(N_slip,3)
         real(8) Msmd(N_slip,3,3)
         real(8) SMsmd(N_slip,3,3)
         real(8) AMsmd(N_slip,3,3)
         real(8) V1smd(N_slip,6)
         real(8) V2smd(N_slip,6)
      contains
c        +----------------------------+
c        +   flow and hardening lows  +
c        +----------------------------+
         subroutine sub_constant_ini_Copper(IB1,IB2)
            implicit none
            integer i,j,is,js,IB1(9),IB2(9)
            real(8) x1,x2,x3
c
            Dvct( 1,:)=[ 0,  1, -1] ; Nvct( 1,:)=[ -1,  1,  1]
            Dvct( 2,:)=[ 1,  0,  1] ; Nvct( 2,:)=[ -1,  1,  1]
            Dvct( 3,:)=[ 1,  1,  0] ; Nvct( 3,:)=[ -1,  1,  1]
            Dvct( 4,:)=[ 0,  1, -1] ; Nvct( 4,:)=[  1,  1,  1]
            Dvct( 5,:)=[ 1,  0, -1] ; Nvct( 5,:)=[  1,  1,  1]
            Dvct( 6,:)=[ 1, -1,  0] ; Nvct( 6,:)=[  1,  1,  1]
            Dvct( 7,:)=[ 0,  1,  1] ; Nvct( 7,:)=[  1,  1, -1]
            Dvct( 8,:)=[ 1,  0,  1] ; Nvct( 8,:)=[  1,  1, -1]
            Dvct( 9,:)=[ 1, -1,  0] ; Nvct( 9,:)=[  1,  1, -1]
            Dvct(10,:)=[ 0,  1,  1] ; Nvct(10,:)=[  1, -1,  1]
            Dvct(11,:)=[ 1,  0, -1] ; Nvct(11,:)=[  1, -1,  1]
            Dvct(12,:)=[ 1,  1,  0] ; Nvct(12,:)=[  1, -1,  1]
c
            Mstiff(:,:)=0
            Mstiff(1,1)=c11
            Mstiff(2,2)=c11
            Mstiff(3,3)=c11
            Mstiff(4,4)=c44*2
            Mstiff(5,5)=c44*2
            Mstiff(6,6)=c44*2
            Mstiff(2,3)=c12
            Mstiff(3,2)=c12
            Mstiff(1,3)=c12
            Mstiff(3,1)=c12
            Mstiff(1,2)=c12
            Mstiff(2,1)=c12
c
            do is=1,N_slip
            Lvct(is,1)=Nvct(is,2)*Dvct(is,3)-Nvct(is,3)*Dvct(is,2)
            Lvct(is,2)=Nvct(is,3)*Dvct(is,1)-Nvct(is,1)*Dvct(is,3)
            Lvct(is,3)=Nvct(is,1)*Dvct(is,2)-Nvct(is,2)*Dvct(is,1)
            x1=dsqrt(sum(Nvct(is,:)**2))
            x2=dsqrt(sum(Dvct(is,:)**2))
            x3=dsqrt(sum(Lvct(is,:)**2))
            Nvct(is,:)=Nvct(is,:)/x1
            Dvct(is,:)=Dvct(is,:)/x2
            Lvct(is,:)=Lvct(is,:)/x3
               do i=1,3
               do j=1,3
                  Msmd(is,i,j)=Dvct(is,i)*Nvct(is,j)
               enddo
               enddo
               SMsmd(is,:,:)=( Msmd(is,:,:)
     &                       +transpose(Msmd(is,:,:)) )/2
               AMsmd(is,:,:)=( Msmd(is,:,:)
     &                       -transpose(Msmd(is,:,:)) )/2
               do i=1,6
                  V1smd(is,i)=
     &            SMsmd(is,ib1(i),ib2(i))
               enddo
               V2smd(is,1:3)=V1smd(is,1:3)
               V2smd(is,4:6)=V1smd(is,4:6)*2
            enddo
c
            do is=1,N_slip
            do js=1,N_slip
               x1=sum(dabs(Nvct(is,:)-Nvct(js,:)))
               if(x1<1.d-10)then
                  HMij(is,js)=c_cpl
               else
                  HMij(is,js)=c_oth
               endif
            enddo
            enddo
c
            return
         endsubroutine

c        c----------------------------c
c        c   get material parameters  c
c        c----------------------------c
         subroutine sub_get_parameter_Copper(
     &                 Nslp_mx,
     &                 Nslp,       
     &                 STFei26,
     &                 smdMi,smdSMi,smdAMi,
     &                 smdVi1,
     &                 smdVi2,
     &                 vd_slp,vl_slp,vn_slp,
     &                 refv_pk2i,refv_IVB,
     &                 IVB_ini)
            implicit none
            integer Nslp_mx,Nslp
            real(8) STFei26(6,6)
            real(8) smdMi(Nslp_mx,3,3) 
            real(8) smdSMi(Nslp_mx,3,3) 
            real(8) smdAMi(Nslp_mx,3,3) 
            real(8) smdVi1(Nslp_mx,6) 
            real(8) smdVi2(Nslp_mx,6) 
            real(8) vd_slp(Nslp_mx,3) 
            real(8) vl_slp(Nslp_mx,3) 
            real(8) vn_slp(Nslp_mx,3) 
            real(8) refv_pk2i,refv_IVB
            real(8) IVB_ini(Nslp_mx)
            Nslp               = N_slip
            STFei26            = Mstiff
            smdMi (1:Nslp,:,:) = Msmd
            smdSMi(1:Nslp,:,:) = SMsmd
            smdAMi(1:Nslp,:,:) = AMsmd
            smdVi1(1:Nslp,:  ) = V1smd
            smdVi2(1:Nslp,:  ) = V2smd
            vd_slp(1:Nslp,:  ) = Dvct
            vl_slp(1:Nslp,:  ) = Lvct
            vn_slp(1:Nslp,:  ) = Nvct
            refv_pk2i          = c44*1.d-6
            refv_IVB           = c44*1.d-6
            IVB_ini(1:Nslp)    = crss0 
            return
         endsubroutine

c        c----------------------------c
c        c   flow and hardening lows  c
c        c----------------------------c
         subroutine sub_flow_harden_Copper(
     &                 Nslp_mx,Iexp_loc,
     &                 IVB_wcp,IVB_bk,
     &                 tau,IVB,
     &                 dgmdt,
     &                 ddgmdt_dtau,
     &                 ddgmdt_dIVB,
     &                 dIVBdt,
     &                 ddIVBdt_ddgmdt,
     &                 ddIVBdt_dIVB,
     &                 ising)
            implicit none
            integer i,j,is,js,ising
            integer Nslp_mx,Iexp_loc
            real(8) tau(Nslp_mx)
            real(8) IVB(Nslp_mx)
            real(8) IVB_wcp(Nslp_mx)
            real(8) IVB_bk(Nslp_mx)
            real(8) IVB_eff(Nslp_mx)
            real(8) dgmdt(Nslp_mx)
            real(8) ddgmdt_dtau(Nslp_mx)
            real(8) ddgmdt_dIVB(Nslp_mx)
            real(8) dIVBdt(Nslp_mx)
            real(8) ddIVBdt_ddgmdt(Nslp_mx,Nslp_mx)
            real(8) ddIVBdt_dIVB(Nslp_mx,Nslp_mx)
            real(8) x1,x2,x3
c
            ising=0
            dgmdt=0
            ddgmdt_dtau=0
            ddgmdt_dIVB=0
            dIVBdt=0
            ddIVBdt_ddgmdt=0
            ddIVBdt_dIVB=0
c-----------resolved shear stress and resistence
            do is=1,N_slip
               x1=crss0*1.d-10
               x2=crsss
               if(IVB(is)<x1 .or. IVB(is)>x2)then
                  ising=112
                  return
               endif
            enddo

c-----------shear rate, derivative of shear rate w.r.t. pk2i,IVB
            do is=1,N_slip
               IVB_eff(is)=IVB(is)+IVB_wcp(is)
               dgmdt(is)=shrt0*(dabs(tau(is)-IVB_bk(is))/IVB_eff(is))
     &                          **pwfl*dsign(1.d0,(tau(is)-IVB_bk(is)))
               if(Iexp_loc/=1)then
                  ddgmdt_dtau(is)=pwfl/IVB_eff(is)*shrt0
     &            *(dabs(tau(is)-IVB_bk(is))/IVB_eff(is))**(pwfl-1)
                  ddgmdt_dIVB(is)=-pwfl*dgmdt(is)/IVB_eff(is)
               endif
            enddo
c--------evolution rate, derivative of evolution rate w.r.t. pk2i,IVB
            do is=1,N_slip
            do js=1,N_slip
               x1=1-IVB(js)/crsss
               dIVBdt(is)=dIVBdt(is) + HMij(is,js)
     &         *hdrt0*dabs(dgmdt(js))*x1**pwhd 
               if(Iexp_loc/=1)then
                  ddIVBdt_ddgmdt(is,js)=HMij(is,js)*hdrt0
     &            *dsign(1.d0,tau(js))*x1**pwhd 
                  ddIVBdt_dIVB(is,js)=HMij(is,js)*hdrt0
     &            *ddgmdt_dIVB(js)*dsign(1.d0,tau(js))*x1**pwhd 
     &            -HMij(is,js)*hdrt0*dabs(dgmdt(js))
     &            *x1**(pwhd-1)*pwhd/crsss
               endif
            enddo
            enddo

            return
         endsubroutine
      endmodule mod_Copper
      
      
      
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c +                                                                         +
c +   Material library: #3==>Ferrite                                        +
c +                                                                         +
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      module mod_Ferrite
         use GlobalValue
         implicit none
         integer,parameter:: N_slip = 12
         real(8),parameter:: c11    = 247.d3
         real(8),parameter:: c12    = 147.d3
         real(8),parameter:: c44    = 125.d3 ! c44=2*(c11-c12) iso
         real(8),parameter:: shrt0  = 1.d-3
         real(8),parameter:: pwfl   = 20.d0
         real(8),parameter:: pwhd   = 2.25d0
         real(8),parameter:: crss0  = 20.d0
         real(8),parameter:: crsss  = 117.d0
         real(8),parameter:: hdrt0  = 180.d0
         real(8),parameter:: c_cpl   = 1.d0
         real(8),parameter:: c_oth   = 1.4d0
         real(8) Mstiff(6,6)
         real(8) HMij(N_slip,N_slip)
         real(8) Dvct(N_slip,3)
         real(8) Lvct(N_slip,3)
         real(8) Nvct(N_slip,3)
         real(8) Msmd(N_slip,3,3)
         real(8) SMsmd(N_slip,3,3)
         real(8) AMsmd(N_slip,3,3)
         real(8) V1smd(N_slip,6)
         real(8) V2smd(N_slip,6)
      contains
c        +----------------------------+
c        +   flow and hardening lows  +
c        +----------------------------+
         subroutine sub_constant_ini_Ferrite(IB1,IB2)
            implicit none
            integer i,j,is,js,IB1(9),IB2(9)
            real(8) x1,x2,x3

            Dvct( 1,:)=[ -1,  1,  1] ; Nvct( 1,:)=[ 0,  1, -1]
            Dvct( 2,:)=[ -1,  1,  1] ; Nvct( 2,:)=[ 1,  0,  1]
            Dvct( 3,:)=[ -1,  1,  1] ; Nvct( 3,:)=[ 1,  1,  0]
            Dvct( 4,:)=[  1,  1,  1] ; Nvct( 4,:)=[ 0,  1, -1]
            Dvct( 5,:)=[  1,  1,  1] ; Nvct( 5,:)=[ 1,  0, -1]
            Dvct( 6,:)=[  1,  1,  1] ; Nvct( 6,:)=[ 1, -1,  0]
            Dvct( 7,:)=[  1,  1, -1] ; Nvct( 7,:)=[ 0,  1,  1]
            Dvct( 8,:)=[  1,  1, -1] ; Nvct( 8,:)=[ 1,  0,  1]
            Dvct( 9,:)=[  1,  1, -1] ; Nvct( 9,:)=[ 1, -1,  0]
            Dvct(10,:)=[  1, -1,  1] ; Nvct(10,:)=[ 0,  1,  1]
            Dvct(11,:)=[  1, -1,  1] ; Nvct(11,:)=[ 1,  0, -1]
            Dvct(12,:)=[  1, -1,  1] ; Nvct(12,:)=[ 1,  1,  0]
c
            Mstiff(:,:)=0
            Mstiff(1,1)=c11
            Mstiff(2,2)=c11
            Mstiff(3,3)=c11
            Mstiff(4,4)=c44*2
            Mstiff(5,5)=c44*2
            Mstiff(6,6)=c44*2
            Mstiff(2,3)=c12
            Mstiff(3,2)=c12
            Mstiff(1,3)=c12
            Mstiff(3,1)=c12
            Mstiff(1,2)=c12
            Mstiff(2,1)=c12
c
            do is=1,N_slip
            Lvct(is,1)=Nvct(is,2)*Dvct(is,3)-Nvct(is,3)*Dvct(is,2)
            Lvct(is,2)=Nvct(is,3)*Dvct(is,1)-Nvct(is,1)*Dvct(is,3)
            Lvct(is,3)=Nvct(is,1)*Dvct(is,2)-Nvct(is,2)*Dvct(is,1)
            x1=dsqrt(sum(Nvct(is,:)**2))
            x2=dsqrt(sum(Dvct(is,:)**2))
            x3=dsqrt(sum(Lvct(is,:)**2))
            Nvct(is,:)=Nvct(is,:)/x1
            Dvct(is,:)=Dvct(is,:)/x2
            Lvct(is,:)=Lvct(is,:)/x3
               do i=1,3
               do j=1,3
                  Msmd(is,i,j)=Dvct(is,i)*Nvct(is,j)
               enddo
               enddo
               SMsmd(is,:,:)=( Msmd(is,:,:)
     &                       +transpose(Msmd(is,:,:)) )/2
               AMsmd(is,:,:)=( Msmd(is,:,:)
     &                       -transpose(Msmd(is,:,:)) )/2
               do i=1,6
                  V1smd(is,i)=
     &            SMsmd(is,ib1(i),ib2(i))
               enddo
               V2smd(is,1:3)=V1smd(is,1:3)
               V2smd(is,4:6)=V1smd(is,4:6)*2
            enddo
c
            do is=1,N_slip
            do js=1,N_slip
               x1=sum(dabs(Nvct(is,:)-Nvct(js,:)))
               if(x1<1.d-10)then
                  HMij(is,js)=c_cpl
               else
                  HMij(is,js)=c_oth
               endif
            enddo
            enddo
c
            return
         endsubroutine

c        c----------------------------c
c        c   get material parameters  c
c        c----------------------------c
         subroutine sub_get_parameter_Ferrite(
     &                 Nslp_mx,
     &                 Nslp,
     &                 STFei26,
     &                 smdMi,smdSMi,smdAMi,
     &                 smdVi1,
     &                 smdVi2,
     &                 vd_slp,vl_slp,vn_slp,
     &                 refv_pk2i,refv_IVB,
     &                 IVB_ini)
            implicit none
            integer Nslp_mx,Nslp
            real(8) STFei26(6,6)
            real(8) smdMi(Nslp_mx,3,3) 
            real(8) smdSMi(Nslp_mx,3,3) 
            real(8) smdAMi(Nslp_mx,3,3) 
            real(8) smdVi1(Nslp_mx,6) 
            real(8) smdVi2(Nslp_mx,6) 
            real(8) vd_slp(Nslp_mx,3) 
            real(8) vl_slp(Nslp_mx,3) 
            real(8) vn_slp(Nslp_mx,3) 
            real(8) refv_pk2i,refv_IVB
            real(8) IVB_ini(Nslp_mx)
            Nslp               = N_slip
            STFei26            = Mstiff
            smdMi (1:Nslp,:,:) = Msmd
            smdSMi(1:Nslp,:,:) = SMsmd
            smdAMi(1:Nslp,:,:) = AMsmd
            smdVi1(1:Nslp,:  ) = V1smd
            smdVi2(1:Nslp,:  ) = V2smd
            vd_slp(1:Nslp,:  ) = Dvct
            vl_slp(1:Nslp,:  ) = Lvct
            vn_slp(1:Nslp,:  ) = Nvct
            refv_pk2i          = c44*1.d-6
            refv_IVB           = c44*1.d-6
            IVB_ini(1:Nslp)    = crss0 
            return
         endsubroutine

c        c----------------------------c
c        c   flow and hardening lows  c
c        c----------------------------c
         subroutine sub_flow_harden_Ferrite(
     &                 Nslp_mx,Iexp_loc,
     &                 IVB_wcp,IVB_bk,
     &                 tau,IVB,
     &                 dgmdt,
     &                 ddgmdt_dtau,
     &                 ddgmdt_dIVB,
     &                 dIVBdt,
     &                 ddIVBdt_ddgmdt,
     &                 ddIVBdt_dIVB,
     &                 ising)
            implicit none
            integer i,j,is,js,ising
            integer Nslp_mx,Iexp_loc
            real(8) tau(Nslp_mx)
            real(8) IVB(Nslp_mx)
            real(8) IVB_wcp(Nslp_mx)
            real(8) IVB_bk(Nslp_mx)
            real(8) IVB_eff(Nslp_mx)
            real(8) dgmdt(Nslp_mx)
            real(8) ddgmdt_dtau(Nslp_mx)
            real(8) ddgmdt_dIVB(Nslp_mx)
            real(8) dIVBdt(Nslp_mx)
            real(8) ddIVBdt_ddgmdt(Nslp_mx,Nslp_mx)
            real(8) ddIVBdt_dIVB(Nslp_mx,Nslp_mx)
            real(8) x1,x2,x3
c
            ising=0
            dgmdt=0
            ddgmdt_dtau=0
            ddgmdt_dIVB=0
            dIVBdt=0
            ddIVBdt_ddgmdt=0
            ddIVBdt_dIVB=0
c-----------resolved shear stress and resistence
            do is=1,N_slip
               x1=crss0*1.d-10
               x2=crsss
               if(IVB(is)<x1 .or. IVB(is)>x2)then
                  ising=112
                  return
               endif
            enddo

c-----------shear rate, derivative of shear rate w.r.t. pk2i,IVB
            do is=1,N_slip
               IVB_eff(is)=IVB(is)+IVB_wcp(is)
               dgmdt(is)=shrt0*(dabs(tau(is)-IVB_bk(is))/IVB_eff(is))
     &                          **pwfl*dsign(1.d0,(tau(is)-IVB_bk(is)))
               if(Iexp_loc/=1)then
                  ddgmdt_dtau(is)=pwfl/IVB_eff(is)*shrt0
     &            *(dabs(tau(is)-IVB_bk(is))/IVB_eff(is))**(pwfl-1)
                  ddgmdt_dIVB(is)=-pwfl*dgmdt(is)/IVB_eff(is)
               endif
            enddo
c--------evolution rate, derivative of evolution rate w.r.t. pk2i,IVB
            do is=1,N_slip
            do js=1,N_slip
               x1=1-IVB(js)/crsss
               dIVBdt(is)=dIVBdt(is) + HMij(is,js)
     &         *hdrt0*dabs(dgmdt(js))*x1**pwhd 
               if(Iexp_loc/=1)then
                  ddIVBdt_ddgmdt(is,js)=HMij(is,js)*hdrt0
     &            *dsign(1.d0,tau(js))*x1**pwhd 
                  ddIVBdt_dIVB(is,js)=HMij(is,js)*hdrt0
     &            *ddgmdt_dIVB(js)*dsign(1.d0,tau(js))*x1**pwhd 
     &            -HMij(is,js)*hdrt0*dabs(dgmdt(js))
     &            *x1**(pwhd-1)*pwhd/crsss
               endif
            enddo
            enddo

            return
         endsubroutine
      endmodule mod_Ferrite      
      
            
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c +                                                                         +
c +   Material library: #4==>Austenite                                      +
c +                                                                         +
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      module mod_Austenite
         use GlobalValue
         implicit none
         
         integer,parameter:: N_slip = 12
         real(8),parameter:: shrt0  = 1.d-3
         real(8),parameter:: Qact   = 200.d3
         real(8),parameter:: pwfl   = 1250
         real(8),parameter:: pwhd   =0
         real(8),parameter:: crsss  =500.d0
         real(8),parameter:: hdrt0  = 0  
         real(8),parameter:: c_cpl  = 1.d0
         real(8),parameter:: c_oth  = 1.4d0
c        parameters for kinematic hardening
         real(8),parameter:: Adyn = 100.d0
         real(8),parameter:: A1 = 65.d4
         real(8),parameter:: B1 = 499.d2
         real(8),parameter:: A2 = 0.d4
         real(8),parameter:: B2 = 0.d3
         real(8),parameter:: A3 = 0.d4
         real(8),parameter:: B3 = 0.d0
         real(8),parameter:: M_OW = 10.0d0
         real(8),parameter:: N_slp = 48
         real(8) Adir  ! temp-dependent
c
         real(8) Mstiff(6,6)
         real(8) HMij(N_slip,N_slip)
         real(8) Dvct(N_slip,3)
         real(8) Lvct(N_slip,3)
         real(8) Nvct(N_slip,3)
         real(8) Msmd(N_slip,3,3)
         real(8) SMsmd(N_slip,3,3)
         real(8) AMsmd(N_slip,3,3)
         real(8) V1smd(N_slip,6)
         real(8) V2smd(N_slip,6)
         real(8) c11                          ! temp-dependent
         real(8) c12                              ! temp-dependent
         real(8) c44                              ! temp-dependent
         real(8) crss0                            ! temp-dependent 
      
      contains
        
c================================================================           

c================================================================
c    define temperature-dependent paramaters for Austenite
c================================================================           
         subroutine param_temp_Austenite(temp)          
            implicit none           
C           define values for elastic constants and CRSS according to
C           Shahmardani&Hartmaier, Metallurgical and Materials Transaction A (2023)
C           https://doi.org/10.1007/s11661-023-06958-5
            real(8) temp
            if(abs(temp-temp_cur)>1.d-5)then
               print*,'Inconsistent temperatures in param_temp_Aust',temp,temp_cur
            end if

            if(temp_cur>1.d0)then
               c11   = 283.16d3 - 102.5d0*temp_cur
               c12   = 118.08d3 - 42.5d0*temp_cur
               c44   = 82.49d3  - 30.0d0*temp_cur
               crss0 = 292.83d0 - 0.21d0*temp_cur  ! NOTE: Error in publication!
               Adir  = 5092.d0 - 4.d0*temp_cur
            else
C              if no temp is defined use room temperature values
               c11   = 256.5d3
               c12   = 111.1d3
               c44   = 77.2d3
               crss0 = 170.d0
               Adir  = 5000.d0
cc             Values for T=923 K
cc               c11 = 196.06d3   
cc               c12 = 81.99d3    
cc               c44 = 57.01d3    
cc               crss0 = 115.d0
cc               Adir = 1.271d3
            end if
            c11_gl = c11
            c12_gl = c12
            c44_gl = c44
            tau0_gl = crss0
            Adir_gl = Adir
            Mstiff(:,:)=0
            Mstiff(1,1)=c11
            Mstiff(2,2)=c11
            Mstiff(3,3)=c11
            Mstiff(4,4)=c44*2
            Mstiff(5,5)=c44*2
            Mstiff(6,6)=c44*2
            Mstiff(2,3)=c12
            Mstiff(3,2)=c12
            Mstiff(1,3)=c12
            Mstiff(3,1)=c12
            Mstiff(1,2)=c12
            Mstiff(2,1)=c12
         return 
         endsubroutine          
            
c================================================================                                                 
c        +----------------------------+
c        +   flow and hardening lows  +
c        +----------------------------+
         subroutine sub_constant_ini_Austenite(IB1,IB2)
            implicit none
            integer i,j,is,js,IB1(9),IB2(9)
            real(8) x1,x2,x3
            call param_temp_Austenite(temp_cur)
            Dvct( 1,:)=[ 0,  1, -1] ; Nvct( 1,:)=[ -1,  1,  1]
            Dvct( 2,:)=[ 1,  0,  1] ; Nvct( 2,:)=[ -1,  1,  1]
            Dvct( 3,:)=[ 1,  1,  0] ; Nvct( 3,:)=[ -1,  1,  1]
            Dvct( 4,:)=[ 0,  1, -1] ; Nvct( 4,:)=[  1,  1,  1]
            Dvct( 5,:)=[ 1,  0, -1] ; Nvct( 5,:)=[  1,  1,  1]
            Dvct( 6,:)=[ 1, -1,  0] ; Nvct( 6,:)=[  1,  1,  1]
            Dvct( 7,:)=[ 0,  1,  1] ; Nvct( 7,:)=[  1,  1, -1]
            Dvct( 8,:)=[ 1,  0,  1] ; Nvct( 8,:)=[  1,  1, -1]
            Dvct( 9,:)=[ 1, -1,  0] ; Nvct( 9,:)=[  1,  1, -1]
            Dvct(10,:)=[ 0,  1,  1] ; Nvct(10,:)=[  1, -1,  1]
            Dvct(11,:)=[ 1,  0, -1] ; Nvct(11,:)=[  1, -1,  1]
            Dvct(12,:)=[ 1,  1,  0] ; Nvct(12,:)=[  1, -1,  1]

cc            Mstiff(:,:)=0
cc            Mstiff(1,1)=c11
cc            Mstiff(2,2)=c11
cc            Mstiff(3,3)=c11
cc            Mstiff(4,4)=c44*2
cc            Mstiff(5,5)=c44*2
cc            Mstiff(6,6)=c44*2
cc            Mstiff(2,3)=c12
cc            Mstiff(3,2)=c12
cc            Mstiff(1,3)=c12
cc            Mstiff(3,1)=c12
cc            Mstiff(1,2)=c12
cc            Mstiff(2,1)=c12

            do is=1,N_slip
            Lvct(is,1)=Nvct(is,2)*Dvct(is,3)-Nvct(is,3)*Dvct(is,2)
            Lvct(is,2)=Nvct(is,3)*Dvct(is,1)-Nvct(is,1)*Dvct(is,3)
            Lvct(is,3)=Nvct(is,1)*Dvct(is,2)-Nvct(is,2)*Dvct(is,1)
            x1=dsqrt(sum(Nvct(is,:)**2))
            x2=dsqrt(sum(Dvct(is,:)**2))
            x3=dsqrt(sum(Lvct(is,:)**2))
            Nvct(is,:)=Nvct(is,:)/x1
            Dvct(is,:)=Dvct(is,:)/x2
            Lvct(is,:)=Lvct(is,:)/x3
               do i=1,3
               do j=1,3
                  Msmd(is,i,j)=Dvct(is,i)*Nvct(is,j)
               enddo
               enddo
               SMsmd(is,:,:)=( Msmd(is,:,:)
     &                       +transpose(Msmd(is,:,:)) )/2
               AMsmd(is,:,:)=( Msmd(is,:,:)
     &                       -transpose(Msmd(is,:,:)) )/2
               do i=1,6
                  V1smd(is,i)=
     &            SMsmd(is,ib1(i),ib2(i))
               enddo
               V2smd(is,1:3)=V1smd(is,1:3)
               V2smd(is,4:6)=V1smd(is,4:6)*2
            enddo
c
            do is=1,N_slip
            do js=1,N_slip
               x1=sum(dabs(Nvct(is,:)-Nvct(js,:)))
               if(x1<1.d-10)then
                  HMij(is,js)=c_cpl
               else
                  HMij(is,js)=c_oth
               endif
            enddo
            enddo
c
            return
         endsubroutine

c        c----------------------------c
c        c   get material parameters  c
c        c----------------------------c
         subroutine sub_get_parameter_Austenite(
     &                 Nslp_mx,
     &                 Nslp,
     &                 STFei26,
     &                 smdMi,smdSMi,smdAMi,
     &                 smdVi1,
     &                 smdVi2,
     &                 vd_slp,vl_slp,vn_slp,
     &                 refv_pk2i,refv_IVB,
     &                 IVB_ini)
            implicit none
            integer Nslp_mx,Nslp
            real(8) STFei26(6,6)
            real(8) smdMi(Nslp_mx,3,3) 
            real(8) smdSMi(Nslp_mx,3,3) 
            real(8) smdAMi(Nslp_mx,3,3) 
            real(8) smdVi1(Nslp_mx,6) 
            real(8) smdVi2(Nslp_mx,6) 
            real(8) vd_slp(Nslp_mx,3) 
            real(8) vl_slp(Nslp_mx,3) 
            real(8) vn_slp(Nslp_mx,3) 
            real(8) refv_pk2i,refv_IVB
            real(8) IVB_ini(Nslp_mx)
c            print*,'sub_get_param: ',crss0, c11, Adir
            call param_temp_Austenite(temp_cur)
            Nslp               = N_slip
            STFei26            = Mstiff
            smdMi (1:Nslp,:,:) = Msmd
            smdSMi(1:Nslp,:,:) = SMsmd
            smdAMi(1:Nslp,:,:) = AMsmd
            smdVi1(1:Nslp,:  ) = V1smd
            smdVi2(1:Nslp,:  ) = V2smd
            vd_slp(1:Nslp,:  ) = Dvct
            vl_slp(1:Nslp,:  ) = Lvct
            vn_slp(1:Nslp,:  ) = Nvct
            refv_pk2i          = c44*1.d-6
            refv_IVB           = c44*1.d-6
            IVB_ini(1:Nslp)    = crss0 

            return
         endsubroutine

c        c----------------------------c
c        c   flow and hardening lows  c
c        c----------------------------c
         subroutine sub_flow_harden_Austenite(
     &                 Nslp_mx,Iexp_loc,
     &                 IVB_wcp,IVB_bk,
     &                 tau,IVB,
     &                 dgmdt,
     &                 ddgmdt_dtau,
     &                 ddgmdt_dIVB,
     &                 dIVBdt,
     &                 ddIVBdt_ddgmdt,
     &                 ddIVBdt_dIVB,
     &                 ising)
            implicit none
            integer i,j,is,js,ising
            integer Nslp_mx,Iexp_loc
            real(8) tau(Nslp_mx)
            real(8) IVB(Nslp_mx)
            real(8) IVB_wcp(Nslp_mx)
            real(8) IVB_bk(Nslp_mx)
            real(8) IVB_eff(Nslp_mx)
            real(8) dgmdt(Nslp_mx)
            real(8) ddgmdt_dtau(Nslp_mx)
            real(8) ddgmdt_dIVB(Nslp_mx)
            real(8) dIVBdt(Nslp_mx)
            real(8) ddIVBdt_ddgmdt(Nslp_mx,Nslp_mx)
            real(8) ddIVBdt_dIVB(Nslp_mx,Nslp_mx)
            real(8) x1,x2,x3
c
            ising=0
            dgmdt=0
            ddgmdt_dtau=0
            ddgmdt_dIVB=0
            dIVBdt=0
            ddIVBdt_ddgmdt=0
            ddIVBdt_dIVB=0
c-----------resolved shear stress and resistence
            do is=1,N_slip
               x1=crss0*1.d-10
               x2=crsss
               if(IVB(is)<x1 .or. IVB(is)>x2)then
                  ising=112
                  return
               endif
            enddo

c-----------shear rate, derivative of shear rate w.r.t. pk2i,IVB
            do is=1,N_slip
               IVB_eff(is)=IVB(is)+IVB_wcp(is)
               dgmdt(is)=shrt0*dexp(-Qact/R_gas/temp_cur)*(dabs(tau(is)-IVB_bk(is))/
     &                   IVB_eff(is))**pwfl*dsign(1.d0,(tau(is)-IVB_bk(is)))
               if(Iexp_loc/=1)then
                  ddgmdt_dtau(is)=pwfl/IVB_eff(is)*shrt0
     &            *(dabs(tau(is)-IVB_bk(is))/IVB_eff(is))**(pwfl-1)*
     &            dexp(-Qact/R_gas/temp_cur)
                  ddgmdt_dIVB(is)=-pwfl*dgmdt(is)/IVB_eff(is)
               endif
            enddo
c--------evolution rate, derivative of evolution rate w.r.t. pk2i,IVB
            do is=1,N_slip
            do js=1,N_slip
               x1=1-IVB(js)/crsss
               dIVBdt(is)=dIVBdt(is) + HMij(is,js)
     &         *hdrt0*dabs(dgmdt(js))*x1**pwhd 
               if(Iexp_loc/=1)then
                  ddIVBdt_ddgmdt(is,js)=HMij(is,js)*hdrt0
     &            *dsign(1.d0,tau(js))*x1**pwhd 
                  ddIVBdt_dIVB(is,js)=HMij(is,js)*hdrt0
     &            *ddgmdt_dIVB(js)*dsign(1.d0,tau(js))*x1**pwhd 
     &            -HMij(is,js)*hdrt0*dabs(dgmdt(js))
     &            *x1**(pwhd-1)*pwhd/crsss
               endif
            enddo
            enddo

            return
         endsubroutine
      endmodule mod_Austenite  
      
      

      
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c +                                                                         +
c +   Material library: #5==>Superalloy                 +
c +                                                                         +
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      module mod_Superalloy
         use GlobalValue
         implicit none
         integer,parameter:: N_slip = 60
         real(8),parameter:: c11    = 193.7d3
         real(8),parameter:: c12    = 135.3d3
         real(8),parameter:: c44    = 93.5d3     
         real(8),parameter:: shrt0m  = 5.0d9
         real(8),parameter:: shrt0p  = 2.0d15 
         real(8),parameter:: shrt0c  = 6.0d12 
         real(8),parameter:: Qactm   = 350.d3
         real(8),parameter:: Qactp   = 500.d3
         real(8),parameter:: Qactc   = 430.d3  
         real(8),parameter:: pwfl   = 10.5d0  !orginal:5.5d0
         real(8),parameter:: pwhd   = 0.05d0         
         real(8),parameter:: crssm0  = 60.d0   !orginal:60.d0
         real(8),parameter:: crssp0  = 500.d0  !orginal:500.d0
         real(8),parameter:: crssc0  = 200.d0
         real(8),parameter:: crss_oro  = 60.d0
         real(8),parameter:: crsss  = 800.d0   !orginal:800.d0
         real(8),parameter:: hdrt0  = 60.d0   !orginal:60.d0
         real(8),parameter:: c_cpl   = 1.d0
         real(8),parameter:: c_oth   = 1.4d0
         real(8),parameter:: c_cplc   = 5.d0
         real(8),parameter:: c_othc   = 10.4d0
         real(8) Mstiff(6,6)
         real(8) HMij(N_slip,N_slip)
         real(8) Dvct(N_slip,3)
         real(8) Lvct(N_slip,3)
         real(8) Nvct(N_slip,3)
         real(8) Msmd(N_slip,3,3)
         real(8) SMsmd(N_slip,3,3)
         real(8) AMsmd(N_slip,3,3)
         real(8) V1smd(N_slip,6)
         real(8) V2smd(N_slip,6)
      contains
c        +----------------------------+
c        +   flow and hardening lows  +
c        +----------------------------+
         subroutine sub_constant_ini_Superalloy(IB1,IB2)
            implicit none
            integer i,j,is,js,IB1(9),IB2(9),ii,ia,ib
            real(8) x1,x2,x3
c
            Dvct( 1,:)=[ 0,  1, -1] ; Nvct( 1,:)=[ -1,  1,  1]
            Dvct( 2,:)=[ 1,  0,  1] ; Nvct( 2,:)=[ -1,  1,  1]
            Dvct( 3,:)=[ 1,  1,  0] ; Nvct( 3,:)=[ -1,  1,  1]
            Dvct( 4,:)=[ 0,  1, -1] ; Nvct( 4,:)=[  1,  1,  1]
            Dvct( 5,:)=[ 1,  0, -1] ; Nvct( 5,:)=[  1,  1,  1]
            Dvct( 6,:)=[ 1, -1,  0] ; Nvct( 6,:)=[  1,  1,  1]
            Dvct( 7,:)=[ 0,  1,  1] ; Nvct( 7,:)=[  1,  1, -1]
            Dvct( 8,:)=[ 1,  0,  1] ; Nvct( 8,:)=[  1,  1, -1]
            Dvct( 9,:)=[ 1, -1,  0] ; Nvct( 9,:)=[  1,  1, -1]
            Dvct(10,:)=[ 0,  1,  1] ; Nvct(10,:)=[  1, -1,  1]
            Dvct(11,:)=[ 1,  0, -1] ; Nvct(11,:)=[  1, -1,  1]
            Dvct(12,:)=[ 1,  1,  0] ; Nvct(12,:)=[  1, -1,  1]
            Dvct(13:24,:)=Dvct( 1:12,:)
            Nvct(13:24,:)=Nvct( 1:12,:)
            Dvct(25:36,:)=Dvct( 1:12,:)
            Nvct(25:36,:)=Nvct( 1:12,:)

            Dvct(37,:)=[ 2,  1,  1]
        Dvct(38,:)=[-1, -2,  1]
        Dvct(39,:)=[-1,  1, -2]
        Dvct(40,:)=[-2,  1,  1]
        Dvct(41,:)=[ 1, -2,  1]
        Dvct(42,:)=[ 1,  1, -2]
        Dvct(43,:)=[-2,  1, -1]
        Dvct(44,:)=[ 1, -2, -1]
        Dvct(45,:)=[ 1,  1,  2]
        Dvct(46,:)=[-2, -1,  1]
        Dvct(47,:)=[ 1,  2,  1]
        Dvct(48,:)=[ 1, -1, -2]
c            Dvct(37:48,:)=Dvct( 1:12,:)
            Nvct(37:48,:)=Nvct( 1:12,:)

            Dvct(49,:)=[-1,  1,  0]; Nvct(49,:)=[0, 0, 1]
            Dvct(50,:)=[ 1,  1,  0]; Nvct(50,:)=[0, 0, 1]
            Dvct(51,:)=[-1,  0,  1]; Nvct(51,:)=[0, 1, 0]
            Dvct(52,:)=[ 1,  0,  1]; Nvct(52,:)=[0, 1, 0]
            Dvct(53,:)=[ 0, -1,  1]; Nvct(53,:)=[1, 0, 0]
            Dvct(54,:)=[ 0,  1,  1]; Nvct(54,:)=[1, 0, 0]
            Dvct(55:60,:)=Dvct(49:54,:)
            Nvct(55:60,:)=Nvct(49:54,:)


            Mstiff(:,:)=0
            Mstiff(1,1)=c11
            Mstiff(2,2)=c11
            Mstiff(3,3)=c11
            Mstiff(4,4)=c44*2
            Mstiff(5,5)=c44*2
            Mstiff(6,6)=c44*2
            Mstiff(2,3)=c12
            Mstiff(3,2)=c12
            Mstiff(1,3)=c12
            Mstiff(3,1)=c12
            Mstiff(1,2)=c12
            Mstiff(2,1)=c12

            do is=1,N_slip
            Lvct(is,1)=Nvct(is,2)*Dvct(is,3)-Nvct(is,3)*Dvct(is,2)
            Lvct(is,2)=Nvct(is,3)*Dvct(is,1)-Nvct(is,1)*Dvct(is,3)
            Lvct(is,3)=Nvct(is,1)*Dvct(is,2)-Nvct(is,2)*Dvct(is,1)
            x1=dsqrt(sum(Nvct(is,:)**2))
            x2=dsqrt(sum(Dvct(is,:)**2))
            x3=dsqrt(sum(Lvct(is,:)**2))
            Nvct(is,:)=Nvct(is,:)/x1
            Dvct(is,:)=Dvct(is,:)/x2
            Lvct(is,:)=Lvct(is,:)/x3
               do i=1,3
               do j=1,3
                  Msmd(is,i,j)=Dvct(is,i)*Nvct(is,j)
               enddo
               enddo
               SMsmd(is,:,:)=( Msmd(is,:,:)
     &                       +transpose(Msmd(is,:,:)) )/2
               AMsmd(is,:,:)=( Msmd(is,:,:)
     &                       -transpose(Msmd(is,:,:)) )/2
               do i=1,6
                  V1smd(is,i)=
     &            SMsmd(is,ib1(i),ib2(i))
               enddo
               V2smd(is,1:3)=V1smd(is,1:3)
               V2smd(is,4:6)=V1smd(is,4:6)*2
            enddo

            do ii=1,4
               ia=1+12*(ii-1)
               ib=12+12*(ii-1)
            do is=ia,ib
            do js=ia,ib
               x1=sum(dabs(Nvct(is,:)-Nvct(js,:)))
               if(x1<1.d-10)then
                  HMij(is,js)=c_cpl
               else
                  HMij(is,js)=c_oth
               endif
            enddo
            enddo
            enddo

            do ii=1,2
               ia=49+6*(ii-1)
               ib=54+6*(ii-1)
            do is=ia,ib
            do js=ia,ib
               x1=sum(dabs(Nvct(is,:)-Nvct(js,:)))
               if(x1<1.d-10)then
                  HMij(is,js)=c_cplc
               else
                  HMij(is,js)=c_othc
               endif
            enddo
            enddo
            enddo

            return
         endsubroutine

c        c----------------------------c
c        c   get material parameters  c
c        c----------------------------c
         subroutine sub_get_parameter_Superalloy(
     &                 Nslp_mx,
     &                 Nslp,
     &                 STFei26,
     &                 smdMi,smdSMi,smdAMi,
     &                 smdVi1,
     &                 smdVi2,
     &                 vd_slp,vl_slp,vn_slp,
     &                 refv_pk2i,refv_IVB,
     &                 IVB_ini)
            implicit none
            integer Nslp_mx,Nslp
            real(8) STFei26(6,6)
            real(8) smdMi(Nslp_mx,3,3) 
            real(8) smdSMi(Nslp_mx,3,3) 
            real(8) smdAMi(Nslp_mx,3,3) 
            real(8) smdVi1(Nslp_mx,6) 
            real(8) smdVi2(Nslp_mx,6) 
            real(8) vd_slp(Nslp_mx,3) 
            real(8) vl_slp(Nslp_mx,3) 
            real(8) vn_slp(Nslp_mx,3) 
            real(8) refv_pk2i,refv_IVB
            real(8) IVB_ini(Nslp_mx) 
            Nslp               = N_slip
            STFei26            = Mstiff
            smdMi (1:Nslp,:,:) = Msmd
            smdSMi(1:Nslp,:,:) = SMsmd
            smdAMi(1:Nslp,:,:) = AMsmd
            smdVi1(1:Nslp,:  ) = V1smd
            smdVi2(1:Nslp,:  ) = V2smd
            vd_slp(1:Nslp,:  ) = Dvct
            vl_slp(1:Nslp,:  ) = Lvct
            vn_slp(1:Nslp,:  ) = Nvct
            refv_pk2i          = c44*1.d-6
            refv_IVB           = c44*1.d-6
            IVB_ini(1:36)      = crssm0 
            IVB_ini(37:48)     = crssp0 
            IVB_ini(49:60)     = crssc0 
            return
         endsubroutine

c        c----------------------------c
c        c   flow and hardening lows  c
c        c----------------------------c
         subroutine sub_flow_harden_Superalloy(
     &                 Nslp_mx,Iexp_loc,
     &                 IVB_cl,IVB_m,IVB_kw,IVB_bk,
     &                 tau,IVB,
     &                 dgmdt,
     &                 ddgmdt_dtau,
     &                 ddgmdt_dIVB,
     &                 dIVBdt,
     &                 ddIVBdt_ddgmdt,
     &                 ddIVBdt_dIVB,
     &                 ising)
            implicit none
            integer i,j,is,js,ising,ii,ia,ib
            integer Nslp_mx,Iexp_loc
            real(8) tau(Nslp_mx)
            real(8) IVB(Nslp_mx)
            real(8) IVB_cl(48)
            real(8) IVB_m(48)            
            real(8) IVB_kw(48) 
            real(8) IVB_bk(Nslp_mx)           
            real(8) IVB_eff(Nslp_mx)
            real(8) dgmdt(Nslp_mx)
            real(8) ddgmdt_dtau(Nslp_mx)
            real(8) ddgmdt_dIVB(Nslp_mx)
            real(8) dIVBdt(Nslp_mx)
            real(8) ddIVBdt_ddgmdt(Nslp_mx,Nslp_mx)
            real(8) ddIVBdt_dIVB(Nslp_mx,Nslp_mx)
            real(8) x1,x2,x3
c
            ising=0
            dgmdt=0
            ddgmdt_dtau=0
            ddgmdt_dIVB=0
            dIVBdt=0
            ddIVBdt_ddgmdt=0
            ddIVBdt_dIVB=0
c-----------resolved shear stress and resistence
            do is=1,36
               x1=crssm0*1.d-10
               x2=crsss
               if(IVB(is)<x1 .or. IVB(is)>x2)then
                  ising=112
                  return
               endif
            enddo

            do is=37,48
               x1=crssp0*1.d-10
               x2=crsss
               if(IVB(is)<x1 .or. IVB(is)>x2)then
                  ising=112
                  return
               endif
            enddo

            do is=49,60
               x1=crssc0*1.d-10
               x2=crsss
               if(IVB(is)<x1 .or. IVB(is)>x2)then
                  ising=112
                  return
               endif
            enddo

c-----------shear rate, derivative of shear rate w.r.t. pk2i,IVB
            do is=1,36
               IVB_eff(is)=IVB(is)+crss_oro-IVB_cl(is)          

               dgmdt(is)=shrt0m*dexp(-Qactm/R_gas/temp_cur)
     &          *(dabs(tau(is)-IVB_bk(is))/IVB_eff(is))
     &          **pwfl*dsign(1.d0,(tau(is)-IVB_bk(is)))
               if(Iexp_loc/=1)then
                  ddgmdt_dtau(is)=pwfl/IVB_eff(is)*shrt0m
     &            *(dabs(tau(is)-IVB_bk(is))/IVB_eff(is))**(pwfl-1)
     &            *dexp(-Qactm/R_gas/temp_cur)
                  ddgmdt_dIVB(is)=-pwfl*dgmdt(is)/IVB_eff(is)
               endif
            enddo

            do is=37,48
               IVB_eff(is)=IVB(is)-IVB_m(is)!+IVB_kw(is)          

               dgmdt(is)=shrt0p*dexp(-Qactp/R_gas/temp_cur)
     &          *(dabs(tau(is)-IVB_bk(is))/IVB_eff(is))
     &          **pwfl*dsign(1.d0,(tau(is)-IVB_bk(is)))
               if(Iexp_loc/=1)then
                  ddgmdt_dtau(is)=pwfl/IVB_eff(is)*shrt0p
     &            *(dabs(tau(is)-IVB_bk(is))/IVB_eff(is))**(pwfl-1)
     &            *dexp(-Qactp/R_gas/temp_cur)
                  ddgmdt_dIVB(is)=-pwfl*dgmdt(is)/IVB_eff(is)
               endif

            enddo

            do is=49,60
               IVB_eff(is)=IVB(is)         
               dgmdt(is)=shrt0c*dexp(-Qactc/R_gas/temp_cur)
     &            *(dabs(tau(is))/IVB_eff(is))**pwfl*dsign(1.d0,tau(is))
               if(Iexp_loc/=1)then
               ddgmdt_dtau(is)=pwfl/IVB_eff(is)*shrt0c
     &       *(dabs(tau(is))/IVB_eff(is))**(pwfl-1)*dexp(-Qactc/R_gas/temp_cur)
               ddgmdt_dIVB(is)=-pwfl*dgmdt(is)/IVB_eff(is)
               endif
            enddo

c--------evolution rate, derivative of evolution rate w.r.t. pk2i,IVB
            do ii=1,4
               ia=1+12*(ii-1)
               ib=12+12*(ii-1)
            do is=ia,ib
            do js=ia,ib
               x1=1-IVB(js)/crsss
               dIVBdt(is)=dIVBdt(is) + HMij(is,js)
     &         *hdrt0*dabs(dgmdt(js))*x1**pwhd 
               if(Iexp_loc/=1)then
                  ddIVBdt_ddgmdt(is,js)=HMij(is,js)*hdrt0
     &            *dsign(1.d0,tau(js))*x1**pwhd 
                  ddIVBdt_dIVB(is,js)=HMij(is,js)*hdrt0
     &            *ddgmdt_dIVB(js)*dsign(1.d0,tau(js))*x1**pwhd 
     &            -HMij(is,js)*hdrt0*dabs(dgmdt(js))
     &            *x1**(pwhd-1)*pwhd/crsss
               endif
            enddo
            enddo
            enddo

            do ii=1,2
               ia=49+6*(ii-1)
               ib=54+6*(ii-1)
            do is=ia,ib
            do js=ia,ib
               x1=1-IVB(js)/crsss
               dIVBdt(is)=dIVBdt(is) + HMij(is,js)
     &         *hdrt0*dabs(dgmdt(js))*x1**pwhd 
               if(Iexp_loc/=1)then
                  ddIVBdt_ddgmdt(is,js)=HMij(is,js)*hdrt0
     &            *dsign(1.d0,tau(js))*x1**pwhd 
                  ddIVBdt_dIVB(is,js)=HMij(is,js)*hdrt0
     &            *ddgmdt_dIVB(js)*dsign(1.d0,tau(js))*x1**pwhd 
     &            -HMij(is,js)*hdrt0*dabs(dgmdt(js))
     &            *x1**(pwhd-1)*pwhd/crsss
               endif
            enddo
            enddo
            enddo

            return
         endsubroutine
      endmodule mod_Superalloy  

c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c +                                                                         +
c +   Material library: #6==>Nickel                                         +
c +                                                                         +
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      module mod_Nickel
         use GlobalValue
         implicit none
         integer,parameter:: N_slip = 12
         real(8),parameter:: c11    = 176.d3
         real(8),parameter:: c12    = 112.d3
         real(8),parameter:: c44    = 64.d3
         real(8),parameter:: shrt0  = 5.d1
         real(8),parameter:: Qact   = 96.d3
         real(8),parameter:: pwfl   = 3.0d0
         real(8),parameter:: pwhd   = 0.05d0
         real(8),parameter:: crss0  = 50.d0
         real(8),parameter:: crsss  = 1500.d0
         real(8),parameter:: hdrt0  = 950.d0
         real(8),parameter:: c_cpl   = 1.d0
         real(8),parameter:: c_oth   = 1.4d0
         real(8) Mstiff(6,6)
         real(8) HMij(N_slip,N_slip)
         real(8) Dvct(N_slip,3)
         real(8) Lvct(N_slip,3)
         real(8) Nvct(N_slip,3)
         real(8) Msmd(N_slip,3,3)
         real(8) SMsmd(N_slip,3,3)
         real(8) AMsmd(N_slip,3,3)
         real(8) V1smd(N_slip,6)
         real(8) V2smd(N_slip,6)
      contains
c        +----------------------------+
c        +   flow and hardening lows  +
c        +----------------------------+
         subroutine sub_constant_ini_Nickel(IB1,IB2)
            implicit none
            integer i,j,is,js,IB1(9),IB2(9)
            real(8) x1,x2,x3
c
            Dvct( 1,:)=[ 0,  1, -1] ; Nvct( 1,:)=[ -1,  1,  1]
            Dvct( 2,:)=[ 1,  0,  1] ; Nvct( 2,:)=[ -1,  1,  1]
            Dvct( 3,:)=[ 1,  1,  0] ; Nvct( 3,:)=[ -1,  1,  1]
            Dvct( 4,:)=[ 0,  1, -1] ; Nvct( 4,:)=[  1,  1,  1]
            Dvct( 5,:)=[ 1,  0, -1] ; Nvct( 5,:)=[  1,  1,  1]
            Dvct( 6,:)=[ 1, -1,  0] ; Nvct( 6,:)=[  1,  1,  1]
            Dvct( 7,:)=[ 0,  1,  1] ; Nvct( 7,:)=[  1,  1, -1]
            Dvct( 8,:)=[ 1,  0,  1] ; Nvct( 8,:)=[  1,  1, -1]
            Dvct( 9,:)=[ 1, -1,  0] ; Nvct( 9,:)=[  1,  1, -1]
            Dvct(10,:)=[ 0,  1,  1] ; Nvct(10,:)=[  1, -1,  1]
            Dvct(11,:)=[ 1,  0, -1] ; Nvct(11,:)=[  1, -1,  1]
            Dvct(12,:)=[ 1,  1,  0] ; Nvct(12,:)=[  1, -1,  1]


            Mstiff(:,:)=0
            Mstiff(1,1)=c11
            Mstiff(2,2)=c11
            Mstiff(3,3)=c11
            Mstiff(4,4)=c44*2
            Mstiff(5,5)=c44*2
            Mstiff(6,6)=c44*2
            Mstiff(2,3)=c12
            Mstiff(3,2)=c12
            Mstiff(1,3)=c12
            Mstiff(3,1)=c12
            Mstiff(1,2)=c12
            Mstiff(2,1)=c12

            do is=1,N_slip
            Lvct(is,1)=Nvct(is,2)*Dvct(is,3)-Nvct(is,3)*Dvct(is,2)
            Lvct(is,2)=Nvct(is,3)*Dvct(is,1)-Nvct(is,1)*Dvct(is,3)
            Lvct(is,3)=Nvct(is,1)*Dvct(is,2)-Nvct(is,2)*Dvct(is,1)
            x1=dsqrt(sum(Nvct(is,:)**2))
            x2=dsqrt(sum(Dvct(is,:)**2))
            x3=dsqrt(sum(Lvct(is,:)**2))
            Nvct(is,:)=Nvct(is,:)/x1
            Dvct(is,:)=Dvct(is,:)/x2
            Lvct(is,:)=Lvct(is,:)/x3
               do i=1,3
               do j=1,3
                  Msmd(is,i,j)=Dvct(is,i)*Nvct(is,j)
               enddo
               enddo
               SMsmd(is,:,:)=( Msmd(is,:,:)
     &                       +transpose(Msmd(is,:,:)) )/2
               AMsmd(is,:,:)=( Msmd(is,:,:)
     &                       -transpose(Msmd(is,:,:)) )/2
               do i=1,6
                  V1smd(is,i)=
     &            SMsmd(is,ib1(i),ib2(i))
               enddo
               V2smd(is,1:3)=V1smd(is,1:3)
               V2smd(is,4:6)=V1smd(is,4:6)*2
            enddo
c
            do is=1,12
            do js=1,12
               x1=sum(dabs(Nvct(is,:)-Nvct(js,:)))
               if(x1<1.d-10)then
                  HMij(is,js)=c_cpl
               else
                  HMij(is,js)=c_oth
               endif
            enddo
            enddo

            return
         endsubroutine

c        c----------------------------c
c        c   get material parameters  c
c        c----------------------------c
         subroutine sub_get_parameter_Nickel(
     &                 Nslp_mx,
     &                 Nslp,
     &                 STFei26,
     &                 smdMi,smdSMi,smdAMi,
     &                 smdVi1,
     &                 smdVi2,
     &                 vd_slp,vl_slp,vn_slp,
     &                 refv_pk2i,refv_IVB,
     &                 IVB_ini)
            implicit none
            integer Nslp_mx,Nslp
            real(8) STFei26(6,6)
            real(8) smdMi(Nslp_mx,3,3) 
            real(8) smdSMi(Nslp_mx,3,3) 
            real(8) smdAMi(Nslp_mx,3,3) 
            real(8) smdVi1(Nslp_mx,6) 
            real(8) smdVi2(Nslp_mx,6) 
            real(8) vd_slp(Nslp_mx,3) 
            real(8) vl_slp(Nslp_mx,3) 
            real(8) vn_slp(Nslp_mx,3) 
            real(8) refv_pk2i,refv_IVB
            real(8) IVB_ini(Nslp_mx)
            Nslp               = N_slip
            STFei26            = Mstiff
            smdMi (1:Nslp,:,:) = Msmd
            smdSMi(1:Nslp,:,:) = SMsmd
            smdAMi(1:Nslp,:,:) = AMsmd
            smdVi1(1:Nslp,:  ) = V1smd
            smdVi2(1:Nslp,:  ) = V2smd
            vd_slp(1:Nslp,:  ) = Dvct
            vl_slp(1:Nslp,:  ) = Lvct
            vn_slp(1:Nslp,:  ) = Nvct
            refv_pk2i          = c44*1.d-6
            refv_IVB           = c44*1.d-6
            IVB_ini            = crss0 
            return
         endsubroutine

c        c----------------------------c
c        c   flow and hardening lows  c
c        c----------------------------c
         subroutine sub_flow_harden_Nickel(
     &                 Nslp_mx,Iexp_loc,
     &                 IVB_wcp,IVB_bk,
     &                 tau,IVB,
     &                 dgmdt,
     &                 ddgmdt_dtau,
     &                 ddgmdt_dIVB,
     &                 dIVBdt,
     &                 ddIVBdt_ddgmdt,
     &                 ddIVBdt_dIVB,
     &                 ising)
            implicit none
            integer i,j,is,js,ising
            integer Nslp_mx,Iexp_loc
            real(8) tau(Nslp_mx)
            real(8) IVB(Nslp_mx)
            real(8) IVB_wcp(Nslp_mx)
            real(8) IVB_eff(Nslp_mx)
            real(8) IVB_bk(Nslp_mx)
            real(8) dgmdt(Nslp_mx)
            real(8) ddgmdt_dtau(Nslp_mx)
            real(8) ddgmdt_dIVB(Nslp_mx)
            real(8) dIVBdt(Nslp_mx)
            real(8) ddIVBdt_ddgmdt(Nslp_mx,Nslp_mx)
            real(8) ddIVBdt_dIVB(Nslp_mx,Nslp_mx)
            real(8) x1,x2,x3
c
            ising=0
            dgmdt=0
            ddgmdt_dtau=0
            ddgmdt_dIVB=0
            dIVBdt=0
            ddIVBdt_ddgmdt=0
            ddIVBdt_dIVB=0
c-----------resolved shear stress and resistence
            do is=1,N_slip
               x1=crss0*1.d-10
               x2=crsss
               if(IVB(is)<x1 .or. IVB(is)>x2)then
                  ising=112
                  return
               endif
            enddo

c-----------shear rate, derivative of shear rate w.r.t. pk2i,IVB
            do is=1,N_slip
               IVB_eff(is)=IVB(is)
               dgmdt(is)=shrt0*exp(-Qact/R_gas/temp_cur)
     &      *(dabs(tau(is))/IVB_eff(is))**pwfl*dsign(1.d0,tau(is))
               if(Iexp_loc/=1)then
             ddgmdt_dtau(is)=pwfl/IVB_eff(is)*shrt0*(dabs(tau(is))
     &                 /IVB_eff(is))**(pwfl-1)*exp(-Qact/R_gas/temp_cur)
             ddgmdt_dIVB(is)=-pwfl*dgmdt(is)/IVB_eff(is)
               endif
            enddo


c--------evolution rate, derivative of evolution rate w.r.t. pk2i,IVB
            do is=1,N_slip
            do js=1,N_slip
               x1=1-IVB(js)/crsss
               dIVBdt(is)=dIVBdt(is) + HMij(is,js)
     &         *hdrt0*dabs(dgmdt(js))*x1**pwhd 
               if(Iexp_loc/=1)then
                  ddIVBdt_ddgmdt(is,js)=HMij(is,js)*hdrt0
     &            *dsign(1.d0,tau(js))*x1**pwhd 
                  ddIVBdt_dIVB(is,js)=HMij(is,js)*hdrt0
     &            *ddgmdt_dIVB(js)*dsign(1.d0,tau(js))*x1**pwhd 
     &            -HMij(is,js)*hdrt0*dabs(dgmdt(js))
     &            *x1**(pwhd-1)*pwhd/crsss
               endif
            enddo
            enddo

            return
         endsubroutine
      endmodule mod_Nickel      

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                               +
! +   Module of total materials                                    +
! +                                                               +
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      module mod_alloys
         use mod_Aluminum
         use mod_Copper
         use mod_Ferrite
         use mod_Austenite
         use mod_Superalloy
         use mod_Nickel
         implicit none
      contains
         ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         ! +   Material parameter initialization                           +
         ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         subroutine mod_alls_ini(IB1,IB2)
            implicit none
            integer IB1(9),IB2(9)

            call sub_constant_ini_Aluminum(IB1,IB2)
            call sub_constant_ini_Copper(IB1,IB2)
            call sub_constant_ini_Ferrite(IB1,IB2)
            call sub_constant_ini_Austenite(IB1,IB2)
            call sub_constant_ini_Superalloy(IB1,IB2)
            call sub_constant_ini_Nickel(IB1,IB2)
            return
         endsubroutine

         ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         ! +   Set temperature-dependent material parameters               +
         ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         subroutine set_param_temp(temp)
            implicit none
            real(8) temp

            call param_temp_Austenite(temp)
            return
         end subroutine set_param_temp

         ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         ! +   Get material parameters for for (ie, ig)                    +
         ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         subroutine get_mat_parameter(ialloy,
     &              Nslp_mx,Nslp,STFei26,
     &              smdMi,smdSMi,smdAMi,smdVi1,smdVi2,
     &              vd_slp,vl_slp,vn_slp,
     &              refv_pk2i,refv_IVB,IVB_ini)
            implicit none
            integer Nslp_mx,Nslp,ialloy
            real(8) STFei26(6,6)
            real(8) smdMi(Nslp_mx,3,3) 
            real(8) smdSMi(Nslp_mx,3,3) 
            real(8) smdAMi(Nslp_mx,3,3) 
            real(8) smdVi1(Nslp_mx,6) 
            real(8) smdVi2(Nslp_mx,6) 
            real(8) vd_slp(Nslp_mx,3) 
            real(8) vl_slp(Nslp_mx,3) 
            real(8) vn_slp(Nslp_mx,3) 
            real(8) refv_pk2i,refv_IVB
            real(8) IVB_ini(Nslp_mx)

            if(ialloy==1)then 
               call sub_get_parameter_Aluminum(
     &              Nslp_mx,Nslp,STFei26,
     &              smdMi,smdSMi,smdAMi,smdVi1,smdVi2,
     &              vd_slp,vl_slp,vn_slp,
     &              refv_pk2i,refv_IVB,IVB_ini)
            elseif(ialloy==2)then 
               call sub_get_parameter_Copper(
     &              Nslp_mx,Nslp,STFei26,
     &              smdMi,smdSMi,smdAMi,smdVi1,smdVi2,
     &              vd_slp,vl_slp,vn_slp,
     &              refv_pk2i,refv_IVB,IVB_ini)
            elseif(ialloy==3)then 
               call sub_get_parameter_Ferrite(
     &              Nslp_mx,Nslp,STFei26,
     &              smdMi,smdSMi,smdAMi,smdVi1,smdVi2,
     &              vd_slp,vl_slp,vn_slp,
     &              refv_pk2i,refv_IVB,IVB_ini)
            elseif(ialloy==4)then 
               call sub_get_parameter_Austenite(
     &              Nslp_mx,Nslp,STFei26,
     &              smdMi,smdSMi,smdAMi,smdVi1,smdVi2,
     &              vd_slp,vl_slp,vn_slp,
     &              refv_pk2i,refv_IVB,IVB_ini)
            elseif(ialloy==5)then 
               call sub_get_parameter_Superalloy(
     &              Nslp_mx,Nslp,STFei26,
     &              smdMi,smdSMi,smdAMi,smdVi1,smdVi2,
     &              vd_slp,vl_slp,vn_slp,
     &              refv_pk2i,refv_IVB,IVB_ini)
            elseif(ialloy==6)then 
               call sub_get_parameter_Nickel(
     &              Nslp_mx,Nslp,STFei26,
     &              smdMi,smdSMi,smdAMi,smdVi1,smdVi2,
     &              vd_slp,vl_slp,vn_slp,
     &              refv_pk2i,refv_IVB,IVB_ini)
            else
               print*,'no this material',ialloy; 
               call xit
            endif         
            return
         endsubroutine
         ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         ! +   Define material behaviour for (ie, ig)                      +
         ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         subroutine cal_flow_harden(ialloy,
     &              Nslp_mx,Iexp_loc,IVB_wcp,IVB_bk,IVB_cl,IVB_m,IVB_kw,tau,
     &              IVB,dgmdt,ddgmdt_dtau,ddgmdt_dIVB,
     &              dIVBdt,ddIVBdt_ddgmdt,ddIVBdt_dIVB,ising)
            implicit none
            integer i,j,is,js,ising,ialloy
            integer Nslp_mx,Iexp_loc
            real(8) tau(Nslp_mx)
            real(8) IVB(Nslp_mx)
            real(8) IVB_wcp(Nslp_mx)
            real(8) IVB_bk(Nslp_mx)
            real(8) IVB_eff(Nslp_mx)
            real(8) IVB_cl(48)
            real(8) IVB_m(48)
            real(8) IVB_kw(48)
            real(8) dgmdt(Nslp_mx)
            real(8) ddgmdt_dtau(Nslp_mx)
            real(8) ddgmdt_dIVB(Nslp_mx)
            real(8) dIVBdt(Nslp_mx)
            real(8) ddIVBdt_ddgmdt(Nslp_mx,Nslp_mx)
            real(8) ddIVBdt_dIVB(Nslp_mx,Nslp_mx)

            if(ialloy==1)then
               call sub_flow_harden_Aluminum(
     &              Nslp_mx,Iexp_loc,IVB_wcp,IVB_bk,tau,IVB,
     &              dgmdt,ddgmdt_dtau,ddgmdt_dIVB,
     &              dIVBdt,ddIVBdt_ddgmdt,ddIVBdt_dIVB,ising)
            elseif(ialloy==2)then
               call sub_flow_harden_Copper(
     &              Nslp_mx,Iexp_loc,IVB_wcp,IVB_bk,tau,IVB,
     &              dgmdt,ddgmdt_dtau,ddgmdt_dIVB,
     &              dIVBdt,ddIVBdt_ddgmdt,ddIVBdt_dIVB,ising)
            elseif(ialloy==3)then
               call sub_flow_harden_Ferrite(
     &              Nslp_mx,Iexp_loc,IVB_wcp,IVB_bk,tau,IVB,
     &              dgmdt,ddgmdt_dtau,ddgmdt_dIVB,
     &              dIVBdt,ddIVBdt_ddgmdt,ddIVBdt_dIVB,ising)
            elseif(ialloy==4)then
               call sub_flow_harden_Austenite(
     &              Nslp_mx,Iexp_loc,IVB_wcp,IVB_bk,tau,IVB,
     &              dgmdt,ddgmdt_dtau,ddgmdt_dIVB,
     &              dIVBdt,ddIVBdt_ddgmdt,ddIVBdt_dIVB,ising)
            elseif(ialloy==5)then
               call sub_flow_harden_Superalloy(
     &              Nslp_mx,Iexp_loc,IVB_cl,IVB_m,IVB_kw,IVB_bk,
     &              tau,IVB,dgmdt,ddgmdt_dtau,ddgmdt_dIVB,
     &              dIVBdt,ddIVBdt_ddgmdt,ddIVBdt_dIVB,ising)
            elseif(ialloy==6)then
               call sub_flow_harden_Nickel(
     &              Nslp_mx,Iexp_loc,IVB_wcp,tau,IVB,IVB_bk,
     &              dgmdt,ddgmdt_dtau,ddgmdt_dIVB,
     &              dIVBdt,ddIVBdt_ddgmdt,ddIVBdt_dIVB,ising)
            else
               print*,'Unknown material ',ialloy
               call xit
            endif         
            return
         endsubroutine
      endmodule mod_alloys

