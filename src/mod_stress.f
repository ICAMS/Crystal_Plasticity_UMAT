c ICAMS CP-UMAT 2026R1
c (c) 2026 by ICAMS, Ruhr University Bochum
c================================================================
c
c    Modules: stress
c    Subroutines to calculate stress, IVB, and material stiffness 
c    at current configuration
c
c================================================================
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                               +
! +   Module "mod_stress" contribute stress calculation           +
! +                                                               +
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      module mod_stress
         use mod_gaussp
         use mod_material
         use mod_wkcoup
         implicit none
         integer Icurrent_dt
         integer :: Nnr_max=200                     
         real(8) :: toler_NRloop=1.d-10         
      contains
         !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         !+                                                                           +
         !+      This subroutine calculates pk2i, IVB by newton raphson algoriths     +
         !+                                                                           +
         !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         recursive subroutine cal_stress_add(ising)
            implicit none
            integer i,j,k,l,m,n,i1,j1,k1,l1,m1,n1,is,js
            integer iNRloop
            integer ising_1,ising
            real(8) x1,x2,x3,x4,y1,z1
            real(8) M1_66(6,6),M1_96(9,6),M1_99(9,9)
            real(8) M2_66(6,6),M2_96(9,6),M2_99(9,9),M3_99(9,9)
            real(8) M1_3333(3,3,3,3),M2_3333(3,3,3,3)
            real(8) MX1(3,3),MX2(3,3),MX3(3,3),MX4(3,3)
            real(8) Qm0(3,3),Qm(3,3),TQm0(3,3),TQm(3,3)
            real(8) vx6_1(6),vx6_2(6)
            real(8) mx43_1(3,3,3,3),mx43_2(3,3,3,3)
            real(8) mxn9_1(48,9),mxn9_2(48,9),mxn33(48,3,3)
            real(8) vx3_1(3),vx3_2(3)
            !1------------------------------------------------------------------1
            !1   calcuate stiffness and schmid matrix at current configuration  1
            !1------------------------------------------------------------------1

            ising=0
            call icams_Eang2Q(eang0(1),eang0(2),eang0(3),Qm0)
            TQm0=transpose(Qm0)
            STFei29(1:6,1:3)=STFei26(1:6,1:3)
            STFei29(1:6,4:6)=STFei26(1:6,4:6)/2
            STFei29(1:6,7:9)=STFei29(1:6,4:6)
            STFei29(7:9,:)=STFei29(4:6,:)
            STFec29=0
            do i=1,9
            do j=1,9
              do i1=1,9
              do j1=1,9
                 STFec29(i,j)=STFec29(i,j)+STFei29(i1,j1)
     &           *Qm0(ib1(i1),ib1(i))*Qm0(ib2(i1),ib2(i))
     &           *Qm0(ib1(j1),ib1(j))*Qm0(ib2(j1),ib2(j))
              enddo
              enddo
              STFec43(ib1(i),ib2(i),ib1(j),ib2(j))=STFec29(i,j)
            enddo
            enddo
            STFec26(1:6,1:3)=STFec29(1:6,1:3)
            STFec26(1:6,4:6)=STFec29(1:6,4:6)*2
            x1=0
            do i=1,3            
            do j=1,3            
               x1=x1+dabs(Fg(i,j)-XI33(i,j))
            enddo
            enddo
            if(x1<1.d-10)then
               cs=0
               matJacb=STFec29(1:6,1:6)
               return
            endif            
            !1------------------------------------------------------------------1
            !1   calcuate Lg, D, W, and so on                                   1
            !1------------------------------------------------------------------1

            ising_1=0
            call icams_determ(Fg,det_Fg)
            call gaussj(Fg,3,IFg,ising_1)
            call gaussj(Fg0,3,IFg0,ising_1)

            if(ising_1/=0)then
               write(6,*) 'non ivertable for Fg'
               ising=3  !Fg is non-invertible,stop current time step
               return
            endif
c            Lg=(XI33-matmul(Fg0,IFg))/dt1
c            Lg=(matmul(Fg,IFg0)-XI33)/dt1
            dLg=matmul(Fg,IFg0)-XI33
c            Lg=Lg_glb

            dDmc=(dLg+transpose(dLg))/2
            dWmc=(dLg-transpose(dLg))/2
            do i=1,6            
               dDvc(i)=dDmc(ib1(i),ib2(i))
            enddo
c            do i=1,9            
c               j=i; if(i>6) j=i-3
c               csM0(ib1(i),ib2(i))=cs0(j)
c            enddo
            !1------------------------------------------------------------------1
            !1   calcuate schmid matrix at current configuration                1
            !1------------------------------------------------------------------1
            do is=1,N_slip
               smdMc(is,:,:)=matmul(matmul(Qm0,smdMi(is,:,:)),TQm0)
               smdSMc(is,:,:)=(smdMc(is,:,:)+smdMc(is,:,:))/2
               smdAMc(is,:,:)=(smdMc(is,:,:)-smdMc(is,:,:))/2
               do i=1,3
                  j=i+3
                  smdVc1(is,i)=smdSMc(is,ib1(i),ib2(i))
                  smdVc1(is,j)=smdSMc(is,ib1(j),ib2(j))
                  smdVc2(is,i)=smdSMc(is,ib1(i),ib2(i))
                  smdVc2(is,j)=smdSMc(is,ib1(j),ib2(j))*2
               enddo
               mx2=matmul(smdAMc(is,:,:),csM0)
               vx6_2=matmul(STFec26,smdVc1(is,:))
c               vx6_2=matmul(STFec26,smdVc2(is,:))
               do i=1,6
                  vx6_1(i)=mx2(ib1(i),ib2(i))+mx2(ib2(i),ib1(i))
               enddo
               STFrlx6(is,:)=vx6_1+vx6_2
            enddo
            mx2=matmul(dWmc,csM0)
            do i=1,6
               vx6_1(i)=mx2(ib1(i),ib2(i))+mx2(ib2(i),ib1(i))
            enddo
            dcs_max=matmul(STFec26,dDvc)-cs0*sum(dDvc(1:3))+vx6_1

            !1------------------------------------------------------1
            !1    Begin newton raphson method to solve cs, IVB      1
            !1------------------------------------------------------1
            do iNRloop=1,Nnr_max
               ising_1=0
               do is=1,N_slip
                  tau(is)=dot_product(cs,smdVc2(is,:))
               enddo
               if (Iwkcoup_sup==0 .and. Iwkcoup_temp==0) then
                  call sub_flow_harden_std(
     &                 Iexp_loc,N_slip,crss0,crsss,hdrt0,pwhd,
     &                 shrt0,pwfl,HMij,IVB_wcp,IVB_bk,tau,
     &                 IVB,dgmdt,ddgmdt_dtau,ddgmdt_dIVB,
     &                 dIVBdt,ddIVBdt_ddgmdt,ddIVBdt_dIVB,ising_1)
               else
                  call sub_flow_harden(
     &                 Iexp_loc,IVB_wcp,IVB_bk,IVB_cl,IVB_m,IVB_kw,
     &                 tau,IVB,dgmdt,ddgmdt_dtau,ddgmdt_dIVB,
     &                 dIVBdt,ddIVBdt_ddgmdt,ddIVBdt_dIVB,ising_1)
               endif
               do is=1,N_slip
                  ddgmdt_dpk2i(is,:)=ddgmdt_dtau(is)*smdVc2(is,:)
               enddo
               ddIVBdt_dpk2i(1:N_slip,:)=
     &         matmul(ddIVBdt_ddgmdt(1:N_slip,1:N_slip), 
     &         ddgmdt_dpk2i(1:N_slip,:))

               if(ising_1/=0)then
                  ising=ising_1  !!! shearrate non-cvg,stop
                  return
               endif
               Gv1=cs-cs0-dcs_max  
     &         +matmul(transpose(STFrlx6(1:N_slip,:)),dgmdt(1:N_slip))*dt1
               Gv2(1:N_slip)=IVB(1:N_slip)-IVB0(1:N_slip)-dIVBdt(1:N_slip)*dt1
               x1=sum(dabs(Gv1))/refv_pk2i
               x2=sum(dabs(Gv2(1:N_slip)))/refv_IVB
      
c               call pm(Fg0,3,3)     
c               call pm(Fg,3,3)     
c               call pm(Fe0,3,3)     
c               call pm(Fe,3,3)     
c               print '(i5,3e14.4)',iNRloop,x1,x2,toler_NRloop

               if(x1+x2<toler_NRloop .or. Iexp_loc==1)then
                  cs=cs0+dcs_max
     &            -matmul(transpose(STFrlx6(1:N_slip,:)),
     &            dgmdt(1:N_slip))*dt1
                  IVB(1:N_slip)=IVB0(1:N_slip)+dIVBdt(1:N_slip)*dt1
                  goto 101 !Converge jump out NR loop and continue
               endif
               ising_1=0
               call cal_NRmatrix(ising_1)
               if(ising_1/=0)then
                  ising=ising_1  !NR matrix cannot be got,stop
                  return
               endif

               dcs=matmul(IeqM66Gv1,
     &         -(Gv1-matmul(eqM6nGv1(:,1:N_slip),Gv2(1:N_slip))))
               dIVB(1:N_slip)=matmul( IeqMnnGv2(1:N_slip,1:N_slip),
     &         -(Gv2(1:N_slip)-matmul(eqMn6Gv2(1:N_slip,:),Gv1)) )

               cs = cs + dcs
               IVB(1:N_slip) = IVB(1:N_slip) + dIVB(1:N_slip) 

            enddo
            ising=5  !Non-coverge for Nnr_max, stop current time step
            return
101         continue
            !1--------------------------------------------------------1
            !1   End of newton raphson method loop                    1
            !1--------------------------------------------------------1
            do i=1,3
            do j=1,3
               Lp(i,j)=dot_product(dgmdt(1:N_slip),smdMc(1:N_slip,i,j))
            enddo
            enddo
            detGM=dgmdt*dt1
            Qm=matmul(XI33+dLg-Lp*dt1,Qm0)
            call caleulang(Qm,eang(1:3),ising_1)
            if(ising_1/=0) eang(1:3)=eang0(1:3)
            call icams_misori(eang00(1:3),eang(1:3),eang(4))

            !1----------------------------------1
            !1    calculate material stiffness  1
            !1----------------------------------1
            if(Iexp_abq==0)then
               ising_1=0
               call cal_NRmatrix(ising_1)
               if(ising_1/=0)then
                  ising=420+ising_1  !NR matrix cannot be got,stop
                  return
               endif
               dIVB_dpk2i=-matmul(IdGv2_dIVB,dGv2_dpk2i)
               do is=1,N_slip
                  vx6_1=ddgmdt_dpk2i(is,:)
     &                 +ddgmdt_dIVB(is)*dIVB_dpk2i(is,:)
                  do i=1,9
                     j=i; if(i>6) j=i-3
                     mxn33(is,ib1(i),ib2(i))=vx6_1(j)
                  enddo                  
               enddo
               mx43_1=0
               mx43_2=0
               do m=1,3
               do n=1,3
               do i=1,3
               do j=1,3
                  mx43_1(m,n,i,j)=XI33(i,m)*XI33(j,n)
                  do is=1,N_slip
                     mx43_1(m,n,i,j)=mx43_1(m,n,i,j)+
     &               mxn33(is,m,n)*STFrlx33(i,j,is)
                  enddo
                  mx43_2(m,n,i,j)=STFec43(m,n,i,j)-csM0(m,n)*XI33(i,j)
               enddo
               enddo
               enddo
               enddo
               do i=1,9
               do j=1,9
                  M1_99(i,j)=mx43_1(ib1(i),ib2(i),ib1(j),ib2(j))
                  M2_99(i,j)=mx43_2(ib1(i),ib2(i),ib1(j),ib2(j))
               enddo
               enddo
               call gaussj(M1_99,9,M3_99,ising)
               if(ising/=0)then
                  matJacb=matJacb0
                  return
               endif
               M1_99=matmul(M3_99,M2_99)
               matJacb=M1_99(1:6,1:6)
            endif

            return
         endsubroutine cal_stress_add

         !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         !+                                                                           +
         !+      This subroutine calculates pk2i, IVB by newton raphson algoriths     +
         !+                                                                           +
         !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         logical function cal_stress_mul_ctx_supported(mp)
            implicit none
            type(mat_param_set), intent(in) :: mp

            cal_stress_mul_ctx_supported =
     &         mp%Iwkcoup_grad == 0 .and.
     &         mp%Iwkcoup_trip == 0 .and.
     &         mp%Iwkcoup_int  == 0 .and.
     &         mp%Iwkcoup_sup  == 0 .and.
     &         mp%Iwkcoup_temp == 0
            return
         endfunction cal_stress_mul_ctx_supported

         subroutine load_gp_context_from_globals(gp)
            implicit none
            type(gp_context), intent(out) :: gp

            gp%dt1 = dt1
            gp%det_Fg = det_Fg
            gp%det_Fe = det_Fe
            gp%det_Fp = det_Fp
            gp%Fg0 = Fg0
            gp%Fg = Fg
            gp%TFg = TFg
            gp%IFg = IFg
            gp%Fe0 = Fe0
            gp%Fe = Fe
            gp%Fp0 = Fp0
            gp%Fp = Fp
            gp%IFp0 = IFp0
            gp%TIFp0 = TIFp0
            gp%IFp = IFp
            gp%TIFp = TIFp
            gp%Qm = Fe0
            gp%CGE = CGE
            gp%CGEe_max = CGEe_max
            gp%Ftrp = Ftrp
            gp%IFtrp = IFtrp
            gp%pk2i = pk2i
            gp%pk2i_max = pk2i_max
            gp%pk2i_gnd = pk2i_gnd
            gp%pk2i_eq = pk2i_eq
            gp%pk2r = pk2r
            gp%cs = cs
            gp%cs0 = cs0
            gp%dpk2i = dpk2i
            gp%pk2i_M = pk2i_M
            gp%pk2r_M = pk2r_M
            gp%csM0 = csM0
            gp%csM = csM
            gp%eang00 = eang00
            gp%eang0 = eang0
            gp%eang = eang
            gp%Lp = Lp
            gp%dFpdt = dFpdt
            gp%MatJacb = MatJacb
            gp%IVB0 = IVB0
            gp%IVB = IVB
            gp%IVB_wcp = IVB_wcp
            gp%IVB_bk = IVB_bk
            gp%tau = tau
            gp%dgmdt = dgmdt
            gp%detGM = detGM
            gp%dIVB = dIVB
            gp%dIVBdt = dIVBdt
            gp%STFrlx6 = STFrlx6
            gp%ddgmdt_dtau = ddgmdt_dtau
            gp%ddgmdt_dIVB = ddgmdt_dIVB
            gp%ddIVBdt_ddgmdt = ddIVBdt_ddgmdt
            gp%ddIVBdt_dIVB = ddIVBdt_dIVB
            gp%ddgmdt_dpk2i = ddgmdt_dpk2i
            gp%ddIVBdt_dpk2i = ddIVBdt_dpk2i
            gp%Gv1 = Gv1
            gp%Gv2 = Gv2
            gp%dGv1_dpk2i = dGv1_dpk2i
            gp%dGv1_dIVB = dGv1_dIVB
            gp%dGv2_dpk2i = dGv2_dpk2i
            gp%dGv2_dIVB = dGv2_dIVB
            gp%IdGv1_dpk2i = IdGv1_dpk2i
            gp%IdGv2_dIVB = IdGv2_dIVB
            gp%eqM66Gv1 = eqM66Gv1
            gp%eqMnnGv2 = eqMnnGv2
            gp%IeqM66Gv1 = IeqM66Gv1
            gp%IeqMnnGv2 = IeqMnnGv2
            gp%eqM6nGv1 = eqM6nGv1
            gp%eqMn6Gv2 = eqMn6Gv2
            gp%dIVB_dpk2i = dIVB_dpk2i
            gp%dTdgmdt_dpk2i = dTdgmdt_dpk2i
            gp%dGv1_dE = dGv1_dE
            gp%dpk2i_dE = dpk2i_dE
            gp%dPk2r_dE = dPk2r_dE
            gp%dIFp_dE = dIFp_dE
            gp%dTIFp_dE = dTIFp_dE
            gp%dSTFrlx6_dE = dSTFrlx6_dE
            gp%dCGEe_mx_dE = dCGEe_mx_dE
            gp%dpk2i_mx_dE = dpk2i_mx_dE
            gp%dFp_dE = dFp_dE
            gp%ddgmdt_dE = ddgmdt_dE
            return
         endsubroutine load_gp_context_from_globals

         subroutine store_gp_context_to_globals(gp)
            implicit none
            type(gp_context), intent(in) :: gp

            dt1 = gp%dt1
            det_Fg = gp%det_Fg
            det_Fe = gp%det_Fe
            det_Fp = gp%det_Fp
            Fg0 = gp%Fg0
            Fg = gp%Fg
            TFg = gp%TFg
            IFg = gp%IFg
            Fe0 = gp%Fe0
            Fe = gp%Fe
            Fp0 = gp%Fp0
            Fp = gp%Fp
            IFp0 = gp%IFp0
            TIFp0 = gp%TIFp0
            IFp = gp%IFp
            TIFp = gp%TIFp
            CGE = gp%CGE
            CGEe_max = gp%CGEe_max
            Ftrp = gp%Ftrp
            IFtrp = gp%IFtrp
            pk2i = gp%pk2i
            pk2i_max = gp%pk2i_max
            pk2i_gnd = gp%pk2i_gnd
            pk2i_eq = gp%pk2i_eq
            pk2r = gp%pk2r
            cs = gp%cs
            cs0 = gp%cs0
            dpk2i = gp%dpk2i
            pk2i_M = gp%pk2i_M
            pk2r_M = gp%pk2r_M
            csM0 = gp%csM0
            csM = gp%csM
            eang00 = gp%eang00
            eang0 = gp%eang0
            eang = gp%eang
            Lp = gp%Lp
            dFpdt = gp%dFpdt
            MatJacb = gp%MatJacb
            IVB0 = gp%IVB0
            IVB = gp%IVB
            IVB_wcp = gp%IVB_wcp
            IVB_bk = gp%IVB_bk
            tau = gp%tau
            dgmdt = gp%dgmdt
            detGM = gp%detGM
            dIVB = gp%dIVB
            dIVBdt = gp%dIVBdt
            STFrlx6 = gp%STFrlx6
            ddgmdt_dtau = gp%ddgmdt_dtau
            ddgmdt_dIVB = gp%ddgmdt_dIVB
            ddIVBdt_ddgmdt = gp%ddIVBdt_ddgmdt
            ddIVBdt_dIVB = gp%ddIVBdt_dIVB
            ddgmdt_dpk2i = gp%ddgmdt_dpk2i
            ddIVBdt_dpk2i = gp%ddIVBdt_dpk2i
            Gv1 = gp%Gv1
            Gv2 = gp%Gv2
            dGv1_dpk2i = gp%dGv1_dpk2i
            dGv1_dIVB = gp%dGv1_dIVB
            dGv2_dpk2i = gp%dGv2_dpk2i
            dGv2_dIVB = gp%dGv2_dIVB
            IdGv1_dpk2i = gp%IdGv1_dpk2i
            IdGv2_dIVB = gp%IdGv2_dIVB
            eqM66Gv1 = gp%eqM66Gv1
            eqMnnGv2 = gp%eqMnnGv2
            IeqM66Gv1 = gp%IeqM66Gv1
            IeqMnnGv2 = gp%IeqMnnGv2
            eqM6nGv1 = gp%eqM6nGv1
            eqMn6Gv2 = gp%eqMn6Gv2
            dIVB_dpk2i = gp%dIVB_dpk2i
            dTdgmdt_dpk2i = gp%dTdgmdt_dpk2i
            dGv1_dE = gp%dGv1_dE
            dpk2i_dE = gp%dpk2i_dE
            dPk2r_dE = gp%dPk2r_dE
            dIFp_dE = gp%dIFp_dE
            dTIFp_dE = gp%dTIFp_dE
            dSTFrlx6_dE = gp%dSTFrlx6_dE
            dCGEe_mx_dE = gp%dCGEe_mx_dE
            dpk2i_mx_dE = gp%dpk2i_mx_dE
            dFp_dE = gp%dFp_dE
            ddgmdt_dE = gp%ddgmdt_dE
            return
         endsubroutine store_gp_context_to_globals

         recursive subroutine cal_stress_mul_legacy(ising,mp)
            implicit none
            integer i,j,k,l,m,n,i1,j1,k1,l1,m1,n1,is,js
            integer iNRloop
            integer ising_1,ising
            type(mat_param_set), intent(in) :: mp
            real(8) x1,x2,x3,x4,y1,z1
            real(8) M1_66(6,6),M1_96(9,6),M1_99(9,9)
            real(8) M2_66(6,6),M2_96(9,6),M2_99(9,9)
            real(8) M1_3333(3,3,3,3),M2_3333(3,3,3,3)
            real(8) MX1(3,3),MX2(3,3),MX3(3,3),MX4(3,3)
            type(gp_context) :: gp

            if(cal_stress_mul_ctx_supported(mp))then
               call load_gp_context_from_globals(gp)
               call cal_stress_mul(gp,mp,ising)
               call store_gp_context_to_globals(gp)
               return
            endif

            !1-----------------------------1
            !1   Whether Fg is distorting  1
            !1-----------------------------1
            ising=0
            call icams_determ(Fg,det_Fg)
            if(det_Fg < 1.d-10)then
               call pm(Fg,3,3)
               write(6,*) 'Strongly distorted, stop'
               ising=20 
               return
            endif
            det_Fe=det_Fg
            TFg=transpose(Fg)
            CGE=matmul(TFg,Fg)
            !1--------------------------------1
            !1   Whether Fp0 is ivertible     1
            !1--------------------------------1
            ising_1=0
            call gaussj(Fp0,3,IFp0,ising_1)
            if(ising_1/=0)then
               write(6,*) 'non ivertable for IFp0'
               ising=1  !Fp0 is non-invertible, stop current time step
               return
            endif
            TIFp0=transpose(IFp0)
            !1-------------------------------------------------1
            !1   Calculate pk2i_max(6)=C*(CGE_max(6)-I(6))/2   1
            !1-------------------------------------------------1
            MX1=matmul( matmul(TIFp0,CGE),IFp0 )  !==> add trip effect
            CGEe_max=matmul( matmul(transpose(IFtrp),MX1), IFtrp ) 
                       
            MX1=(CGEe_max-XI33)/2
            do i=1,6
               pk2i_max(i)=0
               do j=1,6
                  pk2i_max(i)=pk2i_max(i)
     &            +mp%Mstiff(i,j)*MX1(ib1(j),ib2(j))
               enddo
            enddo
            !1------------------------------------------------------------------1
            !1   Calculate STFrlx6=C*( CGE_max*smdMi + (CGE_max*smdMi)^T )/2    1
            !1------------------------------------------------------------------1
            do is=1,mp%N_slip
               MX1=matmul(CGEe_max,mp%Msmd(is,:,:))
               MX2=matmul( matmul(transpose(IFtrp), !==> add trip effect
     &               (MX1+transpose(MX1))/2), IFtrp )
               do i=1,6
                  STFrlx6(is,i)=0
                  do j=1,6
                     STFrlx6(is,i)=STFrlx6(is,i)
     &               +mp%Mstiff(i,j)*MX2(ib1(j),ib2(j))
                  enddo
               enddo
            enddo
            !1------------------------------------------------------1
            !1    Begin newton raphson method to solve pk2i, IVB    1
            !1------------------------------------------------------1
            do iNRloop=1,Nnr_max

            if (mp%Iwkcoup_sup==1) then
c              superalloy model
                pk2i_eqx = pk2i + pk2i_intx   !==> add misfit stress
                pk2i_eqy = pk2i + pk2i_inty
                pk2i_eqz = pk2i + pk2i_intz
                pk2i_eqp = pk2i + pk2i_intp 
 
               do is=1,12
                  tau(is)=dot_product(pk2i_eqx,mp%V2smd(is,:))
               enddo  
               do is=13,24
                  tau(is)=dot_product(pk2i_eqy,mp%V2smd(is,:))
               enddo   
               do is=25,36
                  tau(is)=dot_product(pk2i_eqz,mp%V2smd(is,:))
               enddo   
               do is=37,48
                  tau(is)=dot_product(pk2i_eqp,mp%V2smd(is,:))
               enddo
               
               do is=53,54
                  tau(is)=dot_product(pk2i_eqx,mp%V2smd(is,:))
               enddo 
               do is=51,52
                  tau(is)=dot_product(pk2i_eqy,mp%V2smd(is,:))
               enddo 
               do is=49,50
                  tau(is)=dot_product(pk2i_eqz,mp%V2smd(is,:))
               enddo 
               do is=55,60
                  tau(is)=dot_product(pk2i_eqp,mp%V2smd(is,:))
               enddo

            else   ! not superalloy, normal model
             if (mp%Iwkcoup_grad==0) then
                pk2i_eq = pk2i
             else
                pk2i_eq = pk2i + pk2i_gnd
             endif
             do is=1,mp%N_slip
               tau(is)=dot_product(pk2i_eq,mp%V2smd(is,:))
             enddo 
            endif


               ising_1=0
               if (mp%Iwkcoup_sup==0 .and. mp%Iwkcoup_temp==0) then
                  call sub_flow_harden_std(
     &                 Iexp_loc,mp%N_slip,mp%crss0,mp%crsss,
     &                 mp%hdrt0,mp%pwhd,mp%shrt0,mp%pwfl,
     &                 mp%HMij,IVB_wcp,IVB_bk,tau,
     &                 IVB,dgmdt,ddgmdt_dtau,ddgmdt_dIVB,
     &                 dIVBdt,ddIVBdt_ddgmdt,ddIVBdt_dIVB,ising_1)
               else
                  call sub_flow_harden(
     &                 Iexp_loc,IVB_wcp,IVB_bk,IVB_cl,IVB_m,IVB_kw,
     &                 tau,IVB,dgmdt,ddgmdt_dtau,ddgmdt_dIVB,
     &                 dIVBdt,ddIVBdt_ddgmdt,ddIVBdt_dIVB,ising_1)
               endif
               do is=1,mp%N_slip
                  ddgmdt_dpk2i(is,:)=ddgmdt_dtau(is)*mp%V2smd(is,:)
               enddo
               ddIVBdt_dpk2i(1:mp%N_slip,:)=
     &         matmul(ddIVBdt_ddgmdt(1:mp%N_slip,1:mp%N_slip),
     &         ddgmdt_dpk2i(1:mp%N_slip,:))

               if(ising_1/=0)then
                  ising=ising_1  !!! shearrate non-cvg,stop
                  return
               endif
               Gv1=+pk2i-pk2i_max
     &         +matmul(transpose(STFrlx6(1:mp%N_slip,:)),
     &         dgmdt(1:mp%N_slip))*dt1
               Gv2(1:mp%N_slip)=+IVB(1:mp%N_slip)
     &         -IVB0(1:mp%N_slip)-dIVBdt(1:mp%N_slip)*dt1
               x1=sum(dabs(Gv1))/mp%refv_pk2i
               x2=sum(dabs(Gv2(1:mp%N_slip)))/mp%refv_IVB


               if(x1+x2<toler_NRloop .or. Iexp_loc==1)then
                  pk2i=+pk2i_max
     &                 -matmul(transpose(STFrlx6(1:mp%N_slip,:)),
     &                 dgmdt(1:mp%N_slip))*dt1
                  IVB(1:mp%N_slip)=IVB0(1:mp%N_slip)
     &                 +dIVBdt(1:mp%N_slip)*dt1
                  goto 101 !Converge jump out NR loop and continue
               endif
               ising_1=0
               call cal_NRmatrix(ising_1)
               if(ising_1/=0)then
                  ising=ising_1  !NR matrix cannot be got,stop
                  return
               endif
		
               dpk2i=matmul(IeqM66Gv1,
     &         -(Gv1-matmul(eqM6nGv1(:,1:mp%N_slip),
     &         Gv2(1:mp%N_slip))))
               dIVB(1:mp%N_slip)=matmul(
     &         IeqMnnGv2(1:mp%N_slip,1:mp%N_slip),
     &         -(Gv2(1:mp%N_slip)
     &         -matmul(eqMn6Gv2(1:mp%N_slip,:),Gv1)) )

               pk2i=pk2i+dpk2i
               IVB(1:mp%N_slip) =IVB(1:mp%N_slip)
     &                          + dIVB(1:mp%N_slip)

            enddo
            ising=5  !Non-coverge for Nnr_max, stop current time step
            return
101         continue
            !1--------------------------------------------------------1
            !1   End of newton raphson method, solve cauchy stress    1
            !1--------------------------------------------------------1 

c           if (Iwkcoup_bk/=0) then
c            do is=1,N_slip
c             fem_dgmdt(ie,ig,is)=dgmdt(is)
c            enddo
c           endif
         

            if (mp%Iwkcoup_sup==1) then   ! superalloy model

            do is=1,36
              fem_gamm(ie,ig,is)=fem_gamm0(ie,ig,is)+dabs(dgmdt(is))*dt1
              fem_gamp(ie,ig,is)=0.d0
            enddo

            do is=37,48
              fem_gamm(ie,ig,is)=0.d0
              fem_gamp(ie,ig,is)=fem_gamp0(ie,ig,is)+dabs(dgmdt(is))*dt1
            enddo

            call icams_conv6to33(pk2i,ib1,ib2,pk2i_M)

            do i=1,3
            do j=1,3
               Lpx(i,j)=dot_product(dgmdt(1:12),mp%Msmd(1:12,i,j))
     &               + dot_product(dgmdt(53:54),mp%Msmd(53:54,i,j))
  
               Lpy(i,j)=dot_product(dgmdt(13:24),mp%Msmd(13:24,i,j))
     &               + dot_product(dgmdt(51:52),mp%Msmd(51:52,i,j))

               Lpz(i,j)=dot_product(dgmdt(25:36),mp%Msmd(25:36,i,j))
     &               + dot_product(dgmdt(49:50),mp%Msmd(49:50,i,j))

               Lpp(i,j)=dot_product(dgmdt(37:48),mp%Msmd(37:48,i,j))
     &               + dot_product(dgmdt(55:60),mp%Msmd(55:60,i,j))
            enddo
            enddo

            Fpx=matmul(XI33+Lpx*dt1, fem_Fpx0(ie,ig,:,:))
            Fpy=matmul(XI33+Lpy*dt1, fem_Fpy0(ie,ig,:,:))
            Fpz=matmul(XI33+Lpz*dt1, fem_Fpz0(ie,ig,:,:))
            Fppp=matmul(XI33+Lpp*dt1,fem_Fppp0(ie,ig,:,:))

            fem_Fpx(ie,ig,:,:)=Fpx
            fem_Fpy(ie,ig,:,:)=Fpy
            fem_Fpz(ie,ig,:,:)=Fpz
            fem_Fppp(ie,ig,:,:)=Fppp

            egx=(matmul(transpose(Fpx),Fpx)-XI33)/2
            x1=(egx(1,1)+egx(2,2)+egx(3,3))/3
            egx=egx-x1*XI33

            egy=(matmul(transpose(Fpy),Fpy)-XI33)/2
            x1=(egy(1,1)+egy(2,2)+egy(3,3))/3
            egy=egy-x1*XI33

            egz=(matmul(transpose(Fpz),Fpz)-XI33)/2
            x1=(egz(1,1)+egz(2,2)+egz(3,3))/3
            egz=egz-x1*XI33

            egp=(matmul(transpose(Fppp),Fppp)-XI33)/2
            x1=(egp(1,1)+egp(2,2)+egp(3,3))/3
            egp=egp-x1*XI33

            call icams_conv33to6(egx,ib1,ib2,epsx)
            call icams_conv33to6(egy,ib1,ib2,epsy)
            call icams_conv33to6(egz,ib1,ib2,epsz)
            call icams_conv33to6(egp,ib1,ib2,epsp)

            fem_epsx(ie,ig,:)=epsx
            fem_epsy(ie,ig,:)=epsy
            fem_epsz(ie,ig,:)=epsz
            fem_epsp(ie,ig,:)=epsp
    
            Lp=fx*Lpx+fy*Lpy+fz*Lpz+fpp*Lpp

            else ! not superalloy, normal model

            call icams_conv6to33(pk2i,ib1,ib2,pk2i_M)

            do i=1,3
            do j=1,3
               Lp(i,j)=dot_product(dgmdt(1:mp%N_slip),
     &                              mp%Msmd(1:mp%N_slip,i,j))
            enddo
            enddo

            endif

            Fp=matmul(XI33+Lp*dt1,Fp0)

            call icams_determ(Fp,det_Fp)
            if(det_Fp==0)then
               write(6,*) 'det of Fp is zero'
               ising=22  !Fp is non-invertible,stop current time step
               return
            endif
            Fp=Fp/det_Fp**(1/3.0)
            det_Fp=1

            dFpdt=(Fp-Fp0)/dt1

            ising_1=0
            call gaussj(Fp,3,iFp,ising_1)
            if(ising_1/=0)then
               write(6,*) 'non ivertable for Fp'
               ising=2  !Fp is non-invertible,stop current time step
               return
            endif

            TIFp=transpose(IFp)
            Fe=matmul(Fg,iFp)
            call icams_determ(Fe,det_Fe)
            csM=matmul(matmul(Fe,pk2i_M),transpose(Fe))/det_Fe
            call icams_conv33to6(csM,ib1,ib2,cs)

            ising_1=0
            call gaussj(Fg,3,IFg,ising_1)
            if(ising_1/=0)then
               write(6,*) 'non ivertable for Fg'
               ising=3  !Fp is non-invertible,stop current time step
               return
            endif
            pk2r_M=matmul(matmul(IFg,csM),transpose(IFg))*det_Fe
            call icams_conv33to6(pk2r_M,ib1,ib2,pk2r)

            detGM=dgmdt*dt1
            !1----------------------------------------1
            !1    calculate eang and misorientation   1
            !1----------------------------------------1
            call caleulang(Fe,eang(1:3),ising_1)
            if(ising_1/=0) eang(1:3)=eang0(1:3)


            call icams_misori(eang00(1:3),eang(1:3),eang(4))

            !1----------------------------------1
            !1    calculate material stiffness  1
            !1----------------------------------1
            if(Iexp_abq==0)then
               ising_1=0
               call cal_NRmatrix(ising_1)
               if(ising_1/=0)then
                  ising=420+ising_1  !NR matrix cannot be got,stop
                  return
               endif
               call cal_MatStiffness
            endif

c            print*,'ok',Iexp_abq,ising
c            read*

            return
         endsubroutine cal_stress_mul_legacy

         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         !+                                                                                  +
         !+   THIS ROUTINE calculates tangent matrix for newton raphson algoriths            +
         !+                                                                                  +
         !+   #1: dG1_dpk2i,dG1_dIVB,IdG1_dpk2i                                              +
         !+   #2: dG2_dpk2i,dG2_dIVB,IdG2_dIVB                                               +
         !+   #3: eqM66Gv1,eqMnnGv2,IeqM66Gv1,IeqMnnGv2                                      +
         !+                                                                                  +
         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         recursive subroutine cal_NRmatrix(ising)
            implicit none
            integer i,j,k,l,m,n,i1,j1,k1,l1,m1,n1,is,js
            integer reg_try
            integer ising_1,ising
            real(8) x1,y1,z1,reg_eps,norm_dg1
            real(8) M1_66(6,6),M1_96(9,6),M1_99(9,9)
            real(8) M2_66(6,6),M2_96(9,6),M2_99(9,9)
            real(8) M1_3333(3,3,3,3),M2_3333(3,3,3,3)
            real(8) MX1(3,3),MX2(3,3),MX3(3,3),MX4(3,3)
            real(8) Mreg66(6,6)
            !-----------------------------------------------------------------------------
            ! calculate dG1_dpk2i(6,6),dG1_dIVB(6,:),IdG1_dpk2i(6,6)
            !-----------------------------------------------------------------------------
            ising=0
            dGv1_dpk2i=XI66+matmul(transpose(STFrlx6(1:N_slip,:)),
     &                ddgmdt_dpk2i(1:N_slip,:))*dt1
            do is=1,N_slip
               dGv1_dIVB(:,is)=STFrlx6(is,:)*ddgmdt_dIVB(is)*dt1
            enddo

            ising_1=0
            call gaussj(dGv1_dpk2i,6,idGv1_dpk2i,ising_1)
            if(ising_1/=0)then
               norm_dg1 = maxval(dabs(dGv1_dpk2i))
               reg_eps = max(1.d-14, 1.d-10*norm_dg1)
               do reg_try=1,4
                  Mreg66 = dGv1_dpk2i
                  do i=1,6
                     Mreg66(i,i)=Mreg66(i,i)+reg_eps
                  enddo
                  dGv1_dpk2i = Mreg66
                  ising_1=0
                  call gaussj(dGv1_dpk2i,6,idGv1_dpk2i,ising_1)
                  if(ising_1==0) then
                     write(6,*) 'regularized dGv1_dpk2i, eps=',reg_eps
                     goto 111
                  endif
                  reg_eps = reg_eps*1.d2
               enddo
               write(6,*) 'non ivertable for dGv1_dpk2i'
               ising=41  !dGv1_dpk2i non-invertable,stop
               return
            endif
111         continue
            !-----------------------------------------------------------------------------
            ! calculate dG2_dpk2i(:,6),dG2_dIVB(:,:),IdG2_dIVB(:,:)
            !-----------------------------------------------------------------------------
            dGv2_dpk2i(1:N_slip,:)=-ddIVBdt_dpk2i(1:N_slip,:)*dt1
            dGv2_dIVB(1:N_slip,1:N_slip) =+XInn(1:N_slip,1:N_slip)
     &                            -ddIVBdt_dIVB(1:N_slip,1:N_slip) *dt1

            ising_1=0
            call gaussj( dGv2_dIVB(1:N_slip,1:N_slip),N_slip,
     &                  idGv2_dIVB(1:N_slip,1:N_slip),ising_1)
            if(ising_1/=0)then
               write(6,*) 'non ivertable for dGv2_dIVB'
               ising=42  !dGv2_dIVB non-invertable,stop
               return
            endif
            !-----------------------------------------------------------------------------
            !-----Calculate eqM66Gv1,eqMnnGv2,IeqM66Gv1,IeqMnnGv2
            !-----------------------------------------------------------------------------
            eqM6nGv1=matmul(dGv1_dIVB(:,1:N_slip),
     &                     IdGv2_dIVB(1:N_slip,1:N_slip))
            eqM66Gv1=-matmul(matmul(dGv1_dIVB(:,1:N_slip),
     &         IdGv2_dIVB(1:N_slip,1:N_slip)),dGv2_dpk2i(1:N_slip,:))
     &         +dGv1_dpk2i

            ising_1=0
            call gaussj(eqM66Gv1,6,IeqM66Gv1,ising_1)
            if(ising_1/=0)then
               write(6,*) 'non ivertable for eqM66Gv1'
               ising=43  !eqM66Gv1 non-invertable,stop
               return
            endif
            eqMn6Gv2(1:N_slip,:)=matmul(dGv2_dpk2i(1:N_slip,:),IdGv1_dpk2i)
            eqMnnGv2(1:N_slip,1:N_slip)=-matmul(matmul(dGv2_dpk2i(1:N_slip,:),
     &                          IdGv1_dpk2i),dGv1_dIVB(:,1:N_slip))
     &                         +dGv2_dIVB(1:N_slip,1:N_slip)
            ising_1=0
            call gaussj( eqMnnGv2(1:N_slip,1:N_slip),N_slip,
     &                  IeqMnnGv2(1:N_slip,1:N_slip),ising_1)
            if(ising_1/=0)then
               write(6,*) 'non ivertable for eqMnnGv2'
               ising=44  !eqMnnGv2 non-invertable,stop
               return
            endif
            RETURN
         endsubroutine cal_NRmatrix
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c +                                                                    +
c +   THIS ROUTINE calculates material jacb for umat STFjc_66          +
c +   THIS ROUTINE also can calculates material jacb for uel STFtk_66  +
c +                                                                    +
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         recursive subroutine cal_MatStiffness
         implicit none
         integer i,j,k,l,m,n,i1,j1,k1,l1,m1,n1,is,js
         real(8)    x1,y1,z1
         real(8)    M1_66(6,6),M1_96(9,6),M1_99(9,9)
         real(8)    M2_66(6,6),M2_96(9,6),M2_99(9,9)
         real(8)    M3_66(6,6),M3_96(9,6),M3_99(9,9)
         real(8)    M1_3333(3,3,3,3),M2_3333(3,3,3,3),M3_3333(3,3,3,3)
         real(8)    MX1(3,3),MX2(3,3),MX3(3,3),MX4(3,3)
c-----------------------------------------------------------------------------
c     jacb#1: dpk2i_dE(6,6)
c         #1.1: dIVB_dpk2i(N_slip,6)
c-----------------------------------------------------------------------------
         dIVB_dpk2i=-matmul(IdGv2_dIVB,dGv2_dpk2i)
c-----------------------------------------------------------------------------
c     jacb#1: dpk2i_dE
c         #1.2:   dGv1_dE
c         #1.2.1:   dCGEe_mx_dE(6,6) = d(Fe_mx^T*Fe_mx)_dE
c-----------------------------------------------------------------------------
         do i=1,6
         do j=1,6
            if(j<=3)then
               dCGEe_mx_dE(i,j)=2*TIFp0(ib1(i),ib1(j))
     &                         *IFp0(ib2(j),ib2(i))
            else
               dCGEe_mx_dE(i,j)=
     &         (+2*TIFp0(ib1(i),ib1(j  ))*IFp0(ib2(j  ),ib2(i))
     &          +2*TIFp0(ib1(i),ib1(j+3))*IFp0(ib2(j+3),ib2(i)))/2
            endif
         enddo
         enddo
c-----------------------------------------------------------------------------
c     jacb#1: dpk2i_dE
c         #1.2:   dGv1_dE
c         #1.2.2:   dpk2i_mx_dE(6,6)
c-----------------------------------------------------------------------------
         dpk2i_mx_dE=matmul( STFei26 , dCGEe_mx_dE )/2
c-----------------------------------------------------------------------------
c     jacb#1: dpk2i_dE
c         #1.2:   dGv1_dE
c         #1.2.3:   dSTFrlx6_dE(N_slip,6,6)
c-----------------------------------------------------------------------------
         do is=1,N_slip
            MX1=matmul(IFp0,smdMi(is,:,:))
            MX2=transpose(MX1)
            do i=1,6
            do j=1,6
               M1_66(i,j)=
     &           +2*TIFp0(ib1(i),ib1(j))*MX1 (ib2(j),ib2(j))
     &           +2*MX2  (ib1(i),ib1(j))*IFp0(ib2(j),ib2(j))
            enddo
            enddo
            dSTFrlx6_dE(is,:,:)=matmul(STFei26 , M1_66)/2
         enddo
c-----------------------------------------------------------------------------
c     jacb#1: dpk2i_dE
c         #1.2:   dGv1_dE
c         #1.2.4:   dGv1_dE(6,6)
c-----------------------------------------------------------------------------
         do i=1,6
         do j=1,6
            dGv1_dE(i,j)=
     &      -dpk2i_mx_dE(i,j)
     &      +dot_product(dgmdt,dSTFrlx6_dE(:,i,j))*dt1
         enddo
         enddo
c-----------------------------------------------------------------------------
c     jacb#1: dpk2i_dE
c         #1.3:   dpk2i_dE(6,6)
c-----------------------------------------------------------------------------
         dpk2i_dE=-matmul(IeqM66Gv1,dGv1_dE)
c-----------------------------------------------------------------------------
c     jacb#2: dpk2r_dE
c         #2.1:  dIFp_dE(9,6),dTIFp_dE(9,6)
c-----------------------------------------------------------------------------
         M1_99=0
         M2_99=0
         M3_99=0
         do is=1,N_slip
            MX1=matmul(IFp0,smdMi(is,:,:))
            MX2=transpose(MX1)
            MX4=matmul(smdMi(is,:,:),Fp0)
            do i=1,9
               if(i<=6) i1=i
               if(i >6) i1=i-3
               MX3(ib1(i),ib2(i))=
     &         +ddgmdt_dpk2i(is,i1)
     &         +ddgmdt_dIVB(is)*dIVB_dpk2i(is,i1)
               dTdgmdt_dpk2i(is,i)=+ddgmdt_dpk2i(is,i1)
     &         +ddgmdt_dIVB(is)*dIVB_dpk2i(is,i1)
            enddo
            do i=1,9
            do j=1,9
            M1_99(i,j)=M1_99(i,j)+MX1(ib1(i),ib2(i))*MX3(ib1(j),ib2(j))
            M2_99(i,j)=M2_99(i,j)+MX2(ib1(i),ib2(i))*MX3(ib1(j),ib2(j))
            M3_99(i,j)=M3_99(i,j)+MX4(ib1(i),ib2(i))*MX3(ib1(j),ib2(j))
            enddo
            enddo
         enddo
         dIFp_dE=0
         dTIFp_dE=0
         M1_96(1:6,:)=dpk2i_dE
         M1_96(7,:)=M1_96(4,:)
         M1_96(8,:)=M1_96(5,:)
         M1_96(9,:)=M1_96(6,:)
          dIFp_dE=-matmul(M1_99,M1_96)*dt1
         dTIFp_dE=-matmul(M2_99,M1_96)*dt1

           dFp_dE=-matmul(M3_99,M1_96)*dt1
         ddgmdt_dE=matmul(dTdgmdt_dpk2i,M1_96)*dt1

         dFp_dE(:,4:6)=dFp_dE(:,4:6)*2
         ddgmdt_dE(:,4:6)=ddgmdt_dE(:,4:6)*2
c         do i=1,9
c         do j=1,3
c            dFp_dx0(ib1(i),ib2(i),j)
c     &      =dot_product(dFp_dE(i,:),dEr_dx0(ig,:,j))
c         enddo
c         enddo
c         ddgmdt_dx0=matmul(ddgmdt_dE,dEr_dx0(ig,:,:))

c-----------------------------------------------------------------------------
C     jacb#2: dpk2r_dE
c         #2.2: dpk2r_dE(6,6)
c-----------------------------------------------------------------------------
         dpk2r_dE=0
         MX1=matmul(pk2i_M,TIFp)
         MX2=matmul(IFp,pk2i_M)
         do i=1,6
         do k=1,6
            do m=1,9
               if(m<=6)m1=m
               if(m >6)m1=m-3
               dpk2r_dE(i,k)=dpk2r_dE(i,k)
     &        +XI33(ib1(i),ib1(m))* MX1(ib2(m),ib2(i))* dIFp_dE(m ,k)
     &        + IFp(ib1(i),ib1(m))*TIFp(ib2(m),ib2(i))*dpk2i_dE(m1,k)
     &        + MX2(ib1(i),ib1(m))*XI33(ib2(i),ib2(m))*dTIFp_dE(m ,k)
            enddo
         enddo
         enddo
c-----------------------------------------------------------------------------
C     jacb#3: (STF_JC_3333)_ijkl = + (dpk2r_dE)_mnop*F_im*F_jn*F_ko*F_lp/det(F)
c                                  +  I_ik * CS_lj
c                                  + CS_ik *  I_lj
c                                  - CS_ij *  I_kl
c     MatJacb=(STFjc_66+transpose(STFjc_66))/2
c
C     jacb#3: (STF_TK_3333)_ijkl = + (dpk2r_dE)_mnop*F_im*F_jn*F_ko*F_lp
c     MatJacb=(STF_TK_66+transpose(STF_TK_66))/2/det(F)
c-----------------------------------------------------------------------------
c-------------------------------------------------------------------------------------------
C     Calculate stiffness for Jaummann rate of Cauchy stress STFjc_66 for user material
c-------------------------------------------------------------------------------------------
         do i=1,9
         do j=1,9
            if(i<=6) i1=i
            if(i >6) i1=i-3
            if(j<=6) j1=j
            if(j >6) j1=j-3
            M1_3333(ib1(i),ib2(i),ib1(j),ib2(j))=dpk2r_dE(i1,j1)
         enddo
         enddo
         M2_3333=0
         do i=1,3
         do j=1,3
         do k=1,3
         do l=1,3
            x1=0
            do i1=1,3
            do j1=1,3
            do k1=1,3
            do l1=1,3
               x1=x1+M1_3333(i1,j1,k1,l1)
     &        *Fg(i,i1)*Fg(j,j1)*Fg(k,k1)*Fg(l,l1)
            enddo
            enddo
            enddo
            enddo

            M2_3333(i,j,k,l)=x1/det_Fg
     &                   +XI33(i,k)*csM0(l,j)
     &                   +csM0(i,k)*XI33(l,j)
     &                   +csM0(i,j)*XI33(k,l)
         enddo
         enddo
         enddo
         enddo
         do i=1,6
         do j=1,6
            STFjc_66(I,J)=M2_3333(ib1(i),ib2(i),ib1(j),ib2(j))
         enddo
         enddo
c-----------------------------------------------------------------------------------------
C     Calculate stiffness for Trusdell rate of Kirchhoff stress STFtk_66 for user element
c-----------------------------------------------------------------------------------------
         do i=1,9
         do j=1,9
            if(i<=6) i1=i
            if(i >6) i1=i-3
            if(j<=6) j1=j
            if(j >6) j1=j-3
            M1_3333(ib1(i),ib2(i),ib1(j),ib2(j))=dpk2r_dE(i1,j1)
         enddo
         enddo
         M2_3333=0
         do i=1,3
         do j=1,3
         do k=1,3
         do l=1,3
            do k1=1,3
            do l1=1,3
               M2_3333(i,j,k,l)=M2_3333(i,j,k,l)
     &        +M1_3333(i,j,k1,l1)*Fg(k,k1)*Fg(l,l1)
            enddo
            enddo
         enddo
         enddo
         enddo
         enddo
         do i=1,6
         do j=1,6
            STFtk_66(I,J)=M2_3333(ib1(i),ib2(i),ib1(j),ib2(j))
         enddo
         enddo
c-----------------------------------------------------------------------------
C     jacb#4: Material tangent matrix
c-----------------------------------------------------------------------------
         MatJacb=STFjc_66
         return
         endsubroutine cal_MatStiffness

         recursive subroutine cal_stress_mul(gp,mp,ising)
            implicit none
            type(gp_context), intent(inout) :: gp
            type(mat_param_set), intent(in) :: mp
            integer, intent(out) :: ising
            integer i,j,is,iNRloop,ising_1
            real(8) x1,x2
            real(8) MX1(3,3),MX2(3,3)

            ising=0
            gp%Ftrp=XI33
            gp%IFtrp=XI33
            gp%IVB_wcp=0.d0
            gp%pk2i_gnd=0.d0

            call icams_determ(gp%Fg,gp%det_Fg)
            if(gp%det_Fg < 1.d-10)then
               call pm(gp%Fg,3,3)
               write(6,*) 'Strongly distorted, stop'
               ising=20
               return
            endif
            gp%det_Fe=gp%det_Fg
            gp%TFg=transpose(gp%Fg)
            gp%CGE=matmul(gp%TFg,gp%Fg)

            ising_1=0
            call gaussj(gp%Fp0,3,gp%IFp0,ising_1)
            if(ising_1/=0)then
               write(6,*) 'non ivertable for IFp0'
               ising=1
               return
            endif
            gp%TIFp0=transpose(gp%IFp0)

            MX1=matmul(matmul(gp%TIFp0,gp%CGE),gp%IFp0)
            gp%CGEe_max=MX1
            MX1=(gp%CGEe_max-XI33)/2
            do i=1,6
               gp%pk2i_max(i)=0.d0
               do j=1,6
                  gp%pk2i_max(i)=gp%pk2i_max(i)
     &            +mp%Mstiff(i,j)*MX1(ib1(j),ib2(j))
               enddo
            enddo

            do is=1,mp%N_slip
               MX1=matmul(gp%CGEe_max,mp%Msmd(is,:,:))
               MX2=(MX1+transpose(MX1))/2
               do i=1,6
                  gp%STFrlx6(is,i)=0.d0
                  do j=1,6
                     gp%STFrlx6(is,i)=gp%STFrlx6(is,i)
     &               +mp%Mstiff(i,j)*MX2(ib1(j),ib2(j))
                  enddo
               enddo
            enddo

            do iNRloop=1,Nnr_max
               gp%pk2i_eq=gp%pk2i
               do is=1,mp%N_slip
                  gp%tau(is)=dot_product(gp%pk2i_eq,mp%V2smd(is,:))
               enddo

               ising_1=0
               call sub_flow_harden_std(
     &              Iexp_loc,mp%N_slip,mp%crss0,mp%crsss,
     &              mp%hdrt0,mp%pwhd,mp%shrt0,mp%pwfl,
     &              mp%HMij,gp%IVB_wcp,gp%IVB_bk,gp%tau,
     &              gp%IVB,gp%dgmdt,gp%ddgmdt_dtau,
     &              gp%ddgmdt_dIVB,gp%dIVBdt,
     &              gp%ddIVBdt_ddgmdt,gp%ddIVBdt_dIVB,ising_1)
               if(ising_1/=0)then
                  ising=ising_1
                  return
               endif

               do is=1,mp%N_slip
                  gp%ddgmdt_dpk2i(is,:)=
     &               gp%ddgmdt_dtau(is)*mp%V2smd(is,:)
               enddo
               gp%ddIVBdt_dpk2i(1:mp%N_slip,:)=
     &            matmul(gp%ddIVBdt_ddgmdt(1:mp%N_slip,1:mp%N_slip),
     &            gp%ddgmdt_dpk2i(1:mp%N_slip,:))

               gp%Gv1=gp%pk2i-gp%pk2i_max
     &            +matmul(transpose(gp%STFrlx6(1:mp%N_slip,:)),
     &            gp%dgmdt(1:mp%N_slip))*gp%dt1
               gp%Gv2(1:mp%N_slip)=gp%IVB(1:mp%N_slip)
     &            -gp%IVB0(1:mp%N_slip)
     &            -gp%dIVBdt(1:mp%N_slip)*gp%dt1
               x1=sum(dabs(gp%Gv1))/mp%refv_pk2i
               x2=sum(dabs(gp%Gv2(1:mp%N_slip)))/mp%refv_IVB

               if(x1+x2<toler_NRloop .or. Iexp_loc==1)then
                  gp%pk2i=gp%pk2i_max
     &               -matmul(transpose(gp%STFrlx6(1:mp%N_slip,:)),
     &               gp%dgmdt(1:mp%N_slip))*gp%dt1
                  gp%IVB(1:mp%N_slip)=gp%IVB0(1:mp%N_slip)
     &               +gp%dIVBdt(1:mp%N_slip)*gp%dt1
                  goto 101
               endif

               ising_1=0
               call cal_NRmatrix_ctx(gp,mp,ising_1)
               if(ising_1/=0)then
                  ising=ising_1
                  return
               endif

               gp%dpk2i=matmul(gp%IeqM66Gv1,
     &            -(gp%Gv1-matmul(gp%eqM6nGv1(:,1:mp%N_slip),
     &            gp%Gv2(1:mp%N_slip))))
               gp%dIVB(1:mp%N_slip)=matmul(
     &            gp%IeqMnnGv2(1:mp%N_slip,1:mp%N_slip),
     &            -(gp%Gv2(1:mp%N_slip)
     &            -matmul(gp%eqMn6Gv2(1:mp%N_slip,:),gp%Gv1)))

               gp%pk2i=gp%pk2i+gp%dpk2i
               gp%IVB(1:mp%N_slip)=gp%IVB(1:mp%N_slip)
     &                            +gp%dIVB(1:mp%N_slip)
            enddo

            ising=5
            return
101         continue

            call icams_conv6to33(gp%pk2i,ib1,ib2,gp%pk2i_M)
            do i=1,3
            do j=1,3
               gp%Lp(i,j)=dot_product(gp%dgmdt(1:mp%N_slip),
     &                                 mp%Msmd(1:mp%N_slip,i,j))
            enddo
            enddo

            gp%Fp=matmul(XI33+gp%Lp*gp%dt1,gp%Fp0)
            call icams_determ(gp%Fp,gp%det_Fp)
            if(gp%det_Fp==0.d0)then
               write(6,*) 'det of Fp is zero'
               ising=22
               return
            endif
            gp%Fp=gp%Fp/gp%det_Fp**(1/3.0)
            gp%det_Fp=1.d0
            gp%dFpdt=(gp%Fp-gp%Fp0)/gp%dt1

            ising_1=0
            call gaussj(gp%Fp,3,gp%IFp,ising_1)
            if(ising_1/=0)then
               write(6,*) 'non ivertable for Fp'
               ising=2
               return
            endif

            gp%TIFp=transpose(gp%IFp)
            gp%Fe=matmul(gp%Fg,gp%IFp)
            call icams_determ(gp%Fe,gp%det_Fe)
            gp%csM=matmul(matmul(gp%Fe,gp%pk2i_M),
     &                    transpose(gp%Fe))/gp%det_Fe
            call icams_conv33to6(gp%csM,ib1,ib2,gp%cs)

            ising_1=0
            call gaussj(gp%Fg,3,gp%IFg,ising_1)
            if(ising_1/=0)then
               write(6,*) 'non ivertable for Fg'
               ising=3
               return
            endif
            gp%pk2r_M=matmul(matmul(gp%IFg,gp%csM),
     &                       transpose(gp%IFg))*gp%det_Fe
            call icams_conv33to6(gp%pk2r_M,ib1,ib2,gp%pk2r)

            gp%detGM=gp%dgmdt*gp%dt1
            call caleulang(gp%Fe,gp%eang(1:3),ising_1)
            if(ising_1/=0) gp%eang(1:3)=gp%eang0(1:3)
            call icams_misori(gp%eang00(1:3),gp%eang(1:3),gp%eang(4))

            if(Iexp_abq==0)then
               ising_1=0
               call cal_NRmatrix_ctx(gp,mp,ising_1)
               if(ising_1/=0)then
                  ising=420+ising_1
                  return
               endif
               call cal_MatStiffness_ctx(gp,mp)
            endif

            return
         endsubroutine cal_stress_mul

         recursive subroutine cal_NRmatrix_ctx(gp,mp,ising)
            implicit none
            type(gp_context), intent(inout) :: gp
            type(mat_param_set), intent(in) :: mp
            integer, intent(out) :: ising
            integer i,is,reg_try,ising_1
            real(8) reg_eps,norm_dg1
            real(8) Mreg66(6,6)

            ising=0
            gp%dGv1_dpk2i=XI66
     &         +matmul(transpose(gp%STFrlx6(1:mp%N_slip,:)),
     &         gp%ddgmdt_dpk2i(1:mp%N_slip,:))*gp%dt1
            do is=1,mp%N_slip
               gp%dGv1_dIVB(:,is)=gp%STFrlx6(is,:)
     &            *gp%ddgmdt_dIVB(is)*gp%dt1
            enddo

            ising_1=0
            call gaussj(gp%dGv1_dpk2i,6,gp%IdGv1_dpk2i,ising_1)
            if(ising_1/=0)then
               norm_dg1=maxval(dabs(gp%dGv1_dpk2i))
               reg_eps=max(1.d-14,1.d-10*norm_dg1)
               do reg_try=1,4
                  Mreg66=gp%dGv1_dpk2i
                  do i=1,6
                     Mreg66(i,i)=Mreg66(i,i)+reg_eps
                  enddo
                  gp%dGv1_dpk2i=Mreg66
                  ising_1=0
                  call gaussj(gp%dGv1_dpk2i,6,
     &                        gp%IdGv1_dpk2i,ising_1)
                  if(ising_1==0)then
                     write(6,*) 'regularized dGv1_dpk2i, eps=',reg_eps
                     goto 111
                  endif
                  reg_eps=reg_eps*1.d2
               enddo
               write(6,*) 'non ivertable for dGv1_dpk2i'
               ising=41
               return
            endif
111         continue

            gp%dGv2_dpk2i(1:mp%N_slip,:)=
     &         -gp%ddIVBdt_dpk2i(1:mp%N_slip,:)*gp%dt1
            gp%dGv2_dIVB(1:mp%N_slip,1:mp%N_slip)=
     &         XInn(1:mp%N_slip,1:mp%N_slip)
     &         -gp%ddIVBdt_dIVB(1:mp%N_slip,1:mp%N_slip)*gp%dt1

            ising_1=0
            call gaussj(gp%dGv2_dIVB(1:mp%N_slip,1:mp%N_slip),
     &                  mp%N_slip,
     &                  gp%IdGv2_dIVB(1:mp%N_slip,1:mp%N_slip),ising_1)
            if(ising_1/=0)then
               write(6,*) 'non ivertable for dGv2_dIVB'
               ising=42
               return
            endif

            gp%eqM6nGv1=matmul(gp%dGv1_dIVB(:,1:mp%N_slip),
     &                         gp%IdGv2_dIVB(1:mp%N_slip,1:mp%N_slip))
            gp%eqM66Gv1=-matmul(matmul(
     &         gp%dGv1_dIVB(:,1:mp%N_slip),
     &         gp%IdGv2_dIVB(1:mp%N_slip,1:mp%N_slip)),
     &         gp%dGv2_dpk2i(1:mp%N_slip,:))+gp%dGv1_dpk2i

            ising_1=0
            call gaussj(gp%eqM66Gv1,6,gp%IeqM66Gv1,ising_1)
            if(ising_1/=0)then
               write(6,*) 'non ivertable for eqM66Gv1'
               ising=43
               return
            endif
            gp%eqMn6Gv2(1:mp%N_slip,:)=
     &         matmul(gp%dGv2_dpk2i(1:mp%N_slip,:),gp%IdGv1_dpk2i)
            gp%eqMnnGv2(1:mp%N_slip,1:mp%N_slip)=
     &         -matmul(matmul(gp%dGv2_dpk2i(1:mp%N_slip,:),
     &         gp%IdGv1_dpk2i),gp%dGv1_dIVB(:,1:mp%N_slip))
     &         +gp%dGv2_dIVB(1:mp%N_slip,1:mp%N_slip)
            ising_1=0
            call gaussj(gp%eqMnnGv2(1:mp%N_slip,1:mp%N_slip),
     &                  mp%N_slip,
     &                  gp%IeqMnnGv2(1:mp%N_slip,1:mp%N_slip),ising_1)
            if(ising_1/=0)then
               write(6,*) 'non ivertable for eqMnnGv2'
               ising=44
               return
            endif
            return
         endsubroutine cal_NRmatrix_ctx

         recursive subroutine cal_MatStiffness_ctx(gp,mp)
            implicit none
            type(gp_context), intent(inout) :: gp
            type(mat_param_set), intent(in) :: mp
            integer i,j,k,l,i1,j1,k1,l1,m,m1,is
            real(8) x1
            real(8) M1_66(6,6),M1_96(9,6),M1_99(9,9)
            real(8) M2_99(9,9),M3_99(9,9)
            real(8) M1_3333(3,3,3,3),M2_3333(3,3,3,3)
            real(8) MX1(3,3),MX2(3,3),MX3(3,3),MX4(3,3)
            real(8) STFjc_66_ctx(6,6),STFtk_66_ctx(6,6)

            gp%dIVB_dpk2i=-matmul(gp%IdGv2_dIVB,gp%dGv2_dpk2i)

            do i=1,6
            do j=1,6
               if(j<=3)then
                  gp%dCGEe_mx_dE(i,j)=
     &               2*gp%TIFp0(ib1(i),ib1(j))
     &               *gp%IFp0(ib2(j),ib2(i))
               else
                  gp%dCGEe_mx_dE(i,j)=
     &            (+2*gp%TIFp0(ib1(i),ib1(j  ))
     &               *gp%IFp0(ib2(j  ),ib2(i))
     &             +2*gp%TIFp0(ib1(i),ib1(j+3))
     &               *gp%IFp0(ib2(j+3),ib2(i)))/2
               endif
            enddo
            enddo

            gp%dpk2i_mx_dE=matmul(mp%Mstiff,gp%dCGEe_mx_dE)/2

            do is=1,mp%N_slip
               MX1=matmul(gp%IFp0,mp%Msmd(is,:,:))
               MX2=transpose(MX1)
               do i=1,6
               do j=1,6
                  M1_66(i,j)=
     &              +2*gp%TIFp0(ib1(i),ib1(j))*MX1(ib2(j),ib2(j))
     &              +2*MX2(ib1(i),ib1(j))*gp%IFp0(ib2(j),ib2(j))
               enddo
               enddo
               gp%dSTFrlx6_dE(is,:,:)=matmul(mp%Mstiff,M1_66)/2
            enddo

            do i=1,6
            do j=1,6
               gp%dGv1_dE(i,j)=-gp%dpk2i_mx_dE(i,j)
     &            +dot_product(gp%dgmdt(1:mp%N_slip),
     &            gp%dSTFrlx6_dE(1:mp%N_slip,i,j))*gp%dt1
            enddo
            enddo

            gp%dpk2i_dE=-matmul(gp%IeqM66Gv1,gp%dGv1_dE)

            M1_99=0.d0
            M2_99=0.d0
            M3_99=0.d0
            do is=1,mp%N_slip
               MX1=matmul(gp%IFp0,mp%Msmd(is,:,:))
               MX2=transpose(MX1)
               MX4=matmul(mp%Msmd(is,:,:),gp%Fp0)
               do i=1,9
                  if(i<=6) i1=i
                  if(i >6) i1=i-3
                  MX3(ib1(i),ib2(i))=
     &            +gp%ddgmdt_dpk2i(is,i1)
     &            +gp%ddgmdt_dIVB(is)*gp%dIVB_dpk2i(is,i1)
                  gp%dTdgmdt_dpk2i(is,i)=
     &            +gp%ddgmdt_dpk2i(is,i1)
     &            +gp%ddgmdt_dIVB(is)*gp%dIVB_dpk2i(is,i1)
               enddo
               do i=1,9
               do j=1,9
                  M1_99(i,j)=M1_99(i,j)
     &               +MX1(ib1(i),ib2(i))*MX3(ib1(j),ib2(j))
                  M2_99(i,j)=M2_99(i,j)
     &               +MX2(ib1(i),ib2(i))*MX3(ib1(j),ib2(j))
                  M3_99(i,j)=M3_99(i,j)
     &               +MX4(ib1(i),ib2(i))*MX3(ib1(j),ib2(j))
               enddo
               enddo
            enddo
            gp%dIFp_dE=0.d0
            gp%dTIFp_dE=0.d0
            M1_96(1:6,:)=gp%dpk2i_dE
            M1_96(7,:)=M1_96(4,:)
            M1_96(8,:)=M1_96(5,:)
            M1_96(9,:)=M1_96(6,:)
            gp%dIFp_dE=-matmul(M1_99,M1_96)*gp%dt1
            gp%dTIFp_dE=-matmul(M2_99,M1_96)*gp%dt1
            gp%dFp_dE=-matmul(M3_99,M1_96)*gp%dt1
            gp%ddgmdt_dE=matmul(gp%dTdgmdt_dpk2i,M1_96)*gp%dt1

            gp%dFp_dE(:,4:6)=gp%dFp_dE(:,4:6)*2
            gp%ddgmdt_dE(:,4:6)=gp%ddgmdt_dE(:,4:6)*2

            gp%dpk2r_dE=0.d0
            MX1=matmul(gp%pk2i_M,gp%TIFp)
            MX2=matmul(gp%IFp,gp%pk2i_M)
            do i=1,6
            do k=1,6
               do m=1,9
                  if(m<=6)m1=m
                  if(m >6)m1=m-3
                  gp%dpk2r_dE(i,k)=gp%dpk2r_dE(i,k)
     &           +XI33(ib1(i),ib1(m))*MX1(ib2(m),ib2(i))
     &              *gp%dIFp_dE(m,k)
     &           +gp%IFp(ib1(i),ib1(m))*gp%TIFp(ib2(m),ib2(i))
     &              *gp%dpk2i_dE(m1,k)
     &           +MX2(ib1(i),ib1(m))*XI33(ib2(i),ib2(m))
     &              *gp%dTIFp_dE(m,k)
               enddo
            enddo
            enddo

            do i=1,9
            do j=1,9
               if(i<=6) i1=i
               if(i >6) i1=i-3
               if(j<=6) j1=j
               if(j >6) j1=j-3
               M1_3333(ib1(i),ib2(i),ib1(j),ib2(j))=
     &            gp%dpk2r_dE(i1,j1)
            enddo
            enddo
            M2_3333=0.d0
            do i=1,3
            do j=1,3
            do k=1,3
            do l=1,3
               x1=0.d0
               do i1=1,3
               do j1=1,3
               do k1=1,3
               do l1=1,3
                  x1=x1+M1_3333(i1,j1,k1,l1)
     &           *gp%Fg(i,i1)*gp%Fg(j,j1)*gp%Fg(k,k1)*gp%Fg(l,l1)
               enddo
               enddo
               enddo
               enddo
               M2_3333(i,j,k,l)=x1/gp%det_Fg
     &              +XI33(i,k)*gp%csM0(l,j)
     &              +gp%csM0(i,k)*XI33(l,j)
     &              +gp%csM0(i,j)*XI33(k,l)
            enddo
            enddo
            enddo
            enddo
            do i=1,6
            do j=1,6
               STFjc_66_ctx(i,j)=
     &            M2_3333(ib1(i),ib2(i),ib1(j),ib2(j))
            enddo
            enddo

            do i=1,9
            do j=1,9
               if(i<=6) i1=i
               if(i >6) i1=i-3
               if(j<=6) j1=j
               if(j >6) j1=j-3
               M1_3333(ib1(i),ib2(i),ib1(j),ib2(j))=
     &            gp%dpk2r_dE(i1,j1)
            enddo
            enddo
            M2_3333=0.d0
            do i=1,3
            do j=1,3
            do k=1,3
            do l=1,3
               do k1=1,3
               do l1=1,3
                  M2_3333(i,j,k,l)=M2_3333(i,j,k,l)
     &           +M1_3333(i,j,k1,l1)*gp%Fg(k,k1)*gp%Fg(l,l1)
               enddo
               enddo
            enddo
            enddo
            enddo
            enddo
            do i=1,6
            do j=1,6
               STFtk_66_ctx(i,j)=
     &            M2_3333(ib1(i),ib2(i),ib1(j),ib2(j))
            enddo
            enddo

            gp%MatJacb=STFjc_66_ctx
            return
         endsubroutine cal_MatStiffness_ctx
      endmodule mod_stress
