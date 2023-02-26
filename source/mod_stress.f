! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                               +
! +   Module "mod_stress" contribute stress calculation           +
! +                                                               +
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      module mod_stress
         use mod_gaussp
         use mod_alloys
         use mod_wkcoup
         implicit none
         integer Icurrent_dt,iNRloop
         integer :: Nnr_max=200                     
         real(8) :: toler_NRloop=1.d-10         
      contains
         !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         !+                                                                           +
         !+      This subroutine calculates pk2i, IVB by newton raphson algoriths     +
         !+                                                                           +
         !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         subroutine cal_stress_add(ising)
            implicit none
            integer i,j,k,l,m,n,i1,j1,k1,l1,m1,n1,is,js
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
            do is=1,Nslp
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
               do is=1,Nslp
                  tau(is)=dot_product(cs,smdVc2(is,:))
               enddo
               call cal_flow_harden(ialloy,
     &              Nslp_mx,Iexp_loc,IVB_wcp,IVB_bk,IVB_cl,IVB_m,IVB_kw,
     &              tau,IVB,dgmdt,ddgmdt_dtau,ddgmdt_dIVB,
     &              dIVBdt,ddIVBdt_ddgmdt,ddIVBdt_dIVB,ising_1)
               do is=1,Nslp
                  ddgmdt_dpk2i(is,:)=ddgmdt_dtau(is)*smdVc2(is,:)
               enddo
               ddIVBdt_dpk2i(1:Nslp,:)=
     &         matmul(ddIVBdt_ddgmdt(1:Nslp,1:Nslp), 
     &         ddgmdt_dpk2i(1:Nslp,:))

               if(ising_1/=0)then
                  ising=ising_1  !!! shearrate non-cvg,stop
                  return
               endif
               Gv1=cs-cs0-dcs_max  
     &         +matmul(transpose(STFrlx6(1:Nslp,:)),dgmdt(1:Nslp))*dt1
               Gv2(1:Nslp)=IVB(1:Nslp)-IVB0(1:Nslp)-dIVBdt(1:Nslp)*dt1
               x1=sum(dabs(Gv1))/refv_pk2i
               x2=sum(dabs(Gv2(1:Nslp)))/refv_IVB
      
c               call pm(Fg0,3,3)     
c               call pm(Fg,3,3)     
c               call pm(Fe0,3,3)     
c               call pm(Fe,3,3)     
c               print '(i5,3e14.4)',iNRloop,x1,x2,toler_NRloop

               if(x1+x2<toler_NRloop .or. Iexp_loc==1)then
                  cs=cs0+dcs_max
     &            -matmul(transpose(STFrlx6(1:Nslp,:)),
     &            dgmdt(1:Nslp))*dt1
                  IVB(1:Nslp)=IVB0(1:Nslp)+dIVBdt(1:Nslp)*dt1
                  goto 101 !Converge jump out NR loop and continue
               endif
               ising_1=0
               call cal_NRmatrix(ising_1)
               if(ising_1/=0)then
                  ising=ising_1  !NR matrix cannot be got,stop
                  return
               endif

               dcs=matmul(IeqM66Gv1,
     &         -(Gv1-matmul(eqM6nGv1(:,1:Nslp),Gv2(1:Nslp))))
               dIVB(1:Nslp)=matmul( IeqMnnGv2(1:Nslp,1:Nslp),
     &         -(Gv2(1:Nslp)-matmul(eqMn6Gv2(1:Nslp,:),Gv1)) )

               cs = cs + dcs
               IVB(1:Nslp) = IVB(1:Nslp) + dIVB(1:Nslp) 

            enddo
            ising=5  !Non-coverge for Nnr_max, stop current time step
            return
101         continue
            !1--------------------------------------------------------1
            !1   End of newton raphson method loop                    1
            !1--------------------------------------------------------1
            do i=1,3
            do j=1,3
               Lp(i,j)=dot_product(dgmdt(1:Nslp),smdMc(1:Nslp,i,j))
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
               do is=1,Nslp
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
                  do is=1,Nslp
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
         subroutine cal_stress_mul(ising)
            implicit none
            integer i,j,k,l,m,n,i1,j1,k1,l1,m1,n1,is,js
            integer ising_1,ising
            real(8) x1,x2,x3,x4,y1,z1
            real(8) M1_66(6,6),M1_96(9,6),M1_99(9,9)
            real(8) M2_66(6,6),M2_96(9,6),M2_99(9,9)
            real(8) M1_3333(3,3,3,3),M2_3333(3,3,3,3)
            real(8) MX1(3,3),MX2(3,3),MX3(3,3),MX4(3,3)

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
     &            +STFei26(i,j)*MX1(ib1(j),ib2(j))
               enddo
            enddo
            !1------------------------------------------------------------------1
            !1   Calculate STFrlx6=C*( CGE_max*smdMi + (CGE_max*smdMi)^T )/2    1
            !1------------------------------------------------------------------1
            do is=1,Nslp
               MX1=matmul(CGEe_max,smdMi(is,:,:))
               MX2=matmul( matmul(transpose(IFtrp), !==> add trip effect
     &               (MX1+transpose(MX1))/2), IFtrp )
               do i=1,6
                  STFrlx6(is,i)=0
                  do j=1,6
                     STFrlx6(is,i)=STFrlx6(is,i)
     &               +STFei26(i,j)*MX2(ib1(j),ib2(j))
                  enddo
               enddo
            enddo
            !1------------------------------------------------------1
            !1    Begin newton raphson method to solve pk2i, IVB    1
            !1------------------------------------------------------1
            do iNRloop=1,Nnr_max

            if (Ialloy==5) then
                pk2i_eqx = pk2i + pk2i_intx   !==> add misfit stress
                pk2i_eqy = pk2i + pk2i_inty
                pk2i_eqz = pk2i + pk2i_intz
                pk2i_eqp = pk2i + pk2i_intp 
 
               do is=1,12
                  tau(is)=dot_product(pk2i_eqx,smdVi2(is,:))
               enddo  
               do is=13,24
                  tau(is)=dot_product(pk2i_eqy,smdVi2(is,:))
               enddo   
               do is=25,36
                  tau(is)=dot_product(pk2i_eqz,smdVi2(is,:))
               enddo   
               do is=37,48
                  tau(is)=dot_product(pk2i_eqp,smdVi2(is,:))
               enddo
               
               do is=53,54
                  tau(is)=dot_product(pk2i_eqx,smdVi2(is,:))
               enddo 
               do is=51,52
                  tau(is)=dot_product(pk2i_eqy,smdVi2(is,:))
               enddo 
               do is=49,50
                  tau(is)=dot_product(pk2i_eqz,smdVi2(is,:))
               enddo 
               do is=55,60
                  tau(is)=dot_product(pk2i_eqp,smdVi2(is,:))
               enddo

             endif
             
             if (Ialloy/=5) then
             if (Iwkcoup_grad==0) then
                pk2i_eq = pk2i
             else
                pk2i_eq = pk2i + pk2i_gnd
             endif
               do is=1,Nslp
                  tau(is)=dot_product(pk2i_eq,smdVi2(is,:))
               enddo 
             endif


               ising_1=0
               call cal_flow_harden(ialloy,
     &              Nslp_mx,Iexp_loc,IVB_wcp,IVB_bk,IVB_cl,IVB_m,IVB_kw,
     &              tau,IVB,dgmdt,ddgmdt_dtau,ddgmdt_dIVB,
     &              dIVBdt,ddIVBdt_ddgmdt,ddIVBdt_dIVB,ising_1)
               do is=1,Nslp
                  ddgmdt_dpk2i(is,:)=ddgmdt_dtau(is)*smdVi2(is,:)
               enddo
               ddIVBdt_dpk2i(1:Nslp,:)=
     &         matmul(ddIVBdt_ddgmdt(1:Nslp,1:Nslp), 
     &         ddgmdt_dpk2i(1:Nslp,:))

               if(ising_1/=0)then
                  ising=ising_1  !!! shearrate non-cvg,stop
                  return
               endif
               Gv1=+pk2i-pk2i_max
     &         +matmul(transpose(STFrlx6(1:Nslp,:)),dgmdt(1:Nslp))*dt1
               Gv2(1:Nslp)=+IVB(1:Nslp)-IVB0(1:Nslp)-dIVBdt(1:Nslp)*dt1
               x1=sum(dabs(Gv1))/refv_pk2i
               x2=sum(dabs(Gv2(1:Nslp)))/refv_IVB


               if(x1+x2<toler_NRloop .or. Iexp_loc==1)then
                  pk2i=+pk2i_max-matmul(transpose(STFrlx6(1:Nslp,:)),
     &                 dgmdt(1:Nslp))*dt1
                  IVB(1:Nslp)=IVB0(1:Nslp)+dIVBdt(1:Nslp)*dt1
                  goto 101 !Converge jump out NR loop and continue
               endif
               ising_1=0
               call cal_NRmatrix(ising_1)
               if(ising_1/=0)then
                  ising=ising_1  !NR matrix cannot be got,stop
                  return
               endif
		
               dpk2i=matmul(IeqM66Gv1,
     &         -(Gv1-matmul(eqM6nGv1(:,1:Nslp),Gv2(1:Nslp))))
               dIVB(1:Nslp)=matmul( IeqMnnGv2(1:Nslp,1:Nslp),
     &         -(Gv2(1:Nslp)-matmul(eqMn6Gv2(1:Nslp,:),Gv1)) )

               pk2i=pk2i+dpk2i
               IVB(1:Nslp) =IVB(1:Nslp) + dIVB(1:Nslp) 

            enddo
            ising=5  !Non-coverge for Nnr_max, stop current time step
            return
101         continue
            !1--------------------------------------------------------1
            !1   End of newton raphson method, solve cauchy stress    1
            !1--------------------------------------------------------1 

            if (Iwkcoup_bk/=0) then
             do is=1,Nslp
              fem_dgmdt(ie,ig,is)=dgmdt(is)
             enddo
            endif
         

            if (Ialloy==5) then

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
               Lpx(i,j)=dot_product(dgmdt(1:12),smdMi(1:12,i,j))
     &               + dot_product(dgmdt(53:54),smdMi(53:54,i,j))
  
               Lpy(i,j)=dot_product(dgmdt(13:24),smdMi(13:24,i,j))
     &               + dot_product(dgmdt(51:52),smdMi(51:52,i,j))

               Lpz(i,j)=dot_product(dgmdt(25:36),smdMi(25:36,i,j))
     &               + dot_product(dgmdt(49:50),smdMi(49:50,i,j))

               Lpp(i,j)=dot_product(dgmdt(37:48),smdMi(37:48,i,j))
     &               + dot_product(dgmdt(55:60),smdMi(55:60,i,j))
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

            else

            call icams_conv6to33(pk2i,ib1,ib2,pk2i_M)

            do i=1,3
            do j=1,3
               Lp(i,j)=dot_product(dgmdt(:),smdMi(:,i,j))
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
         endsubroutine cal_stress_mul

         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         !+                                                                                  +
         !+   THIS ROUTINE calculates tangent matrix for newton raphson algoriths            +
         !+                                                                                  +
         !+   #1: dG1_dpk2i,dG1_dIVB,IdG1_dpk2i                                              +
         !+   #2: dG2_dpk2i,dG2_dIVB,IdG2_dIVB                                               +
         !+   #3: eqM66Gv1,eqMnnGv2,IeqM66Gv1,IeqMnnGv2                                      +
         !+                                                                                  +
         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         subroutine cal_NRmatrix(ising)
            implicit none
            integer i,j,k,l,m,n,i1,j1,k1,l1,m1,n1,is,js
            integer ising_1,ising
            real(8) x1,y1,z1
            real(8) M1_66(6,6),M1_96(9,6),M1_99(9,9)
            real(8) M2_66(6,6),M2_96(9,6),M2_99(9,9)
            real(8) M1_3333(3,3,3,3),M2_3333(3,3,3,3)
            real(8) MX1(3,3),MX2(3,3),MX3(3,3),MX4(3,3)
            !-----------------------------------------------------------------------------
            ! calculate dG1_dpk2i(6,6),dG1_dIVB(6,:),IdG1_dpk2i(6,6)
            !-----------------------------------------------------------------------------
            ising=0
            dGv1_dpk2i=XI66+matmul(transpose(STFrlx6(1:Nslp,:)),
     &                ddgmdt_dpk2i(1:Nslp,:))*dt1
            do is=1,Nslp
               dGv1_dIVB(:,is)=STFrlx6(is,:)*ddgmdt_dIVB(is)*dt1
            enddo

            ising_1=0
            call gaussj(dGv1_dpk2i,6,idGv1_dpk2i,ising_1)
            if(ising_1/=0)then
               write(6,*) 'non ivertable for dGv1_dpk2i'
               ising=41  !dGv1_dpk2i non-invertable,stop
               return
            endif
            !-----------------------------------------------------------------------------
            ! calculate dG2_dpk2i(:,6),dG2_dIVB(:,:),IdG2_dIVB(:,:)
            !-----------------------------------------------------------------------------
            dGv2_dpk2i(1:Nslp,:)=-ddIVBdt_dpk2i(1:Nslp,:)*dt1
            dGv2_dIVB(1:Nslp,1:Nslp) =+XInn(1:Nslp,1:Nslp)
     &                            -ddIVBdt_dIVB(1:Nslp,1:Nslp) *dt1

            ising_1=0
            call gaussj( dGv2_dIVB(1:Nslp,1:Nslp),Nslp,
     &                  idGv2_dIVB(1:Nslp,1:Nslp),ising_1)
            if(ising_1/=0)then
               write(6,*) 'non ivertable for dGv2_dIVB'
               ising=42  !dGv2_dIVB non-invertable,stop
               return
            endif
            !-----------------------------------------------------------------------------
            !-----Calculate eqM66Gv1,eqMnnGv2,IeqM66Gv1,IeqMnnGv2
            !-----------------------------------------------------------------------------
            eqM6nGv1=matmul(dGv1_dIVB(:,1:Nslp),
     &                     IdGv2_dIVB(1:Nslp,1:Nslp))
            eqM66Gv1=-matmul(matmul(dGv1_dIVB(:,1:Nslp),
     &         IdGv2_dIVB(1:Nslp,1:Nslp)),dGv2_dpk2i(1:Nslp,:))
     &         +dGv1_dpk2i

            ising_1=0
            call gaussj(eqM66Gv1,6,IeqM66Gv1,ising_1)
            if(ising_1/=0)then
               write(6,*) 'non ivertable for eqM66Gv1'
               ising=43  !eqM66Gv1 non-invertable,stop
               return
            endif
            eqMn6Gv2(1:Nslp,:)=matmul(dGv2_dpk2i(1:Nslp,:),IdGv1_dpk2i)
            eqMnnGv2(1:Nslp,1:Nslp)=-matmul(matmul(dGv2_dpk2i(1:Nslp,:),
     &                          IdGv1_dpk2i),dGv1_dIVB(:,1:Nslp))
     &                         +dGv2_dIVB(1:Nslp,1:Nslp)
            ising_1=0
            call gaussj( eqMnnGv2(1:Nslp,1:Nslp),Nslp,
     &                  IeqMnnGv2(1:Nslp,1:Nslp),ising_1)
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
         subroutine cal_MatStiffness
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
c         #1.1: dIVB_dpk2i(Nslp,6)
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
c         #1.2.3:   dSTFrlx6_dE(Nslp,6,6)
c-----------------------------------------------------------------------------
         do is=1,Nslp
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
         do is=1,Nslp
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
      endmodule mod_stress

