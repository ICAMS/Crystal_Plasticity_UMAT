! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                               +
! +   Module "mod_gaussp" contribute data for one gp              +
! +                                                               +
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      module mod_gaussp
         implicit none
         integer  Nslp
         integer, parameter :: Nslp_mx=120            
c         integer :: Iexp_abq=0  !0: umat,   1: vumat
c         integer :: Iexp_loc=0  !0: impl,   1: expl, for stress cal
c         integer :: Imth_add=0  !0: multidcp, 1: adddcp
         integer :: Iexp_abq=0 
         integer :: Iexp_loc=0 
         integer :: Imth_add=0 
         integer :: first_call=1
         integer ia,iblk,iphase,ialloy
         real(8) dt0,dt1,zeit0,zeit
!---------------------------------------------------------------c
!        general continuum mechanics variables                  c
!---------------------------------------------------------------c
         integer ie,ig
         real(8) det_Fg,det_Fe,det_Fp
         real(8) Fg00(3,3),Fg0(3,3),TFg0(3,3),Fg(3,3),TFg(3,3)
         real(8) Fe00(3,3),Fe0(3,3),TFe0(3,3),Fe(3,3),TFe(3,3)
         real(8) Fp00(3,3),Fp0(3,3),TFp0(3,3),Fp(3,3),TFp(3,3)
         real(8) Fpx(3,3),Fpy(3,3),Fpz(3,3),Fppp(3,3)
         real(8) egx(3,3),egy(3,3),egz(3,3),egp(3,3)
         real(8) IFg0(3,3),TIFg0(3,3),IFg(3,3),TIFg(3,3)
         real(8) IFe0(3,3),TIFe0(3,3),IFe(3,3),TIFe(3,3)
         real(8) IFp0(3,3),TIFp0(3,3),IFp(3,3),TIFp(3,3)
         real(8) pk2i00(6),pk2i0(6),pk2i(6),pk2i_max(6)
         real(8) pk2r(6),cs(6),dpk2i(6)
         real(8) pk2r_M(3,3),pk2i_M(3,3),cs_M(3,3),cs0(6)
         real(8) Lp(3,3),Lpp(3,3),Lpx(3,3),Lpy(3,3),Lpz(3,3)
         real(8) DPv(6),pegrd_1st(6,3),pegrd_2st(6,3,3)
         real(8) dFpdt(3,3)
         real(8) MatJacb0(6,6),MatJacb(6,6)
         real(8) Cstr_max(3,3),Fem(3,3),Femt(3,3)
         real(8) Estr_max(3,3),VEstr_max(6)
         real(8) dpk2i_dt(6),Rvpk2i(6)
         real(8) dRvpk2i_dpk2i(6,6),IdRvpk2i_dpk2i(6,6)
         real(8) STFjc_66(6,6),STFtk_66(6,6)
         real(8) CGE(3,3),CGEe_max(3,3)
         real(8) eang00(4),eang0(4),eang(4),Lg_glb(3,3)
         real(8) Dmc(3,3),Dvc(6),Wmc(3,3),Lg(3,3)
         real(8) dDmc(3,3),dDvc(6),dWmc(3,3),dLg(3,3)
         real(8) smdMc(Nslp_mx,3,3),smdMi(Nslp_mx,3,3)
         real(8) smdSMc(Nslp_mx,3,3),smdAMc(Nslp_mx,3,3)
         real(8) smdSMi(Nslp_mx,3,3),smdAMi(Nslp_mx,3,3)
         real(8) smdVc1(Nslp_mx,6),smdVc2(Nslp_mx,6)
         real(8) smdVi1(Nslp_mx,6),smdVi2(Nslp_mx,6)
         real(8) vd_slp(Nslp_mx,3)
         real(8) vl_slp(Nslp_mx,3)
         real(8) vn_slp(Nslp_mx,3)
         real(8) STFrlx6(Nslp_mx,6),STFrlx33(3,3,Nslp_mx)
         real(8) dSTFrlx6_dE(Nslp_mx,6,6)
         real(8) STFec26(6,6),STFec29(9,9),STFec43(3,3,3,3)
         real(8) STFei26(6,6),STFei29(9,9),STFei43(3,3,3,3)
         real(8) dcsdt(6),dcsdt_max(6),dcs_max(6)
         real(8) csM0(3,3),csM(3,3),dcs(6)
!---------------------------------------------------------------c
!        microstructure(slip system based) relative variables   c
!---------------------------------------------------------------c
         real(8) shgrd00(Nslp_mx,13),shgrd0(Nslp_mx,13)
         real(8) shgrd(Nslp_mx,13),dIVB(Nslp_mx)
         real(8) IVB00(Nslp_mx),IVB0(Nslp_mx),IVB(Nslp_mx)
         real(8) IVB_eff(Nslp_mx),dIVBdt(Nslp_mx)
         real(8) shrt_gnd(Nslp_mx)
         real(8) tau(Nslp_mx)
         real(8) dgmdt(Nslp_mx),detGM(Nslp_mx)
         real(8) vd_ltc(Nslp_mx,3)
         real(8) vl_ltc(Nslp_mx,3)
         real(8) vn_ltc(Nslp_mx,3)
         real(8) IVB_ini(Nslp_mx),refv_IVB,refv_pk2i
!---------------------------------------------------------------c
!        strain gradient effect                                 c
!---------------------------------------------------------------c
         real(8) Rho_gnd(Nslp_mx)
         real(8) IVB_gnd(Nslp_mx)
         real(8) pk2i_gnd(6),pk2i_eq(6)
!---------------------------------------------------------------c
!        phase transformation effect                            c
!---------------------------------------------------------------c
         real(8) Ftrp(3,3)
         real(8) IFtrp(3,3)
         real(8) IVB_trp(Nslp_mx),IVB_wcp(Nslp_mx)
!---------------------------------------------------------------c
!        Kinematic hardening effect                            c
!---------------------------------------------------------------c
         real(8) IVB_bk(Nslp_mx)
!---------------------------------------------------------------c
!        climb, matrix dislocation and KW effect                c
!---------------------------------------------------------------c
         real(8) IVB_cl(48),IVB_m(48),IVB_kw(48)
!---------------------------------------------------------------c
!        misfit stress effect                                   c
!---------------------------------------------------------------c
         real(8) fx,fy,fz,fpp
         real(8) epsx(6),epsy(6),epsz(6),epsp(6)
         real(8) pk2i_intp(6),pk2i_eqp(6)
         real(8) pk2i_intx(6),pk2i_inty(6),pk2i_intz(6)
         real(8) pk2i_eqx(6),pk2i_eqy(6),pk2i_eqz(6)
!---------------------------------------------------------------c
!        vaules for newton-raphson algoriths                    c
!---------------------------------------------------------------c
         real(8) ddgmdt_dtau(Nslp_mx)
         real(8) ddgmdt_dIVB(Nslp_mx)
         real(8) ddIVBdt_ddgmdt(Nslp_mx,Nslp_mx)
         real(8) ddIVBdt_dIVB(Nslp_mx,Nslp_mx)
         real(8) ddgmdt_dpk2i(Nslp_mx,6)
         real(8) ddIVBdt_dpk2i(Nslp_mx,6)
         real(8) GV1(6),GV2(Nslp_mx)
         real(8) dGv1_dpk2i(6,6),dGv1_dIVB( 6,Nslp_mx)
         real(8) dGv2_dpk2i(Nslp_mx,6),dGv2_dIVB(Nslp_mx,Nslp_mx)
         real(8) IdGv1_dpk2i(6,6),IdGv2_dIVB(Nslp_mx,Nslp_mx)
         real(8) eqM66Gv1(6,6),eqMnnGv2(Nslp_mx,Nslp_mx)
         real(8) IeqM66Gv1(6,6),IeqMnnGv2(Nslp_mx,Nslp_mx)
         real(8) eqM6nGv1(6,Nslp_mx),eqMn6Gv2(Nslp_mx,6)
!---------------------------------------------------------------c
!        vaules for material tangent calcuation                 c
!---------------------------------------------------------------c
         real(8) dIVB_dpk2i(Nslp_mx,6)
         real(8) dTdgmdt_dpk2i(Nslp_mx,9)
         real(8) dGv1_dE(6,6)
         real(8) dpk2i_dE(6,6),dPk2r_dE(6,6)
         real(8) dIFp_dE(9,6),dTIFp_dE(9,6)
         real(8) dSTFsh_dE(Nslp_mx,6,6)
         real(8) dCGEe_mx_dE(6,6)
         real(8) dpk2i_mx_dE(6,6)
         real(8) dFp_dE(9,6),ddgmdt_dE(Nslp_mx,6)
         integer ib1(9),ib2(9)
         real(8) XI33(3,3), XI66(6,6), XI99(9,9), XInn(Nslp_mx,Nslp_mx)
         real(8) XI333(3,3,3)
      contains
         !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         !+                                                                           +
         !+      This subroutine initial the constants                                +
         !+                                                                           +
         !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         subroutine mod_gspt_ini
            implicit none
            integer i,j
            real(8) x1,x2
            XI33=0
            do i=1,3
               XI33(i,i)=1                           
            enddo         
            XI333=0
            XI333(1,2,3)=+1
            XI333(1,3,2)=-1
            XI333(2,1,3)=-1
            XI333(2,3,1)=+1
            XI333(3,1,2)=+1
            XI333(3,2,1)=-1
            XI66=0
            do i=1,6
               XI66(i,i)=1                           
            enddo         
            XI99=0
            do i=1,9
               XI99(i,i)=1                           
            enddo         
            XInn=0
            do i=1,Nslp_mx
               XInn(i,i)=1                           
            enddo         
            !-->Marc stress vector sequence
            !IB1=[1,2,3,1,2,1,2,3,3]
            !IB2=[1,2,3,2,3,3,1,2,1]
            !-->Abaqus stress vector sequence
            IB1=[1,2,3,1,1,2,2,3,3]
            IB2=[1,2,3,2,3,3,1,1,2]
         endsubroutine mod_gspt_ini
      endmodule mod_gaussp

