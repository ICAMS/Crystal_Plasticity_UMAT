c ICAMS CP-UMAT 2026R1
c (c) 2026 by ICAMS, Ruhr University Bochum
c================================================================
c
c    Modules: mod_gaussp
c
c================================================================
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                               +
! +   Module "mod_gaussp" contribute data for one gp              +
! +                                                               +
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      module mod_gaussp
         use globalvalue, only: Nslp_mx
         implicit none
         type gp_context
            real(8) :: dt1 = 0.d0
            real(8) :: det_Fg = 0.d0, det_Fe = 0.d0, det_Fp = 0.d0
            real(8) :: Fg0(3,3) = 0.d0, Fg(3,3) = 0.d0
            real(8) :: TFg(3,3) = 0.d0, IFg(3,3) = 0.d0
            real(8) :: Fe0(3,3) = 0.d0, Fe(3,3) = 0.d0
            real(8) :: Fp0(3,3) = 0.d0, Fp(3,3) = 0.d0
            real(8) :: IFp0(3,3) = 0.d0, TIFp0(3,3) = 0.d0
            real(8) :: IFp(3,3) = 0.d0, TIFp(3,3) = 0.d0
            real(8) :: Qm(3,3) = 0.d0
            real(8) :: CGE(3,3) = 0.d0, CGEe_max(3,3) = 0.d0
            real(8) :: Ftrp(3,3) = 0.d0, IFtrp(3,3) = 0.d0
            real(8) :: pk2i(6) = 0.d0, pk2i_max(6) = 0.d0
            real(8) :: pk2i_gnd(6) = 0.d0, pk2i_eq(6) = 0.d0
            real(8) :: pk2r(6) = 0.d0, cs(6) = 0.d0, cs0(6) = 0.d0
            real(8) :: dpk2i(6) = 0.d0
            real(8) :: pk2i_M(3,3) = 0.d0, pk2r_M(3,3) = 0.d0
            real(8) :: csM0(3,3) = 0.d0, csM(3,3) = 0.d0
            real(8) :: eang00(4) = 0.d0, eang0(4) = 0.d0
            real(8) :: eang(4) = 0.d0
            real(8) :: Lp(3,3) = 0.d0, dFpdt(3,3) = 0.d0
            real(8) :: MatJacb(6,6) = 0.d0
            real(8) :: IVB0(Nslp_mx) = 0.d0, IVB(Nslp_mx) = 0.d0
            real(8) :: IVB_wcp(Nslp_mx) = 0.d0
            real(8) :: IVB_bk(Nslp_mx) = 0.d0
            real(8) :: tau(Nslp_mx) = 0.d0, dgmdt(Nslp_mx) = 0.d0
            real(8) :: detGM(Nslp_mx) = 0.d0, dIVB(Nslp_mx) = 0.d0
            real(8) :: dIVBdt(Nslp_mx) = 0.d0
            real(8) :: STFrlx6(Nslp_mx,6) = 0.d0
            real(8) :: ddgmdt_dtau(Nslp_mx) = 0.d0
            real(8) :: ddgmdt_dIVB(Nslp_mx) = 0.d0
            real(8) :: ddIVBdt_ddgmdt(Nslp_mx,Nslp_mx) = 0.d0
            real(8) :: ddIVBdt_dIVB(Nslp_mx,Nslp_mx) = 0.d0
            real(8) :: ddgmdt_dpk2i(Nslp_mx,6) = 0.d0
            real(8) :: ddIVBdt_dpk2i(Nslp_mx,6) = 0.d0
            real(8) :: Gv1(6) = 0.d0, Gv2(Nslp_mx) = 0.d0
            real(8) :: dGv1_dpk2i(6,6) = 0.d0
            real(8) :: dGv1_dIVB(6,Nslp_mx) = 0.d0
            real(8) :: dGv2_dpk2i(Nslp_mx,6) = 0.d0
            real(8) :: dGv2_dIVB(Nslp_mx,Nslp_mx) = 0.d0
            real(8) :: IdGv1_dpk2i(6,6) = 0.d0
            real(8) :: IdGv2_dIVB(Nslp_mx,Nslp_mx) = 0.d0
            real(8) :: eqM66Gv1(6,6) = 0.d0
            real(8) :: eqMnnGv2(Nslp_mx,Nslp_mx) = 0.d0
            real(8) :: IeqM66Gv1(6,6) = 0.d0
            real(8) :: IeqMnnGv2(Nslp_mx,Nslp_mx) = 0.d0
            real(8) :: eqM6nGv1(6,Nslp_mx) = 0.d0
            real(8) :: eqMn6Gv2(Nslp_mx,6) = 0.d0
            real(8) :: dIVB_dpk2i(Nslp_mx,6) = 0.d0
            real(8) :: dTdgmdt_dpk2i(Nslp_mx,9) = 0.d0
            real(8) :: dGv1_dE(6,6) = 0.d0
            real(8) :: dpk2i_dE(6,6) = 0.d0
            real(8) :: dPk2r_dE(6,6) = 0.d0
            real(8) :: dIFp_dE(9,6) = 0.d0
            real(8) :: dTIFp_dE(9,6) = 0.d0
            real(8) :: dSTFrlx6_dE(Nslp_mx,6,6) = 0.d0
            real(8) :: dCGEe_mx_dE(6,6) = 0.d0
            real(8) :: dpk2i_mx_dE(6,6) = 0.d0
            real(8) :: dFp_dE(9,6) = 0.d0
            real(8) :: ddgmdt_dE(Nslp_mx,6) = 0.d0
         end type gp_context
         
         integer :: Iexp_abq=0  !0: umat,   1: vumat
         integer :: Iexp_loc=0  !0: impl,   1: expl, for stress cal
         integer :: Imth_add=0  !0: multidcp, 1: adddcp
         real(8) :: dt1
!---------------------------------------------------------------c
!        general continuum mechanics variables                  c
!---------------------------------------------------------------c
         integer ie,ig  ! should not be global !!!
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
c     Immutable tensor constants: available before any thread enters UMAT.
         integer, parameter :: ib1(9)=[1,2,3,1,1,2,2,3,3]
         integer, parameter :: ib2(9)=[1,2,3,2,3,3,1,1,2]
c     RESHAPE repeats PAD: after the first one, N zeros then one place
c     the remaining diagonal entries N+1 positions apart. Avoid implied-DO
c     indices inside MERGE in parameter expressions (ifort 19 rejects them).
         real(8), parameter :: XI33(3,3)=reshape(
     &      [1.d0],[3,3],pad=[spread(0.d0,1,3),1.d0])
         real(8), parameter :: XI66(6,6)=reshape(
     &      [1.d0],[6,6],pad=[spread(0.d0,1,6),1.d0])
         real(8), parameter :: XI99(9,9)=reshape(
     &      [1.d0],[9,9],pad=[spread(0.d0,1,9),1.d0])
         real(8), parameter :: XInn(Nslp_mx,Nslp_mx)=reshape(
     &      [1.d0],[Nslp_mx,Nslp_mx],pad=[spread(0.d0,1,Nslp_mx),1.d0])
         real(8), parameter :: XI333(3,3,3)=reshape([
     &      0.d0,0.d0,0.d0,0.d0,0.d0,-1.d0,0.d0,1.d0,0.d0,
     &      0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,-1.d0,0.d0,0.d0,
     &      0.d0,-1.d0,0.d0,1.d0,0.d0,0.d0,0.d0,0.d0,0.d0],[3,3,3])
c     This sentinel is enabled only by OpenMP-aware compilation.
         logical, parameter :: cp_openmp_enabled=.false.
!$   &      .or. .true.

!$omp threadprivate(Iexp_abq,Iexp_loc,Imth_add,dt1,ie,ig)
!$omp threadprivate(det_Fg,det_Fe,det_Fp,Fg00,Fg0,TFg0,Fg,TFg)
!$omp threadprivate(Fe00,Fe0,TFe0,Fe,TFe,Fp00,Fp0,TFp0,Fp,TFp)
!$omp threadprivate(Fpx,Fpy,Fpz,Fppp,egx,egy,egz,egp)
!$omp threadprivate(IFg0,TIFg0,IFg,TIFg,IFe0,TIFe0,IFe,TIFe)
!$omp threadprivate(IFp0,TIFp0,IFp,TIFp,pk2i00,pk2i0,pk2i,pk2i_max)
!$omp threadprivate(pk2r,cs,dpk2i,pk2r_M,pk2i_M,cs_M,cs0)
!$omp threadprivate(Lp,Lpp,Lpx,Lpy,Lpz,DPv,pegrd_1st,pegrd_2st,dFpdt)
!$omp threadprivate(MatJacb0,MatJacb,Cstr_max,Fem,Femt,Estr_max,VEstr_max)
!$omp threadprivate(dpk2i_dt,Rvpk2i,dRvpk2i_dpk2i,IdRvpk2i_dpk2i)
!$omp threadprivate(STFjc_66,STFtk_66,CGE,CGEe_max)
!$omp threadprivate(eang00,eang0,eang,Lg_glb,Dmc,Dvc,Wmc,Lg,dDmc,dDvc,dWmc,dLg)
!$omp threadprivate(smdMc,smdMi,smdSMc,smdAMc,smdSMi,smdAMi)
!$omp threadprivate(smdVc1,smdVc2,smdVi1,smdVi2,vd_slp,vl_slp,vn_slp)
!$omp threadprivate(STFrlx6,STFrlx33,dSTFrlx6_dE)
!$omp threadprivate(STFec26,STFec29,STFec43,STFei26,STFei29,STFei43)
!$omp threadprivate(dcsdt,dcsdt_max,dcs_max,csM0,csM,dcs)
!$omp threadprivate(shgrd00,shgrd0,shgrd,dIVB,IVB00,IVB0,IVB)
!$omp threadprivate(IVB_eff,dIVBdt,shrt_gnd,tau,dgmdt,detGM)
!$omp threadprivate(vd_ltc,vl_ltc,vn_ltc,IVB_ini,refv_IVB,refv_pk2i)
!$omp threadprivate(Rho_gnd,IVB_gnd,pk2i_gnd,pk2i_eq,Ftrp,IFtrp,IVB_trp,IVB_wcp)
!$omp threadprivate(IVB_bk,IVB_cl,IVB_m,IVB_kw,fx,fy,fz,fpp)
!$omp threadprivate(epsx,epsy,epsz,epsp,pk2i_intp,pk2i_eqp)
!$omp threadprivate(pk2i_intx,pk2i_inty,pk2i_intz,pk2i_eqx,pk2i_eqy,pk2i_eqz)
!$omp threadprivate(ddgmdt_dtau,ddgmdt_dIVB,ddIVBdt_ddgmdt,ddIVBdt_dIVB)
!$omp threadprivate(ddgmdt_dpk2i,ddIVBdt_dpk2i,GV1,GV2,dGv1_dpk2i,dGv1_dIVB)
!$omp threadprivate(dGv2_dpk2i,dGv2_dIVB,IdGv1_dpk2i,IdGv2_dIVB)
!$omp threadprivate(eqM66Gv1,eqMnnGv2,IeqM66Gv1,IeqMnnGv2,eqM6nGv1,eqMn6Gv2)
!$omp threadprivate(dIVB_dpk2i,dTdgmdt_dpk2i,dGv1_dE,dpk2i_dE,dPk2r_dE)
!$omp threadprivate(dIFp_dE,dTIFp_dE,dCGEe_mx_dE,dpk2i_mx_dE,dFp_dE,ddgmdt_dE)

      contains
         !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         !+                                                                           +
         !+      This subroutine initializes the constants                            +
         !+                                                                           +
         !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         subroutine mod_gspt_ini
            implicit none
c           Compatibility entry point: constants need no runtime writes.
c           Refuse silently shared Euler angles/material caches in bad builds.
            if (.not. cp_openmp_enabled) then
               error stop 'CP-UMAT requires OpenMP: use -fopenmp or -qopenmp /Qopenmp at compile and link'
            endif
         endsubroutine mod_gspt_ini
      endmodule mod_gaussp
