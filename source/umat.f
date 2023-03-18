c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c +                                                                         +
c +   Define global values: Material and model parameters                   +
c +                                                                         +
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      MODULE globalvalue
      IMPLICIT NONE

      CHARACTER(len=12) :: abaqus = 'abq6142'

      ! Parameters and constants                     ! UNIT         Explanation
C     -------------------------------------------------------------------------------                                                                                                                  
      REAL(8), PARAMETER :: R_gas     = 8.314d0      ! Nm/K/mol     Gas constant                        
      REAL(8), PARAMETER :: NAvg      = 6.022141d23  ! 1/mol        Avogardo Number
      REAL(8)            :: temp_cur                 ! K            current temperature
c     temperature-dependent material parameters for austenite
      REAL(8)            :: c11_gl, c12_gl, c44_gl   ! MPa          elast. constants
      REAL(8)            :: tau0_gl                  ! MPa          CRSS
      REAL(8)            :: Adir_gl                  ! MPa          kin. hard. parameter
      
      END MODULE globalvalue
      
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c +                                                                         +
c +   Include other subroutines                                             +
c +                                                                         +
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      include "mod_gaussp.f"            
      include "mod_alloys.f"            
      include "mod_wkcoup.f"
      include "mod_stress.f"            
      include "mod_smthgd.f"
      include "other_code.f"

c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c +                                                                         +
c +   Uexternaldb...                                                        +
c +       lop = 0 beginning of analysis                                     +
c +             1 start of increment                                        +
c +             2 end of increment                                          +
c +             3 end of analysis                                           +
c +             4 beginning of restart                                      +
c +                                                                         +
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine uexternaldb(lop,lrestart,time,dtime,kstep,kinc)
      use globalvalue
      use mod_gaussp
      use mod_wkcoup
      use mod_alloys
      implicit none
      integer lop,lrestart,kstep,kinc
      real(8) time(2),dtime
      real(8) x1,x2,x3,xx
      integer i,j,k

      if(lop==0)then   !ini
         print*, kinc,dtime,time(2), 'First call, initialization...'
         temp_cur= 0.d0  ! Can initial temp be queried in uexternaldb???
         call mod_gspt_ini
         call mod_alls_ini(IB1,IB2)       
         call mod_wkcp_ini(IB1,IB2)
        
      endif
 
      if(lop==1)then
         print*, kinc,dtime,time(2),temp_cur, 'start of time step'
      endif

      if(lop==2)then
         print*, kinc,dtime,time(2),temp_cur, 'call wkcoup_evo subroutine'
         call wkcoup_evolution(dtime)
      endif

      return
      end

c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c +                                                                         +
c +  Abaqus Umat, Abaqus Umat, Abaqus Umat, Abaqus Umat, Abaqus Umat,       +
c +                                                                         +
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine umat(stress,statev,ddsdde,
     &                sse,spd,scd,rpl,
     &                ddsddt,drplde,drpldt,
     &                stran,dstran,
     &                time,dtime,
     &                temp,dtemp,
     &                predef,dpred,
     &                cmname,
     &                ndi,nshr,ntens,nstatv,
     &                props,nprops,coords,drot,
     &                pnewdt,
     &                celent,
     &                dfgrd0,dfgrd1,
     &                noel,npt,layer,kspt,kstep,kinc)
c---------------------------------------------------------------------------
C=======================================================================
C    VARIABLES
C ----------------------------------------------------------------------
C=======================================================================
C
C    STATE VARIABLES
C ----------------------------------------------------------------------
C   SDV01-03   : eang0 (rad)
C   SDV05-10   : cs0
C   SDV11-16   : pk2i0
C   SDV17-22   : pk2i0_gnd
C   SDV23-31   : Fe0
C   SDV32-40   : Fp0
C   SDV41-149  : IVB0, IVB_bk
C   SDV150-155 : elastic strain Ee
C   SDV156-161 : plastic strain in crystal coordinates
C   SDV162-164 : C11, C12, C44 (MPa); temperature-dependent elastic constants
C   SDV165     : tau0 (MPa); temperature-dependent critical resolved shear stress
C   SDV166     : Adir (MPa); temperature-dependent parameter for isotropic hardening
C
C    Material Properties
C ----------------------------------------------------------------------
C    PROPS(1)     : Material identifier -> ialloys
C    PROPS(2...4) : Euler angles (in rad) of grain orientation -> eang00
C 
C ----------------------------------------------------------------------
C=======================================================================
      use globalvalue
      use mod_stress
      use mod_wkcoup
      use mod_gaussp
      
      implicit none

      character*80  cmname !user defined material name
      integer ndi    !number of stress components
      integer nshr   !number of engineering shear stress components
      integer ntens  !size of the stress array (ndi + nshr)
      integer nstatv !number state variables
      integer nprops !number of material constants
      integer layer  !layer number
      integer kspt   !section point number within the current layer
      integer kstep  !step number
      integer noel   !element number
      integer npt    !integration point number
      integer kinc   !increment number
c---------------------------------------------------------------------------
      real(8) drpldt    !jacobian drpl_dt
      real(8) dtime     !time increment dt
      real(8) temp      !temperature
      real(8) dtemp     !increment of temperature.
      real(8) celent    !characteristic element length
      real(8) sse       !specific elastic strain energy
      real(8) spd       !specific plastic dissipation
      real(8) scd       !specific creep dissipation
      real(8) rpl       !volumetric heat generation per unit time
      real(8) pnewdt    !dt_next/dt_now
c---------------------------------------------------------------------------
      real(8) ddsdde(ntens,ntens)  !jacobian ds_de
      real(8) statev(nstatv)       !state variables
      real(8) props (nprops)       !material constants 
      real(8) ddsddt(ntens)        !jacobian ds_dt
      real(8) drplde(ntens)        !jacobian drpl_de
      real(8) stress(ntens)        !stress tensor
      real(8) stran (ntens)        !strains at t0
      real(8) dstran(ntens)        !strain increments
      real(8) dfgrd0(3,3)          !deformation gradient at t0
      real(8) dfgrd1(3,3)          !deformation gradient at t0+dt
      real(8) drot  (3,3)          !rotation increment matrix
      real(8) coords(3)            !coordinates of this point
      real(8) time  (2)            !1:step time; 2:total time, At t0
      real(8) predef(1)            !predefined field variables at t0
      real(8) dpred (1)            !incr of predefined field vrbs
c---------------------------------------------------------------------------
      integer ising,icut
      integer icol,iit_gnd,ix
      integer i,j,k,ii,jj,is,icomp
      real(8) x1,x2,x3,xx
      real(8) mx33_1(3,3)
      real(8) phi1,phi,phi2,Qm(3,3)
      real(8) intLp(3,3),el_st(3,3)
c----------------------------------------------------------------------------
      ialloy = int(props(1))
      eang00(1:3) = props(2:4)
      ie=noel
      ig=npt
      dt1=dtime
      Fg0=dfgrd0
      Fg =dfgrd1
      if (ie.gt.Tnel) then
        print*, 'too many elements, change Tnel'
        stop
      endif

      temp_cur = temp  ! set global variable to current temperature
cc      call set_param_temp(temp)
c
c---------------------------------------------------------c
c     state initialization values for ABQUS statev        c
c---------------------------------------------------------c
      call get_mat_parameter(ialloy,
     &              Nslp_mx,Nslp,STFei26,
     &              smdMi,smdSMi,smdAMi,smdVi1,smdVi2,
     &              vd_slp,vl_slp,vn_slp,
     &              refv_pk2i,refv_IVB,IVB_ini)
      if(time(2) < 1.d-8)then
         call icams_eang2Q(eang00(1),eang00(2),eang00(3),mx33_1)
         statev( 1: 3)=eang00(1:3)               !!-->01-03 eang0
         statev( 5:10)=0                         !!-->05-10 cs0
         statev(11:16)=0                         !!-->11-16 pk2i0
         statev(17:22)=0                         !!-->17-22 pk2i0_gnd
         do i=1,9
            statev(22+i)=mx33_1(ib1(i),ib2(i))   !!-->23-31 Fe0
            statev(31+i)=mx33_1(ib2(i),ib1(i))   !!-->32-40 Fp0
         enddo
         statev(40+0*Nslp+1 : 40+1*Nslp)=IVB_ini(1:Nslp) !!-->41-100 IVB0
      endif

c---------------------------------------------------------c
c     get values from ABQUS statev for this time step     c
c---------------------------------------------------------c
      eang0=statev( 1: 4)
      cs0  =statev( 5:10)
      pk2i =statev(11:16)
      do i=1,9
         j=i; if(i>6) j=i-3
         csM0(ib1(i),ib2(i))=cs0(j)
      enddo
      do i=1,9
         Fe0(ib1(i),ib2(i))=statev(22+i)
         Fe (ib1(i),ib2(i))=statev(22+i)
         Fp0(ib1(i),ib2(i))=statev(31+i)
         Fp (ib1(i),ib2(i))=statev(31+i)
      enddo
      do is=1,Nslp
         IVB0(is)=statev(40+0*Nslp+is)
         IVB (is)=statev(40+0*Nslp+is)
      enddo

c--------------------------------------------------------c
c     considering weak couping effect: trip, gnd, int... c
c--------------------------------------------------------c
      
      Rho_gnd=0
      IVB_gnd=0
      pk2i_gnd=0
      Ftrp=XI33
      IFtrp=XI33
      IVB_trp=0
      pk2i_intp=0
      pk2i_intx=0
      pk2i_inty=0
      pk2i_intz=0
      fx=0
      fy=0
      fz=0
      fpp=0
      IVB_cl=0
      IVB_m=0
      IVB_kw=0
      
      call wkcoup_effect(ie,ig,Ialloy,IB1,IB2,
     &              Nslp_mx,Nslp,STFei26,dt1,
     &              vd_slp,vl_slp,vn_slp,
     &              Rho_gnd, IVB_gnd, pk2i_gnd,
     &              Ftrp, IFtrp, IVB_trp,
     &              fx, fy, fz, fpp, pk2i_intx,  
     &              pk2i_inty, pk2i_intz, pk2i_intp,
     &              IVB_cl,IVB_m,IVB_kw,IVB_bk)
      IVB_wcp=IVB_gnd+IVB_trp

c------------------------------------------------c
c     calculate stress and stiffness: slip        c
c------------------------------------------------c
      if(Imth_add==1)then
         call cal_stress_add(ising)
      else
         call cal_stress_mul(ising)
      endif  

c------------------------------------------------c
c     output to ABAQUS statev                    c
c------------------------------------------------c
      el_st = (matmul(transpose(Fe),Fe))/2
      ! elastic strain Ee-SDV150 to 155
      statev(150) = el_st(1,1)-0.5
      statev(151) = el_st(2,2)-0.5
      statev(152) = el_st(3,3)-0.5
      statev(153) = el_st(1,2)
      statev(154) = el_st(1,3)
      statev(155) = el_st(2,3)
      
      ! plastic strain in crystal coordinates
      intLp=0.5*dt1*(Lp+transpose(Lp))
      statev(156)=statev(156)+intLp(1,1)
      statev(157)=statev(157)+intLp(2,2)
      statev(158)=statev(158)+intLp(3,3)
      statev(159)=statev(159)+intLp(1,2)
      statev(160)=statev(160)+intLp(1,3)
      statev(161)=statev(161)+intLp(2,3) 
      
      ! temperature-dependent material parameters
      statev(162) = c11_gl
      statev(163) = c12_gl
      statev(164) = c44_gl
      statev(165) = tau0_gl
      statev(166) = Adir_gl
      
      if(ising/=0)then
         pnewdt=0.5d0
         return
      else         
         pnewdt=1.1
         stress=cs
         ddsdde=MatJacb
         statev( 1: 4)=eang                   !!-->01-04 eang
         statev( 5:10)=cs                     !!-->05-10 cs
         statev(11:16)=pk2i                   !!-->11-16 pk2i
         statev(17:22)=pk2i_GND               !!-->17-22 internal stress
         do i=1,9                             !!                
            statev(22+i)=Fe(ib1(i),ib2(i))    !!-->23-31 Fe
            statev(31+i)=Fp(ib1(i),ib2(i))    !!-->32-40 Fp
         enddo                                !!
         icomp=40
         do is=1,Nslp                         !!
            icomp=icomp+1                     !! 
            statev(icomp)=IVB(is)             !!-->41-41+NSLP  IVB
         enddo                                !!
         do is=1,48
            icomp=icomp+1
            statev(icomp)=IVB_bk(is)         !! --> 41+NSLP-41+NSLP+48 IVB_bk
         enddo
      	  
c---------------------------------------------------------------------
         if(Iwkcoup_grad/=0)then
            if(Icall_ieig_grad(noel,npt)==0)then
               Icall_ieig_grad(noel,npt)     = 1
               fem_Inf        (noel,npt,1)   = 1
               fem_Fp0        (noel,npt,:,:) = Fp0
               fem_xyz        (noel,npt,:)   = coords*C_unit
            endif
            mx33_1=fem_Fp0(noel,npt,:,:)
            fem_dFp(noel,npt,:,:)=matmul(Fp,transpose(mx33_1))
         endif
c---------------------------------------------------------------------
         if(Iwkcoup_trip/=0)then
            fem_cs(noel,npt,:)=cs
            fem_Fe(noel,npt,:,:)=Fe
            fem_gm (noel,npt,1:NStrp)=
     &      fem_gm0(noel,npt,1:NStrp)+detGM(1:NStrp)
         endif
      endif         

      return
      end
