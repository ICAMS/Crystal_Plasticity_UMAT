c ICAMS CP-UMAT 2026R1
c (c) 2026 by ICAMS, Ruhr University Bochum
c================================================================
c
c    Modules: mod_material
c    Subroutines: init_material, sub_flow_harden (called from mod_stress)
c    Support module flags, used only locally 
c
c================================================================
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c +
c +   Define methods for handling of and integer selection flag (ISF) in form of nibbles (4 bits) 
c +   in one single integer variable (PROPS(10))
c +
c +   Nibble 0: isotropic hardening
c +     0 = none
c +     1 = standard empirical isotropic hardening law with saturation
c +   Nibble 1: kinematic hardening (multiplier: 16)
c +     0 = none
c +     1 = Frederick-Armstrong
c +     2 = Chaboche
c +     3 = Ohno-Wang
c +   Nibble 2: gradient plasticity  (multiplier: 256)
c +     0 = none
c +     1 = gradient plasticity
c +   Nibble 3: internal stresses for superalloys  (multiplier: 4 096)
c +     0 = none
c +     1 = include internal stresses for superalloys
c +   Nibble 4: superalloy-specific mechanisms  (multiplier: 65 536)
c +     0 = none
c +     1 = include superalloy-specific mechanisms
c +   Nibble 5: transformation-induced plasticity (TRIP) (multiplier: 1 048 576)
c +     0 = none
c +     1 = include TRIP effect 
c +   Nibble 6: temperature dependent parameters (multiplier: 16 777 216)
c +     0 = none
c +     1 = include linear temperature dependence of flow behavior and elastic constants
c +   Nibble 7: spare
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      module flags
      use, intrinsic :: iso_fortran_env, only: int32
      implicit none

      ! 8 Nibbles (0–7), Nibble multipier: 16
      integer, parameter :: POS_ISO        = 0
      integer, parameter :: POS_KIN        = 1
      integer, parameter :: POS_GRADIENT   = 2
      integer, parameter :: POS_INT        = 3
      integer, parameter :: POS_SUPERALLOY = 4
      integer, parameter :: POS_TRIP       = 5
      integer, parameter :: POS_THERM      = 6
      integer, parameter :: POS_SPARE      = 7

      contains

      ! --- read nibble ---
      integer function get_nibble(flags, pos)
            implicit none
            integer(int32), intent(in) :: flags
            integer(int32), intent(in) :: pos
            get_nibble = iand(ishft(flags, -4*pos), 15_int32)
      end function get_nibble

      ! --- set nibble ---
      function set_nibble(flags, pos, val) result(out)
            implicit none
            integer(int32), intent(in) :: flags
            integer(int32), intent(in) :: pos, val
            integer(int32) :: out
            integer(int32) :: mask
            integer :: shift

            shift = 4*pos
            mask  = ishft(15_int32, shift)

            out = ior(iand(flags, not(mask)), ishft(iand(val,15), shift))
      end function set_nibble

      end module flags

c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c +                                                                           +
c +   Define material parameters: Constitutive and physical parameters, flags +
c +   Methods for intialization of arrays and evolution of internal variables +
c +                                                                           +
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        MODULE mod_material
        use globalvalue, only : Nslp_mx
        IMPLICIT NONE
c     Indices for parameter extraction from PROPS array
        INTEGER, PARAMETER :: POS_SPCGRP = 9
        INTEGER, PARAMETER :: POS_ISF = 10
        INTEGER, PARAMETER :: POS_C11 = 11
        INTEGER, PARAMETER :: POS_C12 = 12
        INTEGER, PARAMETER :: POS_C44 = 13
        INTEGER, PARAMETER :: POS_NSLP = 14
        INTEGER, PARAMETER :: POS_SHRT0 = 15
        INTEGER, PARAMETER :: POS_PWFL = 16
        INTEGER, PARAMETER :: POS_CRSS0 = 17
        INTEGER, PARAMETER :: LEN_STD = 17
        INTEGER, PARAMETER :: LEN_WH_ISO = 5
        INTEGER, PARAMETER :: LEN_WH_KIN = 8
        INTEGER, PARAMETER :: LEN_SUPER = 12
        INTEGER, PARAMETER :: LEN_GRAD = 8
        INTEGER, parameter :: LEN_THERM = 6
        INTEGER, PARAMETER :: MAX_MATERIAL_CACHE = 16
        INTEGER, PARAMETER :: MAX_PROPS_CACHE = 256

C     Physical parameters and constants              ! UNIT         Explanation
C     -------------------------------------------------------------------------------                                                                                                                  
        REAL(8), PARAMETER :: R_gas     = 8.314d0      ! J/K/mol      Gas constant                        
        REAL(8), PARAMETER :: NAvg      = 6.022141d23  ! 1/mol        Avogardo Number
        REAL(8), PARAMETER :: temp_min  = 294.0d0      ! K            standard temp if not def. in model
        REAL(8)            :: temp_cur                 ! K            current temperature, reference temperature
c     Parameters to activate material-related submodules, values read from PROPS in umat subroutine
c     Block flags: Mandatory flags
        INTEGER :: Isf           !-->superflag with nibbles for activation of submodules (see umat_flags module)
        INTEGER :: Ihard_iso     !-->with(1) isotropic hardening, without(0)
        INTEGER :: Iwkcoup_bk    !-->with(1) FAKH, with(2) CHKH, with(3) OWKH, without(0)
        INTEGER :: Iwkcoup_trip  !-->with(1), without(0): Phase transformation-induced plasticity (TRIP)
        INTEGER :: Iwkcoup_grad  !-->with(1), without(0): Gradient plasticity
        INTEGER :: Iwkcoup_int   !-->with(1), without(0): Superalloy-specific internal stresses 
        INTEGER :: Iwkcoup_sup   !-->with(1), without(0): Superalloy-specific mechanisms
        INTEGER :: Iwkcoup_temp  !-->with(1), without(0): temperature-dependent parameters for austenite
c     Block core: Mandatory constitutive parameters
        INTEGER :: ialloy        !legacy: ICAMS material identifier for alloy, 1 ... 6, other values will be ignored
        INTEGER :: spc_grp       !space group number for slip system generation
        integer :: N_slip        !number of slip systems for this alloy
        real(8) :: c11           !elastic oefficient C_11 
        real(8) :: c12           !elastic oefficient C_12
        real(8) :: c44           !elastic oefficient C_44
        real(8) :: shrt0         !reference shear rate for viscoplastic slip
        real(8) :: pwfl          !stress exponent for power-law viscoplastic slip (= inverse strain rate sensitivity)
        real(8) :: crss0         !initial critical resolved shear stress (CRSS) for slip
c     Block iso: Parameters for isotropic hardening (if Ihard_iso=1)
        real(8) :: crsss = 100.0  ! saturation value of the slip resistance CRSS due to isotropic hardening
        real(8) :: hdrt0 = 0.0    ! reference hardening rate
        real(8) :: c_cpl = 1.0    ! self-hardening parameter for isotropic hardening
        real(8) :: c_oth = 1.4    ! cross-(latent)-hardening parameter for isotropic hardening
        real(8) :: pwhd = 1.0     ! exponent for power-law isotropic hardening
c     Block kin: Parameters for kinematic hardening with Frederick-Armstrong, Chaboche and Ohno_wang models (if Iwkcoup_bk>0)
        integer :: Ns_kin = 48    ! number of effective slip systems considered for kinematic hardening
        real(8) :: Adir = 1.0     ! parameter for kinematic hardening (always first term in FA, OW, and Chab)
        real(8) :: Adyn = 1.0     ! parameter for kinematic hardening (always first term in FA, OW, and Chab)
        real(8) :: M_OW = 0.0     ! Ohno-Wang hardening exponent, only active if Iwkcoup_bk=3
        real(8) :: A2 = 1.0       ! parameter for kinematic hardening, second term in Chaboche, only active for Iwkcoup_bk=2
        real(8) :: B2 = 1.0       ! parameter for kinematic hardening, second term in Chaboche, only active for Iwkcoup_bk=2
        real(8) :: A3 = 1.0       ! parameter for kinematic hardening, third term in Chaboche, only active for Iwkcoup_bk=2
        real(8) :: B3 = 1.0       ! parameter for kinematic hardening, third term in Chaboche, only active for Iwkcoup_bk=2
c     Block superalloy: Parameters for superalloy-specific mechanisms (if Iwkcoup_sup=1)
        real(8) :: shrt0m = 0.0   ! reference shear rate for channel glide in superalloy model
        real(8) :: shrt0p = 0.0   ! reference shear rate for glide in precipitates in superalloy model
        real(8) :: shrt0c = 0.0   ! reference shear rate for interfacial climb in superalloy model
        real(8) :: Qactm = 0.0    ! activation energy for channel glide in superalloy model
        real(8) :: Qactp = 0.0    ! activation energy for glide in precipitates in superalloy model
        real(8) :: Qactc = 0.0    ! activation energy for interfacial climb in superalloy model
        real(8) :: crssm0 = 0.0   ! initial CRSS for channel glide in superalloy model
        real(8) :: crssp0 = 0.0   ! initial CRSS for glide in precipitates in superalloy model
        real(8) :: crssc0 = 0.0   ! initial CRSS for interfacial climb in superalloy model
        real(8) :: crss_oro = 0.0 ! Orowan stress for channel glide in superalloy model
        real(8) :: c_cplc = 1.0   ! self-hardening parameter for channel glide in superalloy model
        real(8) :: c_othc = 1.4   ! cross-(latent)-hardening parameter for channel glide in superalloy model
c     Block grad: Parameters for gradient plasticity (if Iwkcoup_grad=1)
        integer :: Ioct_search = 1 !-->with(1) use octree for neighbor search (only regular meshes!), without(0) use brute force search
        real(8) :: C_taui = 0.01   ! Taylor hardening parameter for gradient plasticity
        real(8) :: L_size = 1.D-9  ! internal length scale parameter for gradient plasticity, multiplied with element size for actual length scale
        integer :: I_gnd_iso = 1   !-->with(1) use isotropic GND-based hardening, without(0) use no hardening for GNDs
        integer :: I_gnd_kin = 1   !-->with(1) use kinematic GND-based hardening, without(0) use no hardening for GNDs
        real(8) :: C_unit = 1.0D-6 ! mesh size parameter for gradient plasticity (unit: m)
        real(8) :: B_lattice = 2.5D-10 ! Burgers vector norm (unit: m)
        real(8) :: CD_smooth = 3.0D-15 ! smoothing parameter for gradient plasticity
c     Block thermal: Parameters for temperature dependence (if Iwkcoup_temp=1)
        real(8) :: Qact = 0.0       ! activation energy for thermally activated glide
        real(8) :: c11_ref = 250.D3 !reference value for c11 at temperature zero, used for linear temperature dependence
        real(8) :: C12_ref = 150.D3 !reference value for c12 at temperature zero, used for linear temperature dependence
        real(8) :: c44_ref = 120.D3 !reference value for c44 at temperature zero, used for linear temperature dependence
        real(8) :: c11_slope = 0.0  !slope for linear temperature dependence of c11
        real(8) :: c12_slope = 0.0  !slope for linear temperature dependence of c12
        real(8) :: c44_slope = 0.0  !slope for linear temperature dependence of c44
        real(8) :: Adir_ref = 1.0   !reference value for Adir at temperature zero, used for linear temperature dependence
        real(8) :: Adir_slope = 0.0 !slope for linear temperature dependence of Adir
        real(8) :: crss0_ref = 100.0 !reference value for crss0 at temperature zero, used for linear temperature dependence
        real(8) :: crss0_slope = 0.0 !slope for linear temperature dependence of crss0
c     Block arrays for slip system geometry and hardening matrices
        real(8) Mstiff(6,6)
        real(8) HMij(Nslp_mx,Nslp_mx)
        real(8) Dvct(Nslp_mx,3)
        real(8) Lvct(Nslp_mx,3)
        real(8) Nvct(Nslp_mx,3)
        real(8) Msmd(Nslp_mx,3,3)
        real(8) SMsmd(Nslp_mx,3,3)
        real(8) AMsmd(Nslp_mx,3,3)
        real(8) V1smd(Nslp_mx,6)
        real(8) V2smd(Nslp_mx,6)

        type mat_param_set
           integer :: Isf = 0
           integer :: ialloy = 0
           integer :: Ihard_iso = 0
           integer :: Iwkcoup_bk = 0
           integer :: Iwkcoup_trip = 0
           integer :: Iwkcoup_grad = 0
           integer :: Iwkcoup_int = 0
           integer :: Iwkcoup_sup = 0
           integer :: Iwkcoup_temp = 0
           integer :: spc_grp = 0
           integer :: N_slip = 0
           real(8) :: c11 = 0.d0
           real(8) :: c12 = 0.d0
           real(8) :: c44 = 0.d0
           real(8) :: shrt0 = 0.d0
           real(8) :: pwfl = 0.d0
           real(8) :: crss0 = 0.d0
           real(8) :: crsss = 100.d0
           real(8) :: hdrt0 = 0.d0
           real(8) :: c_cpl = 1.d0
           real(8) :: c_oth = 1.4d0
           real(8) :: pwhd = 1.d0
           integer :: Ns_kin = 48
           real(8) :: Adir = 1.d0
           real(8) :: Adyn = 1.d0
           real(8) :: M_OW = 0.d0
           real(8) :: A2 = 1.d0
           real(8) :: B2 = 1.d0
           real(8) :: A3 = 1.d0
           real(8) :: B3 = 1.d0
           real(8) :: shrt0m = 0.d0
           real(8) :: shrt0p = 0.d0
           real(8) :: shrt0c = 0.d0
           real(8) :: Qactm = 0.d0
           real(8) :: Qactp = 0.d0
           real(8) :: Qactc = 0.d0
           real(8) :: crssm0 = 0.d0
           real(8) :: crssp0 = 0.d0
           real(8) :: crssc0 = 0.d0
           real(8) :: crss_oro = 0.d0
           real(8) :: c_cplc = 0.d0
           real(8) :: c_othc = 0.d0
           integer :: Ioct_search = 1
           real(8) :: C_taui = 0.01d0
           real(8) :: L_size = 1.d-9
           integer :: I_gnd_iso = 1
           integer :: I_gnd_kin = 1
           real(8) :: C_unit = 1.d-6
           real(8) :: B_lattice = 2.5d-10
           real(8) :: CD_smooth = 3.d-15
           real(8) :: Qact = 0.d0
           real(8) :: c11_ref = 0.d0
           real(8) :: c12_ref = 0.d0
           real(8) :: c44_ref = 0.d0
           real(8) :: c11_slope = 0.d0
           real(8) :: c12_slope = 0.d0
           real(8) :: c44_slope = 0.d0
           real(8) :: Adir_ref = 1.d0
           real(8) :: Adir_slope = 1.d0
           real(8) :: crss0_ref = 100.d0
           real(8) :: crss0_slope = 0.d0
           real(8) :: Mstiff(6,6) = 0.d0
           real(8) :: HMij(Nslp_mx,Nslp_mx) = 0.d0
           real(8) :: Dvct(Nslp_mx,3) = 0.d0
           real(8) :: Lvct(Nslp_mx,3) = 0.d0
           real(8) :: Nvct(Nslp_mx,3) = 0.d0
           real(8) :: Msmd(Nslp_mx,3,3) = 0.d0
           real(8) :: SMsmd(Nslp_mx,3,3) = 0.d0
           real(8) :: AMsmd(Nslp_mx,3,3) = 0.d0
           real(8) :: V1smd(Nslp_mx,6) = 0.d0
           real(8) :: V2smd(Nslp_mx,6) = 0.d0
           real(8) :: IVB_ini(Nslp_mx) = 0.d0
           real(8) :: refv_IVB = 0.d0
           real(8) :: refv_pk2i = 0.d0
        end type mat_param_set

        type material_cache_entry
           logical :: valid = .false.
           integer :: nprops = 0
           real(8) :: props(MAX_PROPS_CACHE) = 0.d0
           type(mat_param_set) :: mp
        end type material_cache_entry

c     Shared reporting registry, protected by material_report_io.
        type material_report_entry
           real(8), allocatable :: props(:)
           type(material_report_entry), pointer :: next => null()
        end type material_report_entry
        type(material_report_entry), pointer :: reported => null()

        type(material_cache_entry) :: material_cache(MAX_MATERIAL_CACHE)
        integer :: material_cache_next = 1
        logical :: material_cache_inited = .false.

!$omp threadprivate(temp_cur)
!$omp threadprivate(Isf,Ihard_iso,Iwkcoup_bk,Iwkcoup_trip)
!$omp threadprivate(Iwkcoup_grad,Iwkcoup_int,Iwkcoup_sup,Iwkcoup_temp)
!$omp threadprivate(ialloy,spc_grp,N_slip,c11,c12,c44,shrt0,pwfl,crss0)
!$omp threadprivate(crsss,hdrt0,c_cpl,c_oth,pwhd)
!$omp threadprivate(Ns_kin,Adir,Adyn,M_OW,A2,B2,A3,B3)
!$omp threadprivate(shrt0m,shrt0p,shrt0c,Qactm,Qactp,Qactc)
!$omp threadprivate(crssm0,crssp0,crssc0,crss_oro,c_cplc,c_othc)
!$omp threadprivate(Ioct_search,C_taui,L_size,I_gnd_iso,I_gnd_kin)
!$omp threadprivate(C_unit,B_lattice,CD_smooth,Qact)
!$omp threadprivate(c11_ref,C12_ref,c44_ref,c11_slope,c12_slope,c44_slope)
!$omp threadprivate(Adir_ref,Adir_slope,crss0_ref,crss0_slope)
!$omp threadprivate(Mstiff,HMij,Dvct,Lvct,Nvct,Msmd,SMsmd,AMsmd,V1smd,V2smd)
!$omp threadprivate(material_cache,material_cache_next)
!$omp threadprivate(material_cache_inited)

      contains
c +--------------------------------------------------------------------------------+
c +   extract a local, per-call parameter snapshot from PROPS                      +
c +--------------------------------------------------------------------------------+
      recursive subroutine extract_material_params(NPROPS, PROPS, mp)
            use flags
            implicit none
            integer, intent(in) :: NPROPS
            real(8), intent(in) :: PROPS(NPROPS)
            type(mat_param_set), intent(out) :: mp
            integer is, isf_loc

c           Four-entry legacy input: identifier and three Euler angles.
c           The identifier does not select an alloy; default is copper.
            if (NPROPS == 4) then
               mp%ialloy = 2
               mp%spc_grp = 225
               mp%Isf = 1
               mp%Ihard_iso = 1
               mp%N_slip = 12
               mp%c11 = 247.d3
               mp%c12 = 147.d3
               mp%c44 = 125.d3
               mp%shrt0 = 1.d-3
               mp%pwfl = 20.d0
               mp%crss0 = 20.d0
               mp%crsss = 117.d0
               mp%hdrt0 = 180.d0
               mp%c_cpl = 1.d0
               mp%c_oth = 1.4d0
               mp%pwhd = 2.25d0
               return
            endif
            if (NPROPS < LEN_STD) then
               write(6,*) 'ERROR: PROPS requires exactly 4 legacy',
     &                    ' entries or at least ', LEN_STD
               stop 1
            endif

            mp%ialloy = int(props(1))
            mp%spc_grp = int(props(POS_SPCGRP))
            isf_loc = int(nint(props(POS_ISF)), int32)
            mp%Isf = isf_loc
            mp%Ihard_iso = get_nibble(isf_loc, POS_ISO)
            mp%Iwkcoup_bk = get_nibble(isf_loc, POS_KIN)
            mp%Iwkcoup_trip = get_nibble(isf_loc, POS_TRIP)
            mp%Iwkcoup_grad = get_nibble(isf_loc, POS_GRADIENT)
            mp%Iwkcoup_int = get_nibble(isf_loc, POS_INT)
            mp%Iwkcoup_sup = get_nibble(isf_loc, POS_SUPERALLOY)
            mp%Iwkcoup_temp = get_nibble(isf_loc, POS_THERM)

c           Validate optional block lengths before reading their values.
            is = LEN_STD
            if (mp%Ihard_iso == 1) is = is + LEN_WH_ISO
            if (mp%Iwkcoup_bk > 0) is = is + LEN_WH_KIN
            if (mp%Iwkcoup_sup == 1) is = is + LEN_SUPER
            if (mp%Iwkcoup_grad == 1) is = is + LEN_GRAD
            if (mp%Iwkcoup_temp == 1) is = is + LEN_THERM
            if (NPROPS < is) then
               write(6,*) 'ERROR: Incomplete PROPS blocks; required:', is,
     &                    ' provided:', NPROPS
               stop 1
            endif

            mp%c11 = props(POS_C11)
            mp%c12 = props(POS_C12)
            mp%c44 = props(POS_C44)
            mp%N_slip = int(props(POS_NSLP))
            mp%shrt0 = props(POS_SHRT0)
            mp%pwfl = props(POS_PWFL)
            mp%crss0 = props(POS_CRSS0)
            is = LEN_STD

            if(mp%Ihard_iso==1)then
               mp%crsss = props(is+1)
               mp%hdrt0 = props(is+2)
               mp%c_cpl = props(is+3)
               mp%c_oth = props(is+4)
               mp%pwhd = props(is+5)
               is = is + LEN_WH_ISO
            endif

            if(mp%Iwkcoup_bk>0)then
               mp%Ns_kin = int(props(is+1))
               mp%Adir = props(is+2)
               mp%Adyn = props(is+3)
               mp%M_OW = props(is+4)
               mp%A2 = props(is+5)
               mp%B2 = props(is+6)
               mp%A3 = props(is+7)
               mp%B3 = props(is+8)
               is = is + LEN_WH_KIN
            endif

            if(mp%Iwkcoup_sup==1)then
               mp%shrt0m = props(is+1)
               mp%shrt0p = props(is+2)
               mp%shrt0c = props(is+3)
               mp%Qactm = props(is+4)
               mp%Qactp = props(is+5)
               mp%Qactc = props(is+6)
               mp%crssm0 = props(is+7)
               mp%crssp0 = props(is+8)
               mp%crssc0 = props(is+9)
               mp%crss_oro = props(is+10)
               mp%c_cplc = props(is+11)
               mp%c_othc = props(is+12)
               is = is + LEN_SUPER
            endif

            if(mp%Iwkcoup_grad==1)then
               mp%Ioct_search = int(props(is+1))
               mp%C_taui = props(is+2)
               mp%L_size = props(is+3)
               mp%I_gnd_iso = int(props(is+4))
               mp%I_gnd_kin = int(props(is+5))
               mp%C_unit = props(is+6)
               mp%B_lattice = props(is+7)
               mp%CD_smooth = props(is+8)
               is = is + LEN_GRAD
            endif

            if(mp%Iwkcoup_temp==1)then
               mp%Qact = props(is+1)
               mp%c11_slope = props(is+2)
               mp%c12_slope = props(is+3)
               mp%c44_slope = props(is+4)
               mp%Adir_slope = props(is+5)
               mp%crss0_slope = props(is+6)
               mp%c11_ref = mp%c11
               mp%c12_ref = mp%c12
               mp%c44_ref = mp%c44
               mp%crss0_ref = mp%crss0
               mp%Adir_ref = mp%Adir
            endif

         return
      endsubroutine extract_material_params

c +--------------------------------------------------------------------------------+
c +   copy legacy module material data into an explicit material snapshot           +
c +--------------------------------------------------------------------------------+
      subroutine snapshot_material_params(mp)
            use mod_gaussp, only : IVB_ini, refv_IVB, refv_pk2i
            implicit none
            type(mat_param_set), intent(inout) :: mp

            mp%Isf = Isf
            mp%ialloy = ialloy
            mp%Ihard_iso = Ihard_iso
            mp%Iwkcoup_bk = Iwkcoup_bk
            mp%Iwkcoup_trip = Iwkcoup_trip
            mp%Iwkcoup_grad = Iwkcoup_grad
            mp%Iwkcoup_int = Iwkcoup_int
            mp%Iwkcoup_sup = Iwkcoup_sup
            mp%Iwkcoup_temp = Iwkcoup_temp
            mp%spc_grp = spc_grp
            mp%N_slip = N_slip
            mp%c11 = c11
            mp%c12 = c12
            mp%c44 = c44
            mp%shrt0 = shrt0
            mp%pwfl = pwfl
            mp%crss0 = crss0
            mp%crsss = crsss
            mp%hdrt0 = hdrt0
            mp%c_cpl = c_cpl
            mp%c_oth = c_oth
            mp%pwhd = pwhd
            mp%Ns_kin = Ns_kin
            mp%Adir = Adir
            mp%Adyn = Adyn
            mp%M_OW = M_OW
            mp%A2 = A2
            mp%B2 = B2
            mp%A3 = A3
            mp%B3 = B3
            mp%shrt0m = shrt0m
            mp%shrt0p = shrt0p
            mp%shrt0c = shrt0c
            mp%Qactm = Qactm
            mp%Qactp = Qactp
            mp%Qactc = Qactc
            mp%crssm0 = crssm0
            mp%crssp0 = crssp0
            mp%crssc0 = crssc0
            mp%crss_oro = crss_oro
            mp%c_cplc = c_cplc
            mp%c_othc = c_othc
            mp%Ioct_search = Ioct_search
            mp%C_taui = C_taui
            mp%L_size = L_size
            mp%I_gnd_iso = I_gnd_iso
            mp%I_gnd_kin = I_gnd_kin
            mp%C_unit = C_unit
            mp%B_lattice = B_lattice
            mp%CD_smooth = CD_smooth
            mp%Qact = Qact
            mp%c11_ref = c11_ref
            mp%c12_ref = c12_ref
            mp%c44_ref = c44_ref
            mp%c11_slope = c11_slope
            mp%c12_slope = c12_slope
            mp%c44_slope = c44_slope
            mp%Adir_ref = Adir_ref
            mp%Adir_slope = Adir_slope
            mp%crss0_ref = crss0_ref
            mp%crss0_slope = crss0_slope
            mp%Mstiff = Mstiff
            mp%HMij = HMij
            mp%Dvct = Dvct
            mp%Lvct = Lvct
            mp%Nvct = Nvct
            mp%Msmd = Msmd
            mp%SMsmd = SMsmd
            mp%AMsmd = AMsmd
            mp%V1smd = V1smd
            mp%V2smd = V2smd
            mp%IVB_ini = IVB_ini
            mp%refv_IVB = refv_IVB
            mp%refv_pk2i = refv_pk2i

            return
      endsubroutine snapshot_material_params

c +--------------------------------------------------------------------------------+
c +   restore legacy module data from an explicit material snapshot                 +
c +--------------------------------------------------------------------------------+
      subroutine activate_material_params(mp)
            use mod_gaussp, only : STFei26,smdMi,smdSMi,smdAMi,
     &                             smdVi1,smdVi2,vd_slp,vl_slp,
     &                             vn_slp,IVB_ini,refv_IVB,refv_pk2i
            implicit none
            type(mat_param_set), intent(in) :: mp

            Isf = mp%Isf
            ialloy = mp%ialloy
            Ihard_iso = mp%Ihard_iso
            Iwkcoup_bk = mp%Iwkcoup_bk
            Iwkcoup_trip = mp%Iwkcoup_trip
            Iwkcoup_grad = mp%Iwkcoup_grad
            Iwkcoup_int = mp%Iwkcoup_int
            Iwkcoup_sup = mp%Iwkcoup_sup
            Iwkcoup_temp = mp%Iwkcoup_temp
            spc_grp = mp%spc_grp
            N_slip = mp%N_slip
            c11 = mp%c11
            c12 = mp%c12
            c44 = mp%c44
            shrt0 = mp%shrt0
            pwfl = mp%pwfl
            crss0 = mp%crss0
            crsss = mp%crsss
            hdrt0 = mp%hdrt0
            c_cpl = mp%c_cpl
            c_oth = mp%c_oth
            pwhd = mp%pwhd
            Ns_kin = mp%Ns_kin
            Adir = mp%Adir
            Adyn = mp%Adyn
            M_OW = mp%M_OW
            A2 = mp%A2
            B2 = mp%B2
            A3 = mp%A3
            B3 = mp%B3
            shrt0m = mp%shrt0m
            shrt0p = mp%shrt0p
            shrt0c = mp%shrt0c
            Qactm = mp%Qactm
            Qactp = mp%Qactp
            Qactc = mp%Qactc
            crssm0 = mp%crssm0
            crssp0 = mp%crssp0
            crssc0 = mp%crssc0
            crss_oro = mp%crss_oro
            c_cplc = mp%c_cplc
            c_othc = mp%c_othc
            Ioct_search = mp%Ioct_search
            C_taui = mp%C_taui
            L_size = mp%L_size
            I_gnd_iso = mp%I_gnd_iso
            I_gnd_kin = mp%I_gnd_kin
            C_unit = mp%C_unit
            B_lattice = mp%B_lattice
            CD_smooth = mp%CD_smooth
            Qact = mp%Qact
            c11_ref = mp%c11_ref
            c12_ref = mp%c12_ref
            c44_ref = mp%c44_ref
            c11_slope = mp%c11_slope
            c12_slope = mp%c12_slope
            c44_slope = mp%c44_slope
            Adir_ref = mp%Adir_ref
            Adir_slope = mp%Adir_slope
            crss0_ref = mp%crss0_ref
            crss0_slope = mp%crss0_slope
            Mstiff = mp%Mstiff
            HMij = mp%HMij
            Dvct = mp%Dvct
            Lvct = mp%Lvct
            Nvct = mp%Nvct
            Msmd = mp%Msmd
            SMsmd = mp%SMsmd
            AMsmd = mp%AMsmd
            V1smd = mp%V1smd
            V2smd = mp%V2smd
            STFei26 = mp%Mstiff
            smdMi = mp%Msmd
            smdSMi = mp%SMsmd
            smdAMi = mp%AMsmd
            smdVi1 = mp%V1smd
            smdVi2 = mp%V2smd
            vd_slp = mp%Dvct
            vl_slp = mp%Lvct
            vn_slp = mp%Nvct
            IVB_ini = mp%IVB_ini
            refv_IVB = mp%refv_IVB
            refv_pk2i = mp%refv_pk2i

            return
      endsubroutine activate_material_params

c +--------------------------------------------------------------------------------+
c +   check whether cached PROPS match the current Abaqus material constants        +
c +--------------------------------------------------------------------------------+
      logical function material_props_match(entry,NPROPS,PROPS)
            implicit none
            type(material_cache_entry), intent(in) :: entry
            integer, intent(in) :: NPROPS
            real(8), intent(in) :: PROPS(NPROPS)

            material_props_match = .false.
            if (.not. entry%valid) return
            if (entry%nprops /= NPROPS) return
            if (NPROPS > MAX_PROPS_CACHE) return
            material_props_match =
     &         all(entry%props(1:NPROPS) == PROPS(1:NPROPS))

            return
      endfunction material_props_match

c +--------------------------------------------------------------------------------+
c +   cache only the standard, local material response covered by mat_param_set     +
c +--------------------------------------------------------------------------------+
      logical function material_cache_supported(mp)
            implicit none
            type(mat_param_set), intent(in) :: mp

            material_cache_supported =
     &         mp%Iwkcoup_grad == 0 .and.
     &         mp%Iwkcoup_trip == 0 .and.
     &         mp%Iwkcoup_int == 0 .and.
     &         mp%Iwkcoup_sup == 0 .and.
     &         mp%Iwkcoup_temp == 0

            return
      endfunction material_cache_supported

c +--------------------------------------------------------------------------------+
c +   initialize the per-thread material cache on first use                         +
c +--------------------------------------------------------------------------------+
      subroutine init_material_cache()
            implicit none
            integer i

            if (material_cache_inited) return

            do i=1,MAX_MATERIAL_CACHE
               material_cache(i)%valid = .false.
               material_cache(i)%nprops = 0
               material_cache(i)%props = 0.d0
            enddo
            material_cache_next = 1
            material_cache_inited = .true.

            return
      endsubroutine init_material_cache

c +--------------------------------------------------------------------------------+
c +   return crystal-frame data; callers explicitly activate legacy module storage +
c +--------------------------------------------------------------------------------+
      recursive subroutine init_material_cached(NPROPS, PROPS, mp)
            implicit none
            integer, intent(in) :: NPROPS
            real(8), intent(in) :: PROPS(NPROPS)
            type(mat_param_set), intent(out) :: mp
            integer i,slot

            call init_material_cache()
            call extract_material_params(NPROPS, PROPS, mp)

            if (material_cache_supported(mp) .and.
     &          NPROPS <= MAX_PROPS_CACHE) then
               do i=1,MAX_MATERIAL_CACHE
                  if (material_props_match(material_cache(i),
     &                                    NPROPS,PROPS)) then
                     mp = material_cache(i)%mp
                     return
                  endif
               enddo
            endif

            call build_material_params(NPROPS, PROPS, mp)

            if (material_cache_supported(mp) .and.
     &          NPROPS <= MAX_PROPS_CACHE) then
               slot = material_cache_next
               material_cache(slot)%valid = .true.
               material_cache(slot)%nprops = NPROPS
               material_cache(slot)%props = 0.d0
               material_cache(slot)%props(1:NPROPS) = PROPS(1:NPROPS)
               material_cache(slot)%mp = mp
               material_cache_next = material_cache_next + 1
               if (material_cache_next > MAX_MATERIAL_CACHE) then
                  material_cache_next = 1
               endif
            endif

            return
      endsubroutine init_material_cached

c +--------------------------------------------------------------------------------+
c +   build constitutive parameters and derived arrays without module side effects  +
c +--------------------------------------------------------------------------------+
      recursive subroutine build_material_params(NPROPS, PROPS, mp)
            use mod_gaussp, only : ib1,ib2
            implicit none
            integer, intent(in) :: NPROPS
            real(8), intent(in) :: PROPS(NPROPS)
            type(mat_param_set), intent(out) :: mp
            integer i,j,is,js, ii,ia,ib
            real(8) x1,x2,x3

            call extract_material_params(NPROPS, PROPS, mp)

c           Initialize stiffness matrix
            mp%Mstiff(:,:)=0.0
            mp%Mstiff(1,1)=mp%c11
            mp%Mstiff(2,2)=mp%c11
            mp%Mstiff(3,3)=mp%c11
            mp%Mstiff(4,4)=mp%c44*2
            mp%Mstiff(5,5)=mp%c44*2
            mp%Mstiff(6,6)=mp%c44*2
            mp%Mstiff(2,3)=mp%c12
            mp%Mstiff(3,2)=mp%c12
            mp%Mstiff(1,3)=mp%c12
            mp%Mstiff(3,1)=mp%c12
            mp%Mstiff(1,2)=mp%c12
            mp%Mstiff(2,1)=mp%c12
c
            if (mp%spc_grp==229) then  ! bcc lattice
               mp%Dvct( 1,:)=[ -1,  1,  1] ; mp%Nvct( 1,:)=[ 0,  1, -1]
               mp%Dvct( 2,:)=[ -1,  1,  1] ; mp%Nvct( 2,:)=[ 1,  0,  1]
               mp%Dvct( 3,:)=[ -1,  1,  1] ; mp%Nvct( 3,:)=[ 1,  1,  0]
               mp%Dvct( 4,:)=[  1,  1,  1] ; mp%Nvct( 4,:)=[ 0,  1, -1]
               mp%Dvct( 5,:)=[  1,  1,  1] ; mp%Nvct( 5,:)=[ 1,  0, -1]
               mp%Dvct( 6,:)=[  1,  1,  1] ; mp%Nvct( 6,:)=[ 1, -1,  0]
               mp%Dvct( 7,:)=[  1,  1, -1] ; mp%Nvct( 7,:)=[ 0,  1,  1]
               mp%Dvct( 8,:)=[  1,  1, -1] ; mp%Nvct( 8,:)=[ 1,  0,  1]
               mp%Dvct( 9,:)=[  1,  1, -1] ; mp%Nvct( 9,:)=[ 1, -1,  0]
               mp%Dvct(10,:)=[  1, -1,  1] ; mp%Nvct(10,:)=[ 0,  1,  1]
               mp%Dvct(11,:)=[  1, -1,  1] ; mp%Nvct(11,:)=[ 1,  0, -1]
               mp%Dvct(12,:)=[  1, -1,  1] ; mp%Nvct(12,:)=[ 1,  1,  0]
            elseif (mp%spc_grp==225) then  ! fcc lattice
               mp%Dvct( 1,:)=[ 0,  1, -1] ; mp%Nvct( 1,:)=[ -1,  1,  1]
               mp%Dvct( 2,:)=[ 1,  0,  1] ; mp%Nvct( 2,:)=[ -1,  1,  1]
               mp%Dvct( 3,:)=[ 1,  1,  0] ; mp%Nvct( 3,:)=[ -1,  1,  1]
               mp%Dvct( 4,:)=[ 0,  1, -1] ; mp%Nvct( 4,:)=[  1,  1,  1]
               mp%Dvct( 5,:)=[ 1,  0, -1] ; mp%Nvct( 5,:)=[  1,  1,  1]
               mp%Dvct( 6,:)=[ 1, -1,  0] ; mp%Nvct( 6,:)=[  1,  1,  1]
               mp%Dvct( 7,:)=[ 0,  1,  1] ; mp%Nvct( 7,:)=[  1,  1, -1]
               mp%Dvct( 8,:)=[ 1,  0,  1] ; mp%Nvct( 8,:)=[  1,  1, -1]
               mp%Dvct( 9,:)=[ 1, -1,  0] ; mp%Nvct( 9,:)=[  1,  1, -1]
               mp%Dvct(10,:)=[ 0,  1,  1] ; mp%Nvct(10,:)=[  1, -1,  1]
               mp%Dvct(11,:)=[ 1,  0, -1] ; mp%Nvct(11,:)=[  1, -1,  1]
               mp%Dvct(12,:)=[ 1,  1,  0] ; mp%Nvct(12,:)=[  1, -1,  1]
            else
               print*, 'space group not implemented for slip system generation: ', mp%spc_grp
               stop
            endif

            if (mp%Iwkcoup_sup==1) then
c              additional definitions for superalloy
               mp%Dvct(13:24,:)=mp%Dvct( 1:12,:)
               mp%Nvct(13:24,:)=mp%Nvct( 1:12,:)
               mp%Dvct(25:36,:)=mp%Dvct( 1:12,:)
               mp%Nvct(25:36,:)=mp%Nvct( 1:12,:)

               mp%Dvct(37,:)=[ 2,  1,  1]
               mp%Dvct(38,:)=[-1, -2,  1]
               mp%Dvct(39,:)=[-1,  1, -2]
               mp%Dvct(40,:)=[-2,  1,  1]
               mp%Dvct(41,:)=[ 1, -2,  1]
               mp%Dvct(42,:)=[ 1,  1, -2]
               mp%Dvct(43,:)=[-2,  1, -1]
               mp%Dvct(44,:)=[ 1, -2, -1]
               mp%Dvct(45,:)=[ 1,  1,  2]
               mp%Dvct(46,:)=[-2, -1,  1]
               mp%Dvct(47,:)=[ 1,  2,  1]
               mp%Dvct(48,:)=[ 1, -1, -2]
c          Dvct(37:48,:)=Dvct( 1:12,:)
               mp%Nvct(37:48,:)=mp%Nvct( 1:12,:)

               mp%Dvct(49,:)=[-1,  1,  0]; mp%Nvct(49,:)=[0, 0, 1]
               mp%Dvct(50,:)=[ 1,  1,  0]; mp%Nvct(50,:)=[0, 0, 1]
               mp%Dvct(51,:)=[-1,  0,  1]; mp%Nvct(51,:)=[0, 1, 0]
               mp%Dvct(52,:)=[ 1,  0,  1]; mp%Nvct(52,:)=[0, 1, 0]
               mp%Dvct(53,:)=[ 0, -1,  1]; mp%Nvct(53,:)=[1, 0, 0]
               mp%Dvct(54,:)=[ 0,  1,  1]; mp%Nvct(54,:)=[1, 0, 0]
               mp%Dvct(55:60,:)=mp%Dvct(49:54,:)
               mp%Nvct(55:60,:)=mp%Nvct(49:54,:)
            endif
c
            do is=1,mp%N_slip
               mp%Lvct(is,1)=mp%Nvct(is,2)*mp%Dvct(is,3)
     &                      -mp%Nvct(is,3)*mp%Dvct(is,2)
               mp%Lvct(is,2)=mp%Nvct(is,3)*mp%Dvct(is,1)
     &                      -mp%Nvct(is,1)*mp%Dvct(is,3)
               mp%Lvct(is,3)=mp%Nvct(is,1)*mp%Dvct(is,2)
     &                      -mp%Nvct(is,2)*mp%Dvct(is,1)
               x1=dsqrt(sum(mp%Nvct(is,:)**2))
               x2=dsqrt(sum(mp%Dvct(is,:)**2))
               x3=dsqrt(sum(mp%Lvct(is,:)**2))
               mp%Nvct(is,:)=mp%Nvct(is,:)/x1
               mp%Dvct(is,:)=mp%Dvct(is,:)/x2
               mp%Lvct(is,:)=mp%Lvct(is,:)/x3
               do i=1,3
               do j=1,3
                  mp%Msmd(is,i,j)=mp%Dvct(is,i)*mp%Nvct(is,j)
               enddo
               enddo
               mp%SMsmd(is,:,:)=( mp%Msmd(is,:,:)
     &                          +transpose(mp%Msmd(is,:,:)) )/2
               mp%AMsmd(is,:,:)=( mp%Msmd(is,:,:)
     &                          -transpose(mp%Msmd(is,:,:)) )/2
               do i=1,6
                  mp%V1smd(is,i)=
     &            mp%SMsmd(is,ib1(i),ib2(i))
               enddo
               mp%V2smd(is,1:3)=mp%V1smd(is,1:3)
               mp%V2smd(is,4:6)=mp%V1smd(is,4:6)*2
            enddo
c
            if (mp%Iwkcoup_sup==0) then
               do is=1,mp%N_slip
               do js=1,mp%N_slip
                  x1=sum(dabs(mp%Nvct(is,:)-mp%Nvct(js,:)))
                  if(x1<1.d-10)then
                     mp%HMij(is,js)=mp%c_cpl
                  else
                     mp%HMij(is,js)=mp%c_oth
                  endif
               enddo
               enddo
            else  ! superalloy
               do ii=1,4
                  ia=1+12*(ii-1)
                  ib=12+12*(ii-1)
               do is=ia,ib
               do js=ia,ib
                  x1=sum(dabs(mp%Nvct(is,:)-mp%Nvct(js,:)))
                  if(x1<1.d-10)then
                     mp%HMij(is,js)=mp%c_cpl
                  else
                     mp%HMij(is,js)=mp%c_oth
                  endif
               enddo
               enddo
               enddo

               do ii=1,2
                  ia=49+6*(ii-1)
                  ib=54+6*(ii-1)
               do is=ia,ib
               do js=ia,ib
                  x1=sum(dabs(mp%Nvct(is,:)-mp%Nvct(js,:)))
                  if(x1<1.d-10)then
                     mp%HMij(is,js)=mp%c_cplc
                  else
                     mp%HMij(is,js)=mp%c_othc
                  endif
               enddo
               enddo
               enddo
            endif
            mp%refv_pk2i = mp%c44*1.d-6
            mp%refv_IVB = mp%c44*1.d-6
            if (mp%Iwkcoup_sup==0) then
               mp%IVB_ini(1:mp%N_slip) = mp%crss0
            else
               mp%IVB_ini(1:36) = mp%crssm0
               mp%IVB_ini(37:48) = mp%crssp0
               mp%IVB_ini(49:60) = mp%crssc0
            endif
c
            call report_material_params(NPROPS, PROPS, mp)
            return
      endsubroutine build_material_params

c +--------------------------------------------------------------------------------+
c +   get constitutive parameters from PROPS and initialize legacy module arrays    +
c +--------------------------------------------------------------------------------+
c     Write each distinct input definition once per process to Abaqus .dat.
c     This registry is shared even though constitutive caches are per-thread.
      subroutine report_material_params(NPROPS, PROPS, mp)
            integer, intent(in) :: NPROPS
            real(8), intent(in) :: PROPS(NPROPS)
            type(mat_param_set), intent(in) :: mp
            type(material_report_entry), pointer :: entry
            logical :: found
            integer :: i
!$omp critical(material_report_io)
            found = .false.
            entry => reported
            do while (associated(entry))
               if (size(entry%props) == NPROPS) then
                  if (all(entry%props == PROPS)) then
                     found = .true.
                     exit
                  endif
               endif
               entry => entry%next
            enddo
            if (.not. found) then
               allocate(entry)
               allocate(entry%props(NPROPS))
               entry%props = PROPS
               entry%next => reported
               reported => entry
               write(6,'(/,A)') 'UMAT material parameter reference'
               if (NPROPS == 4) then
                  write(6,'(A)')
     &              'Legacy fallback: 4 PROPS entries; adapt parameters in subroutine extract_material_params.'
                  write(6,'(A)')
     &              'PROPS(1) ignored for selection; PROPS(2:4) are Euler angles.'
               else
                  write(6,'(A)') 'Material parameters supplied in PROPS.'
               endif
               do i=1,NPROPS
                  write(6,'(A,I0,A,ES24.16)')
     &                  'PROPS(',i,') = ',PROPS(i)
               enddo
               write(6,'(A)')
     &           'Resolved parameters (inactive blocks show defaults):'
               write(6,'(A,I0)') 'Isf = ',mp%Isf
               write(6,'(A,I0)') 'ialloy = ',mp%ialloy
               write(6,'(A,I0)') 'Ihard_iso = ',mp%Ihard_iso
               write(6,'(A,I0)') 'Iwkcoup_bk = ',mp%Iwkcoup_bk
               write(6,'(A,I0)') 'Iwkcoup_trip = ',mp%Iwkcoup_trip
               write(6,'(A,I0)') 'Iwkcoup_grad = ',mp%Iwkcoup_grad
               write(6,'(A,I0)') 'Iwkcoup_int = ',mp%Iwkcoup_int
               write(6,'(A,I0)') 'Iwkcoup_sup = ',mp%Iwkcoup_sup
               write(6,'(A,I0)') 'Iwkcoup_temp = ',mp%Iwkcoup_temp
               write(6,'(A,I0)') 'spc_grp = ',mp%spc_grp
               write(6,'(A,I0)') 'N_slip = ',mp%N_slip
               write(6,'(A,ES24.16)') 'c11 = ',mp%c11
               write(6,'(A,ES24.16)') 'c12 = ',mp%c12
               write(6,'(A,ES24.16)') 'c44 = ',mp%c44
               write(6,'(A,ES24.16)') 'shrt0 = ',mp%shrt0
               write(6,'(A,ES24.16)') 'pwfl = ',mp%pwfl
               write(6,'(A,ES24.16)') 'crss0 = ',mp%crss0
               write(6,'(A,ES24.16)') 'crsss = ',mp%crsss
               write(6,'(A,ES24.16)') 'hdrt0 = ',mp%hdrt0
               write(6,'(A,ES24.16)') 'c_cpl = ',mp%c_cpl
               write(6,'(A,ES24.16)') 'c_oth = ',mp%c_oth
               write(6,'(A,ES24.16)') 'pwhd = ',mp%pwhd
               write(6,'(A,I0)') 'Ns_kin = ',mp%Ns_kin
               write(6,'(A,ES24.16)') 'Adir = ',mp%Adir
               write(6,'(A,ES24.16)') 'Adyn = ',mp%Adyn
               write(6,'(A,ES24.16)') 'M_OW = ',mp%M_OW
               write(6,'(A,ES24.16)') 'A2 = ',mp%A2
               write(6,'(A,ES24.16)') 'B2 = ',mp%B2
               write(6,'(A,ES24.16)') 'A3 = ',mp%A3
               write(6,'(A,ES24.16)') 'B3 = ',mp%B3
               write(6,'(A,ES24.16)') 'shrt0m = ',mp%shrt0m
               write(6,'(A,ES24.16)') 'shrt0p = ',mp%shrt0p
               write(6,'(A,ES24.16)') 'shrt0c = ',mp%shrt0c
               write(6,'(A,ES24.16)') 'Qactm = ',mp%Qactm
               write(6,'(A,ES24.16)') 'Qactp = ',mp%Qactp
               write(6,'(A,ES24.16)') 'Qactc = ',mp%Qactc
               write(6,'(A,ES24.16)') 'crssm0 = ',mp%crssm0
               write(6,'(A,ES24.16)') 'crssp0 = ',mp%crssp0
               write(6,'(A,ES24.16)') 'crssc0 = ',mp%crssc0
               write(6,'(A,ES24.16)') 'crss_oro = ',mp%crss_oro
               write(6,'(A,ES24.16)') 'c_cplc = ',mp%c_cplc
               write(6,'(A,ES24.16)') 'c_othc = ',mp%c_othc
               write(6,'(A,I0)') 'Ioct_search = ',mp%Ioct_search
               write(6,'(A,ES24.16)') 'C_taui = ',mp%C_taui
               write(6,'(A,ES24.16)') 'L_size = ',mp%L_size
               write(6,'(A,I0)') 'I_gnd_iso = ',mp%I_gnd_iso
               write(6,'(A,I0)') 'I_gnd_kin = ',mp%I_gnd_kin
               write(6,'(A,ES24.16)') 'C_unit = ',mp%C_unit
               write(6,'(A,ES24.16)') 'B_lattice = ',mp%B_lattice
               write(6,'(A,ES24.16)') 'CD_smooth = ',mp%CD_smooth
               write(6,'(A,ES24.16)') 'Qact = ',mp%Qact
               write(6,'(A,ES24.16)') 'c11_ref = ',mp%c11_ref
               write(6,'(A,ES24.16)') 'c12_ref = ',mp%c12_ref
               write(6,'(A,ES24.16)') 'c44_ref = ',mp%c44_ref
               write(6,'(A,ES24.16)') 'c11_slope = ',mp%c11_slope
               write(6,'(A,ES24.16)') 'c12_slope = ',mp%c12_slope
               write(6,'(A,ES24.16)') 'c44_slope = ',mp%c44_slope
               write(6,'(A,ES24.16)') 'Adir_ref = ',mp%Adir_ref
               write(6,'(A,ES24.16)') 'Adir_slope = ',mp%Adir_slope
               write(6,'(A,ES24.16)') 'crss0_ref = ',mp%crss0_ref
               write(6,'(A,ES24.16)') 'crss0_slope = ',mp%crss0_slope
               write(6,'(A,ES24.16)') 'R_gas = ',R_gas
               write(6,'(A,ES24.16)') 'NAvg = ',NAvg
               write(6,'(A,ES24.16)') 'temp_min = ',temp_min
               write(6,'(A,/)') 'End UMAT material parameter reference'
            endif
!$omp end critical(material_report_io)
      end subroutine report_material_params

      subroutine init_material(NPROPS, PROPS, mp)
            implicit none
            integer, intent(in) :: NPROPS
            real(8), intent(in) :: PROPS(NPROPS)
            type(mat_param_set), intent(inout), optional :: mp
            type(mat_param_set) :: mp_local

            call build_material_params(NPROPS, PROPS, mp_local)
            call activate_material_params(mp_local)
            if (present(mp)) mp = mp_local

            return
      endsubroutine init_material

c        c----------------------------c
c        c   flow and hardening lows  c
c        c----------------------------c
c        Standard flow and hardening laws in version 2026R1
c        does not support grad, superalloys, temp.dependence, 
c     internal stresses, and TRIP effect
         recursive subroutine sub_flow_harden_std(
     &              Iexp_loc,N_slip_loc,crss0_loc,crsss_loc,
     &              hdrt0_loc,pwhd_loc,shrt0_loc,pwfl_loc,
     &              HMij_loc,IVB_wcp,IVB_bk,tau,
     &              IVB,dgmdt,ddgmdt_dtau,ddgmdt_dIVB,
     &              dIVBdt,ddIVBdt_ddgmdt,ddIVBdt_dIVB,ising)
            use globalvalue, only : Nslp_mx
            implicit none
            integer, intent(in) :: Iexp_loc, N_slip_loc
            real(8), intent(in) :: crss0_loc, crsss_loc
            real(8), intent(in) :: hdrt0_loc, pwhd_loc, shrt0_loc, pwfl_loc
            real(8), intent(in) :: HMij_loc(Nslp_mx,Nslp_mx)
            integer, intent(out) :: ising
            integer is,js
            real(8), intent(in) :: tau(Nslp_mx)
            real(8), intent(in) :: IVB(Nslp_mx)
            real(8), intent(in) :: IVB_wcp(Nslp_mx)
            real(8), intent(in) :: IVB_bk(Nslp_mx)
            real(8) IVB_eff(Nslp_mx)
            real(8), intent(out) :: dgmdt(Nslp_mx)
            real(8), intent(out) :: ddgmdt_dtau(Nslp_mx)
            real(8), intent(out) :: ddgmdt_dIVB(Nslp_mx)
            real(8), intent(out) :: dIVBdt(Nslp_mx)
            real(8), intent(out) :: ddIVBdt_ddgmdt(Nslp_mx,Nslp_mx)
            real(8), intent(out) :: ddIVBdt_dIVB(Nslp_mx,Nslp_mx)
            real(8) x1,x2

            ising=0
            dgmdt=0
            ddgmdt_dtau=0
            ddgmdt_dIVB=0
            dIVBdt=0
            ddIVBdt_ddgmdt=0
            ddIVBdt_dIVB=0

            do is=1,N_slip_loc
               x1=crss0_loc*1.d-10
               x2=crsss_loc
               if(IVB(is)<x1 .or. IVB(is)>x2)then
                  ising=112
                  return
               endif
            enddo

            do is=1,N_slip_loc
               IVB_eff(is)=IVB(is)+IVB_wcp(is)
               dgmdt(is)=shrt0_loc*(dabs(tau(is)-IVB_bk(is))/IVB_eff(is))
     &                        **pwfl_loc*dsign(1.d0,(tau(is)-IVB_bk(is)))
               if(Iexp_loc/=1)then
                  ddgmdt_dtau(is)=pwfl_loc/IVB_eff(is)*shrt0_loc
     &              *(dabs(tau(is)-IVB_bk(is))/IVB_eff(is))**(pwfl_loc-1)
                  ddgmdt_dIVB(is)=-pwfl_loc*dgmdt(is)/IVB_eff(is)
               endif
            enddo

            do is=1,N_slip_loc
            do js=1,N_slip_loc
               x1=1-IVB(js)/crsss_loc
               dIVBdt(is)=dIVBdt(is) + HMij_loc(is,js)
     &         *hdrt0_loc*dabs(dgmdt(js))*x1**pwhd_loc
               if(Iexp_loc/=1)then
                  ddIVBdt_ddgmdt(is,js)=HMij_loc(is,js)*hdrt0_loc
     &            *dsign(1.d0,tau(js))*x1**pwhd_loc
                  ddIVBdt_dIVB(is,js)=HMij_loc(is,js)*hdrt0_loc
     &            *ddgmdt_dIVB(js)*dsign(1.d0,tau(js))*x1**pwhd_loc
     &            -HMij_loc(is,js)*hdrt0_loc*dabs(dgmdt(js))
     &            *x1**(pwhd_loc-1)*pwhd_loc/crsss_loc
               endif
            enddo
            enddo

            return
         endsubroutine

c        Legacy subroutine for flow and hardening laws
c        retains support for grad, superalloys, temp.dependence, 
c        internal stresses, and TRIP effect
         subroutine sub_flow_harden(
     &              Iexp_loc,IVB_wcp,IVB_bk,IVB_cl,IVB_m,IVB_kw,tau,
     &              IVB,dgmdt,ddgmdt_dtau,ddgmdt_dIVB,
     &              dIVBdt,ddIVBdt_ddgmdt,ddIVBdt_dIVB,ising)
            use globalvalue
            use mod_gaussp, only : STFei26
            implicit none
            integer is,js,ising,ii,ia,ib
            integer Iexp_loc
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
            real(8) x1,x2,x3
c
            ising=0
            dgmdt=0
            ddgmdt_dtau=0
            ddgmdt_dIVB=0
            dIVBdt=0
            ddIVBdt_ddgmdt=0
            ddIVBdt_dIVB=0
            if (Iwkcoup_temp==1) then
               x3 = dexp(-Qact/R_gas/temp_cur)
               if(temp_cur>300.d0)then
                  c11   = c11_ref - c11_slope*temp_cur
                  c12   = c12_ref - c12_slope*temp_cur
                  c44   = c44_ref - c44_slope*temp_cur
                  crss0 = crss0_ref - crss0_slope*temp_cur
                  Adir  = Adir_ref - Adir_slope*temp_cur
c                 Update stiffness matrix
                  Mstiff(:,:)=0.0
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
c                 Update arrays in mod_gaussp
                  STFei26 = Mstiff
               endif
            else
               x3 = 1.0
            endif
            if (Iwkcoup_sup==0) then

c--------------resolved shear stress and resistence
               do is=1,N_slip
                  x1=crss0*1.d-10
                  x2=crsss
                  if(IVB(is)<x1 .or. IVB(is)>x2)then
                     ising=112
                     return
                  endif
               enddo

c--------------shear rate, derivative of shear rate w.r.t. pk2i,IVB
               do is=1,N_slip
                  IVB_eff(is)=IVB(is)+IVB_wcp(is)
                  dgmdt(is)=shrt0*x3*(dabs(tau(is)-IVB_bk(is))/IVB_eff(is))
     &                          **pwfl*dsign(1.d0,(tau(is)-IVB_bk(is)))
                  if(Iexp_loc/=1)then
                     ddgmdt_dtau(is)=pwfl/IVB_eff(is)*shrt0*x3
     &                 *(dabs(tau(is)-IVB_bk(is))/IVB_eff(is))**(pwfl-1)
                     ddgmdt_dIVB(is)=-pwfl*dgmdt(is)/IVB_eff(is)
                  endif
               enddo
c--------------evolution rate, derivative of evolution rate w.r.t. pk2i,IVB
               do is=1,N_slip
               do js=1,N_slip
                  x1=1-IVB(js)/crsss
                  dIVBdt(is)=dIVBdt(is) + HMij(is,js)
     &            *hdrt0*dabs(dgmdt(js))*x1**pwhd 
                  if(Iexp_loc/=1)then
                     ddIVBdt_ddgmdt(is,js)=HMij(is,js)*hdrt0
     &               *dsign(1.d0,tau(js))*x1**pwhd 
                     ddIVBdt_dIVB(is,js)=HMij(is,js)*hdrt0
     &               *ddgmdt_dIVB(js)*dsign(1.d0,tau(js))*x1**pwhd 
     &               -HMij(is,js)*hdrt0*dabs(dgmdt(js))
     &               *x1**(pwhd-1)*pwhd/crsss
                  endif
               enddo
               enddo
            else  ! superalloy
c--------------resolved shear stress and resistence
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

c--------------shear rate, derivative of shear rate w.r.t. pk2i,IVB
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
c                 KW deactivated !!!
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

c--------------evolution rate, derivative of evolution rate w.r.t. pk2i,IVB
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
            endif

            return
         endsubroutine

c================================================================
c
c    subroutine kinematic hardening/back stress evolution: Back stress state stored in STATEV
c
c================================================================

c================================================================
      recursive subroutine sub_bk_evolution(dtime,iwkcoup_bk_loc,n_slip_loc,
     &                            adir_loc,adyn_loc,m_ow_loc,
     &                            a2_loc,b2_loc,a3_loc,b3_loc,
     &                            IVB_bk,IVB_bk_ch,dgmdt_loc)
         implicit none
         real(8), intent(in) :: dtime
         integer, intent(in) :: iwkcoup_bk_loc
         integer, intent(in) :: n_slip_loc
         real(8), intent(in) :: adir_loc, adyn_loc, m_ow_loc
         real(8), intent(in) :: a2_loc, b2_loc, a3_loc, b3_loc
         real(8), intent(in) :: dgmdt_loc(Nslp_mx)
         real(8), intent(inout) :: IVB_bk(Nslp_mx)
         real(8), intent(inout) :: IVB_bk_ch(Nslp_mx,3)
         integer :: is
         real(8) :: dIVB_bk(Nslp_mx)
         real(8) :: dIVB_bk_ch(Nslp_mx,3)

         dIVB_bk = 0.d0
         dIVB_bk_ch = 0.d0

         if (iwkcoup_bk_loc == 1) then
            do is=1,n_slip_loc
               dIVB_bk(is) = adir_loc*dgmdt_loc(is) -
     &           adyn_loc*IVB_bk(is)*dabs(dgmdt_loc(is))
               IVB_bk(is) = IVB_bk(is) + dIVB_bk(is)*dtime
            enddo
         else if (iwkcoup_bk_loc == 2) then
            do is=1,n_slip_loc
               dIVB_bk_ch(is,1) = adir_loc*dgmdt_loc(is) -
     &           adyn_loc*IVB_bk_ch(is,1)*dabs(dgmdt_loc(is))
               dIVB_bk_ch(is,2) = a2_loc*dgmdt_loc(is) -
     &           b2_loc*IVB_bk_ch(is,2)*dabs(dgmdt_loc(is))
               dIVB_bk_ch(is,3) = a3_loc*dgmdt_loc(is) -
     &           b3_loc*IVB_bk_ch(is,3)*dabs(dgmdt_loc(is))
               IVB_bk_ch(is,1:3) = IVB_bk_ch(is,1:3) +
     &           dIVB_bk_ch(is,1:3)*dtime
               IVB_bk(is) = IVB_bk_ch(is,1) + IVB_bk_ch(is,2) +
     &           IVB_bk_ch(is,3)
            enddo
         else if (iwkcoup_bk_loc == 3) then
            do is=1,n_slip_loc
               dIVB_bk(is) = adir_loc*dgmdt_loc(is) -
     &           adyn_loc*(dabs(IVB_bk(is))/(adir_loc/adyn_loc))**m_ow_loc
     &           *IVB_bk(is)*dabs(dgmdt_loc(is))
               IVB_bk(is) = IVB_bk(is) + dIVB_bk(is)*dtime
            enddo
         endif

         return
      endsubroutine

c================================================================
      end module mod_material
