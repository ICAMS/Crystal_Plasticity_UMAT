c ICAMS CP-UMAT 2026R1
c (c) 2026 by ICAMS, Ruhr University Bochum
c================================================================
c
c    Modules: globalvalue
c    Subroutines: uexternaldb, umat
c
c================================================================
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c +                                                                         +
c +   Define global values: System parameters                               +
c +                                                                         +
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      MODULE globalvalue
      IMPLICIT NONE

c     Constants for mesh-related arrays (max. mesh size)
      INTEGER, PARAMETER :: Tnel=8000       !-->max. number of elements, change if needed, std: 27000 for 30x30x30 mesh
      INTEGER, PARAMETER :: Tngp=8          !-->max. number of integration points per element, change if needed
      INTEGER, PARAMETER :: Unit_glob=25    !-->unit of global result file, 0: no output 
c     Parameters for FFT mesh grid size used for smoothing and gradient calculation
      INTEGER, PARAMETER :: Tnfx=4 ! standard: 64
      INTEGER, PARAMETER :: Tnfy=4 ! standard: 64
      INTEGER, PARAMETER :: Tnfz=4 ! standard: 64
c     Max. number of slip planes
      integer, parameter :: Nslp_mx = 60     ! maximum number of slip systems for all alloys, change if needed
c     Error margin
      REAL(8), PARAMETER :: EPS = 1.d-8
c     fields for output of homogenized results (average over element and integration points)
      INTEGER :: modNel, modNgp 
      REAL(8), POINTER, DIMENSION(:,:,:) :: Etot => null()
      REAL(8), POINTER, DIMENSION(:,:,:) :: Epl => null()
      REAL(8), POINTER, DIMENSION(:,:,:) :: Sig => null()
      CHARACTER(LEN=255) :: fname_glob

      contains

      subroutine init_global_output_arrays()
      implicit none

      if (.not. associated(Etot)) then
         allocate(Etot(Tnel,Tngp,6))
      endif
      if (.not. associated(Epl)) then
         allocate(Epl(Tnel,Tngp,6))
      endif
      if (.not. associated(Sig)) then
         allocate(Sig(Tnel,Tngp,6))
      endif

      Etot = 0.d0
      Epl  = 0.d0
      Sig  = 0.d0
      modNel = 0
      modNgp = 0

      return
      end subroutine init_global_output_arrays
      
      END MODULE globalvalue

c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c +                                                                         +
c +   Include other subroutines                                             +
c +                                                                         +
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      include "mod_gaussp.f"   
      include "mod_material.f"
      include "mod_octree.f"
      include "other_code.f"        
      include "mod_wkcoup.f"
      include "mod_stress.f"            
      include "mod_smthgd.f"

      module mod_umat_phases
      use globalvalue
      use mod_gaussp
      use mod_material
      use mod_wkcoup
      use mod_stress
      implicit none

      integer, parameter :: isv_ee = 150
      integer, parameter :: isv_ep = 156
      integer, parameter :: isv_temp = 162
      integer, parameter :: isv_spare_beg = 167
      integer, parameter :: n_isv_spare = 10
      integer, parameter :: isv_meta_fixed = isv_spare_beg + n_isv_spare

      contains

      recursive subroutine umat_prepare_phase(nprops,props,nstatv,mp,use_gp_ctx,
     &                              nstate_bk,isv_bk_beg)
      implicit none
      integer, intent(in) :: nprops,nstatv
      real(8), intent(in) :: props(nprops)
      type(mat_param_set), intent(out) :: mp
      logical, intent(out) :: use_gp_ctx
      integer, intent(out) :: nstate_bk,isv_bk_beg
      integer :: nstate_req

      call init_material_cached(nprops, props, mp)
      use_gp_ctx = Imth_add/=1 .and. cal_stress_mul_ctx_supported(mp)
c     Only the legacy solver needs a module material working set.
      if (.not. use_gp_ctx) call activate_material_params(mp)

      nstate_bk = Nslp_mx
      if (mp%Iwkcoup_bk == 2) nstate_bk = 3*Nslp_mx
      isv_bk_beg = isv_meta_fixed
      nstate_req = isv_bk_beg + nstate_bk - 1
      if (nstatv < nstate_req) then
         print*, 'Error: nstatv too small. Need at least ', nstate_req
         stop
      endif

      return
      endsubroutine umat_prepare_phase

      recursive subroutine umat_load_phase(noel,npt,dtime,temp,time,props,
     &                           dfgrd0,dfgrd1,statev,mp,use_gp_ctx,
     &                           nstate_bk,isv_bk_beg,gp,Qm,IVB_bk_ch)
      implicit none
      integer, intent(in) :: noel,npt,nstate_bk,isv_bk_beg
      real(8), intent(in) :: dtime,temp,time(2),props(*)
      real(8), intent(in) :: dfgrd0(3,3),dfgrd1(3,3)
      real(8), intent(inout) :: statev(*)
      type(mat_param_set), intent(in) :: mp
      logical, intent(in) :: use_gp_ctx
      type(gp_context), intent(out) :: gp
      real(8), intent(out) :: Qm(3,3)
      real(8), intent(out) :: IVB_bk_ch(Nslp_mx,3)
      integer i,j,is

c     Build the point context directly from this call, never via module scratch.
      if (noel.gt.Tnel) then
        print*, 'too many elements, change Tnel'
        stop
      endif
      gp%dt1=dtime
      gp%Fg0=dfgrd0
      gp%Fg=dfgrd1
      gp%eang00(1:3)=props(2:4)

      IVB_bk_ch=0.d0

      if(time(2) < EPS)then
         call icams_Eang2Q(props(2),props(3),props(4),Qm)
         statev( 1: 3)=props(2:4)
         statev( 5:10)=0
         statev(11:16)=0
         statev(17:22)=0
         do i=1,9
            statev(22+i)=Qm(ib1(i),ib2(i))
            statev(31+i)=Qm(ib2(i),ib1(i))
         enddo
         statev(40+0*mp%N_slip+1 : 40+1*mp%N_slip)=
     &      mp%IVB_ini(1:mp%N_slip)
         statev(isv_ee:isv_ee+5)=0.d0
         statev(isv_ep:isv_ep+5)=0.d0
         statev(isv_spare_beg)=0.d0
         statev(isv_spare_beg:isv_spare_beg+n_isv_spare-1)=0.d0
         statev(isv_bk_beg:isv_bk_beg+nstate_bk-1)=0.d0
      endif

      gp%eang0=statev( 1: 4)
      gp%cs0  =statev( 5:10)
      gp%pk2i =statev(11:16)
      do i=1,9
         j=i; if(i>6) j=i-3
         gp%csM0(ib1(i),ib2(i))=gp%cs0(j)
      enddo
      do i=1,9
         gp%Fe0(ib1(i),ib2(i))=statev(22+i)
         gp%Fe (ib1(i),ib2(i))=statev(22+i)
         gp%Fp0(ib1(i),ib2(i))=statev(31+i)
         gp%Fp (ib1(i),ib2(i))=statev(31+i)
         Qm (ib1(i),ib2(i))=statev(22+i)
      enddo
      do is=1,mp%N_slip
         gp%IVB0(is)=statev(40+0*mp%N_slip+is)
         gp%IVB (is)=statev(40+0*mp%N_slip+is)
      enddo
      if (mp%Iwkcoup_bk == 2) then
         do is=1,Nslp_mx
            IVB_bk_ch(is,1)=statev(isv_bk_beg+3*(is-1)+0)
            IVB_bk_ch(is,2)=statev(isv_bk_beg+3*(is-1)+1)
            IVB_bk_ch(is,3)=statev(isv_bk_beg+3*(is-1)+2)
            gp%IVB_bk(is)=IVB_bk_ch(is,1)+IVB_bk_ch(is,2)
     &                 +IVB_bk_ch(is,3)
         enddo
      else
         do is=1,Nslp_mx
            gp%IVB_bk(is)=statev(isv_bk_beg+is-1)
         enddo
      endif

      gp%Qm=Qm

c     Explicit compatibility boundary: legacy paths retain thread-private storage.
      if (.not. use_gp_ctx) then
         ie=noel
         ig=npt
         temp_cur=max(temp,temp_min)
         call store_gp_context_to_globals(gp)
      endif

      return
      endsubroutine umat_load_phase

      recursive subroutine umat_solve_phase(dtime,mp,use_gp_ctx,gp,ising,
     &                            IVB_bk_ch)
      implicit none
      real(8), intent(in) :: dtime
      type(mat_param_set), intent(in) :: mp
      logical, intent(in) :: use_gp_ctx
      type(gp_context), intent(inout) :: gp
      integer, intent(out) :: ising
      real(8), intent(inout) :: IVB_bk_ch(Nslp_mx,3)

      if(use_gp_ctx)then
         call cal_stress_mul(gp,mp,ising)
      else
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

         if(mp%Iwkcoup_grad/=0 .or. mp%Iwkcoup_trip/=0 .or.
     &   mp%Iwkcoup_int/=0  .or. mp%Iwkcoup_sup/=0)then
            call wkcoup_effect(ie,ig,mp%Ialloy,IB1,IB2,
     &              Nslp_mx,mp%N_slip,mp%Mstiff,dt1,
     &              mp%Dvct,mp%Lvct,mp%Nvct,
     &              Rho_gnd, IVB_gnd, pk2i_gnd,
     &              Ftrp, IFtrp, IVB_trp,
     &              fx, fy, fz, fpp, pk2i_intx,
     &              pk2i_inty, pk2i_intz, pk2i_intp,
     &              IVB_cl,IVB_m,IVB_kw)
         endif
         IVB_wcp=IVB_gnd+IVB_trp

         if(Imth_add==1)then
            call cal_stress_add(ising)
         else
            call cal_stress_mul_legacy(ising,mp)
         endif

      endif

      if(mp%Iwkcoup_bk/=0)then
         if(use_gp_ctx)then
            call sub_bk_evolution(dtime,mp%Iwkcoup_bk,mp%N_slip,
     &                           mp%Adir,mp%Adyn,mp%M_OW,
     &                           mp%A2,mp%B2,mp%A3,mp%B3,
     &                           gp%IVB_bk,IVB_bk_ch,gp%dgmdt)
         else
            call sub_bk_evolution(dtime,mp%Iwkcoup_bk,mp%N_slip,
     &                           mp%Adir,mp%Adyn,mp%M_OW,
     &                           mp%A2,mp%B2,mp%A3,mp%B3,
     &                           IVB_bk,IVB_bk_ch,dgmdt)
         endif
      endif

      return
      endsubroutine umat_solve_phase

      recursive subroutine umat_store_phase(mp,use_gp_ctx,gp,ising,statev,stress,
     &                            ddsdde,pnewdt,ntens,Qm_loc,
     &                            IVB_bk_ch,ok)
      implicit none
      type(mat_param_set), intent(in) :: mp
      logical, intent(in) :: use_gp_ctx
      type(gp_context), intent(in) :: gp
      integer, intent(in) :: ising,ntens
      real(8), intent(inout) :: statev(*)
      real(8), intent(out) :: stress(ntens)
      real(8), intent(out) :: ddsdde(ntens,ntens)
      real(8), intent(out) :: pnewdt
      real(8), intent(in) :: Qm_loc(3,3)
      real(8), intent(in) :: IVB_bk_ch(Nslp_mx,3)
      logical, intent(out) :: ok
      integer i,is,icomp
      real(8) intLp(3,3),euler_st(3,3),devLp(3,3),trLp,peeq_inc

      if(use_gp_ctx)then
         euler_st = 0.5*(matmul(transpose(gp%Fe),gp%Fe))
      else
         euler_st = 0.5*(matmul(transpose(Fe),Fe))
      endif
      statev(isv_ee+0) = euler_st(1,1)-0.5
      statev(isv_ee+1) = euler_st(2,2)-0.5
      statev(isv_ee+2) = euler_st(3,3)-0.5
      statev(isv_ee+3) = euler_st(1,2)
      statev(isv_ee+4) = euler_st(1,3)
      statev(isv_ee+5) = euler_st(2,3)

      if(use_gp_ctx)then
         intLp=0.5*gp%dt1*(gp%Lp+transpose(gp%Lp))
         intLp = matmul(gp%Qm, matmul(intLp, transpose(gp%Qm)))
      else
         intLp=0.5*dt1*(Lp+transpose(Lp))
         intLp = matmul(Qm_loc, matmul(intLp, transpose(Qm_loc)))
      endif
      statev(isv_ep+0)=statev(isv_ep+0)+intLp(1,1)
      statev(isv_ep+1)=statev(isv_ep+1)+intLp(2,2)
      statev(isv_ep+2)=statev(isv_ep+2)+intLp(3,3)
      statev(isv_ep+3)=statev(isv_ep+3)+intLp(1,2)
      statev(isv_ep+4)=statev(isv_ep+4)+intLp(1,3)
      statev(isv_ep+5)=statev(isv_ep+5)+intLp(2,3)

      trLp=(intLp(1,1)+intLp(2,2)+intLp(3,3))/3.0d0
      devLp=intLp
      devLp(1,1)=devLp(1,1)-trLp
      devLp(2,2)=devLp(2,2)-trLp
      devLp(3,3)=devLp(3,3)-trLp
      peeq_inc=sqrt(2.0d0/3.0d0*
     &              (devLp(1,1)**2+devLp(2,2)**2+devLp(3,3)**2
     &              +2.0d0*(devLp(1,2)**2+devLp(1,3)**2
     &              +devLp(2,3)**2)))
      statev(isv_spare_beg)=statev(isv_spare_beg)+peeq_inc

      statev(isv_temp+0) = mp%c11
      statev(isv_temp+1) = mp%c12
      statev(isv_temp+2) = mp%c44
      statev(isv_temp+3) = mp%crss0
      statev(isv_temp+4) = mp%Adir

      if(ising/=0)then
         pnewdt=0.5d0
         ok=.false.
         return
      endif

      pnewdt=1.1
      ok=.true.
      if(use_gp_ctx)then
         stress=gp%cs
         ddsdde=gp%MatJacb
         statev( 1: 4)=gp%eang
         statev( 5:10)=gp%cs
         statev(11:16)=gp%pk2i
         statev(17:22)=gp%pk2i_gnd
      else
         stress=cs
         ddsdde=MatJacb
         statev( 1: 4)=eang
         statev( 5:10)=cs
         statev(11:16)=pk2i
         statev(17:22)=pk2i_GND
      endif
      do i=1,9
         if(use_gp_ctx)then
            statev(22+i)=gp%Fe(ib1(i),ib2(i))
            statev(31+i)=gp%Fp(ib1(i),ib2(i))
         else
            statev(22+i)=Fe(ib1(i),ib2(i))
            statev(31+i)=Fp(ib1(i),ib2(i))
         endif
      enddo
      icomp=40
      do is=1,mp%N_slip
         icomp=icomp+1
         if(use_gp_ctx)then
            statev(icomp)=gp%IVB(is)
         else
            statev(icomp)=IVB(is)
         endif
      enddo
      if (mp%Iwkcoup_bk == 2) then
         do is=1,Nslp_mx
            statev(isv_meta_fixed+3*(is-1)+0)=IVB_bk_ch(is,1)
            statev(isv_meta_fixed+3*(is-1)+1)=IVB_bk_ch(is,2)
            statev(isv_meta_fixed+3*(is-1)+2)=IVB_bk_ch(is,3)
         enddo
      else
         do is=1,Nslp_mx
            if(use_gp_ctx)then
               statev(isv_meta_fixed+is-1)=gp%IVB_bk(is)
            else
               statev(isv_meta_fixed+is-1)=IVB_bk(is)
            endif
         enddo
      endif

      return
      endsubroutine umat_store_phase

      recursive subroutine umat_store_global_outputs(use_gp_ctx,gp,noel,npt,ntens,
     &                                     statev,stress)
      implicit none
      logical, intent(in) :: use_gp_ctx
      type(gp_context), intent(in) :: gp
      integer, intent(in) :: noel,npt,ntens
      real(8), intent(in) :: statev(*),stress(ntens)
      integer i
      real(8) euler_st(3,3)

      if(use_gp_ctx)then
         euler_st = 0.5*(matmul(transpose(gp%Fg),gp%Fg))
      else
         euler_st = 0.5*(matmul(transpose(Fg),Fg))
      endif
      Etot(NOEL,NPT,1) = euler_st(1,1)-0.5
      Etot(NOEL,NPT,2) = euler_st(2,2)-0.5
      Etot(NOEL,NPT,3) = euler_st(3,3)-0.5
      Etot(NOEL,NPT,4) = euler_st(1,2)
      Etot(NOEL,NPT,5) = euler_st(1,3)
      Etot(NOEL,NPT,6) = euler_st(2,3)
      do i=1,ntens
         Epl(NOEL,NPT,i) = statev(isv_ep+i-1)
         Sig(NOEL,NPT,i) = STRESS(i)
      enddo

      return
      endsubroutine umat_store_global_outputs

      endmodule mod_umat_phases

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
      use mod_material
      implicit none
      integer lop,lrestart,kstep,kinc
      real(8) time(2),dtime
      real(8) sg(6), etg(6), epg(6)
      real(8) :: sg1,sg2,sg3,sg4,sg5,sg6
      real(8) :: ep1,ep2,ep3,ep4,ep5,ep6
      real(8) :: et1,et2,et3,et4,et5,et6
      integer i,j,k,ljname,lxoutdir,io_stat
      character(len=255) xoutdir, job_name

      if(lop==0)then   !ini
         print*, kinc,dtime,time(2), 'First call, initialization...'
         temp_cur = temp_min  ! Can initial temp be queried in uexternaldb???
         call init_global_output_arrays()
         call mod_gspt_ini()
         call mod_wkcp_ini(IB1,IB2)
c        Initialize result file
         if(Unit_glob>0)then
           CALL GETJOBNAME(job_name,ljname)
           CALL GETOUTDIR(xoutdir,lxoutdir)
           fname_glob = xoutdir(1:lxoutdir) // '/' //
     1         job_name(1:ljname) // '_glob.csv'
           print*,'Writing global results to file: ', trim(fname_glob)
           print*,'WARNING: Assuming regular mesh for homogenization!'
           open(unit=Unit_glob, file=trim(fname_glob), status='replace',
     &          action='write', iostat=io_stat)
           if (io_stat /= 0) then
              write(6,*) 'Error opening global result file: ',
     &                   trim(fname_glob), ' iostat=', io_stat
           else
              write(unit=Unit_glob, fmt='("# time, ",3(A42))',
     &              iostat=io_stat)
     &              'Total strain (11, 22, 33, 12, 13, 23), ',
     &              'Plastic strain (11, 22, 33, 12, 13, 23), ',
     &              'Stress (11, 22, 33, 12, 13, 23)'
              if (io_stat /= 0) write(6,*)
     &              'Error writing global header, iostat=', io_stat
              write(unit=Unit_glob, fmt='(19(G0.8,:,","))',
     &              iostat=io_stat)
     &              0., 0., 0., 0., 0., 0., 0.,
     &              0., 0., 0., 0., 0., 0., 0.,
     &              0., 0., 0., 0., 0.
              if (io_stat /= 0) write(6,*)
     &              'Error writing global initial row, iostat=', io_stat
              close(unit=Unit_glob, iostat=io_stat)
              if (io_stat /= 0) write(6,*)
     &              'Error closing global result file, iostat=', io_stat
           endif
         end if
      end if
 
      if(lop==1)then
         print*, kinc,dtime,time(2),temp_cur, 'start of time step'
      endif

      if(lop==2)then
         print*, kinc,dtime,time(2),temp_cur, 'call wkcoup_evo subroutine'

         call wkcoup_evolution(dtime)
         if(Unit_glob>0)then
                sg1=0.d0; sg2=0.d0; sg3=0.d0
                sg4=0.d0; sg5=0.d0; sg6=0.d0
                ep1=0.d0; ep2=0.d0; ep3=0.d0
                ep4=0.d0; ep5=0.d0; ep6=0.d0
                et1=0.d0; et2=0.d0; et3=0.d0
                et4=0.d0; et5=0.d0; et6=0.d0

!$omp parallel do default(shared) private(i,j)
!$omp& reduction(+:sg1,sg2,sg3,sg4,sg5,sg6)
!$omp& reduction(+:ep1,ep2,ep3,ep4,ep5,ep6)
!$omp& reduction(+:et1,et2,et3,et4,et5,et6)
!$omp& collapse(2)
                do i=1,modNel
                   do j=1,modNgp
                      sg1 = sg1 + Sig(i,j,1)
                      sg2 = sg2 + Sig(i,j,2)
                      sg3 = sg3 + Sig(i,j,3)
                      sg4 = sg4 + Sig(i,j,4)
                      sg5 = sg5 + Sig(i,j,5)
                      sg6 = sg6 + Sig(i,j,6)
                      ep1 = ep1 + Epl(i,j,1)
                      ep2 = ep2 + Epl(i,j,2)
                      ep3 = ep3 + Epl(i,j,3)
                      ep4 = ep4 + Epl(i,j,4)
                      ep5 = ep5 + Epl(i,j,5)
                      ep6 = ep6 + Epl(i,j,6)
                      et1 = et1 + Etot(i,j,1)
                      et2 = et2 + Etot(i,j,2)
                      et3 = et3 + Etot(i,j,3)
                      et4 = et4 + Etot(i,j,4)
                      et5 = et5 + Etot(i,j,5)
                      et6 = et6 + Etot(i,j,6)
                   end do
                end do
!$omp end parallel do

                sg(1)=sg1; sg(2)=sg2; sg(3)=sg3
                sg(4)=sg4; sg(5)=sg5; sg(6)=sg6
                epg(1)=ep1; epg(2)=ep2; epg(3)=ep3
                epg(4)=ep4; epg(5)=ep5; epg(6)=ep6
                etg(1)=et1; etg(2)=et2; etg(3)=et3
                etg(4)=et4; etg(5)=et5; etg(6)=et6
           sg  = sg /(modNel*modNgp)
           epg = epg/(modNel*modNgp)
           etg = etg/(modNel*modNgp)
           open(unit=Unit_glob, file=trim(fname_glob), status='old',
     &          position='append', action='write', iostat=io_stat)
           if (io_stat /= 0) then
              write(6,*) 'Error opening global result file: ',
     &                   trim(fname_glob), ' iostat=', io_stat
           else
              write(unit=Unit_glob, fmt='(19(G0.8,:,","))',
     &              iostat=io_stat) time(2), etg, epg, sg
              if (io_stat /= 0) write(6,*)
     &              'Error appending global row, iostat=', io_stat
              close(unit=Unit_glob, iostat=io_stat)
              if (io_stat /= 0) write(6,*)
     &              'Error closing global result file, iostat=', io_stat
           endif
         end if
      endif

      return
      end

c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c +                                                                         +
c +  Abaqus Umat: User Material Subroutine                                  +
c +  Calculate elastic plastic behavior for phenomenological                +
c +  crystal plasticity model with optional submodules                      +
c +  Evaluate stress, plastic strain and consistent tangent stiffness       +
c +  internal variables stored in STATEV array                              +
c +  material parameters read from PROPS array                              +
c +                                                                         +  
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      recursive subroutine umat(stress,statev,ddsdde,
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
C   SDV167     : equivalent plastic strain
C   SDV168-176 : spare fields for future extensions
C   SDV177-... : kinematic hardening back-stress history (size depends on model)
C
C    Material Properties
C ----------------------------------------------------------------------
C    PROPS(1)     : Material identifier -> ialloys (legacy parameter)
C    PROPS(2...4) : Euler angles (in rad) of grain orientation -> eang00
C    PROPS(5...8) : void
C    PROPS(9)     : Integer space group number for slip system generation -> spc_grp
C    PROPS(10)    : Integer selection flag (ISF) with nibbles for activation of material specific options (see umat_flags module)
C    PROPS(11...) : Material parameters according to mapping.csv
C 
C ----------------------------------------------------------------------
C=======================================================================
      use globalvalue
      use mod_material
      use mod_stress
      use mod_wkcoup
      use mod_gaussp
      use mod_umat_phases
      
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
c     Nslp_mx: maximum number of slip systems for all alloys (defined in material module)
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
      integer :: nstate_bk, isv_bk_beg
      type(mat_param_set) :: mp
      type(gp_context) :: gp
      logical :: use_gp_ctx, store_ok
      real(8) mx33_1(3,3)
      real(8) phi1,phi,phi2,Qm(3,3)
      real(8) intLp(3,3),euler_st(3,3),ev(ntens),devLp(3,3),trLp,peeq_inc
      real(8) IVB_bk_ch(Nslp_mx,3)
c----------------------------------------------------------------------------
c     Validate the build before touching any per-thread working data.
      call mod_gspt_ini

!$omp critical (cp_outarr_init)
      if (.not. associated(Etot)) then
         call init_global_output_arrays()
      endif
!$omp end critical (cp_outarr_init)

      IF ((KINC==1).AND.(KSTEP==1)) THEN
!$omp critical (cp_modsize)
        modNel = MAX(modNel, NOEL)
        modNgp = MAX(modNgp, NPT)
!$omp end critical (cp_modsize)
      END IF
c----------------------------------------------------------------------------
      call umat_prepare_phase(nprops,props,nstatv,mp,use_gp_ctx,
     &                        nstate_bk,isv_bk_beg)
      call umat_load_phase(noel,npt,dtime,temp,time,props,dfgrd0,
     &                     dfgrd1,statev,mp,use_gp_ctx,nstate_bk,
     &                     isv_bk_beg,gp,Qm,IVB_bk_ch)
      call umat_solve_phase(dtime,mp,use_gp_ctx,gp,ising,IVB_bk_ch)
      call umat_store_phase(mp,use_gp_ctx,gp,ising,statev,stress,
     &                      ddsdde,pnewdt,ntens,Qm,IVB_bk_ch,store_ok)
      if(.not. store_ok)return

c---------------------------------------------------------------------
         if(mp%Iwkcoup_grad/=0)then
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
         if(mp%Iwkcoup_trip/=0)then
            fem_cs(noel,npt,:)=cs
            fem_Fe(noel,npt,:,:)=Fe
            fem_gm (noel,npt,1:NStrp)=
     &      fem_gm0(noel,npt,1:NStrp)+detGM(1:NStrp)
         endif
      call umat_store_global_outputs(use_gp_ctx,gp,noel,npt,ntens,
     &                               statev,stress)

      return
      end
