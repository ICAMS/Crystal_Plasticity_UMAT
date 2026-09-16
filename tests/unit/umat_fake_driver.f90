program umat_fake_driver
   use globalvalue, only : init_global_output_arrays
   use mod_gaussp, only : mod_gspt_ini, IB1, IB2
   use mod_wkcoup, only : mod_wkcp_ini
   implicit none

   integer, parameter :: ntens = 6
   integer, parameter :: ndi = 3
   integer, parameter :: nshr = 3
   integer, parameter :: nstatv = 240
   integer, parameter :: nprops = 22
   integer, parameter :: max_points = 8000

   integer :: nsteps, npoints, material_id, load_case, p, task, orientation_mode, order_mode, step, unit
   real(8) :: total_strain, dtime, temp
   real(8) :: props(nprops)
   real(8), allocatable :: point_props(:,:), history(:,:,:)
   integer, allocatable :: point_order(:)
   integer, parameter :: record_size = 1 + ntens + ntens*ntens + nstatv
   character(len=1024) :: arg, history_path
   real(8), allocatable :: final_stress(:,:)
   real(8), allocatable :: final_statev(:,:)
   real(8), allocatable :: final_pnewdt(:)
   integer, allocatable :: final_status(:)

   call parse_args(nsteps,npoints,material_id,load_case,total_strain,dtime,temp)
   if (npoints > max_points) then
      write(*,*) 'ERROR,npoints exceeds driver max_points'
      stop 2
   endif

   orientation_mode = 0
   order_mode = 0
   history_path = ''
   if (command_argument_count() >= 8) then
      call get_command_argument(8,arg)
      read(arg,*) orientation_mode
   endif
   if (command_argument_count() >= 9) then
      call get_command_argument(9,arg)
      read(arg,*) order_mode
   endif
   if (command_argument_count() >= 10) call get_command_argument(10,history_path)
   if (npoints < 1 .or. nsteps < 1) stop 2
   allocate(point_props(nprops,npoints),point_order(npoints))
   allocate(history(record_size,nsteps,npoints))
   history = 0.d0
   call fill_standard_props(material_id,props)
   do p=1,npoints
      point_props(:,p) = props
      if (orientation_mode == 1) then
         ! Euler angles in radians; six nonsymmetry-equivalent orientations.
         select case(mod(p-1,6))
         case(0)
            point_props(2:4,p) = [0.d0,0.d0,0.d0]
         case(1)
            point_props(2:4,p) = [0.17d0,0.39d0,0.61d0]
         case(2)
            point_props(2:4,p) = [0.73d0,0.28d0,1.13d0]
         case(3)
            point_props(2:4,p) = [1.21d0,0.87d0,0.43d0]
         case(4)
            point_props(2:4,p) = [2.03d0,1.07d0,0.91d0]
         case(5)
            point_props(2:4,p) = [0.51d0,1.31d0,2.17d0]
         end select
      endif
      point_order(p) = p
   enddo
   if (order_mode == 1) then
      ! Deterministic permutation, retaining stable point IDs and orientations.
      do p=npoints,2,-1
         task = 1 + mod(37*p+11,p)
         unit = point_order(p)
         point_order(p) = point_order(task)
         point_order(task) = unit
      enddo
   endif
   call init_global_output_arrays()
   call mod_gspt_ini()
   call mod_wkcp_ini(IB1,IB2)

   allocate(final_stress(npoints,ntens))
   allocate(final_statev(npoints,nstatv))
   allocate(final_pnewdt(npoints))
   allocate(final_status(npoints))
   final_stress = 0.d0
   final_statev = 0.d0
   final_pnewdt = 0.d0
   final_status = 0

!$omp parallel do default(shared) private(p) schedule(runtime)
   do task=1,npoints
      p = point_order(task)
      call run_material_point(p,nsteps,load_case,total_strain,dtime,temp, &
                              point_props(:,p),final_stress(p,:),final_statev(p,:), &
                              final_pnewdt(p),final_status(p),history(:,:,p))
   enddo
!$omp end parallel do

   write(*,'(A)') 'point,status,pnewdt,s11,s22,s33,s12,s13,s23,ivb1,peeq'
   do p=1,npoints
      write(*,'(I0,",",I0,",",*(G0.16,:,","))') p, final_status(p), &
         final_pnewdt(p), final_stress(p,1:6), final_statev(p,41), &
         final_statev(p,167)
   enddo

   if (len_trim(history_path) > 0) then
      open(newunit=unit,file=trim(history_path),status='replace',action='write')
      do p=1,npoints
         do step=1,nsteps
            write(unit,'(I0,",",I0,",",*(ES25.17E3,:,","))') p,step,history(:,step,p)
         enddo
      enddo
      close(unit)
   endif

contains

   subroutine parse_args(nsteps,npoints,material_id,load_case, &
                         total_strain,dtime,temp)
      implicit none
      integer, intent(out) :: nsteps, npoints, material_id, load_case
      real(8), intent(out) :: total_strain, dtime, temp
      character(len=128) :: arg

      nsteps = 10
      npoints = 1
      material_id = 2
      load_case = 1
      total_strain = 1.d-3
      dtime = 1.d0
      temp = 294.d0

      if (command_argument_count() >= 1) then
         call get_command_argument(1,arg)
         read(arg,*) nsteps
      endif
      if (command_argument_count() >= 2) then
         call get_command_argument(2,arg)
         read(arg,*) npoints
      endif
      if (command_argument_count() >= 3) then
         call get_command_argument(3,arg)
         read(arg,*) material_id
      endif
      if (command_argument_count() >= 4) then
         call get_command_argument(4,arg)
         read(arg,*) load_case
      endif
      if (command_argument_count() >= 5) then
         call get_command_argument(5,arg)
         read(arg,*) total_strain
      endif
      if (command_argument_count() >= 6) then
         call get_command_argument(6,arg)
         read(arg,*) dtime
      endif
      if (command_argument_count() >= 7) then
         call get_command_argument(7,arg)
         read(arg,*) temp
      endif
   end subroutine parse_args

   subroutine fill_standard_props(material_id,props)
      implicit none
      integer, intent(in) :: material_id
      real(8), intent(out) :: props(nprops)

      props = 0.d0
      props(1) = dble(material_id)
      props(2:4) = 0.d0
      props(5:8) = 0.d0

      select case(material_id)
      case(1)
         props(9:22) = (/225.d0,1.d0,247000.d0,147000.d0,125000.d0, &
            12.d0,1.d-6,20.d0,20.d0,1500.d0,60.d0,1.d0,1.4d0,2.25d0/)
      case(2)
         props(9:22) = (/225.d0,1.d0,247000.d0,147000.d0,125000.d0, &
            12.d0,0.001d0,20.d0,20.d0,117.d0,180.d0,1.d0,1.4d0,2.25d0/)
      case(3)
         props(9:22) = (/229.d0,1.d0,247000.d0,147000.d0,125000.d0, &
            12.d0,0.001d0,20.d0,20.d0,117.d0,180.d0,1.d0,1.4d0,2.25d0/)
      case(4)
         props(9:22) = (/225.d0,1.d0,256500.d0,111100.d0,77200.d0, &
            12.d0,3.d-39,1250.d0,170.d0,500.d0,0.d0,1.d0,1.4d0,0.d0/)
      case default
         write(*,*) 'ERROR,unsupported material_id for standard fake driver'
         stop 3
      end select
   end subroutine fill_standard_props

   subroutine run_material_point(point_id,nsteps,load_case,total_strain, &
                                 dtime,temp,props,stress,statev,pnewdt,status,records)
      implicit none
      integer, intent(in) :: point_id, nsteps, load_case
      real(8), intent(in) :: total_strain, dtime, temp
      real(8), intent(in) :: props(nprops)
      real(8), intent(out) :: stress(ntens), statev(nstatv), pnewdt
      real(8), intent(out) :: records(record_size,nsteps)
      integer, intent(out) :: status

      character(len=80) :: cmname
      integer :: step
      real(8) :: ddsdde(ntens,ntens), ddsddt(ntens), drplde(ntens)
      real(8) :: stran(ntens), dstran(ntens), time(2)
      real(8) :: predef(1), dpred(1), coords(3), drot(3,3)
      real(8) :: dfgrd0(3,3), dfgrd1(3,3)
      real(8) :: sse, spd, scd, rpl, drpldt, dtemp, celent
      real(8) :: prev_amount, next_amount

      records = 0.d0
      cmname = 'FAKE-UMAT'
      stress = 0.d0
      statev = 0.d0
      ddsdde = 0.d0
      ddsddt = 0.d0
      drplde = 0.d0
      stran = 0.d0
      dstran = 0.d0
      predef = 0.d0
      dpred = 0.d0
      coords = 0.d0
      coords(1) = dble(point_id)
      drot = 0.d0
      drot(1,1) = 1.d0
      drot(2,2) = 1.d0
      drot(3,3) = 1.d0
      dfgrd0 = 0.d0
      dfgrd1 = 0.d0
      dfgrd0(1,1) = 1.d0
      dfgrd0(2,2) = 1.d0
      dfgrd0(3,3) = 1.d0
      dfgrd1 = dfgrd0
      sse = 0.d0
      spd = 0.d0
      scd = 0.d0
      rpl = 0.d0
      drpldt = 0.d0
      dtemp = 0.d0
      celent = 1.d0
      pnewdt = 1.d0
      status = 0

      do step=1,nsteps
         prev_amount = total_strain*dble(step-1)/dble(nsteps)
         next_amount = total_strain*dble(step)/dble(nsteps)
         call make_deformation(load_case,prev_amount,dfgrd0)
         call make_deformation(load_case,next_amount,dfgrd1)
         time(1) = dtime*dble(step-1)
         time(2) = dtime*dble(step-1)
         call umat(stress,statev,ddsdde,sse,spd,scd,rpl, &
                   ddsddt,drplde,drpldt,stran,dstran,time,dtime, &
                   temp,dtemp,predef,dpred,cmname,ndi,nshr,ntens,nstatv, &
                   props,nprops,coords,drot,pnewdt,celent,dfgrd0,dfgrd1, &
                   point_id,1,1,1,1,step)
         records(:,step) = [pnewdt,stress,reshape(ddsdde,[ntens*ntens]),statev]
         if (pnewdt < 1.d0) then
            status = 1
            exit ! No increment acceptance after a requested cutback.
         endif
      enddo
   end subroutine run_material_point

   subroutine make_deformation(load_case,amount,dfgrd)
      implicit none
      integer, intent(in) :: load_case
      real(8), intent(in) :: amount
      real(8), intent(out) :: dfgrd(3,3)

      dfgrd = 0.d0
      dfgrd(1,1) = 1.d0
      dfgrd(2,2) = 1.d0
      dfgrd(3,3) = 1.d0

      select case(load_case)
      case(1)
         dfgrd(1,1) = 1.d0 + amount
      case(2)
         dfgrd(2,2) = 1.d0 + amount
      case(3)
         dfgrd(3,3) = 1.d0 + amount
      case(4)
         dfgrd(1,2) = amount
      case(5)
         dfgrd(1,3) = amount
      case(6)
         dfgrd(2,3) = amount
      case default
         write(*,*) 'ERROR,unsupported load_case'
         stop 4
      end select
   end subroutine make_deformation

end program umat_fake_driver

subroutine getjobname(job_name,ljname)
   implicit none
   character(len=*), intent(out) :: job_name
   integer, intent(out) :: ljname
   job_name = 'fake_umat'
   ljname = 9
end subroutine getjobname

subroutine getoutdir(xoutdir,lxoutdir)
   implicit none
   character(len=*), intent(out) :: xoutdir
   integer, intent(out) :: lxoutdir
   xoutdir = '/private/tmp'
   lxoutdir = 12
end subroutine getoutdir
