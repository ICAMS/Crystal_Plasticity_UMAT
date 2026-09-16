! Deterministically interleave grain preparation; no timing-dependent race needed.
program orientation_ownership
   use mod_umat_phases
   implicit none
   type(mat_param_set) :: ma, mb
   type(gp_context) :: a, b
   real(8) :: pa(22), pb(22), sa(240), sb(240), q(3,3), expected(3,3)
   real(8) :: bk(Nslp_mx,3), identity(3,3)
   logical :: ca, cb
   integer :: na, nb, ia, ib, i
   call mod_gspt_ini()
   pa=0.d0
   pa(1)=2.d0
   pa(2:4)=[0.17d0,0.39d0,0.61d0]
   pa(9:22)=[225.d0,1.d0,247000.d0,147000.d0,125000.d0, &
             12.d0,0.001d0,20.d0,20.d0,117.d0,180.d0,1.d0,1.4d0,2.25d0]
   pb=pa
   pb(2:4)=[1.21d0,0.87d0,0.43d0]
   identity=0.d0
   do i=1,3
      identity(i,i)=1.d0
   enddo
   sa=0.d0
   sb=0.d0
   call umat_prepare_phase(22,pa,240,ma,ca,na,ia)
   call umat_prepare_phase(22,pb,240,mb,cb,nb,ib)
   if (.not.ca .or. .not.cb) stop 1
   ! Poison the legacy scratch after both prepares: the local path must ignore it.
   eang00=-99.d0
   Fe=-99.d0
   Fp=-99.d0
   csM0=-99.d0
   call umat_load_phase(1,1,1.d0,294.d0,[0.d0,0.d0],pa, &
        identity,identity,sa,ma,ca,na,ia,a,q,bk)
   call icams_Eang2Q(pa(2),pa(3),pa(4),expected)
   if (maxval(abs(a%Fe-expected))>1.d-14) stop 2
   if (maxval(abs(a%Fp-transpose(expected)))>1.d-14) stop 3
   if (maxval(abs(matmul(a%Fe,a%Fp)-identity))>1.d-14) stop 4
   if (any(a%csM0/=0.d0)) stop 5
   call umat_load_phase(2,1,1.d0,294.d0,[0.d0,0.d0],pb, &
        identity,identity,sb,mb,cb,nb,ib,b,q,bk)
   if (maxval(abs(a%Fe-expected))>1.d-14) stop 6
   if (any(Fe/=-99.d0) .or. any(Fp/=-99.d0)) stop 7
   ! Subsequent increments must recover the point's matrices from STATEV.
   sa(23)=sa(23)+0.01d0
   call umat_load_phase(1,1,1.d0,294.d0,[1.d0,1.d0],pa, &
        identity,identity,sa,ma,ca,na,ia,a,q,bk)
   if (a%Fe(1,1)/=sa(23)) stop 8
   if (maxval(abs(a%Fp-transpose(expected)))>1.d-14) stop 9
end program orientation_ownership
