program thread_initialization
   use mod_gaussp
   use omp_lib
   implicit none
   integer :: i, j, k, tid, iteration, failures, team_size, expected
   ! Constants must be valid before mod_gspt_ini, including on a fresh worker.
   failures=0
   team_size=0
   call omp_set_dynamic(.false.)
!$omp parallel num_threads(4) private(tid,i,j,k,iteration,expected) reduction(+:failures)
   tid=omp_get_thread_num()
!$omp single
   team_size=omp_get_num_threads()
!$omp end single
   eang00=dble(tid)+0.25d0
   eang0=dble(tid)+0.5d0
   eang=dble(tid)+0.75d0
!$omp barrier
   do iteration=1,100
      call mod_gspt_ini()
      if (any(eang00/=dble(tid)+0.25d0)) failures=failures+1
      if (any(eang0/=dble(tid)+0.5d0)) failures=failures+1
      if (any(eang/=dble(tid)+0.75d0)) failures=failures+1
      do i=1,Nslp_mx
         do j=1,Nslp_mx
            expected=0
            if (i==j) expected=1
            if (XInn(i,j)/=expected) failures=failures+1
            if (i<=9.and.j<=9) then
               if (XI99(i,j)/=expected) failures=failures+1
            endif
            if (i<=6.and.j<=6) then
               if (XI66(i,j)/=expected) failures=failures+1
            endif
            if (i<=3.and.j<=3) then
               if (XI33(i,j)/=expected) failures=failures+1
            endif
         enddo
      enddo
      do i=1,3
         do j=1,3
            do k=1,3
               expected=(j-i)*(k-i)*(k-j)/2
               if (XI333(i,j,k)/=expected) failures=failures+1
            enddo
         enddo
      enddo
      if (any(ib1/=[1,2,3,1,1,2,2,3,3])) failures=failures+1
      if (any(ib2/=[1,2,3,2,3,3,1,1,2])) failures=failures+1
   enddo
!$omp end parallel
   if (team_size/=4) error stop 'test requires four active threads'
   if (failures/=0) error stop 'thread-local angles or tensor constants corrupted'
end program thread_initialization
