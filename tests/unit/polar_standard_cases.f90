program polar_standard_cases
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   real(8) :: identity(3,3), q(3,3), stretch(3,3), f(3,3), u(3,3), r(3,3)
   real(8) :: rotation(3,3), gamma, scale
   integer :: i, case_id, status
   identity=0.d0
   do i=1,3
      identity(i,i)=1.d0
   enddo
   ! Independent prescribed rotation about z, 0.7 radians.
   q=identity
   q(1,1)=cos(0.7d0); q(1,2)=-sin(0.7d0)
   q(2,1)=sin(0.7d0); q(2,2)=cos(0.7d0)
   do case_id=1,7
      stretch=identity
      rotation=identity
      select case(case_id)
      case(1) ! identity
      case(2) ! rigid rotation
         rotation=q
      case(3) ! uniform dilation
         stretch=2.d0*identity
      case(4) ! uniaxial stretch, repeated transverse stretches
         stretch(1,1)=2.d0
      case(5) ! three distinct stretches
         stretch(2,2)=2.d0; stretch(3,3)=3.d0
      case(6) ! rotated principal stretch axes, plus rigid rotation
         stretch(1,1)=0.8d0; stretch(2,2)=1.1d0; stretch(3,3)=1.4d0
         stretch=matmul(q,matmul(stretch,transpose(q)))
         rotation=q
      case(7) ! simple shear with analytic right polar factors
         gamma=0.5d0
         scale=sqrt(4.d0+gamma**2)
         rotation(1,1)=2.d0/scale; rotation(2,2)=2.d0/scale
         rotation(1,2)=gamma/scale; rotation(2,1)=-gamma/scale
         stretch(1,1)=2.d0/scale
         stretch(1,2)=gamma/scale; stretch(2,1)=gamma/scale
         stretch(2,2)=(2.d0+gamma**2)/scale
      end select
      f=matmul(rotation,stretch)
      status=0
      call polar_decomp(f,u,r,status)
      write(*,'(I0,",",I0,",",L1,",",*(ES24.16E3,:,","))') &
         case_id,status,all(ieee_is_finite(u)).and.all(ieee_is_finite(r)), &
         maxval(abs(u-stretch)),maxval(abs(r-rotation)), &
         maxval(abs(matmul(transpose(r),r)-identity)), &
         maxval(abs(matmul(r,u)-f)),maxval(abs(u-transpose(u)))
   enddo
end program polar_standard_cases
