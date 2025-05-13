c================================================================
c
c    octree-search code
c
c================================================================
C This code is based on implementation from:
C
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C Copyright (c) 2019, Henrique Miranda
C All rights reserved.
C
C Redistribution and use in source and binary forms, with or without
C modification, are permitted provided that the following conditions 
C are met:
C* Redistributions of source code must retain the above copyright
C  notice, this list of conditions and the following disclaimer.
C* Redistributions in binary form must reproduce the above copyright
C  notice, this list of conditions and the following disclaimer in the
C  documentation and/or other materials provided with the distribution.
C* Neither the name of the fortran_snippets project nor the
C  names of its contributors may be used to endorse or promote products
C  derived from this software without specific prior written permission.
C
C THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
C "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
C LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS 
C FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE 
C COPYRIGHT HOLDER BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
C SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
C LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
C LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) 
C HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, 
C STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
C ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED 
C OF THE POSSIBILITY OF SUCH DAMAGE.
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
C https://github.com/henriquemiranda/fortran_snippets
C 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C    editor:    Rene Zandomeni
C    date:      13.03.2025
C    changelog:
C    13.03.2025 No deallocation of octree implemented, not necessary 
C               in current framework. Subroutine is only called once 
C               during runtime
C    12.05.2025 For larger models, deallocation may make sense. 
C               Next Step deallocation octrees
C    
C    
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      module mod_octree
c
      implicit none
c
      real(8), parameter :: zero=0.0d0
      real(8), parameter :: half=0.5d0 
      real(8), parameter :: one=1.d0
      real(8), parameter :: two=2.d0 
      real(8), protected :: shift(3,2,8)
c
      type :: oct_leaf_t
      type(oct_leaf_t), pointer :: childs(:)
      integer, allocatable :: ids(:)
      end type oct_leaf_t
      type :: oct_t
      real(8) :: hi(3)
      real(8) :: lo(3)
      integer :: max_points
      real(8), pointer :: points(:,:)
      type(oct_leaf_t) :: first
      end type oct_t
c
      contains
c
      !-----------------------------------------
      ! Initialize octree with all points
      !-----------------------------------------
      type(oct_t)function octree_init(points,max_points,lo,hi) 
     &  result (new)
        integer,intent(in) :: max_points
        real(8), target, intent(in) :: points(:,:)
        real(8), intent(in) :: lo(3), hi(3)
        real,parameter :: ieps = 0.1D0
        integer,allocatable :: ids(:)
        integer :: ii, jj, kk, npoints, ioctant
        ! determine shift
        do ii=0,1
          do jj=0,1
            do kk=0,1
              ioctant = ii*4+jj*2+kk+1
              shift(:,1,ioctant) = [half*ii,half*jj,half*kk]
              shift(:,2,ioctant)=shift(:,1,ioctant)+[half,half,half]
            end do
          end do
        end do
        ! box dimensions
        new%lo = lo
        new%hi = hi
        ! first octree contains all the points
        npoints = size(points,2)
        new%points => points
        new%max_points = max_points
        ids = [(ii,ii=1,npoints)]
        new%first = octree_node_build(new,new%lo,new%hi,npoints,ids)
c
      end function octree_init
      !-----------------------------------------
      ! Recursive function to build octree octants with points
      !-----------------------------------------
      type(oct_leaf_t)
     & recursive function octree_node_build(octree,lo,hi,nids,ids) 
     &  result (new)
        type(oct_t),intent(in) :: octree
        integer,intent(in) :: nids
        integer,intent(in) :: ids(nids)
        integer :: id, counter, ioctant, ipoint
        integer :: octants(nids)
        integer :: new_ids(nids)
        real(8) :: lo(3), hi(3), new_lo(3), new_hi(3)
c
        if (nids<octree%max_points) then
          allocate(new%ids(nids))
          new%ids = ids(:nids)
          return
        end if
      ! call get octants
        call get_octants(lo,hi,nids,ids,octree%points,octants)
c
        allocate(new%childs(8))
        ! Get number of points in this octant
        do ioctant=1,8
          counter = 0
          do id=1,nids
            ipoint = ids(id)
            if (octants(id) /= ioctant) cycle
              counter = counter + 1
            new_ids(counter) = ipoint
          end do
      ! Build this octant
          call get_lo_hi(lo,hi,new_lo,new_hi,ioctant)
          new%childs(ioctant) = 
     &    octree_node_build(octree,new_lo,new_hi,counter,new_ids)
        end do
c
      end function octree_node_build
      !-----------------------------------------
      ! Function to find  nearest neighbor in tree of given point
      !-----------------------------------------
      integer function octree_find(octree,point,dist) 
     & result(closest_id)
        ! find the closest point in the box that contains it
        type(oct_t),target,intent(in) :: octree
        real(8),intent(in) :: point(3)
        real(8),intent(out) :: dist
        type(oct_leaf_t),pointer :: octn
        integer :: id, ipoint, ioctant
        real(8) :: hi(3),lo(3),hi_out(3),lo_out(3)
        real(8) :: trial_dist
        closest_id = 0
        dist = huge(dist)
        octn => octree%first
        lo = octree%lo
        hi = octree%hi
        ! check if the point is inside the initial box
        trial_dist = box_dist(lo,hi,point)
        if (trial_dist>0) then
          closest_id = -1
          return
        end if
        do
          ! if leaf node
          if (allocated(octn%ids)) then
            do id=1,size(octn%ids)
              ipoint = octn%ids(id)
              trial_dist = dist_points(octree%points(:,ipoint),point)
              if (trial_dist > dist) cycle
                dist = trial_dist
                closest_id = ipoint
            end do
            return
          end if
          ! get octant of this point
          ioctant = get_octant_lohi(lo,hi,point)
          ! point to this node
          octn => octn%childs(ioctant)
          ! get lo and hi
          call get_lo_hi(lo,hi,lo_out,hi_out,ioctant)
          lo = lo_out; hi = hi_out
        end do
c
      end function octree_find
      !-----------------------------------------
      ! Calculate squared distance between two points
      !-----------------------------------------
      pure real(8) function dist_points(p1,p2) result(dist)
        real(8),intent(in) :: p1(3),p2(3)
        dist = pow2(p1(1)-p2(1))+
     &         pow2(p1(2)-p2(2))+
     &         pow2(p1(3)-p2(3))
c
      end function dist_points
      !-----------------------------------------
      ! Calculate squared distance point to box
      !-----------------------------------------
      pure real(8) function box_dist(lo,hi,po) result(dist)
        real(8),intent(in) :: lo(3), hi(3), po(3)
        dist = zero
        if (po(1)<lo(1)) dist = dist + pow2(po(1)-lo(1))
        if (po(1)>hi(1)) dist = dist + pow2(po(1)-hi(1))
        if (po(2)<lo(2)) dist = dist + pow2(po(2)-lo(2))
        if (po(2)>hi(2)) dist = dist + pow2(po(2)-hi(2))
        if (po(3)<lo(3)) dist = dist + pow2(po(3)-lo(3))
        if (po(3)>hi(3)) dist = dist + pow2(po(3)-hi(3))
c
      end function box_dist
      !-----------------------------------------
      ! Calculate power 2
      !-----------------------------------------
      pure real(8) function pow2(x) result(x2)
        real(8),intent(in) :: x
        x2 = x*x
c
      end function pow2
      !-----------------------------------------
      ! Get octant
      !-----------------------------------------
      pure integer function get_octant(mi,po) result(ioctant)
        real(8),intent(in) :: po(3), mi(3)
        integer :: ii,jj,kk
        ii = 0; if (po(1)>=mi(1)) ii = 1
        jj = 0; if (po(2)>=mi(2)) jj = 1
        kk = 0; if (po(3)>=mi(3)) kk = 1
        ioctant = ii*4+jj*2+kk+1
c
      end function get_octant
      !-----------------------------------------
      ! Get min max of octant
      !-----------------------------------------
      pure integer function get_octant_lohi(lo,hi,po) result(ioctant)
        real(8),intent(in) :: lo(3),hi(3),po(3)
        real(8) :: mi(3)
        mi = half*(hi+lo)
        ioctant = get_octant(mi,po)
c
      end function get_octant_lohi
      !-----------------------------------------
      ! From the list of points get the corrsponding octant
      !-----------------------------------------
      pure subroutine get_octants(lo,hi,nids,ids,points,octants)
        real(8),intent(in) :: lo(3), hi(3)
        real(8),intent(in) :: points(:,:)
        real(8) :: mi(3)
        integer,intent(in) :: nids
        integer,intent(in) :: ids(nids)
        integer,intent(out) :: octants(nids)
        integer :: id, ipoint
        ! calculate midpoint
        mi = half*(hi+lo)
        do id=1,nids
          ipoint = ids(id)
          octants(id) = get_octant(mi,points(:,ipoint))
        end do
c
      end subroutine get_octants
      !-----------------------------------------
      ! Split octant into new octants
      !-----------------------------------------
      pure subroutine get_lo_hi(lo_in,hi_in,lo_out,hi_out,ioctant)
        integer,intent(in) :: ioctant
        real(8),intent(in) :: lo_in(3), hi_in(3)
        real(8),intent(out) :: lo_out(3), hi_out(3)
        real(8) :: de(3)
        de = hi_in-lo_in
        lo_out = lo_in + shift(:,1,ioctant)*de
        hi_out = lo_in + shift(:,2,ioctant)*de
c
      end subroutine get_lo_hi
c
c-----End of module octree
c
      end module mod_octree