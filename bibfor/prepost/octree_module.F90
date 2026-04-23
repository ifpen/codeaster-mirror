! --------------------------------------------------------------------
! Copyright (C) 1991 - 2026 - EDF - www.code-aster.org
! This file is part of code_aster.
!
! code_aster is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! code_aster is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with code_aster.  If not, see <http://www.gnu.org/licenses/>.
! --------------------------------------------------------------------
!
module octree_module
!
    implicit none
!
    private
!
#include "asterf_types.h"
#include "asterc/r8maem.h"
!
    integer(kind=8), parameter, private :: nb_div = 2
    aster_logical, parameter, private :: debug = ASTER_TRUE
!
    type BBox
        real(kind=8) :: center(3) = 0.d0, h(3) = 0.d0
        integer(kind=8) :: dim = 3
! ----- member functions
    contains
        procedure, public, pass :: has_intersection
        procedure, public, pass :: pts_is_inside
        procedure, public, pass :: is_strictly_inside => is_strictly_inside_node
    end type
!
    type OctreeNode
        real(kind=8) :: xmin(3) = 0.d0, dx(3) = 0.0
        type(OctreeNode), pointer :: children(:, :, :) => null()
        integer(kind=8), allocatable :: list_pts(:)
        real(kind=8), allocatable :: coordinates(:, :)
        integer(kind=8) :: nb_pts = 0, dim = 3, level = 0
! ----- member functions
    contains
        procedure, public, pass :: split => split_node
        procedure, public, pass :: free => free_node
        procedure, public, pass :: is_strictly_inside => is_strictly_inside_pts
        procedure, public, pass :: find_children
        procedure, public, pass :: is_leaf
        procedure, private, pass :: get_pts_around => get_pts_around_node
        procedure, public, pass :: get_number_of_level => get_number_of_level_node
    end type OctreeNode
!
    type Octree
        real(kind=8) :: xmin(3) = 0.d0, xmax(3) = 0.d0
        integer(kind=8) :: max_pts_by_node = 2, max_level = 50, dim = 3, nb_pt = 0
        integer(kind=8), allocatable :: list_pts(:)
        type(OctreeNode) :: node
! ----- member functions
    contains
        procedure, public, pass :: init => init_octree
        procedure, public, pass :: free => free_octree
        procedure, public, pass :: get_pts_around => get_pts_around_tree
        procedure, public, pass :: get_alt_1_pts_around => get_al1_pts_around_tree
        procedure, public, pass :: get_number_of_level => get_number_of_level_tree
    end type
!
!===================================================================================================
!
!===================================================================================================
!
    public :: Octree, OctreeNode
    private :: split_node, free_node, init_octree, free_octree
!
contains
!
    function has_intersection(this, onode) result(inter)
        class(BBox), intent(in) :: this
        class(OctreeNode), intent(in) :: onode
!
        aster_logical :: inter
        real(kind=8) :: center_on(3)
        inter = ASTER_FALSE
!
        center_on = onode%xmin+0.5d0*onode%dx
!
        if (abs(center_on(1)-this%center(1)) <= 0.5d0*(onode%dx(1)+this%h(1))) then
            if (abs(center_on(2)-this%center(2)) <= 0.5d0*(onode%dx(2)+this%h(2))) then
                if (this%dim == 3) then
                    if (abs(center_on(3)-this%center(3)) <= 0.5d0*(onode%dx(3)+this%h(3))) then
                        inter = ASTER_TRUE
                    end if
                else
                    inter = ASTER_TRUE
                end if
            end if
        end if
!
    end function
!
! ==================================================================================================
!
    function pts_is_inside(this, coor) result(inter)
        class(BBox), intent(in) :: this
        real(kind=8), intent(in) :: coor(3)
!
        aster_logical :: inter
        inter = ASTER_FALSE
!
        if (abs(coor(1)-this%center(1)) <= 0.5d0*this%h(1)) then
            if (abs(coor(2)-this%center(2)) <= 0.5d0*this%h(2)) then
                if (this%dim == 3) then
                    if (abs(coor(3)-this%center(3)) <= 0.5d0*this%h(3)) then
                        inter = ASTER_TRUE
                    end if
                else
                    inter = ASTER_TRUE
                end if
            end if
        end if
!
    end function
!
! ==================================================================================================
!
!
    function is_strictly_inside_node(this, onode) result(inside)
        class(BBox), intent(in) :: this
        class(OctreeNode), intent(in) :: onode
!
        aster_logical :: inside
        inside = ASTER_FALSE
!
        if (onode%xmin(1) <= (this%center(1)-0.5d0*this%h(1)) .and. &
            (this%center(1)+0.5d0*this%h(1)) <= (onode%xmin(1)+onode%dx(1))) then
            if (onode%xmin(2) <= (this%center(2)-0.5d0*this%h(2)) .and. &
                (this%center(2)+0.5d0*this%h(2)) <= (onode%xmin(2)+onode%dx(2))) then
                if (this%dim == 3) then
                    if (onode%xmin(3) <= (this%center(3)-0.5d0*this%h(3)) .and. &
                        (this%center(3)+0.5d0*this%h(3)) <= (onode%xmin(3)+onode%dx(3))) then
                        inside = ASTER_TRUE
                    end if
                else
                    inside = ASTER_TRUE
                end if
            end if
        end if
!
    end function
!
! ==================================================================================================
!
    recursive subroutine free_node(this)
        class(OctreeNode), intent(inout) :: this
!
        integer(kind=8) :: i, j, k
!
        if (allocated(this%list_pts)) then
            this%nb_pts = 0
            deallocate (this%list_pts)
            deallocate (this%coordinates)
        else
            if (associated(this%children)) then
                do i = 1, size(this%children, 1)
                    do j = 1, size(this%children, 2)
                        do k = 1, size(this%children, 3)
                            call this%children(i, j, k)%free()
                        end do
                    end do
                end do
                deallocate (this%children)
            end if
        end if
!
    end subroutine free_node
!
! ==================================================================================================
!
!
    function is_leaf(this) result(leaf)
        class(OctreeNode), intent(in) :: this
!
        aster_logical :: leaf
!
        leaf = this%nb_pts == 0 .or. allocated(this%coordinates)
!
    end function
!
! ==================================================================================================
!
!
    function is_strictly_inside_pts(this, coor, distance) result(inside)
        class(OctreeNode), intent(in) :: this
        real(kind=8), intent(in) :: distance, coor(3)
!
        aster_logical :: inside
        inside = ASTER_FALSE
!
        if (this%xmin(1) <= (coor(1)-distance) .and. &
            (coor(1)+distance) <= (this%xmin(1)+this%dx(1))) then
            if (this%xmin(2) <= (coor(2)-distance) .and. &
                (coor(2)+distance) <= (this%xmin(2)+this%dx(2))) then
                if (this%dim == 3) then
                    if (this%xmin(3) <= (coor(3)-distance) .and. &
                        (coor(3)+distance) <= (this%xmin(3)+this%dx(3))) then
                        inside = ASTER_TRUE
                    end if
                else
                    inside = ASTER_TRUE
                end if
            end if
        end if
!
    end function
!
! ==================================================================================================
!
    subroutine find_children(this, coor, child_id)
        class(OctreeNode), intent(in) :: this
        real(kind=8), intent(in) :: coor(3)
        integer(kind=8), intent(inout) :: child_id(3)
!
        integer(kind=8) :: i_dim
        real(kind=8) :: rid, dx_child(3)
!
        dx_child = this%dx/real(nb_div, kind=8)
        child_id = 1
        do i_dim = 1, this%dim
            rid = (coor(i_dim)-this%xmin(i_dim))/dx_child(i_dim)
            child_id(i_dim) = max(1, min(nb_div, ceiling(rid, kind=8)))
        end do
!
    end subroutine
!
! ==================================================================================================
!
    subroutine free_octree(this)
        class(Octree), intent(inout) :: this
!
        call this%node%free()
!
    end subroutine free_octree
!
! ==================================================================================================
!
    subroutine init_octree(this, nb_pts, coordinates, max_pts_by_node, max_level)
        class(Octree), intent(inout) :: this
        integer(kind=8), intent(in) :: nb_pts, max_pts_by_node, max_level
        real(kind=8), intent(in) :: coordinates(3, nb_pts)
!
        integer(kind=8), allocatable :: sub_pts(:)
        integer(kind=8) :: i_pt, i_dim
        real(kind=8) :: x_offset(3)
!
        this%max_pts_by_node = max_pts_by_node
        this%max_level = max_level
!
        this%xmin = R8MAEM()
        this%xmax = -R8MAEM()
!
        if (nb_pts == 0) then
            goto 999
        end if
!
        do i_pt = 1, nb_pts
            do i_dim = 1, 3
                this%xmin(i_dim) = min(this%xmin(i_dim), coordinates(i_dim, i_pt))
                this%xmax(i_dim) = max(this%xmax(i_dim), coordinates(i_dim, i_pt))
            end do
        end do
!
        this%dim = 3
        if (abs(this%xmax(3)-this%xmin(3)) < 1d-12) then
            this%dim = 2
        end if
        ! increase slithly the first bounding box
        x_offset = abs(this%xmax-this%xmin)
        this%xmin(1:this%dim) = this%xmin(1:this%dim)-0.05*x_offset(1:this%dim)
        this%xmax(1:this%dim) = this%xmax(1:this%dim)+0.05*x_offset(1:this%dim)
!
        this%node%xmin = this%xmin
        this%node%dx = 0.d0
        this%node%level = 0
        do i_dim = 1, this%dim
            this%node%dx(i_dim) = (this%xmax(i_dim)-this%node%xmin(i_dim))
        end do
!
        allocate (sub_pts(nb_pts))
!
        do i_pt = 1, nb_pts
            sub_pts(i_pt) = i_pt
        end do
!
        call this%node%split(nb_pts, coordinates, this%dim, max_pts_by_node, max_level, &
                             nb_pts, sub_pts)
        deallocate (sub_pts)
!
999     continue
!
    end subroutine init_octree
!
! ==================================================================================================
!
    subroutine split_node(this, nb_pts, coordinates, dim, max_nb_pts, max_level, &
                          nb_sub_pt, list_sub_pt)
        class(OctreeNode), intent(inout) :: this
        integer(kind=8), intent(in) :: nb_pts, dim, nb_sub_pt, list_sub_pt(nb_sub_pt), max_nb_pts
        integer(kind=8), intent(in) :: max_level
        real(kind=8), intent(in) :: coordinates(3, nb_pts)
!
        integer(kind=8), allocatable :: sub_pts(:, :, :, :)
        integer(kind=8) :: i_pt, oi(3), pt_id, i, j, k, nb_sub_pts(3, 3, 3), kmax
        real(kind=8) :: dx_child(3)
!
        this%nb_pts = nb_sub_pt
        this%dim = dim
        if (nb_sub_pt <= max_nb_pts .or. this%level == max_level) then
            if (nb_sub_pt > 0) then
                allocate (this%list_pts(nb_sub_pt))
                allocate (this%coordinates(3, nb_sub_pt))
!
                this%list_pts(1:nb_sub_pt) = list_sub_pt(1:nb_sub_pt)
                do i_pt = 1, nb_sub_pt
                    this%coordinates(1:3, i_pt) = coordinates(1:3, list_sub_pt(i_pt))
                end do
            end if
        else
!
            kmax = nb_div
            if (this%dim < 3) then
                kmax = 1
            end if
!
            allocate (this%children(nb_div, nb_div, kmax))
!
            dx_child = this%dx/real(nb_div, kind=8)
!
            do i = 1, nb_div
                do j = 1, nb_div
                    do k = 1, kmax
                        this%children%level = this%level+1
                        this%children(i, j, k)%dx = dx_child
                        this%children(i, j, k)%xmin(1) = this%xmin(1)+ &
                                                         (i-1)*this%children(i, j, k)%dx(1)
                        this%children(i, j, k)%xmin(2) = this%xmin(2)+ &
                                                         (j-1)*this%children(i, j, k)%dx(2)
                        this%children(i, j, k)%xmin(3) = this%xmin(3)+ &
                                                         (k-1)*this%children(i, j, k)%dx(3)
                    end do
                end do
            end do
!
            allocate (sub_pts(nb_div, nb_div, kmax, nb_sub_pt))
            nb_sub_pts = 0
!
            do i_pt = 1, nb_sub_pt
                pt_id = list_sub_pt(i_pt)
                call this%find_children(coordinates(:, pt_id), oi)
                nb_sub_pts(oi(1), oi(2), oi(3)) = nb_sub_pts(oi(1), oi(2), oi(3))+1
                sub_pts(oi(1), oi(2), oi(3), nb_sub_pts(oi(1), oi(2), oi(3))) = pt_id
            end do
!
            do i = 1, nb_div
                do j = 1, nb_div
                    do k = 1, kmax
                        call this%children(i, j, k)%split(nb_pts, coordinates, dim, max_nb_pts, &
                                                          max_level, &
                                                          nb_sub_pts(i, j, k), sub_pts(i, j, k, :))
                    end do
                end do
            end do
!
            deallocate (sub_pts)
        end if
!
    end subroutine split_node
!
! ==================================================================================================
!
    subroutine get_pts_around_tree(this, coor, distance, nb_pts, list_pts)
!
!       Find at points inside the boundign box around the given point.
!
        class(Octree), intent(in) :: this
        integer(kind=8), intent(inout) :: nb_pts, list_pts(*)
        real(kind=8), intent(in) :: coor(3), distance
!
        type(BBox) :: box
!
        box%center = coor
        box%h = distance
        box%dim = this%dim
!
        nb_pts = 0
        call this%node%get_pts_around(box, nb_pts, list_pts)
!
    end subroutine
!
!
! ==================================================================================================
!
    subroutine get_al1_pts_around_tree(this, coor, distance, nb_pts, list_pts)
!
!       Find at least one point around the given point.
!       Double size of bounding box until necessary
!
        class(Octree), intent(in) :: this
        integer(kind=8), intent(inout) :: nb_pts, list_pts(*)
        real(kind=8), intent(in) :: coor(3), distance
!
        type(BBox) :: box
!
        box%center = coor
        box%h = distance
        box%dim = this%dim
!
        nb_pts = 0
        do while (nb_pts == 0)
            call this%node%get_pts_around(box, nb_pts, list_pts)
            box%h = 2.d0*box%h
        end do
!
    end subroutine
!
! ==================================================================================================
!
    recursive subroutine get_pts_around_node(this, box, nb_pts, list_pts)
        class(OctreeNode), intent(in) :: this
        integer(kind=8), intent(inout) :: nb_pts, list_pts(*)
        type(BBox), intent(in) :: box
!
        integer(kind=8) :: ocx, ocy, ocz, ocz_max, i_pt
!
        if (this%nb_pts > 0) then
            if (box%has_intersection(this)) then
                if (this%is_leaf()) then
                    do i_pt = 1, this%nb_pts
                        if (box%pts_is_inside(this%coordinates(:, i_pt))) then
                            nb_pts = nb_pts+1
                            list_pts(nb_pts) = this%list_pts(i_pt)
                        end if
                    end do
                else
                    ocz_max = nb_div
                    if (this%dim < 3) then
                        ocz_max = 1
                    end if
!
                    do ocx = 1, nb_div
                        do ocy = 1, nb_div
                            do ocz = 1, ocz_max
                                call this%children(ocx, ocy, ocz)%get_pts_around(box, &
                                                                                 nb_pts, list_pts)
                            end do
                        end do
                    end do
                end if
            end if
        end if
!
    end subroutine
!
! ==================================================================================================
!
    integer(kind=8) function get_number_of_level_tree(this)
!
!       Find at points inside the boundign box around the given point.
!
        class(Octree), intent(in) :: this
!
        get_number_of_level_tree = 0
        get_number_of_level_tree = this%node%get_number_of_level()
!
    end function
!
! ==================================================================================================
!
    recursive integer(kind=8) function get_number_of_level_node(this)
        class(OctreeNode), intent(in) :: this
!
        integer(kind=8) :: ocz_max, ocx, ocy, ocz, child_level
!
        get_number_of_level_node = this%level
        if (.not. this%is_leaf()) then
            ocz_max = nb_div
            if (this%dim < 3) then
                ocz_max = 1
            end if
!
            do ocx = 1, nb_div
                do ocy = 1, nb_div
                    do ocz = 1, ocz_max
                        child_level = this%children(ocx, ocy, ocz)%get_number_of_level()
                        get_number_of_level_node = max(get_number_of_level_node, &
                                                       child_level)
                    end do
                end do
            end do
        end if
!
    end function
!
end module
