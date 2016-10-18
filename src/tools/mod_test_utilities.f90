module mod_test_utilities
#include <messenger.h>
    use mod_kinds,                  only: rk,ik
    use mod_constants,              only: ZERO, ONE, TWO, THREE, FOUR, FIVE, SIX
    use mod_chidg_mpi,              only: IRANK
    use mod_plot3d_utilities,       only: get_block_points_plot3d, &
                                          get_block_elements_plot3d, &
                                          get_block_boundary_faces_plot3d
    use mod_hdf_utilities,          only: initialize_file_hdf, add_domain_hdf, &
                                          open_domain_hdf, close_domain_hdf, &
                                          set_bc_patch_hdf, add_bc_state_hdf, &
                                          set_contains_grid_hdf, close_file_hdf, close_hdf
    use mod_bc,                     only: create_bc
    use mod_gridgen_cylinder,       only: create_mesh_file__cylinder_abutting

    use type_point,                 only: point_t
    use type_bc_state,              only: bc_state_t
    use type_domain_connectivity,   only: domain_connectivity_t
    use hdf5
    implicit none


contains



    !>  Create an actual ChiDG-formatted grid file that could be
    !!  read in by a test. Also with initialized boundary conditions.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   10/17/2016
    !!
    !!
    !---------------------------------------------------------------------------
    subroutine create_mesh_file(selector, filename)
        character(*),   intent(in)  :: selector
        character(*),   intent(in)  :: filename

        integer(ik) :: ierr


        ! Generate grid file base on selector case.
        select case (trim(selector))
            case("D2_E8_M1 : Overlapping : Matching")
                call create_mesh_file_D2E8M1_overlapping_matching(filename)
            case("Cylinder : Diagonal : Matching")
                call create_mesh_file__cylinder_abutting(filename)

            case default
                call chidg_signal(FATAL,"create_mesh_file: There was no valid case that matched the incoming string")

        end select


    end subroutine create_mesh_file
    !***************************************************************************







    !>  Generate a set of points for a mesh. String input calls specialized
    !!  procedure for generating the points
    !!
    !!  @author Nathan A. Wukie
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   5/24/2016
    !!
    !!  @param[in]      string          Character string used to select a specialized 
    !!                                  meshgen call
    !!  @param[inout]   nodes           Array of node coordinates for the grid
    !!  @param[inout]   connectivity    Connectivity data for the grid
    !--------------------------------------------------------------------
    subroutine create_mesh(string,nodes,connectivity)
        character(*),                   intent(in)      :: string
        type(point_t),  allocatable,    intent(inout)   :: nodes(:)
        type(domain_connectivity_t),    intent(inout)   :: connectivity

        integer(ik)                                     :: idomain, mapping, ielem
        integer(ik),    allocatable                     :: elements(:,:)
        real(rk),       allocatable, dimension(:,:,:)   :: xcoords,ycoords,zcoords

        select case (trim(string))
            case ('1x1x1','111')
                call meshgen_1x1x1_linear(xcoords,ycoords,zcoords)

            case ('1x1x1_unit','111u')
                call meshgen_1x1x1_unit_linear(xcoords,ycoords,zcoords)

            case ('3x3x3','333')
                call meshgen_3x3x3_linear(xcoords,ycoords,zcoords)

            case ('3x3x3_unit','333u')
                call meshgen_3x3x3_unit_linear(xcoords,ycoords,zcoords)

            case ('2x2x2','222')
                call meshgen_2x2x2_linear(xcoords,ycoords,zcoords)

            case ('2x2x1','221')
                call meshgen_2x2x1_linear(xcoords,ycoords,zcoords)

            case ('3x3x1','331')
                call meshgen_3x3x1_linear(xcoords,ycoords,zcoords)

            case ('4x1x1','411')
                call meshgen_4x1x1_linear(xcoords,ycoords,zcoords)

            case ('3x1x1','311')
                call meshgen_3x1x1_linear(xcoords,ycoords,zcoords)

            case ('2x1x1','211')
                call meshgen_2x1x1_linear(xcoords,ycoords,zcoords)

            case ('40x15x1')
                call meshgen_40x15x1_linear(xcoords,ycoords,zcoords)

            case ('15x15x1')
                call meshgen_15x15x1_linear(xcoords,ycoords,zcoords)

            case ('15x15x2')
                call meshgen_15x15x2_linear(xcoords,ycoords,zcoords)

            case ('15x15x3')
                call meshgen_15x15x3_linear(xcoords,ycoords,zcoords)


            case default
                call chidg_signal(FATAL,'String identifying mesh generation routine was not recognized')
        end select


        !
        ! Generate nodes, connectivity
        !
        mapping = 1
        idomain = 1
        nodes    = get_block_points_plot3d(xcoords,ycoords,zcoords)
        elements = get_block_elements_plot3d(xcoords,ycoords,zcoords,mapping,idomain)

        call connectivity%init(size(elements,1),size(nodes))
        do ielem = 1,size(elements,1)
            call connectivity%data(ielem)%init(1)
            call connectivity%data(ielem)%set_element_partition(IRANK)
            connectivity%data(ielem)%data = elements(ielem,:)
        end do

    end subroutine create_mesh
    !****************************************************************************











    !> Generate a set of points defining a 1x1x1 element mesh
    !!
    !!  @author Nathan A. Wukie
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   5/24/2016
    !!
    !!  @param[inout]   pts     points_t array of rank-3 that gets allocated, 
    !!                          filled, and returned
    !--------------------------------------------------------------------------------------
    subroutine meshgen_1x1x1_linear(xcoords,ycoords,zcoords)
        real(rk),   allocatable :: xcoords(:,:,:)
        real(rk),   allocatable :: ycoords(:,:,:)
        real(rk),   allocatable :: zcoords(:,:,:)

        integer(ik) :: ipt_xi, ipt_eta, ipt_zeta, ierr, npts_xi, &
                       npts_eta, npts_zeta
        real(rk)    :: x,y,z, dx, dy, dz

        ! elements (1x1x1) - linear

        npts_xi   = 2
        npts_eta  = 2
        npts_zeta = 2

        dx = 1._rk
        dy = 1._rk
        dz = 1._rk

        allocate(xcoords(npts_xi,npts_eta,npts_zeta), ycoords(npts_xi,npts_eta,npts_zeta), &
                 zcoords(npts_xi,npts_eta,npts_zeta), stat=ierr)
        if (ierr /= 0) call AllocationError

        z = ZERO
        do ipt_zeta = 1,npts_zeta
            y = ZERO
            do ipt_eta = 1,npts_eta
                x = ZERO
                do ipt_xi = 1,npts_xi


                    xcoords(ipt_xi,ipt_eta,ipt_zeta) = x
                    ycoords(ipt_xi,ipt_eta,ipt_zeta) = y
                    zcoords(ipt_xi,ipt_eta,ipt_zeta) = z


                    x = x + dx
                end do
                y = y + dy
            end do
            z = z + dz
        end do



    end subroutine meshgen_1x1x1_linear
    !**************************************************************************************














    !> Generate a set of points defining a 1x1x1 element mesh
    !!
    !!  @author Nathan A. Wukie
    !!  @param[inout]   pts     points_t array of rank-3 that gets allocated, 
    !!                          filled, and returned
    !--------------------------------------------------------------------------------------
    subroutine meshgen_1x1x1_unit_linear(xcoords,ycoords,zcoords)
        real(rk),   allocatable :: xcoords(:,:,:)
        real(rk),   allocatable :: ycoords(:,:,:)
        real(rk),   allocatable :: zcoords(:,:,:)

        integer(ik) :: ipt_xi, ipt_eta, ipt_zeta, ierr, npts_xi, &
                       npts_eta, npts_zeta
        real(rk)    :: x,y,z, dx, dy, dz

        ! elements (1x1x1) - linear

        npts_xi   = 2
        npts_eta  = 2
        npts_zeta = 2

        dx = 2._rk
        dy = 2._rk
        dz = 2._rk

        allocate(xcoords(npts_xi,npts_eta,npts_zeta), ycoords(npts_xi,npts_eta,npts_zeta), &
                 zcoords(npts_xi,npts_eta,npts_zeta), stat=ierr)
        if (ierr /= 0) call AllocationError

        z = -ONE
        do ipt_zeta = 1,npts_zeta
            y = -ONE
            do ipt_eta = 1,npts_eta
                x = -ONE
                do ipt_xi = 1,npts_xi

                    xcoords(ipt_xi,ipt_eta,ipt_zeta) = x
                    ycoords(ipt_xi,ipt_eta,ipt_zeta) = y
                    zcoords(ipt_xi,ipt_eta,ipt_zeta) = z

                    x = x + dx
                end do
                y = y + dy
            end do
            z = z + dz
        end do



    end subroutine meshgen_1x1x1_unit_linear
    !**************************************************************************************












    !> Generate a set of points defining a 2x2x2 element mesh
    !!
    !!  @author Nathan A. Wukie
    !!  @param[inout]   pts     points_t array of rank-3 that gets allocated, 
    !!                          filled, and returned
    !--------------------------------------------------------------------------------------
    subroutine meshgen_2x2x2_linear(xcoords,ycoords,zcoords)
        real(rk),   allocatable :: xcoords(:,:,:)
        real(rk),   allocatable :: ycoords(:,:,:)
        real(rk),   allocatable :: zcoords(:,:,:)

        integer(ik) :: ipt_xi, ipt_eta, ipt_zeta, ipt, ierr, &
                       npts_xi, npts_eta, npts_zeta

        integer(ik), parameter      :: npt = 27
        real(rk), dimension(npt)    :: x,y,z

        ! elements (2x2x2) - linear
        !
        !          *-------*-------*
        !         /       /       /|
        !        *-------*-------* |
        !       /       /       /| *
        !      *-------*-------* |/|
        !      |       |       | * |
        !      |       |       |/| *
        !      *-------*-------* |/
        !      |       |       | *
        !      |       |       |/
        !      *-------*-------*
        !
        !
        x = [ZERO, ONE, TWO, ZERO, ONE, TWO, ZERO, ONE, TWO, &
             ZERO, ONE, TWO, ZERO, ONE, TWO, ZERO, ONE, TWO, &
             ZERO, ONE, TWO, ZERO, ONE, TWO, ZERO, ONE, TWO]

        y = [ZERO, ZERO, ZERO, ONE, ONE, ONE, TWO, TWO, TWO, &
             ZERO, ZERO, ZERO, ONE, ONE, ONE, TWO, TWO, TWO, &
             ZERO, ZERO, ZERO, ONE, ONE, ONE, TWO, TWO, TWO]

        z = [ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, &
             ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, &
             TWO, TWO, TWO, TWO, TWO, TWO, TWO, TWO, TWO]

        npts_xi   = 3
        npts_eta  = 3
        npts_zeta = 3



        allocate(xcoords(npts_xi,npts_eta,npts_zeta), ycoords(npts_xi,npts_eta,npts_zeta), &
                 zcoords(npts_xi,npts_eta,npts_zeta), stat=ierr)
        if (ierr /= 0) call AllocationError


        ipt = 1
        do ipt_zeta = 1,npts_zeta
            do ipt_eta = 1,npts_eta
                do ipt_xi = 1,npts_xi
                    xcoords(ipt_xi,ipt_eta,ipt_zeta) = x(ipt)
                    ycoords(ipt_xi,ipt_eta,ipt_zeta) = y(ipt)
                    zcoords(ipt_xi,ipt_eta,ipt_zeta) = z(ipt)
                    ipt = ipt + 1
                end do
            end do
        end do



    end subroutine meshgen_2x2x2_linear
    !**************************************************************************************













    !> Generate a set of points defining a 2x2x1 element mesh
    !!
    !!  @author Nathan A. Wukie
    !!  @param[inout]   pts     points_t array of rank-3 that gets allocated, 
    !!                          filled, and returned
    !--------------------------------------------------------------------------------------
    subroutine meshgen_2x2x1_linear(xcoords,ycoords,zcoords)
        real(rk),   allocatable :: xcoords(:,:,:)
        real(rk),   allocatable :: ycoords(:,:,:)
        real(rk),   allocatable :: zcoords(:,:,:)

        integer(ik) :: ipt_xi, ipt_eta, ipt_zeta, ipt, ierr, &
                       npts_xi, npts_eta, npts_zeta

        integer(ik), parameter      :: npt = 18
        real(rk), dimension(npt)    :: x,y,z

        ! elements (2x2x1) - linear
        !
        !          *-------*
        !         /       /|
        !        *-------* |
        !       /       /| * 
        !      *-------* |/|
        !      |       | * |
        !      |       |/| * 
        !      *-------* |/
        !      |       | * 
        !      |       |/ 
        !      *-------*
        !
        !
        x = [ZERO, ONE, TWO, ZERO, ONE, TWO, ZERO, ONE, TWO, &
             ZERO, ONE, TWO, ZERO, ONE, TWO, ZERO, ONE, TWO]
             

        y = [ZERO, ZERO, ZERO, ONE, ONE, ONE, TWO, TWO, TWO, &
             ZERO, ZERO, ZERO, ONE, ONE, ONE, TWO, TWO, TWO]

        z = [ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, &
             ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE]

        npts_xi   = 3
        npts_eta  = 3
        npts_zeta = 2


        allocate(xcoords(npts_xi,npts_eta,npts_zeta), ycoords(npts_xi,npts_eta,npts_zeta), &
                 zcoords(npts_xi,npts_eta,npts_zeta), stat=ierr)
        if (ierr /= 0) call AllocationError

        ipt = 1
        do ipt_zeta = 1,npts_zeta
            do ipt_eta = 1,npts_eta
                do ipt_xi = 1,npts_xi
                    xcoords(ipt_xi,ipt_eta,ipt_zeta) = x(ipt)
                    ycoords(ipt_xi,ipt_eta,ipt_zeta) = y(ipt)
                    zcoords(ipt_xi,ipt_eta,ipt_zeta) = z(ipt)
                    ipt = ipt + 1
                end do
            end do
        end do



    end subroutine meshgen_2x2x1_linear
    !***************************************************************************************















    !> Generate a set of points defining a 3x3x3 element mesh
    !!
    !!  @author Nathan A. Wukie
    !!  @param[inout]   pts     points_t array of rank-3 that gets allocated, 
    !!                          filled, and returned
    !---------------------------------------------------------------------------------------
    subroutine meshgen_3x3x3_linear(xcoords,ycoords,zcoords)
        real(rk),   allocatable :: xcoords(:,:,:)
        real(rk),   allocatable :: ycoords(:,:,:)
        real(rk),   allocatable :: zcoords(:,:,:)

        integer(ik) :: ipt_xi, ipt_eta, ipt_zeta, ipt, ierr, &
                       npts_xi, npts_eta, npts_zeta

        integer(ik), parameter      :: npt = 64
        real(rk), dimension(npt)    :: x,y,z

        ! elements (3x3x3) - linear
        !
        !            *-------*-------*-------*
        !           /       /       /       /|
        !          *-------*-------*-------* |
        !         /       /       /       /| *
        !        *-------*-------*-------* |/|
        !       /       /       /       /| * |
        !      *-------*-------*-------* |/| *
        !      |       |       |       | * |/|
        !      |       |       |       |/| * |
        !      *-------*-------*-------* |/| *
        !      |       |       |       | * |/
        !      |       |       |       |/| *
        !      *-------*-------*-------* |/
        !      |       |       |       | *
        !      |       |       |       |/
        !      *-------*-------*-------*
        !
        !
        x = [ZERO, ONE, TWO, THREE, ZERO, ONE, TWO, THREE, ZERO, ONE, TWO, THREE, &
             ZERO, ONE, TWO, THREE, ZERO, ONE, TWO, THREE, ZERO, ONE, TWO, THREE, &
             ZERO, ONE, TWO, THREE, ZERO, ONE, TWO, THREE, ZERO, ONE, TWO, THREE, &
             ZERO, ONE, TWO, THREE, ZERO, ONE, TWO, THREE, ZERO, ONE, TWO, THREE, &
             ZERO, ONE, TWO, THREE, ZERO, ONE, TWO, THREE, ZERO, ONE, TWO, THREE, &
             ZERO, ONE, TWO, THREE]

        y = [ZERO, ZERO, ZERO, ZERO, ONE, ONE, ONE, ONE, TWO, TWO, TWO, TWO, &
             THREE, THREE, THREE, THREE, &
             ZERO, ZERO, ZERO, ZERO, ONE, ONE, ONE, ONE, TWO, TWO, TWO, TWO, &
             THREE, THREE, THREE, THREE, &
             ZERO, ZERO, ZERO, ZERO, ONE, ONE, ONE, ONE, TWO, TWO, TWO, TWO, &
             THREE, THREE, THREE, THREE, &
             ZERO, ZERO, ZERO, ZERO, ONE, ONE, ONE, ONE, TWO, TWO, TWO, TWO, &
             THREE, THREE, THREE, THREE]

        z = [ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, &
             ZERO, ZERO, ZERO, ZERO, ZERO, &
             ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, & 
             ONE, ONE, ONE, &
             TWO, TWO, TWO, TWO, TWO, TWO, TWO, TWO, TWO, TWO, TWO, TWO, TWO, &
             TWO, TWO, TWO, &
             THREE, THREE, THREE, THREE, THREE, THREE, THREE, THREE, THREE, THREE, &
             THREE, THREE, THREE, THREE, THREE, THREE]

        npts_xi   = 4
        npts_eta  = 4
        npts_zeta = 4


        allocate(xcoords(npts_xi,npts_eta,npts_zeta), ycoords(npts_xi,npts_eta,npts_zeta), &
                 zcoords(npts_xi,npts_eta,npts_zeta), stat=ierr)
        if (ierr /= 0) call AllocationError

        ipt = 1
        do ipt_zeta = 1,npts_zeta
            do ipt_eta = 1,npts_eta
                do ipt_xi = 1,npts_xi
                    xcoords(ipt_xi,ipt_eta,ipt_zeta) = x(ipt)
                    ycoords(ipt_xi,ipt_eta,ipt_zeta) = y(ipt)
                    zcoords(ipt_xi,ipt_eta,ipt_zeta) = z(ipt)
                    ipt = ipt + 1
                end do
            end do
        end do


    end subroutine meshgen_3x3x3_linear
    !***************************************************************************************
















    !> Generate a set of points defining a 3x3x3 unit-element mesh
    !!
    !!  @author Nathan A. Wukie
    !!  @param[inout]   pts     points_t array of rank-3 that gets allocated, 
    !!                          filled, and returned
    !---------------------------------------------------------------------
    subroutine meshgen_3x3x3_unit_linear(xcoords,ycoords,zcoords)
        real(rk),   allocatable :: xcoords(:,:,:)
        real(rk),   allocatable :: ycoords(:,:,:)
        real(rk),   allocatable :: zcoords(:,:,:)

        integer(ik) :: ipt_xi, ipt_eta, ipt_zeta, ipt, ierr, &
                       npts_xi, npts_eta, npts_zeta

        integer(ik), parameter      :: npt = 64
        real(rk), dimension(npt)    :: x,y,z

        ! elements (3x3x3) - unit linear
        !
        !            *-------*-------*-------*
        !           /       /       /       /|
        !          *-------*-------*-------* |
        !         /       /       /       /| *
        !        *-------*-------*-------* |/|
        !       /       /       /       /| * |
        !      *-------*-------*-------* |/| *
        !      |       |       |       | * |/|
        !      |       |       |       |/| * |
        !      *-------*-------*-------* |/| *
        !      |       |       |       | * |/
        !      |       |       |       |/| *
        !      *-------*-------*-------* |/
        !      |       |       |       | *
        !      |       |       |       |/
        !      *-------*-------*-------*
        !
        !
        x = [ZERO, TWO, FOUR, SIX, ZERO, TWO, FOUR, SIX, ZERO, TWO, FOUR, SIX, ZERO, TWO, FOUR, SIX, &
             ZERO, TWO, FOUR, SIX, ZERO, TWO, FOUR, SIX, ZERO, TWO, FOUR, SIX, ZERO, TWO, FOUR, SIX, &
             ZERO, TWO, FOUR, SIX, ZERO, TWO, FOUR, SIX, ZERO, TWO, FOUR, SIX, ZERO, TWO, FOUR, SIX, &
             ZERO, TWO, FOUR, SIX, ZERO, TWO, FOUR, SIX, ZERO, TWO, FOUR, SIX, ZERO, TWO, FOUR, SIX]

        y = [ZERO, ZERO, ZERO, ZERO, TWO, TWO, TWO, TWO, FOUR, FOUR, FOUR, FOUR, SIX, SIX, SIX, SIX, &
             ZERO, ZERO, ZERO, ZERO, TWO, TWO, TWO, TWO, FOUR, FOUR, FOUR, FOUR, SIX, SIX, SIX, SIX, &
             ZERO, ZERO, ZERO, ZERO, TWO, TWO, TWO, TWO, FOUR, FOUR, FOUR, FOUR, SIX, SIX, SIX, SIX, &
             ZERO, ZERO, ZERO, ZERO, TWO, TWO, TWO, TWO, FOUR, FOUR, FOUR, FOUR, SIX, SIX, SIX, SIX]

        z = [ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, &
             TWO, TWO, TWO, TWO, TWO, TWO, TWO, TWO, TWO, TWO, TWO, TWO, TWO, TWO, TWO, TWO, &
             FOUR, FOUR, FOUR, FOUR, FOUR, FOUR, FOUR, FOUR, FOUR, FOUR, FOUR, FOUR, FOUR, FOUR, FOUR, FOUR, &
             SIX, SIX, SIX, SIX, SIX, SIX, SIX, SIX, SIX, SIX, SIX, SIX, SIX, SIX, SIX, SIX]

        npts_xi   = 4
        npts_eta  = 4
        npts_zeta = 4


        allocate(xcoords(npts_xi,npts_eta,npts_zeta), ycoords(npts_xi,npts_eta,npts_zeta), &
                 zcoords(npts_xi,npts_eta,npts_zeta), stat=ierr)
        if (ierr /= 0) call AllocationError


        ipt = 1
        do ipt_zeta = 1,npts_zeta
            do ipt_eta = 1,npts_eta
                do ipt_xi = 1,npts_xi
                    xcoords(ipt_xi,ipt_eta,ipt_zeta) = x(ipt)
                    ycoords(ipt_xi,ipt_eta,ipt_zeta) = y(ipt)
                    zcoords(ipt_xi,ipt_eta,ipt_zeta) = z(ipt)
                    ipt = ipt + 1
                end do
            end do
        end do



    end subroutine meshgen_3x3x3_unit_linear
    !***************************************************************************************












    !> Generate a set of points defining a 3x3x1 element mesh
    !!
    !!  @author Nathan A. Wukie
    !!  @param[inout]   pts     points_t array of rank-3 that gets allocated, 
    !!                          filled, and returned
    !---------------------------------------------------------------------------------------
    subroutine meshgen_3x3x1_linear(xcoords,ycoords,zcoords)
        real(rk),   allocatable :: xcoords(:,:,:)
        real(rk),   allocatable :: ycoords(:,:,:)
        real(rk),   allocatable :: zcoords(:,:,:)

        integer(ik) :: ipt_xi, ipt_eta, ipt_zeta, ipt, ierr, &
                       npts_xi, npts_eta, npts_zeta

        integer(ik), parameter      :: npt = 32
        real(rk), dimension(npt)    :: x,y,z

        ! elements (3x3x1) - linear
        !
        !            *-------*
        !           /       /| 
        !          *-------* |
        !         /       /| *   
        !        *-------* |/|
        !       /       /| * | 
        !      *-------* |/| *
        !      |       | * |/| 
        !      |       |/| * |
        !      *-------* |/| *
        !      |       | * |/  
        !      |       |/| *  
        !      *-------* |/
        !      |       | *
        !      |       |/ 
        !      *-------*
        !
        !
        x = [ZERO, ONE, TWO, THREE, ZERO, ONE, TWO, THREE, ZERO, ONE, TWO, THREE, ZERO, ONE, TWO, THREE, &
             ZERO, ONE, TWO, THREE, ZERO, ONE, TWO, THREE, ZERO, ONE, TWO, THREE, ZERO, ONE, TWO, THREE]

        y = [ZERO, ZERO, ZERO, ZERO, ONE, ONE, ONE, ONE, TWO, TWO, TWO, TWO, THREE, THREE, THREE, THREE, &
             ZERO, ZERO, ZERO, ZERO, ONE, ONE, ONE, ONE, TWO, TWO, TWO, TWO, THREE, THREE, THREE, THREE]

        z = [ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, &
             ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE]

        npts_xi   = 4
        npts_eta  = 4
        npts_zeta = 2


        allocate(xcoords(npts_xi,npts_eta,npts_zeta), ycoords(npts_xi,npts_eta,npts_zeta), &
                 zcoords(npts_xi,npts_eta,npts_zeta), stat=ierr)
        if (ierr /= 0) call AllocationError

        ipt = 1
        do ipt_zeta = 1,npts_zeta
            do ipt_eta = 1,npts_eta
                do ipt_xi = 1,npts_xi
                    xcoords(ipt_xi,ipt_eta,ipt_zeta) = x(ipt)
                    ycoords(ipt_xi,ipt_eta,ipt_zeta) = y(ipt)
                    zcoords(ipt_xi,ipt_eta,ipt_zeta) = z(ipt)
                    ipt = ipt + 1
                end do
            end do
        end do


    end subroutine meshgen_3x3x1_linear
    !***************************************************************************************












    !> Generate a set of points defining a 4x1x1 element mesh
    !!
    !!  @author Nathan A. Wukie
    !!  @param[inout]   pts     points_t array of rank-3 that gets allocated, 
    !!                          filled, and returned
    !---------------------------------------------------------------------------------------
    subroutine meshgen_4x1x1_linear(xcoords,ycoords,zcoords)
        real(rk),   allocatable :: xcoords(:,:,:)
        real(rk),   allocatable :: ycoords(:,:,:)
        real(rk),   allocatable :: zcoords(:,:,:)

        integer(ik) :: ipt_xi, ipt_eta, ipt_zeta, ipt, ierr, &
                       npts_xi, npts_eta, npts_zeta

        integer(ik), parameter      :: npt = 20
        real(rk), dimension(npt)    :: x,y,z

        ! elements (4x1x1) - linear
        !
        !      *------*------*------*------*
        !      |      |      |      |      | 
        !      |      |      |      |      | 
        !      *------*------*------*------*
        !



        x = [ZERO, ONE, TWO, THREE, FOUR, &       
             ZERO, ONE, TWO, THREE, FOUR, &       
             ZERO, ONE, TWO, THREE, FOUR, &       
             ZERO, ONE, TWO, THREE, FOUR]
             

        y = [ZERO, ZERO, ZERO, ZERO, ZERO, &
             ONE, ONE, ONE, ONE, ONE, &
             ZERO, ZERO, ZERO, ZERO, ZERO, &
             ONE, ONE, ONE, ONE, ONE]

        z = [ZERO, ZERO, ZERO, ZERO, ZERO, &
             ZERO, ZERO, ZERO, ZERO, ZERO, &
             ONE, ONE, ONE, ONE, ONE, &
             ONE, ONE, ONE, ONE, ONE]

        npts_xi   = 5
        npts_eta  = 2
        npts_zeta = 2


        allocate(xcoords(npts_xi,npts_eta,npts_zeta), ycoords(npts_xi,npts_eta,npts_zeta), &
                 zcoords(npts_xi,npts_eta,npts_zeta), stat=ierr)
        if (ierr /= 0) call AllocationError


        ipt = 1
        do ipt_zeta = 1,npts_zeta
            do ipt_eta = 1,npts_eta
                do ipt_xi = 1,npts_xi
                    xcoords(ipt_xi,ipt_eta,ipt_zeta) = x(ipt)
                    ycoords(ipt_xi,ipt_eta,ipt_zeta) = y(ipt)
                    zcoords(ipt_xi,ipt_eta,ipt_zeta) = z(ipt)
                    ipt = ipt + 1
                end do
            end do
        end do



    end subroutine meshgen_4x1x1_linear
    !***************************************************************************************












    !> Generate a set of points defining a 2x1x1 element mesh
    !!
    !!  @author Nathan A. Wukie
    !!  @param[inout]   pts     points_t array of rank-3 that gets allocated, 
    !!                          filled, and returned
    !---------------------------------------------------------------------------------------
    subroutine meshgen_2x1x1_linear(xcoords,ycoords,zcoords)
        real(rk),   allocatable :: xcoords(:,:,:)
        real(rk),   allocatable :: ycoords(:,:,:)
        real(rk),   allocatable :: zcoords(:,:,:)

        integer(ik) :: ipt_xi, ipt_eta, ipt_zeta, ipt, ierr, &
                       npts_xi, npts_eta, npts_zeta

        integer(ik), parameter      :: npt = 12
        real(rk), dimension(npt)    :: x,y,z

        ! elements (3x1x1) - linear
        !
        !      *------*------*
        !      |      |      | 
        !      |      |      | 
        !      *------*------*
        !



        x = [ZERO, ONE, TWO, &       
             ZERO, ONE, TWO, &       
             ZERO, ONE, TWO, &       
             ZERO, ONE, TWO]
             

        y = [ZERO, ZERO, ZERO, &
             ONE, ONE, ONE, &
             ZERO, ZERO, ZERO, &
             ONE, ONE, ONE]

        z = [ZERO, ZERO, ZERO, &
             ZERO, ZERO, ZERO, &
             ONE, ONE, ONE, &
             ONE, ONE, ONE]


        npts_xi   = 3
        npts_eta  = 2
        npts_zeta = 2

        allocate(xcoords(npts_xi,npts_eta,npts_zeta), ycoords(npts_xi,npts_eta,npts_zeta), &
                 zcoords(npts_xi,npts_eta,npts_zeta), stat=ierr)
        if (ierr /= 0) call AllocationError


        ipt = 1
        do ipt_zeta = 1,npts_zeta
            do ipt_eta = 1,npts_eta
                do ipt_xi = 1,npts_xi
                    xcoords(ipt_xi,ipt_eta,ipt_zeta) = x(ipt)
                    ycoords(ipt_xi,ipt_eta,ipt_zeta) = y(ipt)
                    zcoords(ipt_xi,ipt_eta,ipt_zeta) = z(ipt)
                    ipt = ipt + 1
                end do
            end do
        end do


    end subroutine meshgen_2x1x1_linear
    !**************************************************************************************












    !> Generate a set of points defining a 3x1x1 element mesh
    !!
    !!  @author Nathan A. Wukie
    !!  @param[inout]   pts     points_t array of rank-3 that gets allocated, 
    !!                          filled, and returned
    !---------------------------------------------------------------------------------------
    subroutine meshgen_3x1x1_linear(xcoords,ycoords,zcoords)
        real(rk),   allocatable :: xcoords(:,:,:)
        real(rk),   allocatable :: ycoords(:,:,:)
        real(rk),   allocatable :: zcoords(:,:,:)

        integer(ik) :: ipt_xi, ipt_eta, ipt_zeta, ipt, ierr, &
                       npts_xi, npts_eta, npts_zeta

        integer(ik), parameter      :: npt = 16
        real(rk), dimension(npt)    :: x,y,z

        ! elements (3x1x1) - linear
        !
        !      *------*------*------*
        !      |      |      |      | 
        !      |      |      |      | 
        !      *------*------*------*
        !



        x = [ZERO, ONE, TWO, THREE, &       
             ZERO, ONE, TWO, THREE, &       
             ZERO, ONE, TWO, THREE, &       
             ZERO, ONE, TWO, THREE]
             

        y = [ZERO, ZERO, ZERO, ZERO, &
             ONE, ONE, ONE, ONE, &
             ZERO, ZERO, ZERO, ZERO, &
             ONE, ONE, ONE, ONE]

        z = [ZERO, ZERO, ZERO, ZERO, &
             ZERO, ZERO, ZERO, ZERO, &
             ONE, ONE, ONE, ONE, &
             ONE, ONE, ONE, ONE]


        npts_xi   = 4
        npts_eta  = 2
        npts_zeta = 2


        allocate(xcoords(npts_xi,npts_eta,npts_zeta), ycoords(npts_xi,npts_eta,npts_zeta), &
                 zcoords(npts_xi,npts_eta,npts_zeta), stat=ierr)
        if (ierr /= 0) call AllocationError

        ipt = 1
        do ipt_zeta = 1,npts_zeta
            do ipt_eta = 1,npts_eta
                do ipt_xi = 1,npts_xi
                    xcoords(ipt_xi,ipt_eta,ipt_zeta) = x(ipt)
                    ycoords(ipt_xi,ipt_eta,ipt_zeta) = y(ipt)
                    zcoords(ipt_xi,ipt_eta,ipt_zeta) = z(ipt)
                    ipt = ipt + 1
                end do
            end do
        end do


    end subroutine meshgen_3x1x1_linear
    !**************************************************************************************



















    !> Generate a set of points defining a 4x1x1 element mesh
    !!
    !!  @author Nathan A. Wukie
    !!  @param[inout]   pts     points_t array of rank-3 that gets allocated, 
    !!                          filled, and returned
    !--------------------------------------------------------------------------------------
    subroutine meshgen_40x15x1_linear(xcoords,ycoords,zcoords)
        real(rk),   allocatable :: xcoords(:,:,:)
        real(rk),   allocatable :: ycoords(:,:,:)
        real(rk),   allocatable :: zcoords(:,:,:)

        integer(ik) :: ipt_xi, ipt_eta, ipt_zeta, ipt, ierr, &
                       npts_xi, npts_eta, npts_zeta

        real(rk)    :: x,y,z, dx, dy, dz

        ! elements (40x15x1) - linear

        npts_xi   = 41
        npts_eta  = 16
        npts_zeta = 2

        dx = 0.5_rk
        dy = 0.5_rk
        dz = 1._rk

        allocate(xcoords(npts_xi,npts_eta,npts_zeta), ycoords(npts_xi,npts_eta,npts_zeta), &
                 zcoords(npts_xi,npts_eta,npts_zeta), stat=ierr)
        if (ierr /= 0) call AllocationError

        z = ZERO
        do ipt_zeta = 1,npts_zeta
            y = ZERO
            do ipt_eta = 1,npts_eta
                x = ZERO
                do ipt_xi = 1,npts_xi

                    xcoords(ipt_xi,ipt_eta,ipt_zeta) = x
                    ycoords(ipt_xi,ipt_eta,ipt_zeta) = y
                    zcoords(ipt_xi,ipt_eta,ipt_zeta) = z

                    x = x + dx
                end do
                y = y + dy
            end do
            z = z + dz
        end do


    end subroutine meshgen_40x15x1_linear
    !**************************************************************************************













    !> Generate a set of points defining a 4x1x1 element mesh
    !!
    !!  @author Nathan A. Wukie
    !!  @param[inout]   pts     points_t array of rank-3 that gets allocated, 
    !!                          filled, and returned
    !---------------------------------------------------------------------------------------
    subroutine meshgen_15x15x1_linear(xcoords,ycoords,zcoords)
        real(rk),   allocatable :: xcoords(:,:,:)
        real(rk),   allocatable :: ycoords(:,:,:)
        real(rk),   allocatable :: zcoords(:,:,:)

        integer(ik) :: ipt_xi, ipt_eta, ipt_zeta, ipt, ierr, &
                       npts_xi, npts_eta, npts_zeta

        real(rk)    :: x,y,z,dx,dy,dz

        ! elements (15x15x1) - linear

        npts_xi   = 16
        npts_eta  = 16
        npts_zeta = 2

        dx = 0.5_rk
        dy = 0.5_rk
        dz = 1.0_rk


        allocate(xcoords(npts_xi,npts_eta,npts_zeta), ycoords(npts_xi,npts_eta,npts_zeta), &
                 zcoords(npts_xi,npts_eta,npts_zeta), stat=ierr)
        if (ierr /= 0) call AllocationError

        z = ZERO
        do ipt_zeta = 1,npts_zeta
            y = ZERO
            do ipt_eta = 1,npts_eta
                x = ZERO
                do ipt_xi = 1,npts_xi

                    xcoords(ipt_xi,ipt_eta,ipt_zeta) = x
                    ycoords(ipt_xi,ipt_eta,ipt_zeta) = y
                    zcoords(ipt_xi,ipt_eta,ipt_zeta) = z
                    ipt = ipt + 1

                    x = x + dx
                end do
                y = y + dy
            end do
            z = z + dz
        end do


    end subroutine meshgen_15x15x1_linear
    !***************************************************************************************












    !> Generate a set of points defining a 15x15x2 element mesh
    !!
    !!  @author Nathan A. Wukie
    !!  @param[inout]   pts     points_t array of rank-3 that gets allocated, 
    !!                          filled, and returned
    !---------------------------------------------------------------------------------------
    subroutine meshgen_15x15x2_linear(xcoords,ycoords,zcoords)
        real(rk),   allocatable :: xcoords(:,:,:)
        real(rk),   allocatable :: ycoords(:,:,:)
        real(rk),   allocatable :: zcoords(:,:,:)

        integer(ik) :: ipt_xi, ipt_eta, ipt_zeta, ipt, ierr, &
                       npts_xi,npts_eta,npts_zeta
        real(rk)    :: x,y,z,dx,dy,dz

        ! elements (15x15x2) - linear

        npts_xi   = 16
        npts_eta  = 16
        npts_zeta = 3

        dx = 0.5_rk
        dy = 0.5_rk
        dz = 0.5_rk

        allocate(xcoords(npts_xi,npts_eta,npts_zeta), ycoords(npts_xi,npts_eta,npts_zeta), &
                 zcoords(npts_xi,npts_eta,npts_zeta), stat=ierr)
        if (ierr /= 0) call AllocationError

        z = ZERO
        do ipt_zeta = 1,npts_zeta
            y = ZERO
            do ipt_eta = 1,npts_eta
                x = ZERO
                do ipt_xi = 1,npts_xi

                    xcoords(ipt_xi,ipt_eta,ipt_zeta) = x
                    ycoords(ipt_xi,ipt_eta,ipt_zeta) = y
                    zcoords(ipt_xi,ipt_eta,ipt_zeta) = z

                    x = x + dx
                end do
                y = y + dy
            end do
            z = z + dz
        end do




    end subroutine meshgen_15x15x2_linear
    !***************************************************************************************












    !> Generate a set of points defining a 15x15x3 element mesh
    !!
    !!  @author Nathan A. Wukie
    !!  @param[inout]   pts     points_t array of rank-3 that gets allocated, 
    !!                          filled, and returned
    !---------------------------------------------------------------------------------------
    subroutine meshgen_15x15x3_linear(xcoords,ycoords,zcoords)
        real(rk),   allocatable :: xcoords(:,:,:)
        real(rk),   allocatable :: ycoords(:,:,:)
        real(rk),   allocatable :: zcoords(:,:,:)

        integer(ik) :: ipt_xi, ipt_eta, ipt_zeta, ipt, ierr, &
                       npts_xi, npts_eta, npts_zeta

        real(rk)    :: x,y,z, dx, dy, dz

        ! elements (15x15x3) - linear

        npts_xi   = 16
        npts_eta  = 16
        npts_zeta = 4

        dx = 0.5_rk
        dy = 0.5_rk
        dz = 0.5_rk

        allocate(xcoords(npts_xi,npts_eta,npts_zeta), ycoords(npts_xi,npts_eta,npts_zeta), &
                 zcoords(npts_xi,npts_eta,npts_zeta), stat=ierr)
        if (ierr /= 0) call AllocationError

        z = ZERO
        do ipt_zeta = 1,npts_zeta
            y = ZERO
            do ipt_eta = 1,npts_eta
                x = ZERO
                do ipt_xi = 1,npts_xi

                    xcoords(ipt_xi,ipt_eta,ipt_zeta) = x
                    ycoords(ipt_xi,ipt_eta,ipt_zeta) = y
                    zcoords(ipt_xi,ipt_eta,ipt_zeta) = z

                    x = x + dx
                end do
                y = y + dy
            end do
            z = z + dz
        end do



    end subroutine meshgen_15x15x3_linear
    !**************************************************************************************







    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   10/17/2016
    !!
    !!
    !!
    !---------------------------------------------------------------------------------------
    subroutine create_mesh_file_D2E8M1_overlapping_matching(filename)
        character(*),   intent(in)  :: filename

        class(bc_state_t),  allocatable                 :: bc_state
        character(8)                                    :: faces(5)
        integer(HID_T)                                  :: file_id, dom1_id, dom2_id, bcface_id
        integer(ik)                                     :: spacedim, mapping, bcface, ierr
        type(point_t),  allocatable                     :: nodes1(:), nodes2(:)
        integer(ik),    allocatable                     :: elements1(:,:), elements2(:,:) 
        integer(ik),    allocatable                     :: faces1(:,:), faces2(:,:)
        real(rk),       allocatable, dimension(:,:,:)   :: xcoords1, xcoords2, ycoords, zcoords
        real(rk)                                        :: xmax


        ! Create/initialize file
        file_id = initialize_file_hdf(filename)
        

        ! Generate coordinates for first block
        call meshgen_2x2x2_linear(xcoords1,ycoords,zcoords)


        !
        ! Create second block by copying and translating first block.
        !
        xmax = maxval(xcoords1)
        xcoords2 = xcoords1 + 0.95*xmax


        !
        ! Get nodes/elements
        !
        mapping = 1
        nodes1    = get_block_points_plot3d(xcoords1,ycoords,zcoords)
        nodes2    = get_block_points_plot3d(xcoords2,ycoords,zcoords)
        elements1 = get_block_elements_plot3d(xcoords1,ycoords,zcoords,mapping,idomain=1)
        elements2 = get_block_elements_plot3d(xcoords2,ycoords,zcoords,mapping,idomain=2)


        !
        ! Add domains
        !
        spacedim = 3
        call add_domain_hdf(file_id,"01",nodes1,elements1,"Scalar Advection",spacedim)
        call add_domain_hdf(file_id,"02",nodes2,elements2,"Scalar Advection",spacedim)


        !
        ! Set boundary conditions patch connectivities
        !
        dom1_id = open_domain_hdf(file_id,"01")
        dom2_id = open_domain_hdf(file_id,"02")

        do bcface = 1,6
            ! Get face node indices for boundary 'bcface'
            faces1 = get_block_boundary_faces_plot3d(xcoords1,ycoords,zcoords,mapping,bcface)
            faces2 = get_block_boundary_faces_plot3d(xcoords2,ycoords,zcoords,mapping,bcface)

            ! Set bc patch face indices
            call set_bc_patch_hdf(dom1_id,faces1,bcface)
            call set_bc_patch_hdf(dom2_id,faces2,bcface)
        end do !bcface


        !
        ! Create bc_state, "Scalar Extrapolate"
        !
        call create_bc("Scalar Extrapolate", bc_state)


        !
        ! Set boundary condition states for Domain 1: Leave XI_MAX empty
        !
        faces = ["XI_MIN  ", "ETA_MIN ", "ETA_MAX ", "ZETA_MIN", "ZETA_MAX"]
        do bcface = 1,size(faces)
            call h5gopen_f(dom1_id,"BoundaryConditions/"//trim(adjustl(faces(bcface))),bcface_id,ierr)
            call add_bc_state_hdf(bcface_id,bc_state)
            call h5gclose_f(bcface_id,ierr)
        end do

        !
        ! Set boundary condition states for Domain 1: Leave XI_MIN empty
        !
        faces = ["XI_MAX  ", "ETA_MIN ", "ETA_MAX ", "ZETA_MIN", "ZETA_MAX"]
        do bcface = 1,size(faces)
            call h5gopen_f(dom2_id,"BoundaryConditions/"//trim(adjustl(faces(bcface))),bcface_id,ierr)
            call add_bc_state_hdf(bcface_id,bc_state)
            call h5gclose_f(bcface_id,ierr)
        end do


        ! Set 'Contains Grid'
        call set_contains_grid_hdf(file_id,"True")

        ! Close file
        call close_domain_hdf(dom1_id)
        call close_domain_hdf(dom2_id)
        call close_file_hdf(file_id)
        call close_hdf()

    end subroutine create_mesh_file_D2E8M1_overlapping_matching
    !*************************************************************************************










end module mod_test_utilities
