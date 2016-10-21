module mod_gridgen_blocks
#include <messenger.h>
    use mod_kinds,              only: rk, ik
    use mod_constants,          only: ZERO, ONE, TWO, THREE, FOUR, SIX, PI
    use mod_bc,                 only: create_bc
    use mod_plot3d_utilities,   only: get_block_points_plot3d, &
                                      get_block_elements_plot3d, &
                                      get_block_boundary_faces_plot3d
    use mod_hdf_utilities,      only: initialize_file_hdf, add_domain_hdf, &
                                      open_file_hdf, close_file_hdf, &
                                      open_domain_hdf, close_domain_hdf, &
                                      set_bc_patch_hdf, add_bc_state_hdf, &
                                      set_contains_grid_hdf, close_hdf, open_hdf

    use type_point,             only: point_t
    use type_bc_state,          only: bc_state_t
    use type_bc_state_wrapper,  only: bc_state_wrapper_t
    use hdf5
    implicit none





contains

    !-------------------------------------------------------------------------------------
    !!
    !!
    !!  Create File: Grid + BC's
    !!  -----------------------------
    !!  create_mesh_file__singleblock
    !!  create_mesh_file__multiblock
    !!  create_mesh_file__D2E8M1
    !!
    !!
    !!  Generate grid: point arrays
    !!  --------------------------
    !!  meshgen_1x1x1_linear
    !!  meshgen_1x1x1_unit_linear
    !!  meshgen_2x2x2_linear
    !!  meshgen_2x2x1_linear
    !!  meshgen_3x3x3_linear
    !!  meshgen_3x3x3_unit_linear
    !!  meshgen_3x3x1_linear
    !!  meshgen_4x1x1_linear
    !!  meshgen_4x2x2_linear
    !!  meshgen_3x1x1_linear
    !!  meshgen_2x1x1_linear
    !!  meshgen_40x15x1_linear
    !!  meshgen_15x15x1_linear
    !!  meshgen_15x15x2_linear
    !!  meshgen_15x15x3_linear
    !!
    !!
    !!
    !**************************************************************************************



    !>  Write a ChiDG-formatted grid file consisting of:
    !!
    !!      - one block domain, D1
    !!      - Linear element mapping, M1
    !!      - boundary conditions initialized to Scalar Extrapolate.
    !!
    !!  Particular block can be specified by the input string 'grid':
    !!      'grid' = "D1 E1 M1"
    !!      'grid' = "D1 E8 M1"
    !!      'grid' = "D1 E16 M1"
    !!      'grid' = "D1 E27 M1"
    !!
    !!  @author Nathan A. Wukie
    !!  @date   10/17/2016
    !!
    !!
    !!
    !---------------------------------------------------------------------------------------
    subroutine create_mesh_file__singleblock(filename,grid,equation_set1,bc_states1,nelem_xi,nelem_eta,nelem_zeta,clusterx)
        character(*),               intent(in)              :: filename
        character(*),               intent(in)              :: grid
        character(*),               intent(in), optional    :: equation_set1
        type(bc_state_wrapper_t),   intent(in), optional    :: bc_states1(:)
        integer(ik),                intent(in), optional    :: nelem_xi
        integer(ik),                intent(in), optional    :: nelem_eta
        integer(ik),                intent(in), optional    :: nelem_zeta
        integer(ik),                intent(in), optional    :: clusterx

        character(:),       allocatable                 :: user_msg
        class(bc_state_t),  allocatable                 :: bc_state
        character(8)                                    :: face_strings(6)
        integer(HID_T)                                  :: file_id, dom_id, bcface_id
        integer(ik)                                     :: spacedim, mapping, bcface, ierr
        type(point_t),  allocatable                     :: nodes(:)
        integer(ik),    allocatable                     :: elements(:,:) 
        integer(ik),    allocatable                     :: faces(:,:)
        real(rk),       allocatable, dimension(:,:,:)   :: xcoords, ycoords, zcoords


        ! Create/initialize file
        file_id = initialize_file_hdf(filename)
        


        ! Generate coordinates for first block
        select case (trim(grid))
            case("D1 E1 M1")
                call meshgen_1x1x1_linear(xcoords,ycoords,zcoords)
            case("D1 E4 M1")
                call meshgen_4x1x1_linear(xcoords,ycoords,zcoords)
            case("D1 E16 M1")
                call meshgen_4x2x2_linear(xcoords,ycoords,zcoords)
            case("D1 E27 M1")
                call meshgen_3x3x3_linear(xcoords,ycoords,zcoords)
            case("D1 NxNxN")
                if ( present(nelem_xi) .and. present(nelem_eta) .and. present(nelem_zeta) ) then
                    call meshgen_NxNxN_linear(nelem_xi,nelem_eta,nelem_zeta,xcoords,ycoords,zcoords,clusterx)
                else
                    user_msg = "create_mesh_file__singleblock: For 'D1 NxNxN', need to specify &
                                the optional inputs 'nelem_xi', 'nelem_eta', 'nelem_zeta'."
                    call chidg_signal(FATAL,user_msg)
                end if
            case default
                call chidg_signal(FATAL,"create_mesh_file__singleblock: Invalid string to select grid block")
        end select



        !
        ! Get nodes/elements
        !
        mapping = 1
        nodes    = get_block_points_plot3d(xcoords,ycoords,zcoords)
        elements = get_block_elements_plot3d(xcoords,ycoords,zcoords,mapping,idomain=1)


        !
        ! Add domains
        !
        spacedim = 3

        if ( present(equation_set1) ) then
            call add_domain_hdf(file_id,"01",nodes,elements,equation_set1,spacedim)
        else
            call add_domain_hdf(file_id,"01",nodes,elements,"Scalar Advection",spacedim)
        end if


        !
        ! Set boundary conditions patch connectivities
        !
        dom_id = open_domain_hdf(file_id,"01")

        do bcface = 1,6
            ! Get face node indices for boundary 'bcface'
            faces = get_block_boundary_faces_plot3d(xcoords,ycoords,zcoords,mapping,bcface)

            ! Set bc patch face indices
            call set_bc_patch_hdf(dom_id,faces,bcface)
        end do !bcface


        !
        ! Create bc_state, "Scalar Extrapolate"
        !
        call create_bc("Scalar Extrapolate", bc_state)


        !
        ! Set boundary condition states for Domain 1: Leave XI_MAX empty
        !
        face_strings = ["XI_MIN  ","XI_MAX  ", "ETA_MIN ", "ETA_MAX ", "ZETA_MIN", "ZETA_MAX"]
        do bcface = 1,size(face_strings)
            call h5gopen_f(dom_id,"BoundaryConditions/"//trim(adjustl(face_strings(bcface))),bcface_id,ierr)

            if (present(bc_states1)) then
                call add_bc_state_hdf(bcface_id,bc_states1(bcface)%state)
            else
                call add_bc_state_hdf(bcface_id,bc_state)
            end if

            call h5gclose_f(bcface_id,ierr)
        end do



        ! Set 'Contains Grid'
        call set_contains_grid_hdf(file_id,"True")

        ! Close file
        call close_domain_hdf(dom_id)
        call close_file_hdf(file_id)
        call close_hdf()

    end subroutine create_mesh_file__singleblock
    !*************************************************************************************










    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   10/19/2016
    !!
    !!
    !!
    !-------------------------------------------------------------------------------------
    subroutine create_mesh_file__multiblock(filename,block1,block2,equation_set1,equation_set2,bc_states1,bc_states2)
        character(*),               intent(in)              :: filename
        character(*),               intent(in)              :: block1
        character(*),               intent(in)              :: block2
        character(*),               intent(in), optional    :: equation_set1
        character(*),               intent(in), optional    :: equation_set2
        type(bc_state_wrapper_t),   intent(in), optional    :: bc_states1(:)
        type(bc_state_wrapper_t),   intent(in), optional    :: bc_states2(:)

        class(bc_state_t),  allocatable                 :: bc_state
        character(8)                                    :: faces(6)
        integer(HID_T)                                  :: file_id, dom1_id, dom2_id, bcface1_id, bcface2_id
        integer(ik)                                     :: spacedim, mapping, bcface, ierr
        type(point_t),  allocatable                     :: nodes1(:), nodes2(:)
        integer(ik),    allocatable                     :: elements1(:,:), elements2(:,:) 
        integer(ik),    allocatable                     :: faces1(:,:), faces2(:,:)
        real(rk),       allocatable, dimension(:,:,:)   :: xcoords1, ycoords1, zcoords1, &
                                                           xcoords2, ycoords2, zcoords2
        real(rk)                                        :: xmax_block1,xmin_block2


        ! Create/initialize file
        file_id = initialize_file_hdf(filename)
        

        select case (trim(block1))
            case("D1 E1 M1")
                call meshgen_1x1x1_linear(xcoords1,ycoords1,zcoords1)
            case("D1 E2 M1")
                call meshgen_2x1x1_linear(xcoords1,ycoords1,zcoords1)
            case("D1 E27 M1")
                call meshgen_3x3x3_linear(xcoords1,ycoords1,zcoords1)
            case default
                call chidg_signal(FATAL,"create_mesh_file__multiblock: Invalid block1 string")
        end select




        select case (trim(block2))
            case("D1 E1 M1")
                !call meshgen_1x1x1_linear(xcoords2,ycoords2,zcoords2)
                call meshgen_NxNxN_linear(1,1,1,xcoords2,ycoords2,zcoords2)
            case("D1 E2 M1")
                !call meshgen_2x1x1_linear(xcoords2,ycoords2,zcoords2)
                call meshgen_NxNxN_linear(2,1,1,xcoords2,ycoords2,zcoords2)
            case("D1 E27 M1")
                !call meshgen_3x3x3_linear(xcoords2,ycoords2,zcoords2)
                call meshgen_NxNxN_linear(3,3,3,xcoords2,ycoords2,zcoords2)
            case default
                call chidg_signal(FATAL,"create_mesh_file__multiblock: Invalid block2 string")
        end select


        !
        ! Translate block2 to end of block1
        !
        xmax_block1 = maxval(xcoords1)
        xmin_block2 = minval(xcoords2)
        xcoords2 = xcoords2 + (xmax_block1-xmin_block2)



        !
        ! Get nodes/elements
        !
        mapping = 1
        nodes1    = get_block_points_plot3d(xcoords1,ycoords1,zcoords1)
        nodes2    = get_block_points_plot3d(xcoords2,ycoords2,zcoords2)
        elements1 = get_block_elements_plot3d(xcoords1,ycoords1,zcoords1,mapping,idomain=1)
        elements2 = get_block_elements_plot3d(xcoords2,ycoords2,zcoords2,mapping,idomain=2)


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
            faces1 = get_block_boundary_faces_plot3d(xcoords1,ycoords1,zcoords1,mapping,bcface)
            faces2 = get_block_boundary_faces_plot3d(xcoords2,ycoords2,zcoords2,mapping,bcface)

            ! Set bc patch face indices
            call set_bc_patch_hdf(dom1_id,faces1,bcface)
            call set_bc_patch_hdf(dom2_id,faces2,bcface)
        end do !bcface


        !
        ! Create bc_state, "Scalar Extrapolate"
        !
        call create_bc("Scalar Extrapolate", bc_state)


        !
        ! Set boundary condition states
        !
        faces = ["XI_MIN  ","XI_MAX  ", "ETA_MIN ", "ETA_MAX ", "ZETA_MIN", "ZETA_MAX"]
        do bcface = 1,size(faces)
            call h5gopen_f(dom1_id,"BoundaryConditions/"//trim(adjustl(faces(bcface))),bcface1_id,ierr)
            call h5gopen_f(dom2_id,"BoundaryConditions/"//trim(adjustl(faces(bcface))),bcface2_id,ierr)

            call add_bc_state_hdf(bcface1_id,bc_state)
            call add_bc_state_hdf(bcface2_id,bc_state)

            call h5gclose_f(bcface1_id,ierr)
            call h5gclose_f(bcface2_id,ierr)
        end do



        ! Set 'Contains Grid'
        call set_contains_grid_hdf(file_id,"True")

        ! Close file
        call close_domain_hdf(dom1_id)
        call close_domain_hdf(dom2_id)
        call close_file_hdf(file_id)
        call close_hdf()

    end subroutine create_mesh_file__multiblock
    !*************************************************************************************

























    !>  Write a ChiDG-formatted grid file consisting of:
    !!
    !!      - two block domains
    !!      - overlapping slightly
    !!      - boundary conditions initialized to Scalar Extrapolate. Interior boundaries
    !!        not set so they are detected as Chimera.
    !!
    !!  Incoming Parameter, 'matching' specifies if the elements should overlap with
    !!  a single, or potentially multiple elements
    !!
    !!     Block 1           Block 2 : matching=True        Block 2 : matching=False
    !!  .-----.-----.             .-----.-----.                  .-----.-----.
    !!  |     |     |             |     |     |                  |     |     |
    !!  |     |     |             |     |     |                  |     |     |
    !!  .-----.-----.             .-----.-----.                  |     |     |
    !!  |     |     |             |     |     |                  .-----.-----.
    !!  |     |     |             |     |     |                  |     |     |
    !!  .-----.-----.             .-----.-----.                  .-----.-----.
    !!
    !!  Abutting
    !!
    !!       abutting = .true.        abutting = .false.
    !!           ----.----               ----.-.----
    !!               |                       | | 
    !!               |                       | | 
    !!           ----.----               ----.-.----
    !!               |                       | |   
    !!               |                       | | 
    !!           ----.----               ----.-.----
    !!
    !!  Overlap
    !!
    !!       matching = .true.        matching = .false.
    !!          ----.-.----              ----.-.----
    !!              | |                      | | 
    !!              | |                      | | 
    !!          ----.-.----              ----|-.   
    !!              | |                      .-|----
    !!              | |                      | | 
    !!          ----.-.----              ----.-.----
    !!
    !!
    !!  @author Nathan A. Wukie
    !!  @date   10/17/2016
    !!
    !!
    !!
    !---------------------------------------------------------------------------------------
    subroutine create_mesh_file__D2E8M1(filename,abutting,matching,equation_set1,equation_set2,bc_states1,bc_states2)
        character(*),               intent(in)              :: filename
        logical,                    intent(in)              :: abutting
        logical,                    intent(in)              :: matching
        character(*),               intent(in), optional    :: equation_set1
        character(*),               intent(in), optional    :: equation_set2
        type(bc_state_wrapper_t),   intent(in), optional    :: bc_states1(:)
        type(bc_state_wrapper_t),   intent(in), optional    :: bc_states2(:)

        class(bc_state_t),  allocatable                 :: bc_state
        character(8)                                    :: faces(5)
        integer(HID_T)                                  :: file_id, dom1_id, dom2_id, bcface_id
        integer(ik)                                     :: spacedim, mapping, bcface, ierr
        type(point_t),  allocatable                     :: nodes1(:), nodes2(:)
        integer(ik),    allocatable                     :: elements1(:,:), elements2(:,:) 
        integer(ik),    allocatable                     :: faces1(:,:), faces2(:,:)
        real(rk),       allocatable, dimension(:,:,:)   :: xcoords1, xcoords2, ycoords1, ycoords2, zcoords
        real(rk)                                        :: xmax,ymax


        ! Create/initialize file
        file_id = initialize_file_hdf(filename)
        

        ! Generate coordinates for first block
        call meshgen_2x2x2_linear(xcoords1,ycoords1,zcoords)


        !
        ! Create second block by copying and translating first block.
        !
        xmax = maxval(xcoords1)
        ymax = maxval(ycoords1)



        !
        ! If abutting=.true., Create block2 by copying block1 and translating
        ! it by xmax of block1.
        !
        ! If abutting=.false., only translate by a fraction of xmax so there is overlap
        !
        if (abutting) then
            xcoords2 = xcoords1 + xmax
        else
            xcoords2 = xcoords1 + 0.90_rk*xmax
        end if
        ycoords2 = ycoords1


        !
        ! If matching=false, shift center plane of points so that the overlapping
        ! faces between blocks to not match exactly, and might have multiple Chimera donors
        !
        if (.not. matching) then
            ycoords2(:,2,:) = ycoords2(:,2,:) - 0.2*ymax
        end if


        !
        ! Get nodes/elements
        !
        mapping = 1
        nodes1    = get_block_points_plot3d(xcoords1,ycoords1,zcoords)
        nodes2    = get_block_points_plot3d(xcoords2,ycoords2,zcoords)
        elements1 = get_block_elements_plot3d(xcoords1,ycoords1,zcoords,mapping,idomain=1)
        elements2 = get_block_elements_plot3d(xcoords2,ycoords2,zcoords,mapping,idomain=2)


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
            faces1 = get_block_boundary_faces_plot3d(xcoords1,ycoords1,zcoords,mapping,bcface)
            faces2 = get_block_boundary_faces_plot3d(xcoords2,ycoords2,zcoords,mapping,bcface)

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
        ! Set boundary condition states for Domain 2: Leave XI_MIN empty
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

    end subroutine create_mesh_file__D2E8M1
    !*************************************************************************************










    !>  Generate a set of points defining a:
    !!      - nelem_xi by nelem_eta by nelem_zeta element, single-block, mesh
    !!
    !!  @author Nathan A. Wukie
    !!  @date   10/20/2016
    !!
    !!  @param[inout]   pts     points_t array of rank-3 that gets allocated, 
    !!                          filled, and returned
    !!
    !--------------------------------------------------------------------------------------
    subroutine meshgen_NxNxN_linear(nelem_xi,nelem_eta,nelem_zeta,xcoords,ycoords,zcoords,clusterx)
        integer(ik)             :: nelem_xi
        integer(ik)             :: nelem_eta
        integer(ik)             :: nelem_zeta
        real(rk),   allocatable :: xcoords(:,:,:)
        real(rk),   allocatable :: ycoords(:,:,:)
        real(rk),   allocatable :: zcoords(:,:,:)
        integer(ik), optional   :: clusterx

        integer(ik) :: ipt_xi, ipt_eta, ipt_zeta, ierr, &
                       npts_xi, npts_eta, npts_zeta
        real(rk)    :: x,y,z


        npts_xi   = nelem_xi   + 1
        npts_eta  = nelem_eta  + 1
        npts_zeta = nelem_zeta + 1


        allocate(xcoords(npts_xi,npts_eta,npts_zeta), ycoords(npts_xi,npts_eta,npts_zeta), &
                 zcoords(npts_xi,npts_eta,npts_zeta), stat=ierr)
        if (ierr /= 0) call AllocationError

        do ipt_zeta = 1,npts_zeta
            do ipt_eta = 1,npts_eta
                do ipt_xi = 1,npts_xi

                    if (present(clusterx)) then
                        if ( clusterx == -1 ) then
                            x = ONE - tanh( (PI/TWO)*(ONE - real(ipt_xi-1,rk)/real(npts_xi-1,rk) ) )/tanh(PI/TWO)
                        else if ( clusterx == 1 ) then
                            call chidg_signal(FATAL,"meshgen_NxNxN_linear: 'clusterx'=1 not yet implemented.")
                        else
                            call chidg_signal(FATAL,"meshgen_NxNxN_linear: Invalid value for 'clusterx'. -1,0,1.")
                        end if
                    else
                        x = real(ipt_xi-1,rk)/real(npts_xi-1,rk)
                    end if

                    if (ipt_xi == npts_xi) then
                        x = ONE
                    end if

                    y = real(ipt_eta-1,rk)/real(npts_eta-1,rk)
                    z = real(ipt_zeta-1,rk)/real(npts_zeta-1,rk)

                    xcoords(ipt_xi,ipt_eta,ipt_zeta) = x
                    ycoords(ipt_xi,ipt_eta,ipt_zeta) = y
                    zcoords(ipt_xi,ipt_eta,ipt_zeta) = z

                end do
            end do
        end do


    end subroutine meshgen_NxNxN_linear
    !**************************************************************************************











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
        real(rk)                    :: x,y,z

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
                    x = ONE*real(ipt_xi  -1,rk)/real(npts_xi  -1,rk)
                    y = ONE*real(ipt_eta -1,rk)/real(npts_eta -1,rk)
                    z = ONE*real(ipt_zeta-1,rk)/real(npts_zeta-1,rk)

                    xcoords(ipt_xi,ipt_eta,ipt_zeta) = x
                    ycoords(ipt_xi,ipt_eta,ipt_zeta) = y
                    zcoords(ipt_xi,ipt_eta,ipt_zeta) = z
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
        !real(rk), dimension(npt)    :: x,y,z
        real(rk)                    :: x,y,z

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
!        x = [ZERO, ONE, TWO, THREE, ZERO, ONE, TWO, THREE, ZERO, ONE, TWO, THREE, &
!             ZERO, ONE, TWO, THREE, ZERO, ONE, TWO, THREE, ZERO, ONE, TWO, THREE, &
!             ZERO, ONE, TWO, THREE, ZERO, ONE, TWO, THREE, ZERO, ONE, TWO, THREE, &
!             ZERO, ONE, TWO, THREE, ZERO, ONE, TWO, THREE, ZERO, ONE, TWO, THREE, &
!             ZERO, ONE, TWO, THREE, ZERO, ONE, TWO, THREE, ZERO, ONE, TWO, THREE, &
!             ZERO, ONE, TWO, THREE]
!
!        y = [ZERO, ZERO, ZERO, ZERO, ONE, ONE, ONE, ONE, TWO, TWO, TWO, TWO, &
!             THREE, THREE, THREE, THREE, &
!             ZERO, ZERO, ZERO, ZERO, ONE, ONE, ONE, ONE, TWO, TWO, TWO, TWO, &
!             THREE, THREE, THREE, THREE, &
!             ZERO, ZERO, ZERO, ZERO, ONE, ONE, ONE, ONE, TWO, TWO, TWO, TWO, &
!             THREE, THREE, THREE, THREE, &
!             ZERO, ZERO, ZERO, ZERO, ONE, ONE, ONE, ONE, TWO, TWO, TWO, TWO, &
!             THREE, THREE, THREE, THREE]
!
!        z = [ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, &
!             ZERO, ZERO, ZERO, ZERO, ZERO, &
!             ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, & 
!             ONE, ONE, ONE, &
!             TWO, TWO, TWO, TWO, TWO, TWO, TWO, TWO, TWO, TWO, TWO, TWO, TWO, &
!             TWO, TWO, TWO, &
!             THREE, THREE, THREE, THREE, THREE, THREE, THREE, THREE, THREE, THREE, &
!             THREE, THREE, THREE, THREE, THREE, THREE]

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
                    x = (real(ipt_xi  -1,rk))/real(npts_xi  -1,rk)
                    y = (real(ipt_eta -1,rk))/real(npts_eta -1,rk)
                    z = (real(ipt_zeta-1,rk))/real(npts_zeta-1,rk)

                    !xcoords(ipt_xi,ipt_eta,ipt_zeta) = x(ipt)
                    !ycoords(ipt_xi,ipt_eta,ipt_zeta) = y(ipt)
                    !zcoords(ipt_xi,ipt_eta,ipt_zeta) = z(ipt)
                    xcoords(ipt_xi,ipt_eta,ipt_zeta) = x
                    ycoords(ipt_xi,ipt_eta,ipt_zeta) = y
                    zcoords(ipt_xi,ipt_eta,ipt_zeta) = z
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










    !> Generate a set of points defining a 4x2x2 element mesh
    !!
    !!  @author Nathan A. Wukie
    !!  @param[inout]   pts     points_t array of rank-3 that gets allocated, 
    !!                          filled, and returned
    !---------------------------------------------------------------------------------------
    subroutine meshgen_4x2x2_linear(xcoords,ycoords,zcoords)
        real(rk),   allocatable :: xcoords(:,:,:)
        real(rk),   allocatable :: ycoords(:,:,:)
        real(rk),   allocatable :: zcoords(:,:,:)

        real(rk)    :: x,y,z
        integer(ik) :: ipt_xi, ipt_eta, ipt_zeta, ierr, &
                       npts_xi, npts_eta, npts_zeta


        ! elements (4x2x2) - linear
        !  *---*---*---*---*
        !  |\   \   \   \   \
        !  * *---*---*---*---*
        !  |\|\   \   \   \   \
        !  * * *---*---*---*---*
        !   \|\|   |   |   |   | 
        !    * *---*---*---*---*
        !     \|   |   |   |   | 
        !      *---*---*---*---*
        !


        npts_xi   = 5
        npts_eta  = 3
        npts_zeta = 3


        allocate(xcoords(npts_xi,npts_eta,npts_zeta), ycoords(npts_xi,npts_eta,npts_zeta), &
                 zcoords(npts_xi,npts_eta,npts_zeta), stat=ierr)
        if (ierr /= 0) call AllocationError


        do ipt_zeta = 1,npts_zeta
            do ipt_eta = 1,npts_eta
                do ipt_xi = 1,npts_xi
                    x = TWO*real(ipt_xi  -1,rk)/real(npts_xi  -1,rk)
                    y = ONE*real(ipt_eta -1,rk)/real(npts_eta -1,rk)
                    z = ONE*real(ipt_zeta-1,rk)/real(npts_zeta-1,rk)

                    xcoords(ipt_xi,ipt_eta,ipt_zeta) = x
                    ycoords(ipt_xi,ipt_eta,ipt_zeta) = y
                    zcoords(ipt_xi,ipt_eta,ipt_zeta) = z
                end do
            end do
        end do



    end subroutine meshgen_4x2x2_linear
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






end module mod_gridgen_blocks
