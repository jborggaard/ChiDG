module mod_prescribed_mesh_motion_function
#include <messenger.h>
    use mod_kinds,                              only: rk, ik
    use type_prescribed_mesh_motion_function,   only: prescribed_mesh_motion_function_t
    use type_pmmfvector,                        only: pmmfvector_t

    !
    ! Import prescribed_mesh_motion_functions
    !
    use pmmf_static,                            only: static_pmmf
    use pmmf_sinusoidal,                        only: sinusoidal_pmmf
    implicit none



    !
    ! Global vector of registered prescribed_mesh_motion_functions
    !
    type(pmmfvector_t)          :: registered_pmmfs
    logical                     :: initialized = .false.

contains


    !>  Register prescribed_mesh_motion_functions in a module vector.
    !!
    !!  This allows the available prescribed_mesh_motion_functions to be queried in the same way that they 
    !!  are registered for allocation. Adapted from mod_functions
    !!
    !!  @author Eric Wolf
    !!  @date   3/15/2017
    !!
    !!
    !--------------------------------------------------------------------------------------------
    subroutine register_prescribed_mesh_motion_functions()
        integer :: npmmfs, ipmmf

        !
        ! Instantiate prescribed_mesh_motion_functions
        !
        type(static_pmmf)                               :: static
        type(sinusoidal_pmmf)                           :: sinusoidal
        

        if ( .not. initialized ) then
            !
            ! Register in global vector
            !
            call registered_pmmfs%push_back(static)
            call registered_pmmfs%push_back(sinusoidal)
       
            !
            ! Initialize each boundary condition in set. Doesn't need modified.
            !
            npmmfs = registered_pmmfs%size()
            do ipmmf = 1,npmmfs
                call registered_pmmfs%data(ipmmf)%pmmf%init()
            end do

            !
            ! Confirm initialization
            !
            initialized = .true.

        end if

    end subroutine register_prescribed_mesh_motion_functions
    !********************************************************************************************











    !> Factory method for allocating concrete prescribed_mesh_motion_functions
    !!
    !!      - Allocate a concrete prescribed_mesh_motion_function_t based on the incoming string specification.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/9/2016
    !!
    !!  @param[in]      string  Character string used to select the appropriate boundary condition
    !!  @param[inout]   pmmf      Allocatable boundary condition
    !!
    !--------------------------------------------------------------
    subroutine create_prescribed_mesh_motion_function(pmmf,string)
        class(prescribed_mesh_motion_function_t),  allocatable,    intent(inout)   :: pmmf
        character(*),                       intent(in)      :: string

        integer(ik) :: ierr, pmmfindex


        if ( allocated(pmmf) ) then
            deallocate(pmmf)
        end if



        !
        ! Find equation set in 'registered_bcs' vector
        !
        pmmfindex = registered_pmmfs%index_by_name(trim(string))



        !
        ! Check prescribed_mesh_motion_function was found in 'registered_pmmfs'
        !
        if (pmmfindex == 0) call chidg_signal_one(FATAL,"create_prescribed_mesh_motion_function: prescribed_mesh_motion_function not recognized", trim(string))



        !
        ! Allocate conrete prescribed_mesh_motion_function_t instance
        !
        allocate(pmmf, source=registered_pmmfs%data(pmmfindex)%pmmf, stat=ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"create_prescribed_mesh_motion_function: error allocating prescribed_mesh_motion_function from global vector.")



        !
        ! Check boundary condition was allocated
        !
        if ( .not. allocated(pmmf) ) call chidg_signal(FATAL,"create_prescribed_mesh_motion_function: error allocating concrete prescribed_mesh_motion_function.")



    end subroutine create_prescribed_mesh_motion_function
    !*****************************************************************************************













    !>  Print a list of the registered prescribed_mesh_motion_functions. Really just a utility for 'chidg edit' to 
    !!  dynamically list the available 'prescribed_mesh_motion_function_t's.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/8/2016
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine list_prescribed_mesh_motion_functions()
        integer                         :: npmmfs, ipmmf
        character(len=:),   allocatable :: pmmf_name

        npmmfs = registered_pmmfs%size()


        do ipmmf = 1,npmmfs

            pmmf_name = registered_pmmfs%data(ipmmf)%pmmf%get_name()
            call write_line(trim(pmmf_name))

        end do ! ipmmf


    end subroutine list_prescribed_mesh_motion_functions
    !******************************************************************************************









end module mod_prescribed_mesh_motion_function
