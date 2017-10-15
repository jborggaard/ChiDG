module tm_bc_operator
    use mod_constants,      only: HALF, ONE, TWO, ZERO
    use mod_kinds,          only: ik, rk
    use type_operator,      only: operator_t
    use type_chidg_worker,  only: chidg_worker_t
    use type_properties,    only: properties_t
    use DNAD_D
    implicit none




    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/12/2016
    !!
    !!
    !------------------------------------------------------------------------------------
    type, public, extends(operator_t)   :: tm_bc_operator_t


    contains

        procedure   :: init
        procedure   :: compute

    end type tm_bc_operator_t
    !*************************************************************************************






contains


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(tm_bc_operator_t),   intent(inout) :: self
        
        !
        ! Set operator name
        !
        call self%set_name("TM BC Flux")

        !
        ! Set operator type
        !
        call self%set_operator_type("BC Advective Flux")

        !
        ! Set operator equations
        !
        call self%add_primary_field("Pressure")

    end subroutine init
    !********************************************************************************





    !> Specialized compute routine for Extrapolation Boundary Condition
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !!  @param[in]      mesh    Mesh data containing elements and faces for the domain
    !!  @param[inout]   sdata   Solver data containing solution vector, rhs, linearization, etc.
    !!  @param[inout]   prop    properties_t object containing equations and material_t objects
    !!  @param[in]      face    face_info_t containing indices on location and misc information
    !!  @param[in]      flux    function_into_t containing info on the function being computed
    !!
    !-------------------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(tm_bc_operator_t),    intent(inout)   :: self
        type(chidg_worker_t),       intent(inout)   :: worker
        class(properties_t),        intent(inout)   :: prop


        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:)   ::  &
            p, density_bc, u_bc, v_bc,              &
            mom1_bc, mom2_bc,                       &
            flux_1,  flux_2,  flux_3

        real(rk),   allocatable, dimension(:)   :: r


        !
        ! Interpolate boundary condition state to face quadrature nodes
        !
        p = worker%get_field("Pressure", 'value', 'boundary')


        !
        ! Get model fields
        !
        density_bc = worker%get_field('Density',    'value', 'boundary')
        mom1_bc    = worker%get_field('Momentum-1', 'value', 'boundary')
        mom2_bc    = worker%get_field('Momentum-2', 'value', 'boundary')

        if (worker%coordinate_system() == 'Cylindrical') then
            r = worker%coordinate('1','face interior')
            mom2_bc = mom2_bc/r
        end if



        u_bc = mom1_bc/density_bc
        v_bc = mom2_bc/density_bc



        !=================================================
        !                   Momentum-1
        !=================================================
        flux_1 = density_bc
        flux_2 = density_bc*v_bc*v_bc  +  p
        flux_3 = density_bc
        flux_1 = ZERO
        flux_3 = ZERO

        call worker%integrate_boundary_condition('Pressure','Advection',flux_1,flux_2,flux_3)


    end subroutine compute
    !**********************************************************************************************







end module tm_bc_operator
