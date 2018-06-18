module bc_giles_HB_base
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO, ONE, TWO, FOUR, HALF, ME, CYLINDRICAL,    &
                                      XI_MIN, XI_MAX, ETA_MIN, ETA_MAX, ZETA_MIN, ZETA_MAX, PI
    use mod_fluid,              only: Rgas, cp, gam
    use mod_interpolation,      only: interpolate_linear, interpolate_linear_ad
    use mod_gridspace,          only: linspace
    use mod_dft,                only: dft, idft_eval
    use mod_chimera,            only: find_gq_donor, find_gq_donor_parallel

    use type_point,             only: point_t
    use type_mesh,              only: mesh_t
    use type_bc_state,          only: bc_state_t
    use type_bc_patch,          only: bc_patch_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use type_face_info,         only: face_info_t, face_info_constructor
    use type_element_info,      only: element_info_t
    use mod_chidg_mpi,          only: IRANK
    use mod_interpolate,        only: interpolate_face_autodiff
    use mpi_f08,                only: MPI_REAL8, MPI_AllReduce, mpi_comm, MPI_INTEGER, MPI_BCast, MPI_MIN, MPI_MAX
    use ieee_arithmetic,        only: ieee_is_nan
    use DNAD_D
    implicit none



    !>  Name: Outlet - 3D Giles
    !!
    !!  Options:
    !!      : Average Pressure
    !!
    !!  Behavior:
    !!      
    !!  References:
    !!              
    !!  @author Nathan A. Wukie
    !!  @date   2/8/2018
    !!
    !---------------------------------------------------------------------------------
    type, public, abstract, extends(bc_state_t) :: giles_HB_base_t

        integer(ik) :: nr = 10
        !integer(ik) :: nfourier_space = 14
        !integer(ik) :: nfourier_space = 20
        integer(ik) :: nfourier_space = 8

        real(rk),   allocatable :: r(:)
        real(rk),   allocatable :: theta(:,:)   ! (nr,ntheta)
        real(rk)                :: theta_ref    

        type(element_info_t),   allocatable :: donor(:,:)
        real(rk),               allocatable :: donor_node(:,:,:)

    contains

        procedure   :: init                 ! Set-up bc state with options/name etc.
        procedure   :: init_bc_coupling     ! Initialize global coupling
        procedure   :: init_bc_postcomm     ! Implement specialized initialization

        procedure   :: get_q_interior

        procedure   :: compute_temporal_dft
        procedure   :: compute_spatial_dft
        procedure   :: compute_spatial_idft_gq
        procedure   :: compute_temporal_idft_gq
        procedure   :: compute_absorbing_inlet
        procedure   :: compute_absorbing_outlet
        procedure   :: interpolate_raux_to_rgq
        procedure   :: primitive_to_eigenmodes
        procedure   :: eigenmodes_to_primitive
        procedure   :: compute_eigenvalues
        procedure   :: compute_boundary_average
        procedure   :: analyze_bc_geometry
        procedure   :: initialize_fourier_discretization

        procedure   :: primitive_to_characteristics
        procedure   :: characteristics_to_primitive

    end type giles_HB_base_t
    !*********************************************************************************




contains



    !>
    !!
    !!  @author Nathan A. average_pressure 
    !!  @date   2/8/2017
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(giles_HB_base_t),   intent(inout) :: self

        ! Set name, family
        call self%set_name('Outlet - Giles Quasi3D Unsteady HB')
        call self%set_family('Outlet')

        ! Add functions
        call self%bcproperties%add('Average Pressure',    'Required')
        call self%bcproperties%add('Pitch',               'Required')
        call self%bcproperties%add('Spatial Periodicity', 'Required')

    end subroutine init
    !********************************************************************************




    !>  Initialize boundary group coupling.
    !!
    !!  Call global coupling routine to initialize implicit coupling between each
    !!  element with every other element on the boundary, a result of averaging
    !!  operations.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/18/2018
    !!
    !--------------------------------------------------------------------------------
    subroutine init_bc_coupling(self,mesh,group_ID,bc_comm)
        class(giles_HB_base_t), intent(inout)   :: self
        type(mesh_t),           intent(inout)   :: mesh
        integer(ik),            intent(in)      :: group_ID
        type(mpi_comm),         intent(in)      :: bc_comm

        call self%init_bc_coupling_global(mesh,group_ID,bc_comm)

    end subroutine init_bc_coupling
    !********************************************************************************




    !>  Default specialized initialization procedure. This is called from the base 
    !!  bc%init procedure and can be overwritten by derived types to implement 
    !!  specialized initiailization details.
    !!
    !!  By default, this routine does nothing. However, a particular bc_state_t could 
    !!  reimplement this routine to perform some specialized initialization calculations 
    !!  during initialization.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/27/2018
    !!
    !-------------------------------------------------------------------------------------
    subroutine init_bc_postcomm(self,mesh,group_ID,bc_comm)
        class(giles_HB_base_t),   intent(inout)   :: self
        type(mesh_t),                                intent(inout)   :: mesh
        integer(ik),                                 intent(in)      :: group_ID
        type(mpi_comm),                              intent(in)      :: bc_comm

        call self%analyze_bc_geometry(mesh,group_ID,bc_comm)

        call self%initialize_fourier_discretization(mesh,group_ID,bc_comm)

    end subroutine init_bc_postcomm
    !*************************************************************************************



    
    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/25/2018
    !!
    !------------------------------------------------------------------------------------
    subroutine get_q_interior(self,worker,bc_comm,    &
                              density,                &
                              vel1,                   &
                              vel2,                   &
                              vel3,                   &
                              pressure)
        class(giles_HB_base_t),         intent(inout)   :: self
        type(chidg_worker_t),           intent(inout)   :: worker
        type(mpi_comm),                 intent(in)      :: bc_comm
        type(AD_D),     allocatable,    intent(inout)   :: density(:,:,:)
        type(AD_D),     allocatable,    intent(inout)   :: vel1(:,:,:)
        type(AD_D),     allocatable,    intent(inout)   :: vel2(:,:,:)
        type(AD_D),     allocatable,    intent(inout)   :: vel3(:,:,:)
        type(AD_D),     allocatable,    intent(inout)   :: pressure(:,:,:)

        type(AD_D), allocatable, dimension(:,:,:) ::  &
            mom1, mom2, mom3, energy

        integer(ik) :: iradius, itime, nradius, ntheta, ntime, ierr

        ! Define Fourier space discretization to determine
        ! number of theta-samples are being taken
        nradius = size(self%r)
        ntheta  = size(self%theta,2)
        ntime   = worker%time_manager%ntime

        ! Allocate storage for discrete time instances
        allocate(density( nradius,ntheta,ntime),  &
                 mom1(    nradius,ntheta,ntime),  &
                 mom2(    nradius,ntheta,ntime),  &
                 mom3(    nradius,ntheta,ntime),  &
                 vel1(    nradius,ntheta,ntime),  &
                 vel2(    nradius,ntheta,ntime),  &
                 vel3(    nradius,ntheta,ntime),  &
                 energy(  nradius,ntheta,ntime),  &
                 pressure(nradius,ntheta,ntime), stat=ierr)
        if (ierr /= 0) call AllocationError


        ! Perform Fourier decomposition at each radial station.
        do iradius = 1,nradius
            do itime = 1,ntime
                ! Interpolate solution to physical_nodes at current radial station: [ntheta]
                density(iradius,:,itime) = worker%interpolate_field_general('Density',    donors=self%donor(iradius,:), donor_nodes=self%donor_node(iradius,:,:), itime=itime)
                mom1(iradius,:,itime)    = worker%interpolate_field_general('Momentum-1', donors=self%donor(iradius,:), donor_nodes=self%donor_node(iradius,:,:), itime=itime)
                mom2(iradius,:,itime)    = worker%interpolate_field_general('Momentum-2', donors=self%donor(iradius,:), donor_nodes=self%donor_node(iradius,:,:), itime=itime)
                mom3(iradius,:,itime)    = worker%interpolate_field_general('Momentum-3', donors=self%donor(iradius,:), donor_nodes=self%donor_node(iradius,:,:), itime=itime)
                energy(iradius,:,itime)  = worker%interpolate_field_general('Energy',     donors=self%donor(iradius,:), donor_nodes=self%donor_node(iradius,:,:), itime=itime)

                if (worker%coordinate_system() == 'Cylindrical') then
                    mom2(iradius,:,itime) = mom2(iradius,:,itime)/self%r(iradius)  ! convert to tangential momentum
                end if
            end do
        end do

        ! Compute velocities and pressure at each time
        vel1 = mom1/density
        vel2 = mom2/density
        vel3 = mom3/density
        pressure = (gam-ONE)*(energy - HALF*(mom1*mom1 + mom2*mom2 + mom3*mom3)/density)


    end subroutine get_q_interior
    !************************************************************************************






    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/12/2018
    !!
    !------------------------------------------------------------------------------------
    subroutine compute_temporal_dft(self,worker,bc_comm,                    &
                                    density_grid, vel1_grid, vel2_grid, vel3_grid, pressure_grid, c_grid, &
                                    density_Ft_real,  density_Ft_imag,      &
                                    vel1_Ft_real,     vel1_Ft_imag,         &
                                    vel2_Ft_real,     vel2_Ft_imag,         &
                                    vel3_Ft_real,     vel3_Ft_imag,         &
                                    pressure_Ft_real, pressure_Ft_imag,     &
                                    c_Ft_real,        c_Ft_imag)
        class(giles_HB_base_t),  intent(inout)   :: self
        type(chidg_worker_t),                       intent(inout)   :: worker
        type(mpi_comm),                             intent(in)      :: bc_comm
        type(AD_D),     allocatable,                intent(inout)   :: density_grid(:,:,:)
        type(AD_D),     allocatable,                intent(inout)   :: vel1_grid(:,:,:)
        type(AD_D),     allocatable,                intent(inout)   :: vel2_grid(:,:,:)
        type(AD_D),     allocatable,                intent(inout)   :: vel3_grid(:,:,:)
        type(AD_D),     allocatable,                intent(inout)   :: pressure_grid(:,:,:)
        type(AD_D),     allocatable,                intent(inout)   :: c_grid(:,:,:)
        type(AD_D),     allocatable,                intent(inout)   :: density_Ft_real(:,:,:)
        type(AD_D),     allocatable,                intent(inout)   :: density_Ft_imag(:,:,:)
        type(AD_D),     allocatable,                intent(inout)   :: vel1_Ft_real(:,:,:)
        type(AD_D),     allocatable,                intent(inout)   :: vel1_Ft_imag(:,:,:)
        type(AD_D),     allocatable,                intent(inout)   :: vel2_Ft_real(:,:,:)
        type(AD_D),     allocatable,                intent(inout)   :: vel2_Ft_imag(:,:,:)
        type(AD_D),     allocatable,                intent(inout)   :: vel3_Ft_real(:,:,:)
        type(AD_D),     allocatable,                intent(inout)   :: vel3_Ft_imag(:,:,:)
        type(AD_D),     allocatable,                intent(inout)   :: pressure_Ft_real(:,:,:)
        type(AD_D),     allocatable,                intent(inout)   :: pressure_Ft_imag(:,:,:)
        type(AD_D),     allocatable,                intent(inout)   :: c_Ft_real(:,:,:)
        type(AD_D),     allocatable,                intent(inout)   :: c_Ft_imag(:,:,:)


        type(AD_D), allocatable,    dimension(:)    ::                                          &
            density_real_tmp, vel1_real_tmp, vel2_real_tmp, vel3_real_tmp, pressure_real_tmp, c_real_tmp,   &
            density_imag_tmp, vel1_imag_tmp, vel2_imag_tmp, vel3_imag_tmp, pressure_imag_tmp, c_imag_tmp

        type(AD_D)  :: density_bar, vel1_bar, vel2_bar, vel3_bar, pressure_bar, c_bar

        integer(ik) :: nradius, ntheta, iradius, itheta, imode, itime, ntime, ierr


        ! Define Fourier space discretization to determine
        ! number of theta-samples are being taken
        ntheta  = size(self%theta,2)
        nradius = size(self%r)
        ntime   = worker%mesh%domain(worker%element_info%idomain_l)%elems(worker%element_info%ielement_l)%ntime

        ! Allocate storage for temporal dft
        allocate(density_Ft_real( nradius,ntheta,ntime), density_Ft_imag( nradius,ntheta,ntime),  &
                 vel1_Ft_real(    nradius,ntheta,ntime), vel1_Ft_imag(    nradius,ntheta,ntime),  &
                 vel2_Ft_real(    nradius,ntheta,ntime), vel2_Ft_imag(    nradius,ntheta,ntime),  &
                 vel3_Ft_real(    nradius,ntheta,ntime), vel3_Ft_imag(    nradius,ntheta,ntime),  &
                 pressure_Ft_real(nradius,ntheta,ntime), pressure_Ft_imag(nradius,ntheta,ntime),  &
                 c_Ft_real(       nradius,ntheta,ntime), c_Ft_imag(       nradius,ntheta,ntime), stat=ierr)
        if (ierr /= 0) call AllocationError

        ! Perform Fourier decomposition at each radial station.
        do iradius = 1,nradius
            do itheta = 1,ntheta
                call dft(density_grid(iradius,itheta,:),  ZERO*density_grid(iradius,itheta,:),  density_real_tmp,  density_imag_tmp )
                call dft(vel1_grid(iradius,itheta,:),     ZERO*vel1_grid(iradius,itheta,:),     vel1_real_tmp,     vel1_imag_tmp    )
                call dft(vel2_grid(iradius,itheta,:),     ZERO*vel2_grid(iradius,itheta,:),     vel2_real_tmp,     vel2_imag_tmp    )
                call dft(vel3_grid(iradius,itheta,:),     ZERO*vel3_grid(iradius,itheta,:),     vel3_real_tmp,     vel3_imag_tmp    )
                call dft(pressure_grid(iradius,itheta,:), ZERO*pressure_grid(iradius,itheta,:), pressure_real_tmp, pressure_imag_tmp)
                call dft(c_grid(iradius,itheta,:),        ZERO*c_grid(iradius,itheta,:),        c_real_tmp,        c_imag_tmp       )

                density_Ft_real( iradius,itheta,:) = density_real_tmp
                density_Ft_imag( iradius,itheta,:) = density_imag_tmp
                vel1_Ft_real(    iradius,itheta,:) = vel1_real_tmp
                vel1_Ft_imag(    iradius,itheta,:) = vel1_imag_tmp
                vel2_Ft_real(    iradius,itheta,:) = vel2_real_tmp
                vel2_Ft_imag(    iradius,itheta,:) = vel2_imag_tmp
                vel3_Ft_real(    iradius,itheta,:) = vel3_real_tmp
                vel3_Ft_imag(    iradius,itheta,:) = vel3_imag_tmp
                pressure_Ft_real(iradius,itheta,:) = pressure_real_tmp
                pressure_Ft_imag(iradius,itheta,:) = pressure_imag_tmp
                c_Ft_real(       iradius,itheta,:) = c_real_tmp
                c_Ft_imag(       iradius,itheta,:) = c_imag_tmp

            end do !itheta
        end do !iradius

    end subroutine compute_temporal_dft
    !********************************************************************************



    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/25/2018
    !!
    !--------------------------------------------------------------------------------
    subroutine compute_spatial_dft(self,worker,bc_comm,                   &
                                   density_Ft_real,  density_Ft_imag,     &
                                   vel1_Ft_real,     vel1_Ft_imag,        &
                                   vel2_Ft_real,     vel2_Ft_imag,        &
                                   vel3_Ft_real,     vel3_Ft_imag,        &
                                   pressure_Ft_real, pressure_Ft_imag,    &
                                   c_Ft_real,        c_Ft_imag,           &
                                   density_Fts_real,  density_Fts_imag,   &
                                   vel1_Fts_real,     vel1_Fts_imag,      &
                                   vel2_Fts_real,     vel2_Fts_imag,      &
                                   vel3_Fts_real,     vel3_Fts_imag,      &
                                   pressure_Fts_real, pressure_Fts_imag,  &
                                   c_Fts_real,        c_Fts_imag)
        class(giles_HB_base_t),  intent(inout)   :: self
        type(chidg_worker_t),                       intent(inout)   :: worker
        type(mpi_comm),                             intent(in)      :: bc_comm
        type(AD_D),     allocatable,                intent(inout)   :: density_Ft_real(:,:,:)
        type(AD_D),     allocatable,                intent(inout)   :: density_Ft_imag(:,:,:)
        type(AD_D),     allocatable,                intent(inout)   :: vel1_Ft_real(:,:,:)
        type(AD_D),     allocatable,                intent(inout)   :: vel1_Ft_imag(:,:,:)
        type(AD_D),     allocatable,                intent(inout)   :: vel2_Ft_real(:,:,:)
        type(AD_D),     allocatable,                intent(inout)   :: vel2_Ft_imag(:,:,:)
        type(AD_D),     allocatable,                intent(inout)   :: vel3_Ft_real(:,:,:)
        type(AD_D),     allocatable,                intent(inout)   :: vel3_Ft_imag(:,:,:)
        type(AD_D),     allocatable,                intent(inout)   :: pressure_Ft_real(:,:,:)
        type(AD_D),     allocatable,                intent(inout)   :: pressure_Ft_imag(:,:,:)
        type(AD_D),     allocatable,                intent(inout)   :: c_Ft_real(:,:,:)
        type(AD_D),     allocatable,                intent(inout)   :: c_Ft_imag(:,:,:)
        type(AD_D),     allocatable,                intent(inout)   :: density_Fts_real(:,:,:)
        type(AD_D),     allocatable,                intent(inout)   :: density_Fts_imag(:,:,:)
        type(AD_D),     allocatable,                intent(inout)   :: vel1_Fts_real(:,:,:)
        type(AD_D),     allocatable,                intent(inout)   :: vel1_Fts_imag(:,:,:)
        type(AD_D),     allocatable,                intent(inout)   :: vel2_Fts_real(:,:,:)
        type(AD_D),     allocatable,                intent(inout)   :: vel2_Fts_imag(:,:,:)
        type(AD_D),     allocatable,                intent(inout)   :: vel3_Fts_real(:,:,:)
        type(AD_D),     allocatable,                intent(inout)   :: vel3_Fts_imag(:,:,:)
        type(AD_D),     allocatable,                intent(inout)   :: pressure_Fts_real(:,:,:)
        type(AD_D),     allocatable,                intent(inout)   :: pressure_Fts_imag(:,:,:)
        type(AD_D),     allocatable,                intent(inout)   :: c_Fts_real(:,:,:)
        type(AD_D),     allocatable,                intent(inout)   :: c_Fts_imag(:,:,:)

        type(AD_D), allocatable,    dimension(:)    ::                                          &
            density_real_tmp, vel1_real_tmp, vel2_real_tmp, vel3_real_tmp, pressure_real_tmp, c_real_tmp, &
            density_imag_tmp, vel1_imag_tmp, vel2_imag_tmp, vel3_imag_tmp, pressure_imag_tmp, c_imag_tmp

        integer(ik)             :: nradius, ntheta, iradius, itheta, imode, itime, ntime, ierr
        real(rk)                :: shift_r, shift_i
        real(rk),   allocatable :: spatial_periodicity(:)

        ! Define Fourier space discretization to determine
        ! number of theta-samples being taken
        ntheta  = size(self%theta,2)
        nradius = size(self%r)
        ntime   = worker%mesh%domain(worker%element_info%idomain_l)%elems(worker%element_info%ielement_l)%ntime
        spatial_periodicity = self%bcproperties%compute('Spatial Periodicity', time=ZERO,coord=[point_t(ZERO,ZERO,ZERO)])
        
        ! Allocate storage in result
        allocate(density_Fts_real( nradius,ntheta,ntime), density_Fts_imag( nradius,ntheta,ntime),  &
                 vel1_Fts_real(    nradius,ntheta,ntime), vel1_Fts_imag(    nradius,ntheta,ntime),  &
                 vel2_Fts_real(    nradius,ntheta,ntime), vel2_Fts_imag(    nradius,ntheta,ntime),  &
                 vel3_Fts_real(    nradius,ntheta,ntime), vel3_Fts_imag(    nradius,ntheta,ntime),  &
                 pressure_Fts_real(nradius,ntheta,ntime), pressure_Fts_imag(nradius,ntheta,ntime),  &
                 c_Fts_real(       nradius,ntheta,ntime), c_Fts_imag(       nradius,ntheta,ntime), stat=ierr)
        if (ierr /= 0) call AllocationError

        ! Perform Fourier decomposition at each radial station.
        do iradius = 1,nradius
            do itime = 1,ntime

                ! DFT in space
                call dft(density_Ft_real( iradius,:,itime), density_Ft_imag( iradius,:,itime), density_real_tmp,  density_imag_tmp,  negate=.true.)
                call dft(vel1_Ft_real(    iradius,:,itime), vel1_Ft_imag(    iradius,:,itime), vel1_real_tmp,     vel1_imag_tmp,     negate=.true.)
                call dft(vel2_Ft_real(    iradius,:,itime), vel2_Ft_imag(    iradius,:,itime), vel2_real_tmp,     vel2_imag_tmp,     negate=.true.)
                call dft(vel3_Ft_real(    iradius,:,itime), vel3_Ft_imag(    iradius,:,itime), vel3_real_tmp,     vel3_imag_tmp,     negate=.true.)
                call dft(pressure_Ft_real(iradius,:,itime), pressure_Ft_imag(iradius,:,itime), pressure_real_tmp, pressure_imag_tmp, negate=.true.)
                call dft(c_Ft_real(       iradius,:,itime), c_Ft_imag(       iradius,:,itime), c_real_tmp,        c_imag_tmp,        negate=.true.)

                ! Adjust Fourier coefficients so their phase is relative to self%theta_ref
                ! instead of the minimum theta of the transform.
                !
                !       q(relative to theta_ref) = q(relative to theta_min) * e^(j 2pi imode delta_theta/spatial_periodicity)
                !
                ! NOTE: self%theta(:,1) are defined to be the DFT-theta_min at each radius
                !
                do imode = 1,size(density_real_tmp)
                    !shift_r = realpart(exp(cmplx(ZERO,ONE)*real(imode-1,rk)*TWO*PI*(self%theta_ref-self%theta(iradius,1))/spatial_periodicity(1)))
                    !shift_i = imagpart(exp(cmplx(ZERO,ONE)*real(imode-1,rk)*TWO*PI*(self%theta_ref-self%theta(iradius,1))/spatial_periodicity(1)))

                    !density_Fts_real( iradius,imode,itime) = density_real_tmp(imode)*shift_r  - density_imag_tmp(imode)*shift_i
                    !vel1_Fts_real(    iradius,imode,itime) = vel1_real_tmp(imode)*shift_r     - vel1_imag_tmp(imode)*shift_i
                    !vel2_Fts_real(    iradius,imode,itime) = vel2_real_tmp(imode)*shift_r     - vel2_imag_tmp(imode)*shift_i
                    !vel3_Fts_real(    iradius,imode,itime) = vel3_real_tmp(imode)*shift_r     - vel3_imag_tmp(imode)*shift_i
                    !pressure_Fts_real(iradius,imode,itime) = pressure_real_tmp(imode)*shift_r - pressure_imag_tmp(imode)*shift_i

                    !density_Fts_imag( iradius,imode,itime) = density_imag_tmp(imode)*shift_r  + density_real_tmp(imode)*shift_i
                    !vel1_Fts_imag(    iradius,imode,itime) = vel1_imag_tmp(imode)*shift_r     + vel1_real_tmp(imode)*shift_i
                    !vel2_Fts_imag(    iradius,imode,itime) = vel2_imag_tmp(imode)*shift_r     + vel2_real_tmp(imode)*shift_i
                    !vel3_Fts_imag(    iradius,imode,itime) = vel3_imag_tmp(imode)*shift_r     + vel3_real_tmp(imode)*shift_i
                    !pressure_Fts_imag(iradius,imode,itime) = pressure_imag_tmp(imode)*shift_r + pressure_real_tmp(imode)*shift_i

                    density_Fts_real( iradius,imode,itime) = density_real_tmp(imode) 
                    vel1_Fts_real(    iradius,imode,itime) = vel1_real_tmp(imode)    
                    vel2_Fts_real(    iradius,imode,itime) = vel2_real_tmp(imode)    
                    vel3_Fts_real(    iradius,imode,itime) = vel3_real_tmp(imode)    
                    pressure_Fts_real(iradius,imode,itime) = pressure_real_tmp(imode)
                    c_Fts_real(       iradius,imode,itime) = c_real_tmp(imode)

                    density_Fts_imag( iradius,imode,itime) = density_imag_tmp(imode) 
                    vel1_Fts_imag(    iradius,imode,itime) = vel1_imag_tmp(imode)    
                    vel2_Fts_imag(    iradius,imode,itime) = vel2_imag_tmp(imode)    
                    vel3_Fts_imag(    iradius,imode,itime) = vel3_imag_tmp(imode)    
                    pressure_Fts_imag(iradius,imode,itime) = pressure_imag_tmp(imode)
                    c_Fts_imag(       iradius,imode,itime) = c_imag_tmp(imode)
                end do !imode

            end do !itime
        end do !iradius

        print*, 'WARNING: need scaling since dft is only over a single passage.'
        print*, 'WARNING: check correct pitch in phase shift.'

    end subroutine compute_spatial_dft
    !*********************************************************************************



    !>  Compute inverse Fourier transform of spatio-temporal Fourier coefficients
    !!  to quadrature nodes for the current face.
    !!
    !!  Input: q_hat(r_gq)
    !!  Output: q_check(r_gq,theta_gq)
    !!
    !!  @author Nathan A. Wukie
    !!  @date   6/18/2018
    !!
    !----------------------------------------------------------------------------
    subroutine compute_spatial_idft_gq(self,worker,bc_comm, &
                                       density_hat_real,    &
                                       density_hat_imag,    &
                                       vel1_hat_real,       &
                                       vel1_hat_imag,       &
                                       vel2_hat_real,       &
                                       vel2_hat_imag,       &
                                       vel3_hat_real,       &
                                       vel3_hat_imag,       &
                                       pressure_hat_real,   &
                                       pressure_hat_imag,   &
                                       density_check_real,  &
                                       density_check_imag,  &
                                       vel1_check_real,     &
                                       vel1_check_imag,     &
                                       vel2_check_real,     &
                                       vel2_check_imag,     &
                                       vel3_check_real,     &
                                       vel3_check_imag,     &
                                       pressure_check_real, &
                                       pressure_check_imag)
        class(giles_HB_base_t),     intent(inout)   :: self
        type(chidg_worker_t),       intent(inout)   :: worker
        type(mpi_comm),             intent(in)      :: bc_comm
        type(AD_D), allocatable,    intent(inout)   :: density_hat_real(:,:,:)
        type(AD_D), allocatable,    intent(inout)   :: density_hat_imag(:,:,:)
        type(AD_D), allocatable,    intent(inout)   :: vel1_hat_real(:,:,:)
        type(AD_D), allocatable,    intent(inout)   :: vel1_hat_imag(:,:,:)
        type(AD_D), allocatable,    intent(inout)   :: vel2_hat_real(:,:,:)
        type(AD_D), allocatable,    intent(inout)   :: vel2_hat_imag(:,:,:)
        type(AD_D), allocatable,    intent(inout)   :: vel3_hat_real(:,:,:)
        type(AD_D), allocatable,    intent(inout)   :: vel3_hat_imag(:,:,:)
        type(AD_D), allocatable,    intent(inout)   :: pressure_hat_real(:,:,:)
        type(AD_D), allocatable,    intent(inout)   :: pressure_hat_imag(:,:,:)
        type(AD_D), allocatable,    intent(inout)   :: density_check_real(:,:)
        type(AD_D), allocatable,    intent(inout)   :: density_check_imag(:,:)
        type(AD_D), allocatable,    intent(inout)   :: vel1_check_real(:,:)
        type(AD_D), allocatable,    intent(inout)   :: vel1_check_imag(:,:)
        type(AD_D), allocatable,    intent(inout)   :: vel2_check_real(:,:)
        type(AD_D), allocatable,    intent(inout)   :: vel2_check_imag(:,:)
        type(AD_D), allocatable,    intent(inout)   :: vel3_check_real(:,:)
        type(AD_D), allocatable,    intent(inout)   :: vel3_check_imag(:,:)
        type(AD_D), allocatable,    intent(inout)   :: pressure_check_real(:,:)
        type(AD_D), allocatable,    intent(inout)   :: pressure_check_imag(:,:)

        type(point_t),  allocatable, dimension(:)   :: coords
        real(rk),       allocatable, dimension(:)   :: pitch
        real(rk)                                    :: theta_offset
        integer(ik) :: igq, itime

        ! Get BC properties
        pitch = self%bcproperties%compute('Pitch', worker%time(),worker%coords())

        ! Get gq coordinates 
        coords = worker%coords()

        ! Inverse DFT of spatio-temporal Fourier modes to give temporal Fourier
        ! modes at quadrature nodes.
        density_check_real  = ZERO*density_hat_real(:,1,:)
        vel1_check_real     = ZERO*density_hat_real(:,1,:)
        vel2_check_real     = ZERO*density_hat_real(:,1,:)
        vel3_check_real     = ZERO*density_hat_real(:,1,:)
        pressure_check_real = ZERO*density_hat_real(:,1,:)
        density_check_imag  = ZERO*density_hat_real(:,1,:)
        vel1_check_imag     = ZERO*density_hat_real(:,1,:)
        vel2_check_imag     = ZERO*density_hat_real(:,1,:)
        vel3_check_imag     = ZERO*density_hat_real(:,1,:)
        pressure_check_imag = ZERO*density_hat_real(:,1,:)

        do igq = 1,size(coords)
            do itime = 1,size(density_hat_real,3)
                !theta_offset = coords(igq)%c2_ - self%theta_ref
                theta_offset = coords(igq)%c2_ - self%theta(1,1)
                ! **** WARNING: probably want ipdft_eval here ****
                call idft_eval(density_hat_real(igq,:,itime),       &
                               density_hat_imag(igq,:,itime),       &
                               [theta_offset]/pitch(1),             &
                               density_check_real(igq:igq,itime),   &
                               density_check_imag(igq:igq,itime),negate=.true.)

                call idft_eval(vel1_hat_real(igq,:,itime),          &
                               vel1_hat_imag(igq,:,itime),          &
                               [theta_offset]/pitch(1),             &
                               vel1_check_real(igq:igq,itime),      &
                               vel1_check_imag(igq:igq,itime),negate=.true.)

                call idft_eval(vel2_hat_real(igq,:,itime),          &
                               vel2_hat_imag(igq,:,itime),          &
                               [theta_offset]/pitch(1),             &
                               vel2_check_real(igq:igq,itime),      &
                               vel2_check_imag(igq:igq,itime),negate=.true.)

                call idft_eval(vel3_hat_real(igq,:,itime),          &
                               vel3_hat_imag(igq,:,itime),          &
                               [theta_offset]/pitch(1),             &
                               vel3_check_real(igq:igq,itime),      &
                               vel3_check_imag(igq:igq,itime),negate=.true.)

                call idft_eval(pressure_hat_real(igq,:,itime),      &
                               pressure_hat_imag(igq,:,itime),      &
                               [theta_offset]/pitch(1),             &
                               pressure_check_real(igq:igq,itime),  &
                               pressure_check_imag(igq:igq,itime),negate=.true.)
            end do !itime
        end do !igq

    end subroutine compute_spatial_idft_gq
    !***************************************************************************








    !>  Compute inverse Fourier transform of spatio-temporal Fourier coefficients
    !!  to quadrature nodes for the current face.
    !!
    !!  Input: q_check(r_gq,theta_gq)
    !!  Output: q(r_gq,theta_gq,itime)
    !!
    !!  @author Nathan A. Wukie
    !!  @date   6/18/2018
    !!
    !----------------------------------------------------------------------------
    subroutine compute_temporal_idft_gq(self,worker,bc_comm,    &
                                        density_check_real,     &
                                        density_check_imag,     &
                                        vel1_check_real,        &
                                        vel1_check_imag,        &
                                        vel2_check_real,        &
                                        vel2_check_imag,        &
                                        vel3_check_real,        &
                                        vel3_check_imag,        &
                                        pressure_check_real,    &
                                        pressure_check_imag,    &
                                        density, vel1, vel2, vel3, pressure)
        class(giles_HB_base_t),     intent(inout)   :: self
        type(chidg_worker_t),       intent(inout)   :: worker
        type(mpi_comm),             intent(in)      :: bc_comm
        type(AD_D), allocatable,    intent(inout)   :: density_check_real(:,:)
        type(AD_D), allocatable,    intent(inout)   :: density_check_imag(:,:)
        type(AD_D), allocatable,    intent(inout)   :: vel1_check_real(:,:)
        type(AD_D), allocatable,    intent(inout)   :: vel1_check_imag(:,:)
        type(AD_D), allocatable,    intent(inout)   :: vel2_check_real(:,:)
        type(AD_D), allocatable,    intent(inout)   :: vel2_check_imag(:,:)
        type(AD_D), allocatable,    intent(inout)   :: vel3_check_real(:,:)
        type(AD_D), allocatable,    intent(inout)   :: vel3_check_imag(:,:)
        type(AD_D), allocatable,    intent(inout)   :: pressure_check_real(:,:)
        type(AD_D), allocatable,    intent(inout)   :: pressure_check_imag(:,:)
        type(AD_D), allocatable,    intent(inout)   :: density(:)
        type(AD_D), allocatable,    intent(inout)   :: vel1(:)
        type(AD_D), allocatable,    intent(inout)   :: vel2(:)
        type(AD_D), allocatable,    intent(inout)   :: vel3(:)
        type(AD_D), allocatable,    intent(inout)   :: pressure(:)

        type(AD_D),     allocatable, dimension(:)   :: expect_zero, density_tmp, vel1_tmp, vel2_tmp, vel3_tmp, pressure_tmp
        integer(ik) :: igq

        ! Allocate storage for boundary state
        expect_zero = [AD_D(1)]
        density  = ZERO*density_check_real(:,1)
        vel1     = ZERO*density_check_real(:,1)
        vel2     = ZERO*density_check_real(:,1)
        vel3     = ZERO*density_check_real(:,1)
        pressure = ZERO*density_check_real(:,1)

        ! Inverse DFT of temporal Fourier modes to give primitive variables
        ! at quarature nodes for the current time instance.
        density_tmp  = [ZERO*density_check_real(1,1)]
        vel1_tmp     = [ZERO*density_check_real(1,1)]
        vel2_tmp     = [ZERO*density_check_real(1,1)]
        vel3_tmp     = [ZERO*density_check_real(1,1)]
        pressure_tmp = [ZERO*density_check_real(1,1)]
        do igq = 1,size(density_check_real,1)
            ! **** WARNING: probably want ipdft_eval here ****
            call idft_eval(density_check_real(igq,:),                                       &
                           density_check_imag(igq,:),                                       &
                           [real(worker%itime-1,rk)/real(worker%time_manager%ntime,rk)],    &
                           density_tmp,                                                     &
                           expect_zero,symmetric=.true.)
            if (abs(expect_zero(1)) > 0.0000001) print*, 'WARNING: inverse transform returning complex values.'

            call idft_eval(vel1_check_real(igq,:),                                          &
                           vel1_check_imag(igq,:),                                          &
                           [real(worker%itime-1,rk)/real(worker%time_manager%ntime,rk)],    &
                           vel1_tmp,                                                        &
                           expect_zero,symmetric=.true.)
            if (abs(expect_zero(1)) > 0.0000001) print*, 'WARNING: inverse transform returning complex values.'

            call idft_eval(vel2_check_real(igq,:),                                          &
                           vel2_check_imag(igq,:),                                          &
                           [real(worker%itime-1,rk)/real(worker%time_manager%ntime,rk)],    &
                           vel2_tmp,                                                        &
                           expect_zero,symmetric=.true.)
            if (abs(expect_zero(1)) > 0.0000001) print*, 'WARNING: inverse transform returning complex values.'

            call idft_eval(vel3_check_real(igq,:),                                          &
                           vel3_check_imag(igq,:),                                          &
                           [real(worker%itime-1,rk)/real(worker%time_manager%ntime,rk)],    &
                           vel3_tmp,                                                        &
                           expect_zero,symmetric=.true.)
            if (abs(expect_zero(1)) > 0.0000001) print*, 'WARNING: inverse transform returning complex values.'

            call idft_eval(pressure_check_real(igq,:),                                      &
                           pressure_check_imag(igq,:),                                      &
                           [real(worker%itime-1,rk)/real(worker%time_manager%ntime,rk)],    &
                           pressure_tmp,                                                    &
                           expect_zero,symmetric=.true.)
            if (abs(expect_zero(1)) > 0.0000001) print*, 'WARNING: inverse transform returning complex values.'

            ! Accumulate total contribution from unsteady modes
            density(igq)  = density(igq)  + density_tmp(1)
            vel1(igq)     = vel1(igq)     + vel1_tmp(1)
            vel2(igq)     = vel2(igq)     + vel2_tmp(1)
            vel3(igq)     = vel3(igq)     + vel3_tmp(1)
            pressure(igq) = pressure(igq) + pressure_tmp(1)

        end do

    end subroutine compute_temporal_idft_gq
    !****************************************************************************








    !>  Interpolate Fourier coefficients from auxiliary grid to quadrature grid.
    !!
    !!  q_hat(r_gq) = I(q_hat(r_aux))
    !!
    !!  @author Nathan A. Wukie
    !!  @date   6/18/2018
    !!
    !-------------------------------------------------------------------------------------
    subroutine interpolate_raux_to_rgq(self,worker,bc_comm,                             &
                                       density_hat_real,        density_hat_imag,       &
                                       vel1_hat_real,           vel1_hat_imag,          &
                                       vel2_hat_real,           vel2_hat_imag,          &
                                       vel3_hat_real,           vel3_hat_imag,          &
                                       pressure_hat_real,       pressure_hat_imag,      &
                                       density_hat_real_gq,     density_hat_imag_gq,    &
                                       vel1_hat_real_gq,        vel1_hat_imag_gq,       &
                                       vel2_hat_real_gq,        vel2_hat_imag_gq,       &
                                       vel3_hat_real_gq,        vel3_hat_imag_gq,       &
                                       pressure_hat_real_gq,    pressure_hat_imag_gq)
        class(giles_HB_base_t),     intent(inout)   :: self
        type(chidg_worker_t),       intent(inout)   :: worker
        type(mpi_comm),             intent(in)      :: bc_comm
        type(AD_D), allocatable,    intent(in)      :: density_hat_real(:,:,:)
        type(AD_D), allocatable,    intent(in)      :: density_hat_imag(:,:,:)
        type(AD_D), allocatable,    intent(in)      :: vel1_hat_real(:,:,:)
        type(AD_D), allocatable,    intent(in)      :: vel1_hat_imag(:,:,:)
        type(AD_D), allocatable,    intent(in)      :: vel2_hat_real(:,:,:)
        type(AD_D), allocatable,    intent(in)      :: vel2_hat_imag(:,:,:)
        type(AD_D), allocatable,    intent(in)      :: vel3_hat_real(:,:,:)
        type(AD_D), allocatable,    intent(in)      :: vel3_hat_imag(:,:,:)
        type(AD_D), allocatable,    intent(in)      :: pressure_hat_real(:,:,:)
        type(AD_D), allocatable,    intent(in)      :: pressure_hat_imag(:,:,:)
        type(AD_D), allocatable,    intent(inout)   :: density_hat_real_gq(:,:,:)
        type(AD_D), allocatable,    intent(inout)   :: density_hat_imag_gq(:,:,:)
        type(AD_D), allocatable,    intent(inout)   :: vel1_hat_real_gq(:,:,:)
        type(AD_D), allocatable,    intent(inout)   :: vel1_hat_imag_gq(:,:,:)
        type(AD_D), allocatable,    intent(inout)   :: vel2_hat_real_gq(:,:,:)
        type(AD_D), allocatable,    intent(inout)   :: vel2_hat_imag_gq(:,:,:)
        type(AD_D), allocatable,    intent(inout)   :: vel3_hat_real_gq(:,:,:)
        type(AD_D), allocatable,    intent(inout)   :: vel3_hat_imag_gq(:,:,:)
        type(AD_D), allocatable,    intent(inout)   :: pressure_hat_real_gq(:,:,:)
        type(AD_D), allocatable,    intent(inout)   :: pressure_hat_imag_gq(:,:,:)
        
        type(point_t),  allocatable, dimension(:)   :: coords
        integer(ik) :: igq, itheta, itime, ierr


        ! Interpolate spatio-temporal Fourier coefficients to quadrature nodes
        ! linear interpolation between radial coordinates.
        coords = worker%coords()
        allocate(density_hat_real_gq( size(coords),size(density_hat_real,2),size(density_hat_real,3)), &
                 density_hat_imag_gq( size(coords),size(density_hat_real,2),size(density_hat_real,3)), &
                 vel1_hat_real_gq(    size(coords),size(density_hat_real,2),size(density_hat_real,3)), &
                 vel1_hat_imag_gq(    size(coords),size(density_hat_real,2),size(density_hat_real,3)), &
                 vel2_hat_real_gq(    size(coords),size(density_hat_real,2),size(density_hat_real,3)), &
                 vel2_hat_imag_gq(    size(coords),size(density_hat_real,2),size(density_hat_real,3)), &
                 vel3_hat_real_gq(    size(coords),size(density_hat_real,2),size(density_hat_real,3)), &
                 vel3_hat_imag_gq(    size(coords),size(density_hat_real,2),size(density_hat_real,3)), &
                 pressure_hat_real_gq(size(coords),size(density_hat_real,2),size(density_hat_real,3)), &
                 pressure_hat_imag_gq(size(coords),size(density_hat_real,2),size(density_hat_real,3)), &
                 stat=ierr)
        if (ierr /= 0) call AllocationError
        density_hat_real_gq  = ZERO*density_hat_real(1,1,1)
        density_hat_imag_gq  = ZERO*density_hat_real(1,1,1)
        vel1_hat_real_gq     = ZERO*density_hat_real(1,1,1)
        vel1_hat_imag_gq     = ZERO*density_hat_real(1,1,1)
        vel2_hat_real_gq     = ZERO*density_hat_real(1,1,1)
        vel2_hat_imag_gq     = ZERO*density_hat_real(1,1,1)
        vel3_hat_real_gq     = ZERO*density_hat_real(1,1,1)
        vel3_hat_imag_gq     = ZERO*density_hat_real(1,1,1)
        pressure_hat_real_gq = ZERO*density_hat_real(1,1,1)
        pressure_hat_imag_gq = ZERO*density_hat_real(1,1,1)
        do igq = 1,size(coords)
            do itheta = 1,size(density_hat_real,2)
                do itime = 1,size(density_hat_real,3)
                    density_hat_real_gq( igq,itheta,itime) = interpolate_linear_ad(self%r,density_hat_real(:,itheta,itime), coords(igq)%c1_)
                    density_hat_imag_gq( igq,itheta,itime) = interpolate_linear_ad(self%r,density_hat_imag(:,itheta,itime), coords(igq)%c1_)
                    vel1_hat_real_gq(    igq,itheta,itime) = interpolate_linear_ad(self%r,vel1_hat_real(:,itheta,itime),    coords(igq)%c1_)
                    vel1_hat_imag_gq(    igq,itheta,itime) = interpolate_linear_ad(self%r,vel1_hat_imag(:,itheta,itime),    coords(igq)%c1_)
                    vel2_hat_real_gq(    igq,itheta,itime) = interpolate_linear_ad(self%r,vel2_hat_real(:,itheta,itime),    coords(igq)%c1_)
                    vel2_hat_imag_gq(    igq,itheta,itime) = interpolate_linear_ad(self%r,vel2_hat_imag(:,itheta,itime),    coords(igq)%c1_)
                    vel3_hat_real_gq(    igq,itheta,itime) = interpolate_linear_ad(self%r,vel3_hat_real(:,itheta,itime),    coords(igq)%c1_)
                    vel3_hat_imag_gq(    igq,itheta,itime) = interpolate_linear_ad(self%r,vel3_hat_imag(:,itheta,itime),    coords(igq)%c1_)
                    pressure_hat_real_gq(igq,itheta,itime) = interpolate_linear_ad(self%r,pressure_hat_real(:,itheta,itime),coords(igq)%c1_)
                    pressure_hat_imag_gq(igq,itheta,itime) = interpolate_linear_ad(self%r,pressure_hat_imag(:,itheta,itime),coords(igq)%c1_)
                end do !itime
            end do !itheta
        end do !igq

    end subroutine interpolate_raux_to_rgq
    !********************************************************************************





    !>  DANIEL'S FORMULATION
    !!
    !!
    !!
    !!
    !--------------------------------------------------------------------------------
    subroutine compute_absorbing_inlet(self,worker,bc_comm,                     &
                                       density_real_m,    density_imag_m,       &
                                       vel1_real_m,       vel1_imag_m,          &
                                       vel2_real_m,       vel2_imag_m,          &
                                       vel3_real_m,       vel3_imag_m,          &
                                       pressure_real_m,   pressure_imag_m,      &
                                       c_real_m,          c_imag_m,             &
                                       density_real_p,    density_imag_p,       &
                                       vel1_real_p,       vel1_imag_p,          &
                                       vel2_real_p,       vel2_imag_p,          &
                                       vel3_real_p,       vel3_imag_p,          &
                                       pressure_real_p,   pressure_imag_p,      &
                                       c_real_p,          c_imag_p,             &
                                       density_real_abs,  density_imag_abs,     &
                                       vel1_real_abs,     vel1_imag_abs,        &
                                       vel2_real_abs,     vel2_imag_abs,        &
                                       vel3_real_abs,     vel3_imag_abs,        &
                                       pressure_real_abs, pressure_imag_abs)
        class(giles_HB_base_t),  intent(inout)   :: self
        type(chidg_worker_t),    intent(inout)   :: worker
        type(mpi_comm),          intent(in)      :: bc_comm
        type(AD_D),              intent(inout)   :: density_real_m(:,:,:)
        type(AD_D),              intent(inout)   :: density_imag_m(:,:,:)
        type(AD_D),              intent(inout)   :: vel1_real_m(:,:,:)
        type(AD_D),              intent(inout)   :: vel1_imag_m(:,:,:)
        type(AD_D),              intent(inout)   :: vel2_real_m(:,:,:)
        type(AD_D),              intent(inout)   :: vel2_imag_m(:,:,:)
        type(AD_D),              intent(inout)   :: vel3_real_m(:,:,:)
        type(AD_D),              intent(inout)   :: vel3_imag_m(:,:,:)
        type(AD_D),              intent(inout)   :: pressure_real_m(:,:,:)
        type(AD_D),              intent(inout)   :: pressure_imag_m(:,:,:)
        type(AD_D),              intent(inout)   :: c_real_m(:,:,:)
        type(AD_D),              intent(inout)   :: c_imag_m(:,:,:)
        type(AD_D),              intent(inout)   :: density_real_p(:,:,:)
        type(AD_D),              intent(inout)   :: density_imag_p(:,:,:)
        type(AD_D),              intent(inout)   :: vel1_real_p(:,:,:)
        type(AD_D),              intent(inout)   :: vel1_imag_p(:,:,:)
        type(AD_D),              intent(inout)   :: vel2_real_p(:,:,:)
        type(AD_D),              intent(inout)   :: vel2_imag_p(:,:,:)
        type(AD_D),              intent(inout)   :: vel3_real_p(:,:,:)
        type(AD_D),              intent(inout)   :: vel3_imag_p(:,:,:)
        type(AD_D),              intent(inout)   :: pressure_real_p(:,:,:)
        type(AD_D),              intent(inout)   :: pressure_imag_p(:,:,:)
        type(AD_D),              intent(inout)   :: c_real_p(:,:,:)
        type(AD_D),              intent(inout)   :: c_imag_p(:,:,:)
        type(AD_D),              intent(inout), allocatable   :: density_real_abs(:,:,:)
        type(AD_D),              intent(inout), allocatable   :: density_imag_abs(:,:,:)
        type(AD_D),              intent(inout), allocatable   :: vel1_real_abs(:,:,:)
        type(AD_D),              intent(inout), allocatable   :: vel1_imag_abs(:,:,:)
        type(AD_D),              intent(inout), allocatable   :: vel2_real_abs(:,:,:)
        type(AD_D),              intent(inout), allocatable   :: vel2_imag_abs(:,:,:)
        type(AD_D),              intent(inout), allocatable   :: vel3_real_abs(:,:,:)
        type(AD_D),              intent(inout), allocatable   :: vel3_imag_abs(:,:,:)
        type(AD_D),              intent(inout), allocatable   :: pressure_real_abs(:,:,:)
        type(AD_D),              intent(inout), allocatable   :: pressure_imag_abs(:,:,:)

        type(AD_D), allocatable, dimension(:,:,:) ::                &
            a1_real, a2_real, a3_real, a4_real, a5_real,            &
            a1_imag, a2_imag, a3_imag, a4_imag, a5_imag,            &
            a1_real_m, a2_real_m, a3_real_m, a4_real_m, a5_real_m,  &
            a1_imag_m, a2_imag_m, a3_imag_m, a4_imag_m, a5_imag_m,  &
            a1_real_p, a2_real_p, a3_real_p, a4_real_p, a5_real_p,  &
            a1_imag_p, a2_imag_p, a3_imag_p, a4_imag_p, a5_imag_p

        ! Project to eigenmodes
        call self%primitive_to_eigenmodes(worker,bc_comm,                   &
                                          density_real_p(:,1,1),            &
                                          vel1_real_p(:,1,1),               &
                                          vel2_real_p(:,1,1),               &
                                          vel3_real_p(:,1,1),               &
                                          pressure_real_p(:,1,1),           &
                                          c_real_p(:,1,1),                  &
                                          density_real_m,  density_imag_m,  &
                                          vel1_real_m,     vel1_imag_m,     &
                                          vel2_real_m,     vel2_imag_m,     &
                                          vel3_real_m,     vel3_imag_m,     &
                                          pressure_real_m, pressure_imag_m, &
                                          a1_real_m,       a1_imag_m,       &
                                          a2_real_m,       a2_imag_m,       &
                                          a3_real_m,       a3_imag_m,       &
                                          a4_real_m,       a4_imag_m,       &
                                          a5_real_m,       a5_imag_m)



        call self%primitive_to_eigenmodes(worker,bc_comm,                   &
                                          density_real_p(:,1,1),            &
                                          vel1_real_p(:,1,1),               &
                                          vel2_real_p(:,1,1),               &
                                          vel3_real_p(:,1,1),               &
                                          pressure_real_p(:,1,1),           &
                                          c_real_p(:,1,1),                  &
                                          density_real_p,  density_imag_p,  &
                                          vel1_real_p,     vel1_imag_p,     &
                                          vel2_real_p,     vel2_imag_p,     &
                                          vel3_real_p,     vel3_imag_p,     &
                                          pressure_real_p, pressure_imag_p, &
                                          a1_real_p,       a1_imag_p,       &
                                          a2_real_p,       a2_imag_p,       &
                                          a3_real_p,       a3_imag_p,       &
                                          a4_real_p,       a4_imag_p,       &
                                          a5_real_p,       a5_imag_p)


!        ! Zero out incoming amplitudes
!        a1_real(:,:,:) = ZERO
!        a1_imag(:,:,:) = ZERO
!        a2_real(:,:,:) = ZERO
!        a2_imag(:,:,:) = ZERO
!        a3_real(:,:,:) = ZERO
!        a3_imag(:,:,:) = ZERO
!        a5_real(:,:,:) = ZERO
!        a5_imag(:,:,:) = ZERO
!
!        ! User-specified amplitude
!        a1_real(:,2,2)  = 0.001_rk  !entropy
!        a1_real(:,27,3) = 0.001_rk  !entropy
!
!        !a3_real(:,2,2)  = 0.001_rk  !vorticity
!        !a5_real(:,27,2) = 0.001_rk  !downstream pressure


        ! Incoming perturbations from exterior state
        a1_real = a1_real_p
        a1_imag = a1_imag_p
        a2_real = a2_real_p
        a2_imag = a2_imag_p
        a3_real = a3_real_p
        a3_imag = a3_imag_p
        a5_real = a5_real_p
        a5_imag = a5_imag_p


        ! Zero out incoming amplitudes
        a1_real(:,:,:) = ZERO
        a1_imag(:,:,:) = ZERO
        a2_real(:,:,:) = ZERO
        a2_imag(:,:,:) = ZERO
        a3_real(:,:,:) = ZERO
        a3_imag(:,:,:) = ZERO
        a5_real(:,:,:) = ZERO
        a5_imag(:,:,:) = ZERO


!        a1_real(:,2,2) = 0.001_rk   ! entropy
!        a3_real(:,2,2) = 0.001_rk   ! vorticity
        a5_real(:,2,2) = 0.001_rk   ! downstream pressure




        ! Outgoing perturbations from interior state
        a4_real = a4_real_m
        a4_imag = a4_imag_m


        ! To initialize average
        density_real_abs  = density_real_p
        density_imag_abs  = density_imag_p
        vel1_real_abs     = vel1_real_p
        vel1_imag_abs     = vel1_imag_p
        vel2_real_abs     = vel2_real_p
        vel2_imag_abs     = vel2_imag_p
        vel3_real_abs     = vel3_real_p
        vel3_imag_abs     = vel3_imag_p
        pressure_real_abs = pressure_real_p
        pressure_imag_abs = pressure_imag_p


        call self%eigenmodes_to_primitive(worker,bc_comm,                       &
                                          density_real_p(:,1,1),                &
                                          vel1_real_p(:,1,1),                   &
                                          vel2_real_p(:,1,1),                   &
                                          vel3_real_p(:,1,1),                   &
                                          pressure_real_p(:,1,1),               &
                                          c_real_p(:,1,1),                      &
                                          a1_real,           a1_imag,           &
                                          a2_real,           a2_imag,           &
                                          a3_real,           a3_imag,           &
                                          a4_real,           a4_imag,           &
                                          a5_real,           a5_imag,           &
                                          density_real_abs,  density_imag_abs,  &
                                          vel1_real_abs,     vel1_imag_abs,     &
                                          vel2_real_abs,     vel2_imag_abs,     &
                                          vel3_real_abs,     vel3_imag_abs,     &
                                          pressure_real_abs, pressure_imag_abs)


!        ! Compute 1-4 characteristics from extrapolation: difference in radius-local mean and boundary average
!        call self%compute_boundary_average(worker,bc_comm,density_bar,vel1_bar,vel2_bar,vel3_bar,pressure_bar, &
!                                                          density_avg,vel1_avg,vel2_avg,vel3_avg,pressure_avg)
!        c_avg = sqrt(gam*pressure_avg/density_avg)
!
!
!
!        ! Get boundary condition Total Temperature, Total Pressure, and normal vector
!        PT   = self%bcproperties%compute('Total Pressure',   worker%time(),worker%coords())
!        TT   = self%bcproperties%compute('Total Temperature',worker%time(),worker%coords())
!
!        ! Get user-input normal vector and normalize
!        n1 = self%bcproperties%compute('Normal-1', worker%time(), worker%coords())
!        n2 = self%bcproperties%compute('Normal-2', worker%time(), worker%coords())
!        n3 = self%bcproperties%compute('Normal-3', worker%time(), worker%coords())
!
!        !   Explicit allocation to handle GCC bug:
!        !       GCC/GFortran Bugzilla Bug 52162 
!        allocate(nmag(size(n1)), stat=ierr)
!        if (ierr /= 0) call AllocationError
!
!        nmag = sqrt(n1*n1 + n2*n2 + n3*n3)
!        n1 = n1/nmag
!        n2 = n2/nmag
!        n3 = n3/nmag
!
!
!        ! Override spatio-temporal mean according to specified total conditions
!        T_avg = TT(1)*(pressure_avg/PT(1))**((gam-ONE)/gam)
!        density_avg = pressure_avg/(T_avg*Rgas)
!        vmag = sqrt(TWO*cp*(TT(1)-T_avg))
!        vel1_avg = n1(1)*vmag
!        vel2_avg = n2(1)*vmag
!        vel3_avg = n3(1)*vmag
!
!        density_real_abs(:,1,1)  = density_avg
!        vel1_real_abs(:,1,1)     = vel1_avg
!        vel2_real_abs(:,1,1)     = vel2_avg
!        vel3_real_abs(:,1,1)     = vel3_avg
!        pressure_real_abs(:,1,1) = pressure_avg



    end subroutine compute_absorbing_inlet
    !********************************************************************************







    !> DANIEL'S FORMULATION
    !!
    !!
    !!
    !!
    !--------------------------------------------------------------------------------
    subroutine compute_absorbing_outlet(self,worker,bc_comm,                    &
                                        density_real_m,    density_imag_m,      &
                                        vel1_real_m,       vel1_imag_m,         &
                                        vel2_real_m,       vel2_imag_m,         &
                                        vel3_real_m,       vel3_imag_m,         &
                                        pressure_real_m,   pressure_imag_m,     &
                                        c_real_m,          c_imag_m,            &
                                        density_real_p,    density_imag_p,      &
                                        vel1_real_p,       vel1_imag_p,         &
                                        vel2_real_p,       vel2_imag_p,         &
                                        vel3_real_p,       vel3_imag_p,         &
                                        pressure_real_p,   pressure_imag_p,     &
                                        c_real_p,          c_imag_p,            &
                                        density_real_abs,  density_imag_abs,    &
                                        vel1_real_abs,     vel1_imag_abs,       &
                                        vel2_real_abs,     vel2_imag_abs,       &
                                        vel3_real_abs,     vel3_imag_abs,       &
                                        pressure_real_abs, pressure_imag_abs)
        class(giles_HB_base_t),     intent(inout)   :: self
        type(chidg_worker_t),       intent(inout)   :: worker
        type(mpi_comm),             intent(in)      :: bc_comm
        type(AD_D),                 intent(inout)   :: density_real_m(:,:,:)
        type(AD_D),                 intent(inout)   :: density_imag_m(:,:,:)
        type(AD_D),                 intent(inout)   :: vel1_real_m(:,:,:)
        type(AD_D),                 intent(inout)   :: vel1_imag_m(:,:,:)
        type(AD_D),                 intent(inout)   :: vel2_real_m(:,:,:)
        type(AD_D),                 intent(inout)   :: vel2_imag_m(:,:,:)
        type(AD_D),                 intent(inout)   :: vel3_real_m(:,:,:)
        type(AD_D),                 intent(inout)   :: vel3_imag_m(:,:,:)
        type(AD_D),                 intent(inout)   :: pressure_real_m(:,:,:)
        type(AD_D),                 intent(inout)   :: pressure_imag_m(:,:,:)
        type(AD_D),                 intent(inout)   :: c_real_m(:,:,:)
        type(AD_D),                 intent(inout)   :: c_imag_m(:,:,:)
        type(AD_D),                 intent(inout)   :: density_real_p(:,:,:)
        type(AD_D),                 intent(inout)   :: density_imag_p(:,:,:)
        type(AD_D),                 intent(inout)   :: vel1_real_p(:,:,:)
        type(AD_D),                 intent(inout)   :: vel1_imag_p(:,:,:)
        type(AD_D),                 intent(inout)   :: vel2_real_p(:,:,:)
        type(AD_D),                 intent(inout)   :: vel2_imag_p(:,:,:)
        type(AD_D),                 intent(inout)   :: vel3_real_p(:,:,:)
        type(AD_D),                 intent(inout)   :: vel3_imag_p(:,:,:)
        type(AD_D),                 intent(inout)   :: pressure_real_p(:,:,:)
        type(AD_D),                 intent(inout)   :: pressure_imag_p(:,:,:)
        type(AD_D),                 intent(inout)   :: c_real_p(:,:,:)
        type(AD_D),                 intent(inout)   :: c_imag_p(:,:,:)
        type(AD_D),                 intent(inout)   :: density_real_abs(:,:,:)
        type(AD_D),                 intent(inout)   :: density_imag_abs(:,:,:)
        type(AD_D),                 intent(inout)   :: vel1_real_abs(:,:,:)
        type(AD_D),                 intent(inout)   :: vel1_imag_abs(:,:,:)
        type(AD_D),                 intent(inout)   :: vel2_real_abs(:,:,:)
        type(AD_D),                 intent(inout)   :: vel2_imag_abs(:,:,:)
        type(AD_D),                 intent(inout)   :: vel3_real_abs(:,:,:)
        type(AD_D),                 intent(inout)   :: vel3_imag_abs(:,:,:)
        type(AD_D),                 intent(inout)   :: pressure_real_abs(:,:,:)
        type(AD_D),                 intent(inout)   :: pressure_imag_abs(:,:,:)

        type(AD_D), allocatable, dimension(:,:,:) ::                &
            a1_real, a2_real, a3_real, a4_real, a5_real,            &
            a1_imag, a2_imag, a3_imag, a4_imag, a5_imag,            &
            a1_real_m, a2_real_m, a3_real_m, a4_real_m, a5_real_m,  &
            a1_imag_m, a2_imag_m, a3_imag_m, a4_imag_m, a5_imag_m,  &
            a1_real_p, a2_real_p, a3_real_p, a4_real_p, a5_real_p,  &
            a1_imag_p, a2_imag_p, a3_imag_p, a4_imag_p, a5_imag_p

        real(rk),       allocatable, dimension(:)   :: p_user, pitch

        ! Retrieve target average pressure
        p_user = self%bcproperties%compute('Average Pressure',worker%time(),worker%coords())
        pitch  = self%bcproperties%compute('Pitch',worker%time(),worker%coords())

        ! Project interior to eigenmodes
        call self%primitive_to_eigenmodes(worker,bc_comm,                   &
                                          density_real_m(:,1,1),            &
                                          vel1_real_m(:,1,1),               &
                                          vel2_real_m(:,1,1),               &
                                          vel3_real_m(:,1,1),               &
                                          pressure_real_m(:,1,1),           &
                                          c_real_m(:,1,1),                  &
                                          density_real_m,  density_imag_m,  &
                                          vel1_real_m,     vel1_imag_m,     &
                                          vel2_real_m,     vel2_imag_m,     &
                                          vel3_real_m,     vel3_imag_m,     &
                                          pressure_real_m, pressure_imag_m, &
                                          a1_real_m,       a1_imag_m,       &
                                          a2_real_m,       a2_imag_m,       &
                                          a3_real_m,       a3_imag_m,       &
                                          a4_real_m,       a4_imag_m,       &
                                          a5_real_m,       a5_imag_m)

        ! Project exterior to eigenmodes
        call self%primitive_to_eigenmodes(worker,bc_comm,                   &
                                          density_real_p(:,1,1),            &
                                          vel1_real_p(:,1,1),               &
                                          vel2_real_p(:,1,1),               &
                                          vel3_real_p(:,1,1),               &
                                          pressure_real_p(:,1,1),           &
                                          c_real_p(:,1,1),                  &
                                          density_real_p,  density_imag_p,  &
                                          vel1_real_p,     vel1_imag_p,     &
                                          vel2_real_p,     vel2_imag_p,     &
                                          vel3_real_p,     vel3_imag_p,     &
                                          pressure_real_p, pressure_imag_p, &
                                          a1_real_p,       a1_imag_p,       &
                                          a2_real_p,       a2_imag_p,       &
                                          a3_real_p,       a3_imag_p,       &
                                          a4_real_p,       a4_imag_p,       &
                                          a5_real_p,       a5_imag_p)


        ! Outoing amplitudes from interior
        a1_real = a1_real_m
        a1_imag = a1_imag_m
        a2_real = a2_real_m
        a2_imag = a2_imag_m
        a3_real = a3_real_m
        a3_imag = a3_imag_m
        a4_real = a4_real_m
        a4_imag = a4_imag_m

        ! Incoming amplitudes from exterior
        a5_real = a5_real_p
        a5_imag = a5_imag_p

!        ! Zero out incoming amplitudes
!        a5_real(:,:,:) = ZERO
!        a5_imag(:,:,:) = ZERO


        ! Convert back to primitive variables
        call self%eigenmodes_to_primitive(worker,bc_comm,                       &
                                          density_real_p(:,1,1),                &
                                          vel1_real_p(:,1,1),                   &
                                          vel2_real_p(:,1,1),                   &
                                          vel3_real_p(:,1,1),                   &
                                          pressure_real_p(:,1,1),               &
                                          c_real_p(:,1,1),                      &
                                          a1_real,          a1_imag,            &
                                          a2_real,          a2_imag,            &
                                          a3_real,          a3_imag,            &
                                          a4_real,          a4_imag,            &
                                          a5_real,          a5_imag,            &
                                          density_real_abs, density_imag_abs,   &
                                          vel1_real_abs,    vel1_imag_abs,      &
                                          vel2_real_abs,    vel2_imag_abs,      &
                                          vel3_real_abs,    vel3_imag_abs,      &
                                          pressure_real_abs,pressure_imag_abs)

!        pressure_real(:,1,1) = p_user(1)

    end subroutine compute_absorbing_outlet
    !********************************************************************************





    !>  DANIEL'S FORMULATION
    !!
    !!
    !!
    !!
    !--------------------------------------------------------------------------------
    subroutine primitive_to_eigenmodes(self,worker,bc_comm,             &
                                       density_bar_r,                   &
                                       vel1_bar_r,                      &
                                       vel2_bar_r,                      &
                                       vel3_bar_r,                      &
                                       pressure_bar_r,                  &
                                       c_bar_r,                         &
                                       density_real,  density_imag,     &
                                       vel1_real,     vel1_imag,        &
                                       vel2_real,     vel2_imag,        &
                                       vel3_real,     vel3_imag,        &
                                       pressure_real, pressure_imag,    &
                                       a1_real,       a1_imag,          &
                                       a2_real,       a2_imag,          &
                                       a3_real,       a3_imag,          &
                                       a4_real,       a4_imag,          &
                                       a5_real,       a5_imag)
        class(giles_HB_base_t),     intent(inout)   :: self
        type(chidg_worker_t),       intent(inout)   :: worker
        type(mpi_comm),             intent(in)      :: bc_comm
        type(AD_D),                 intent(in)      :: density_bar_r(:)
        type(AD_D),                 intent(in)      :: vel1_bar_r(:)
        type(AD_D),                 intent(in)      :: vel2_bar_r(:)
        type(AD_D),                 intent(in)      :: vel3_bar_r(:)
        type(AD_D),                 intent(in)      :: pressure_bar_r(:)
        type(AD_D),                 intent(in)      :: c_bar_r(:)
        type(AD_D),                 intent(inout)   :: density_real(:,:,:)
        type(AD_D),                 intent(inout)   :: density_imag(:,:,:)
        type(AD_D),                 intent(inout)   :: vel1_real(:,:,:)
        type(AD_D),                 intent(inout)   :: vel1_imag(:,:,:)
        type(AD_D),                 intent(inout)   :: vel2_real(:,:,:)
        type(AD_D),                 intent(inout)   :: vel2_imag(:,:,:)
        type(AD_D),                 intent(inout)   :: vel3_real(:,:,:)
        type(AD_D),                 intent(inout)   :: vel3_imag(:,:,:)
        type(AD_D),                 intent(inout)   :: pressure_real(:,:,:)
        type(AD_D),                 intent(inout)   :: pressure_imag(:,:,:)
        type(AD_D), allocatable,    intent(inout)   :: a1_real(:,:,:)
        type(AD_D), allocatable,    intent(inout)   :: a1_imag(:,:,:)
        type(AD_D), allocatable,    intent(inout)   :: a2_real(:,:,:)
        type(AD_D), allocatable,    intent(inout)   :: a2_imag(:,:,:)
        type(AD_D), allocatable,    intent(inout)   :: a3_real(:,:,:)
        type(AD_D), allocatable,    intent(inout)   :: a3_imag(:,:,:)
        type(AD_D), allocatable,    intent(inout)   :: a4_real(:,:,:)
        type(AD_D), allocatable,    intent(inout)   :: a4_imag(:,:,:)
        type(AD_D), allocatable,    intent(inout)   :: a5_real(:,:,:)
        type(AD_D), allocatable,    intent(inout)   :: a5_imag(:,:,:)

        type(AD_D)  :: k1, k2, k3, k4_real, k4_imag, k5_real, k5_imag, &
                       density_bar, vel1_bar, vel2_bar, vel3_bar, pressure_bar, c_bar, &
                       c_real, c_imag, denom, Tinv_real(5,5), Tinv_imag(5,5)

        real(rk),       allocatable, dimension(:)   :: pitch, unorm3
        real(rk)                                    :: theta_offset, omega, kz, lm
        integer(ik)                                 :: iradius, igq, ierr, itheta, ntheta, itime, ntime, nr

        pitch  = self%bcproperties%compute('Pitch',worker%time(),worker%coords())
        nr     = size(density_real,1)
        ntheta = size(density_real,2)
        ntime  = size(density_real,3)

        ! Allocate storage
        a1_real = ZERO*density_real
        a2_real = ZERO*density_real
        a3_real = ZERO*density_real
        a4_real = ZERO*density_real
        a5_real = ZERO*density_real
        a1_imag = ZERO*density_real
        a2_imag = ZERO*density_real
        a3_imag = ZERO*density_real
        a4_imag = ZERO*density_real
        a5_imag = ZERO*density_real

        ! Project
        print*, 'WARNING! CHECK DEFINITION OF lm PITCH!'
        do iradius = 1,nr
            ! Get radius-local average
            density_bar  = density_bar_r(iradius)
            vel1_bar     = vel1_bar_r(iradius)
            vel2_bar     = vel2_bar_r(iradius)
            vel3_bar     = vel3_bar_r(iradius)
            pressure_bar = pressure_bar_r(iradius)
            c_bar        = c_bar_r(iradius)

            ! starting with 2 here because the first mode is treated with 1D characteristics
            do itheta = 1,ntheta
                do itime = 1,ntime
                    
                    ! Space-time average handled at the bottom
                    if (itime == 1 .and. itheta == 1) then
                    
                    else

                    ! Get temporal/spatial frequencies
                    omega = get_omega(worker,itime)
                    lm    = get_lm(itheta,ntheta)

                    ! Compute wavenumbers
                    call self%compute_eigenvalues(worker,lm,omega,vel1_bar,vel2_bar,vel3_bar,c_bar, &
                                                  kz, k1, k2, k3, k4_real, k4_imag, k5_real, k5_imag)


                    ! Zero Tinv
                    Tinv_real = ZERO*density_real(1,1,1)
                    Tinv_imag = ZERO*density_real(1,1,1)


                    ! Assemble (1,:)
                    Tinv_real(1,1) = ONE/density_bar
                    Tinv_real(1,5) = -ONE/(density_bar*c_bar*c_bar)

                    ! Assemble (2,:)
                    Tinv_real(2,3) = ONE/c_bar

                    ! Assemble (3,:)
                    Tinv_real(3,2) = (-kz/(c_bar*(k1*k1+kz*kz)))
                    Tinv_real(3,4) = ( k1/(c_bar*(k1*k1+kz*kz)))
                    Tinv_real(3,5) = (-kz/(density_bar*c_bar*vel3_bar*(k1*k1+kz*kz)))


                    ! Assemble (4,:)
                    ! Contribution from vel3:
                    denom  = TWO*c_bar*c_bar*((k4_real*k1 + kz*kz)**TWO + k4_imag*k4_imag*k1*k1)
                    c_real = -vel3_bar*( (k4_real*k1 - k1*k1)*(k4_real*k1 + kz*kz) + (k4_imag*k4_imag*k1*k1) )/denom
                    c_imag = -vel3_bar*( (k4_real*k1 + kz*kz)*k4_imag*k1 - (k4_real*k1 - k1*k1)*k4_imag*k1 )/denom
                    Tinv_real(4,2) = c_real
                    Tinv_imag(4,2) = c_imag

                    ! Contribution from vel2: 
                    c_real = -vel3_bar*( (k4_real*kz - k1*kz)*(k4_real*k1 + kz*kz) + (k4_imag*k4_imag*kz*k1) )/denom
                    c_imag = -vel3_bar*( (k4_real*k1 + kz*kz)*k4_imag*kz - (k4_real*kz - k1*kz)*k4_imag*k1 )/denom
                    Tinv_real(4,4) = c_real
                    Tinv_imag(4,4) = c_imag

                    ! Contribution from pressure:
                    Tinv_real(4,5) = ONE/(TWO*density_bar*c_bar*c_bar)


                    ! Assemble (5,:)
                    ! Contribution from vel3:
                    denom  = TWO*c_bar*c_bar*((k5_real*k1 + kz*kz)**TWO + k5_imag*k5_imag*k1*k1)
                    c_real = -vel3_bar*( (k5_real*k1 - k1*k1)*(k5_real*k1 + kz*kz) + (k5_imag*k5_imag*k1*k1) )/denom
                    c_imag = -vel3_bar*( (k5_real*k1 + kz*kz)*k5_imag*k1 - (k5_real*k1 - k1*k1)*k5_imag*k1 )/denom
                    Tinv_real(5,2) = c_real
                    Tinv_imag(5,2) = c_imag

                    ! Contribution from vel2: 
                    c_real = -vel3_bar*( (k5_real*kz - k1*kz)*(k5_real*k1 + kz*kz) + (k5_imag*k5_imag*kz*k1) )/denom
                    c_imag = -vel3_bar*( (k5_real*k1 + kz*kz)*k5_imag*kz - (k5_real*kz - k1*kz)*k5_imag*k1 )/denom
                    Tinv_real(5,4) = c_real
                    Tinv_imag(5,4) = c_imag

                    ! Contribution from pressure:
                    Tinv_real(5,5) = ONE/(TWO*density_bar*c_bar*c_bar)


                    ! Project
                    a1_real(iradius,itheta,itime) = Tinv_real(1,1)*density_real(iradius,itheta,itime)  - Tinv_imag(1,1)*density_imag(iradius,itheta,itime) + &
                                                    Tinv_real(1,2)*vel3_real(iradius,itheta,itime)     - Tinv_imag(1,2)*vel3_imag(iradius,itheta,itime)    + &
                                                    Tinv_real(1,3)*vel1_real(iradius,itheta,itime)     - Tinv_imag(1,3)*vel1_imag(iradius,itheta,itime)    + &
                                                    Tinv_real(1,4)*vel2_real(iradius,itheta,itime)     - Tinv_imag(1,4)*vel2_imag(iradius,itheta,itime)    + &
                                                    Tinv_real(1,5)*pressure_real(iradius,itheta,itime) - Tinv_imag(1,5)*pressure_imag(iradius,itheta,itime)

                    a1_imag(iradius,itheta,itime) = Tinv_real(1,1)*density_imag(iradius,itheta,itime)  + Tinv_imag(1,1)*density_real(iradius,itheta,itime) + &
                                                    Tinv_real(1,2)*vel3_imag(iradius,itheta,itime)     + Tinv_imag(1,2)*vel3_real(iradius,itheta,itime)    + &
                                                    Tinv_real(1,3)*vel1_imag(iradius,itheta,itime)     + Tinv_imag(1,3)*vel1_real(iradius,itheta,itime)    + &
                                                    Tinv_real(1,4)*vel2_imag(iradius,itheta,itime)     + Tinv_imag(1,4)*vel2_real(iradius,itheta,itime)    + &
                                                    Tinv_real(1,5)*pressure_imag(iradius,itheta,itime) + Tinv_imag(1,5)*pressure_real(iradius,itheta,itime)


                    a2_real(iradius,itheta,itime) = Tinv_real(2,1)*density_real(iradius,itheta,itime)  - Tinv_imag(2,1)*density_imag(iradius,itheta,itime) + &
                                                    Tinv_real(2,2)*vel3_real(iradius,itheta,itime)     - Tinv_imag(2,2)*vel3_imag(iradius,itheta,itime)    + &
                                                    Tinv_real(2,3)*vel1_real(iradius,itheta,itime)     - Tinv_imag(2,3)*vel1_imag(iradius,itheta,itime)    + &
                                                    Tinv_real(2,4)*vel2_real(iradius,itheta,itime)     - Tinv_imag(2,4)*vel2_imag(iradius,itheta,itime)    + &
                                                    Tinv_real(2,5)*pressure_real(iradius,itheta,itime) - Tinv_imag(2,5)*pressure_imag(iradius,itheta,itime)

                    a2_imag(iradius,itheta,itime) = Tinv_real(2,1)*density_imag(iradius,itheta,itime)  + Tinv_imag(2,1)*density_real(iradius,itheta,itime) + &
                                                    Tinv_real(2,2)*vel3_imag(iradius,itheta,itime)     + Tinv_imag(2,2)*vel3_real(iradius,itheta,itime)    + &
                                                    Tinv_real(2,3)*vel1_imag(iradius,itheta,itime)     + Tinv_imag(2,3)*vel1_real(iradius,itheta,itime)    + &
                                                    Tinv_real(2,4)*vel2_imag(iradius,itheta,itime)     + Tinv_imag(2,4)*vel2_real(iradius,itheta,itime)    + &
                                                    Tinv_real(2,5)*pressure_imag(iradius,itheta,itime) + Tinv_imag(2,5)*pressure_real(iradius,itheta,itime)


                    a3_real(iradius,itheta,itime) = Tinv_real(3,1)*density_real(iradius,itheta,itime)  - Tinv_imag(3,1)*density_imag(iradius,itheta,itime) + &
                                                    Tinv_real(3,2)*vel3_real(iradius,itheta,itime)     - Tinv_imag(3,2)*vel3_imag(iradius,itheta,itime)    + &
                                                    Tinv_real(3,3)*vel1_real(iradius,itheta,itime)     - Tinv_imag(3,3)*vel1_imag(iradius,itheta,itime)    + &
                                                    Tinv_real(3,4)*vel2_real(iradius,itheta,itime)     - Tinv_imag(3,4)*vel2_imag(iradius,itheta,itime)    + &
                                                    Tinv_real(3,5)*pressure_real(iradius,itheta,itime) - Tinv_imag(3,5)*pressure_imag(iradius,itheta,itime)

                    a3_imag(iradius,itheta,itime) = Tinv_real(3,1)*density_imag(iradius,itheta,itime)  + Tinv_imag(3,1)*density_real(iradius,itheta,itime) + &
                                                    Tinv_real(3,2)*vel3_imag(iradius,itheta,itime)     + Tinv_imag(3,2)*vel3_real(iradius,itheta,itime)    + &
                                                    Tinv_real(3,3)*vel1_imag(iradius,itheta,itime)     + Tinv_imag(3,3)*vel1_real(iradius,itheta,itime)    + &
                                                    Tinv_real(3,4)*vel2_imag(iradius,itheta,itime)     + Tinv_imag(3,4)*vel2_real(iradius,itheta,itime)    + &
                                                    Tinv_real(3,5)*pressure_imag(iradius,itheta,itime) + Tinv_imag(3,5)*pressure_real(iradius,itheta,itime)


                    a4_real(iradius,itheta,itime) = Tinv_real(4,1)*density_real(iradius,itheta,itime)  - Tinv_imag(4,1)*density_imag(iradius,itheta,itime) + &
                                                    Tinv_real(4,2)*vel3_real(iradius,itheta,itime)     - Tinv_imag(4,2)*vel3_imag(iradius,itheta,itime)    + &
                                                    Tinv_real(4,3)*vel1_real(iradius,itheta,itime)     - Tinv_imag(4,3)*vel1_imag(iradius,itheta,itime)    + &
                                                    Tinv_real(4,4)*vel2_real(iradius,itheta,itime)     - Tinv_imag(4,4)*vel2_imag(iradius,itheta,itime)    + &
                                                    Tinv_real(4,5)*pressure_real(iradius,itheta,itime) - Tinv_imag(4,5)*pressure_imag(iradius,itheta,itime)

                    a4_imag(iradius,itheta,itime) = Tinv_real(4,1)*density_imag(iradius,itheta,itime)  + Tinv_imag(4,1)*density_real(iradius,itheta,itime) + &
                                                    Tinv_real(4,2)*vel3_imag(iradius,itheta,itime)     + Tinv_imag(4,2)*vel3_real(iradius,itheta,itime)    + &
                                                    Tinv_real(4,3)*vel1_imag(iradius,itheta,itime)     + Tinv_imag(4,3)*vel1_real(iradius,itheta,itime)    + &
                                                    Tinv_real(4,4)*vel2_imag(iradius,itheta,itime)     + Tinv_imag(4,4)*vel2_real(iradius,itheta,itime)    + &
                                                    Tinv_real(4,5)*pressure_imag(iradius,itheta,itime) + Tinv_imag(4,5)*pressure_real(iradius,itheta,itime)


                    a5_real(iradius,itheta,itime) = Tinv_real(5,1)*density_real(iradius,itheta,itime)  - Tinv_imag(5,1)*density_imag(iradius,itheta,itime) + &
                                                    Tinv_real(5,2)*vel3_real(iradius,itheta,itime)     - Tinv_imag(5,2)*vel3_imag(iradius,itheta,itime)    + &
                                                    Tinv_real(5,3)*vel1_real(iradius,itheta,itime)     - Tinv_imag(5,3)*vel1_imag(iradius,itheta,itime)    + &
                                                    Tinv_real(5,4)*vel2_real(iradius,itheta,itime)     - Tinv_imag(5,4)*vel2_imag(iradius,itheta,itime)    + &
                                                    Tinv_real(5,5)*pressure_real(iradius,itheta,itime) - Tinv_imag(5,5)*pressure_imag(iradius,itheta,itime)

                    a5_imag(iradius,itheta,itime) = Tinv_real(5,1)*density_imag(iradius,itheta,itime)  + Tinv_imag(5,1)*density_real(iradius,itheta,itime) + &
                                                    Tinv_real(5,2)*vel3_imag(iradius,itheta,itime)     + Tinv_imag(5,2)*vel3_real(iradius,itheta,itime)    + &
                                                    Tinv_real(5,3)*vel1_imag(iradius,itheta,itime)     + Tinv_imag(5,3)*vel1_real(iradius,itheta,itime)    + &
                                                    Tinv_real(5,4)*vel2_imag(iradius,itheta,itime)     + Tinv_imag(5,4)*vel2_real(iradius,itheta,itime)    + &
                                                    Tinv_real(5,5)*pressure_imag(iradius,itheta,itime) + Tinv_imag(5,5)*pressure_real(iradius,itheta,itime)

                    end if
                end do !itime
            end do !itheta
        end do !iradius



    end subroutine primitive_to_eigenmodes
    !********************************************************************************





    !>  DANIEL'S FORMULATION
    !!
    !!
    !!
    !!
    !--------------------------------------------------------------------------------
    subroutine eigenmodes_to_primitive(self,worker,bc_comm,             &
                                       density_bar_r, vel1_bar_r, vel2_bar_r, vel3_bar_r, pressure_bar_r, c_bar_r, &
                                       a1_real,       a1_imag,          &
                                       a2_real,       a2_imag,          &
                                       a3_real,       a3_imag,          &
                                       a4_real,       a4_imag,          &
                                       a5_real,       a5_imag,          &
                                       density_real,  density_imag,     &
                                       vel1_real,     vel1_imag,        &
                                       vel2_real,     vel2_imag,        &
                                       vel3_real,     vel3_imag,        &
                                       pressure_real, pressure_imag)
        class(giles_HB_base_t),     intent(inout)   :: self
        type(chidg_worker_t),       intent(inout)   :: worker
        type(mpi_comm),             intent(in)      :: bc_comm
        type(AD_D),                 intent(in)      :: density_bar_r(:)
        type(AD_D),                 intent(in)      :: vel1_bar_r(:)
        type(AD_D),                 intent(in)      :: vel2_bar_r(:)
        type(AD_D),                 intent(in)      :: vel3_bar_r(:)
        type(AD_D),                 intent(in)      :: pressure_bar_r(:)
        type(AD_D),                 intent(in)      :: c_bar_r(:)
        type(AD_D),                 intent(in)      :: a1_real(:,:,:)
        type(AD_D),                 intent(in)      :: a1_imag(:,:,:)
        type(AD_D),                 intent(in)      :: a2_real(:,:,:)
        type(AD_D),                 intent(in)      :: a2_imag(:,:,:)
        type(AD_D),                 intent(in)      :: a3_real(:,:,:)
        type(AD_D),                 intent(in)      :: a3_imag(:,:,:)
        type(AD_D),                 intent(in)      :: a4_real(:,:,:)
        type(AD_D),                 intent(in)      :: a4_imag(:,:,:)
        type(AD_D),                 intent(in)      :: a5_real(:,:,:)
        type(AD_D),                 intent(in)      :: a5_imag(:,:,:)
        type(AD_D),                 intent(inout)   :: density_real(:,:,:)
        type(AD_D),                 intent(inout)   :: density_imag(:,:,:)
        type(AD_D),                 intent(inout)   :: vel1_real(:,:,:)
        type(AD_D),                 intent(inout)   :: vel1_imag(:,:,:)
        type(AD_D),                 intent(inout)   :: vel2_real(:,:,:)
        type(AD_D),                 intent(inout)   :: vel2_imag(:,:,:)
        type(AD_D),                 intent(inout)   :: vel3_real(:,:,:)
        type(AD_D),                 intent(inout)   :: vel3_imag(:,:,:)
        type(AD_D),                 intent(inout)   :: pressure_real(:,:,:)
        type(AD_D),                 intent(inout)   :: pressure_imag(:,:,:)


        type(AD_D)  :: k1, k2, k3, k4_real, k4_imag, k5_real, k5_imag, &
                       density_bar, vel1_bar, vel2_bar, vel3_bar, pressure_bar, &
                       c_bar,  denom, c_real, c_imag, T_real(5,5), T_imag(5,5)
        
        real(rk),       allocatable, dimension(:)   :: pitch, unorm3
        real(rk)                                    :: theta_offset, omega, kz, lm
        integer(ik)                                 :: iradius, igq, ierr, itheta, ntheta, itime, ntime, nr

        pitch  = self%bcproperties%compute('Pitch',worker%time(),worker%coords())
        nr     = size(density_real,1)
        ntheta = size(density_real,2)
        ntime  = size(density_real,3)

        ! Project
        print*, 'WARNING! CHECK DEFINITION OF lm PITCH!'
        do iradius = 1,nr
            ! Get radius-local average
            density_bar  = density_bar_r(iradius)
            vel1_bar     = vel1_bar_r(iradius)
            vel2_bar     = vel2_bar_r(iradius)
            vel3_bar     = vel3_bar_r(iradius)
            pressure_bar = pressure_bar_r(iradius)
            c_bar        = c_bar_r(iradius)

            ! starting with 2 here because the first mode is treated with 1D characteristics
            do itheta = 1,ntheta
                do itime = 1,ntime
                    
                    ! Space-time average handled at the bottom
                    if (itime == 1 .and. itheta == 1) then
                    
                    else

                    ! Get temporal/spatial frequencies
                    omega = get_omega(worker,itime)
                    lm    = get_lm(itheta,ntheta)

                    ! Compute wavenumbers
                    call self%compute_eigenvalues(worker,lm,omega,vel1_bar,vel2_bar,vel3_bar,c_bar, &
                                                  kz, k1, k2, k3, k4_real, k4_imag, k5_real, k5_imag)

                    ! First zero fields
                    density_real(iradius,itheta,itime)  = ZERO
                    density_imag(iradius,itheta,itime)  = ZERO
                    vel1_real(iradius,itheta,itime)     = ZERO
                    vel1_imag(iradius,itheta,itime)     = ZERO
                    vel2_real(iradius,itheta,itime)     = ZERO
                    vel2_imag(iradius,itheta,itime)     = ZERO
                    vel3_real(iradius,itheta,itime)     = ZERO
                    vel3_imag(iradius,itheta,itime)     = ZERO
                    pressure_real(iradius,itheta,itime) = ZERO
                    pressure_imag(iradius,itheta,itime) = ZERO
                    T_real = ZERO*density_real(1,1,1)
                    T_imag = ZERO*density_real(1,1,1)


                    ! Assemble (1,:)
                    T_real(1,1) = density_bar
                    T_real(1,4) = density_bar
                    T_real(1,5) = density_bar

                    
                    ! Assemble (2,:)
                    T_real(2,3) = -c_bar*kz

                    ! Contribution from a4
                    denom = vel3_bar*( (k4_real-k1)**TWO + k4_imag*k4_imag)
                    c_real = (-c_bar*c_bar*(k4_real*(k4_real-k1) + k4_imag*k4_imag)/denom)
                    c_imag = (c_bar*c_bar*(k4_imag*k1)/denom)
                    T_real(2,4) = c_real
                    T_imag(2,4) = c_imag

                    ! Contribution from a5
                    denom = vel3_bar*( (k5_real-k1)**TWO + k5_imag*k5_imag)
                    c_real = (-c_bar*c_bar*(k5_real*(k5_real-k1) + k5_imag*k5_imag)/denom)
                    c_imag = (c_bar*c_bar*(k5_imag*k1)/denom)
                    T_real(2,5) = c_real
                    T_imag(2,5) = c_imag


                    ! Assemble (3,:)
                    T_real(3,2) = c_bar


                    ! Assemble (4,:)
                    ! Contribution from a3
                    T_real(4,3) = c_bar*k1

                    ! Contribution from a4
                    denom = vel3_bar*((k4_real-k1)**TWO + k4_imag*k4_imag)
                    c_real = -c_bar*c_bar*kz*(k4_real-k1)/denom
                    c_imag =  c_bar*c_bar*kz*k4_imag/denom
                    T_real(4,4) = c_real
                    T_imag(4,4) = c_imag

                    ! Contribution from a5
                    denom = vel3_bar*((k5_real-k1)**TWO + k5_imag*k5_imag)
                    c_real = -c_bar*c_bar*kz*(k5_real-k1)/denom
                    c_imag =  c_bar*c_bar*kz*k5_imag/denom
                    T_real(4,5) = c_real
                    T_imag(4,5) = c_imag


                    ! Assemble (5,:)
                    T_real(5,4) = density_bar*c_bar*c_bar
                    T_real(5,5) = density_bar*c_bar*c_bar



                    ! Density
                    density_real(iradius,itheta,itime) = T_real(1,1)*a1_real(iradius,itheta,itime) - T_imag(1,1)*a1_imag(iradius,itheta,itime) + &
                                                         T_real(1,2)*a2_real(iradius,itheta,itime) - T_imag(1,2)*a2_imag(iradius,itheta,itime) + &
                                                         T_real(1,3)*a3_real(iradius,itheta,itime) - T_imag(1,3)*a3_imag(iradius,itheta,itime) + &
                                                         T_real(1,4)*a4_real(iradius,itheta,itime) - T_imag(1,4)*a4_imag(iradius,itheta,itime) + &
                                                         T_real(1,5)*a5_real(iradius,itheta,itime) - T_imag(1,5)*a5_imag(iradius,itheta,itime)

                    density_imag(iradius,itheta,itime) = T_real(1,1)*a1_imag(iradius,itheta,itime) + T_imag(1,1)*a1_real(iradius,itheta,itime) + &
                                                         T_real(1,2)*a2_imag(iradius,itheta,itime) + T_imag(1,2)*a2_real(iradius,itheta,itime) + &
                                                         T_real(1,3)*a3_imag(iradius,itheta,itime) + T_imag(1,3)*a3_real(iradius,itheta,itime) + &
                                                         T_real(1,4)*a4_imag(iradius,itheta,itime) + T_imag(1,4)*a4_real(iradius,itheta,itime) + &
                                                         T_real(1,5)*a5_imag(iradius,itheta,itime) + T_imag(1,5)*a5_real(iradius,itheta,itime)

                    ! Vel3
                    vel3_real(iradius,itheta,itime) = T_real(2,1)*a1_real(iradius,itheta,itime) - T_imag(2,1)*a1_imag(iradius,itheta,itime) + &
                                                      T_real(2,2)*a2_real(iradius,itheta,itime) - T_imag(2,2)*a2_imag(iradius,itheta,itime) + &
                                                      T_real(2,3)*a3_real(iradius,itheta,itime) - T_imag(2,3)*a3_imag(iradius,itheta,itime) + &
                                                      T_real(2,4)*a4_real(iradius,itheta,itime) - T_imag(2,4)*a4_imag(iradius,itheta,itime) + &
                                                      T_real(2,5)*a5_real(iradius,itheta,itime) - T_imag(2,5)*a5_imag(iradius,itheta,itime)

                    vel3_imag(iradius,itheta,itime) = T_real(2,1)*a1_imag(iradius,itheta,itime) + T_imag(2,1)*a1_real(iradius,itheta,itime) + &
                                                      T_real(2,2)*a2_imag(iradius,itheta,itime) + T_imag(2,2)*a2_real(iradius,itheta,itime) + &
                                                      T_real(2,3)*a3_imag(iradius,itheta,itime) + T_imag(2,3)*a3_real(iradius,itheta,itime) + &
                                                      T_real(2,4)*a4_imag(iradius,itheta,itime) + T_imag(2,4)*a4_real(iradius,itheta,itime) + &
                                                      T_real(2,5)*a5_imag(iradius,itheta,itime) + T_imag(2,5)*a5_real(iradius,itheta,itime)

                    ! Vel1
                    vel1_real(iradius,itheta,itime) = T_real(3,1)*a1_real(iradius,itheta,itime) - T_imag(3,1)*a1_imag(iradius,itheta,itime) + &
                                                      T_real(3,2)*a2_real(iradius,itheta,itime) - T_imag(3,2)*a2_imag(iradius,itheta,itime) + &
                                                      T_real(3,3)*a3_real(iradius,itheta,itime) - T_imag(3,3)*a3_imag(iradius,itheta,itime) + &
                                                      T_real(3,4)*a4_real(iradius,itheta,itime) - T_imag(3,4)*a4_imag(iradius,itheta,itime) + &
                                                      T_real(3,5)*a5_real(iradius,itheta,itime) - T_imag(3,5)*a5_imag(iradius,itheta,itime)

                    vel1_imag(iradius,itheta,itime) = T_real(3,1)*a1_imag(iradius,itheta,itime) + T_imag(3,1)*a1_real(iradius,itheta,itime) + &
                                                      T_real(3,2)*a2_imag(iradius,itheta,itime) + T_imag(3,2)*a2_real(iradius,itheta,itime) + &
                                                      T_real(3,3)*a3_imag(iradius,itheta,itime) + T_imag(3,3)*a3_real(iradius,itheta,itime) + &
                                                      T_real(3,4)*a4_imag(iradius,itheta,itime) + T_imag(3,4)*a4_real(iradius,itheta,itime) + &
                                                      T_real(3,5)*a5_imag(iradius,itheta,itime) + T_imag(3,5)*a5_real(iradius,itheta,itime)

                    ! Vel2
                    vel2_real(iradius,itheta,itime) = T_real(4,1)*a1_real(iradius,itheta,itime) - T_imag(4,1)*a1_imag(iradius,itheta,itime) + &
                                                      T_real(4,2)*a2_real(iradius,itheta,itime) - T_imag(4,2)*a2_imag(iradius,itheta,itime) + &
                                                      T_real(4,3)*a3_real(iradius,itheta,itime) - T_imag(4,3)*a3_imag(iradius,itheta,itime) + &
                                                      T_real(4,4)*a4_real(iradius,itheta,itime) - T_imag(4,4)*a4_imag(iradius,itheta,itime) + &
                                                      T_real(4,5)*a5_real(iradius,itheta,itime) - T_imag(4,5)*a5_imag(iradius,itheta,itime)

                    vel2_imag(iradius,itheta,itime) = T_real(4,1)*a1_imag(iradius,itheta,itime) + T_imag(4,1)*a1_real(iradius,itheta,itime) + &
                                                      T_real(4,2)*a2_imag(iradius,itheta,itime) + T_imag(4,2)*a2_real(iradius,itheta,itime) + &
                                                      T_real(4,3)*a3_imag(iradius,itheta,itime) + T_imag(4,3)*a3_real(iradius,itheta,itime) + &
                                                      T_real(4,4)*a4_imag(iradius,itheta,itime) + T_imag(4,4)*a4_real(iradius,itheta,itime) + &
                                                      T_real(4,5)*a5_imag(iradius,itheta,itime) + T_imag(4,5)*a5_real(iradius,itheta,itime)

                    ! Pressure
                    pressure_real(iradius,itheta,itime) = T_real(5,1)*a1_real(iradius,itheta,itime) - T_imag(5,1)*a1_imag(iradius,itheta,itime) + &
                                                          T_real(5,2)*a2_real(iradius,itheta,itime) - T_imag(5,2)*a2_imag(iradius,itheta,itime) + &
                                                          T_real(5,3)*a3_real(iradius,itheta,itime) - T_imag(5,3)*a3_imag(iradius,itheta,itime) + &
                                                          T_real(5,4)*a4_real(iradius,itheta,itime) - T_imag(5,4)*a4_imag(iradius,itheta,itime) + &
                                                          T_real(5,5)*a5_real(iradius,itheta,itime) - T_imag(5,5)*a5_imag(iradius,itheta,itime)

                    pressure_imag(iradius,itheta,itime) = T_real(5,1)*a1_imag(iradius,itheta,itime) + T_imag(5,1)*a1_real(iradius,itheta,itime) + &
                                                          T_real(5,2)*a2_imag(iradius,itheta,itime) + T_imag(5,2)*a2_real(iradius,itheta,itime) + &
                                                          T_real(5,3)*a3_imag(iradius,itheta,itime) + T_imag(5,3)*a3_real(iradius,itheta,itime) + &
                                                          T_real(5,4)*a4_imag(iradius,itheta,itime) + T_imag(5,4)*a4_real(iradius,itheta,itime) + &
                                                          T_real(5,5)*a5_imag(iradius,itheta,itime) + T_imag(5,5)*a5_real(iradius,itheta,itime)

                    end if

                end do !itime
            end do !itheta
        end do !iradius


    end subroutine eigenmodes_to_primitive
    !********************************************************************************
    



    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   5/8/2018
    !!
    !--------------------------------------------------------------------------------
    subroutine compute_eigenvalues(self, worker, lm, omega, vel1_bar, vel2_bar, vel3_bar, c_bar, &
                                   kz, k1, k2, k3, k4_real, k4_imag, k5_real, k5_imag)
        class(giles_HB_base_t), intent(inout)   :: self
        type(chidg_worker_t),   intent(inout)   :: worker
        real(rk),               intent(in)      :: lm
        real(rk),               intent(in)      :: omega
        type(AD_D),             intent(in)      :: vel1_bar
        type(AD_D),             intent(in)      :: vel2_bar
        type(AD_D),             intent(in)      :: vel3_bar
        type(AD_D),             intent(in)      :: c_bar
        real(rk),               intent(inout)   :: kz
        type(AD_D),             intent(inout)   :: k1
        type(AD_D),             intent(inout)   :: k2
        type(AD_D),             intent(inout)   :: k3
        type(AD_D),             intent(inout)   :: k4_real
        type(AD_D),             intent(inout)   :: k4_imag
        type(AD_D),             intent(inout)   :: k5_real
        type(AD_D),             intent(inout)   :: k5_imag

        real(rk),   allocatable, dimension(:)   :: pitch, unorm3
        type(AD_D)  :: pyra, ktmp_real, ktmp_imag, vg4, vg5

        complex(rk) :: omega_c, pyra_c, k4_c, k5_c
        real(rk)    :: vel2_bar_r, vel3_bar_r, c_bar_r


        pitch  = self%bcproperties%compute('Pitch',worker%time(),worker%coords())

        kz = TWO*PI*lm/pitch(1)
        pyra = (omega - kz*vel2_bar)**TWO - kz*kz*(c_bar**TWO - vel3_bar**TWO)

        ! Compute k1,k2,k3
        k1 = (omega - kz*vel2_bar)/vel3_bar
        k2 = (omega - kz*vel2_bar)/vel3_bar
        k3 = (omega - kz*vel2_bar)/vel3_bar

        if (pyra >= ZERO) then
            k4_real = (-vel3_bar*(omega - kz*vel2_bar) + c_bar*sqrt(pyra))/(c_bar**TWO - vel3_bar**TWO)
            k5_real = (-vel3_bar*(omega - kz*vel2_bar) - c_bar*sqrt(pyra))/(c_bar**TWO - vel3_bar**TWO)
            k4_imag = ZERO*k4_real
            k5_imag = ZERO*k4_real

!            ! test direction of propagation by perturbing omega
!            vel3_bar_r = vel3_bar%x_ad_
!            vel2_bar_r = vel2_bar%x_ad_
!            c_bar_r = c_bar%x_ad_
!            omega_c = cmplx(omega,-0.00001_rk)
!            pyra_c = (omega_c - kz*vel2_bar_r)**TWO - kz*kz*(c_bar_r**TWO - vel3_bar_r**TWO)
!            k4_c = (-vel3_bar_r*(omega_c - kz*vel2_bar_r) + c_bar_r*sqrt(pyra_c))/(c_bar_r**TWO - vel3_bar_r**TWO)
!            k5_c = (-vel3_bar_r*(omega_c - kz*vel2_bar_r) - c_bar_r*sqrt(pyra_c))/(c_bar_r**TWO - vel3_bar_r**TWO)
!
!            if (imagpart(k4_c) < ZERO) then
!                k4_imag = -0.0000001_rk
!            else
!                k4_imag = 0.0000001_rk
!            end if
!
!            if (imagpart(k5_c) < ZERO) then
!                k5_imag = -0.0000001_rk
!            else
!                k5_imag = 0.0000001_rk
!            end if


            vg4 = -(c_bar*c_bar - vel3_bar*vel3_bar)/(vel3_bar - (c_bar*(omega-kz*vel2_bar)/sqrt(pyra)))
            vg5 = -(c_bar*c_bar - vel3_bar*vel3_bar)/(vel3_bar + (c_bar*(omega-kz*vel2_bar)/sqrt(pyra)))


        else
            k4_real = (-vel3_bar*(omega - kz*vel2_bar))/(c_bar**TWO - vel3_bar**TWO)
            k4_imag =  c_bar*sqrt(-pyra)/(c_bar**TWO - vel3_bar**TWO)
            k5_real = (-vel3_bar*(omega - kz*vel2_bar))/(c_bar**TWO - vel3_bar**TWO)
            k5_imag = -c_bar*sqrt(-pyra)/(c_bar**TWO - vel3_bar**TWO)

!            ! test direction of propagation by perturbing omega
!            vel3_bar_r = vel3_bar%x_ad_
!            vel2_bar_r = vel2_bar%x_ad_
!            c_bar_r = c_bar%x_ad_
!            omega_c = cmplx(omega,-0.00001_rk)
!            pyra_c = (omega_c - kz*vel2_bar_r)**TWO - kz*kz*(c_bar_r**TWO - vel3_bar_r**TWO)
!            k4_c = (-vel3_bar_r*(omega_c - kz*vel2_bar_r) + c_bar_r*sqrt(pyra_c))/(c_bar_r**TWO - vel3_bar_r**TWO)
!            k5_c = (-vel3_bar_r*(omega_c - kz*vel2_bar_r) - c_bar_r*sqrt(pyra_c))/(c_bar_r**TWO - vel3_bar_r**TWO)
!
!            if ( ((imagpart(k4_c) < ZERO) .and. (k4_imag > ZERO)) .or. &
!                 ((imagpart(k4_c) > ZERO) .and. (k4_imag < ZERO)) ) then
!                k4_imag = -k4_imag
!            end if
!
!            if ( ((imagpart(k5_c) < ZERO) .and. (k5_imag > ZERO)) .or. &
!                 ((imagpart(k5_c) > ZERO) .and. (k5_imag < ZERO)) ) then
!                k5_imag = -k5_imag
!            end if

            vg4 = ZERO*k4_real
            vg5 = ZERO*k4_real
            if (k4_imag < ZERO) then
                vg4 = ONE
            else
                vg4 = -ONE
            end if

            if (k5_imag < ZERO) then
                vg5 = ONE
            else
                vg5 = -ONE
            end if

        end if

!        unorm3 = worker%unit_normal(3)
!        if (k5_imag*unorm3(1) < ZERO) then
!            ktmp_real = k5_real
!            ktmp_imag = k5_imag
!
!            k5_real = k4_real
!            k5_imag = k4_imag
!            k4_real = ktmp_real
!            k4_imag = ktmp_imag
!        end if


        unorm3 = worker%unit_normal(3)
        if (unorm3(1) < ZERO) then
            if (vg5 < ZERO) then
                ktmp_real = k5_real
                ktmp_imag = k5_imag

                k5_real = k4_real
                k5_imag = k4_imag
                k4_real = ktmp_real
                k4_imag = ktmp_imag
            end if
        end if

        
        if (unorm3(1) > ZERO) then
            if (vg5 > ZERO) then
                ktmp_real = k5_real
                ktmp_imag = k5_imag

                k5_real = k4_real
                k5_imag = k4_imag
                k4_real = ktmp_real
                k4_imag = ktmp_imag
            end if
        end if

    end subroutine compute_eigenvalues
    !********************************************************************************




    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/25/2018
    !!
    !--------------------------------------------------------------------------------
    subroutine primitive_to_characteristics(self,worker,bc_comm,                    &
                                            density_Fts_real,  density_Fts_imag,    &
                                            vel1_Fts_real,     vel1_Fts_imag,       &
                                            vel2_Fts_real,     vel2_Fts_imag,       &
                                            vel3_Fts_real,     vel3_Fts_imag,       &
                                            pressure_Fts_real, pressure_Fts_imag,   &
                                            c1_hat_real,       c1_hat_imag,         &
                                            c2_hat_real,       c2_hat_imag,         &
                                            c3_hat_real,       c3_hat_imag,         &
                                            c4_hat_real,       c4_hat_imag,         &
                                            c5_hat_real,       c5_hat_imag)
        class(giles_HB_base_t),         intent(inout)   :: self
        type(chidg_worker_t),           intent(inout)   :: worker
        type(mpi_comm),                 intent(in)      :: bc_comm
        type(AD_D),                     intent(in)      :: density_Fts_real(:,:,:)
        type(AD_D),                     intent(in)      :: density_Fts_imag(:,:,:)
        type(AD_D),                     intent(in)      :: vel1_Fts_real(:,:,:)
        type(AD_D),                     intent(in)      :: vel1_Fts_imag(:,:,:)
        type(AD_D),                     intent(in)      :: vel2_Fts_real(:,:,:)
        type(AD_D),                     intent(in)      :: vel2_Fts_imag(:,:,:)
        type(AD_D),                     intent(in)      :: vel3_Fts_real(:,:,:)
        type(AD_D),                     intent(in)      :: vel3_Fts_imag(:,:,:)
        type(AD_D),                     intent(in)      :: pressure_Fts_real(:,:,:)
        type(AD_D),                     intent(in)      :: pressure_Fts_imag(:,:,:)
        type(AD_D),     allocatable,    intent(inout)   :: c1_hat_real(:,:,:)
        type(AD_D),     allocatable,    intent(inout)   :: c1_hat_imag(:,:,:)
        type(AD_D),     allocatable,    intent(inout)   :: c2_hat_real(:,:,:)
        type(AD_D),     allocatable,    intent(inout)   :: c2_hat_imag(:,:,:)
        type(AD_D),     allocatable,    intent(inout)   :: c3_hat_real(:,:,:)
        type(AD_D),     allocatable,    intent(inout)   :: c3_hat_imag(:,:,:)
        type(AD_D),     allocatable,    intent(inout)   :: c4_hat_real(:,:,:)
        type(AD_D),     allocatable,    intent(inout)   :: c4_hat_imag(:,:,:)
        type(AD_D),     allocatable,    intent(inout)   :: c5_hat_real(:,:,:)
        type(AD_D),     allocatable,    intent(inout)   :: c5_hat_imag(:,:,:)

        type(AD_D), allocatable,    dimension(:)    ::                                          &
            density, mom1, mom2, mom3, energy, vel1, vel2, vel3, pressure,                      &
            c1_real_tmp,      c2_real_tmp,   c3_real_tmp,   c4_real_tmp,   c5_real_tmp,         &
            c1_imag_tmp,      c2_imag_tmp,   c3_imag_tmp,   c4_imag_tmp,   c5_imag_tmp,         &
            density_bar, vel1_bar, vel2_bar, vel3_bar, pressure_bar

        type(AD_D)  :: c_bar

        integer(ik)             :: nmodes, nradius, iradius, itheta, ntheta, imode, ierr, ntime, itime
        real(rk)                :: shift_r, shift_i
        real(rk),   allocatable :: pitch(:)

        ! Get spatio-temporal average at radial stations
        density_bar  = density_Fts_real(:,1,1)
        vel1_bar     = vel1_Fts_real(:,1,1)
        vel2_bar     = vel2_Fts_real(:,1,1)
        vel3_bar     = vel3_Fts_real(:,1,1)
        pressure_bar = pressure_Fts_real(:,1,1)

        ! Define Fourier discretization
        nmodes  = self%nfourier_space
        ntheta  = 1 + (nmodes-1)*2
        nradius = size(self%r)
        ntime = worker%time_manager%ntime

        pitch  = self%bcproperties%compute('Pitch',time=ZERO,coord=[point_t(ZERO,ZERO,ZERO)])

        ! Allocate storage in result
        allocate(c1_hat_real(nradius,ntheta,ntime), c1_hat_imag(nradius,ntheta,ntime),  &
                 c2_hat_real(nradius,ntheta,ntime), c2_hat_imag(nradius,ntheta,ntime),  &
                 c3_hat_real(nradius,ntheta,ntime), c3_hat_imag(nradius,ntheta,ntime),  &
                 c4_hat_real(nradius,ntheta,ntime), c4_hat_imag(nradius,ntheta,ntime),  &
                 c5_hat_real(nradius,ntheta,ntime), c5_hat_imag(nradius,ntheta,ntime), stat=ierr)
        if (ierr /= 0) call AllocationError
        c1_hat_real = ZERO*density_bar(1)
        c2_hat_real = ZERO*density_bar(1)
        c3_hat_real = ZERO*density_bar(1)
        c4_hat_real = ZERO*density_bar(1)
        c5_hat_real = ZERO*density_bar(1)
        c1_hat_imag = ZERO*density_bar(1)
        c2_hat_imag = ZERO*density_bar(1)
        c3_hat_imag = ZERO*density_bar(1)
        c4_hat_imag = ZERO*density_bar(1)
        c5_hat_imag = ZERO*density_bar(1)

        ! Convert Fourier modes of primitive varibles to 1D characteristics
        do iradius = 1,nradius
            c_bar = sqrt(gam*pressure_bar(iradius)/density_bar(iradius))
            do itheta = 1,ntheta
                do itime = 1,ntime
                    c1_hat_real(iradius,itheta,itime) = -(c_bar*c_bar)*density_Fts_real(iradius,itheta,itime)             +  (ONE)*pressure_Fts_real(iradius,itheta,itime)
                    c2_hat_real(iradius,itheta,itime) = (density_bar(iradius)*c_bar)*vel1_Fts_real(iradius,itheta,itime)
                    c3_hat_real(iradius,itheta,itime) = (density_bar(iradius)*c_bar)*vel2_Fts_real(iradius,itheta,itime)
                    c4_hat_real(iradius,itheta,itime) = (density_bar(iradius)*c_bar)*vel3_Fts_real(iradius,itheta,itime)  +  (ONE)*pressure_Fts_real(iradius,itheta,itime)
                    c5_hat_real(iradius,itheta,itime) = -(density_bar(iradius)*c_bar)*vel3_Fts_real(iradius,itheta,itime) +  (ONE)*pressure_Fts_real(iradius,itheta,itime)
                                               
                    c1_hat_imag(iradius,itheta,itime) = -(c_bar*c_bar)*density_Fts_imag(iradius,itheta,itime)             +  (ONE)*pressure_Fts_imag(iradius,itheta,itime)
                    c2_hat_imag(iradius,itheta,itime) = (density_bar(iradius)*c_bar)*vel1_Fts_imag(iradius,itheta,itime)
                    c3_hat_imag(iradius,itheta,itime) = (density_bar(iradius)*c_bar)*vel2_Fts_imag(iradius,itheta,itime)
                    c4_hat_imag(iradius,itheta,itime) = (density_bar(iradius)*c_bar)*vel3_Fts_imag(iradius,itheta,itime)  +  (ONE)*pressure_Fts_imag(iradius,itheta,itime)
                    c5_hat_imag(iradius,itheta,itime) = -(density_bar(iradius)*c_bar)*vel3_Fts_imag(iradius,itheta,itime) +  (ONE)*pressure_Fts_imag(iradius,itheta,itime)
                end do !itime
            end do !itheta
        end do !iradius

    end subroutine primitive_to_characteristics
    !*********************************************************************************





    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/25/2018
    !!
    !--------------------------------------------------------------------------------
    subroutine characteristics_to_primitive(self,worker,bc_comm,                  &
                                            c1_hat_real,       c1_hat_imag,       &
                                            c2_hat_real,       c2_hat_imag,       &
                                            c3_hat_real,       c3_hat_imag,       &
                                            c4_hat_real,       c4_hat_imag,       &
                                            c5_hat_real,       c5_hat_imag,       &
                                            density_Fts_real,  density_Fts_imag,  &
                                            vel1_Fts_real,     vel1_Fts_imag,     &
                                            vel2_Fts_real,     vel2_Fts_imag,     &
                                            vel3_Fts_real,     vel3_Fts_imag,     &
                                            pressure_Fts_real, pressure_Fts_imag)
        class(giles_HB_base_t),         intent(inout)   :: self
        type(chidg_worker_t),           intent(inout)   :: worker
        type(mpi_comm),                 intent(in)      :: bc_comm
        type(AD_D),                     intent(in)      :: c1_hat_real(:,:,:)
        type(AD_D),                     intent(in)      :: c1_hat_imag(:,:,:)
        type(AD_D),                     intent(in)      :: c2_hat_real(:,:,:)
        type(AD_D),                     intent(in)      :: c2_hat_imag(:,:,:)
        type(AD_D),                     intent(in)      :: c3_hat_real(:,:,:)
        type(AD_D),                     intent(in)      :: c3_hat_imag(:,:,:)
        type(AD_D),                     intent(in)      :: c4_hat_real(:,:,:)
        type(AD_D),                     intent(in)      :: c4_hat_imag(:,:,:)
        type(AD_D),                     intent(in)      :: c5_hat_real(:,:,:)
        type(AD_D),                     intent(in)      :: c5_hat_imag(:,:,:)
        type(AD_D),                     intent(inout)   :: density_Fts_real(:,:,:)
        type(AD_D),                     intent(inout)   :: density_Fts_imag(:,:,:)
        type(AD_D),                     intent(inout)   :: vel1_Fts_real(:,:,:)
        type(AD_D),                     intent(inout)   :: vel1_Fts_imag(:,:,:)
        type(AD_D),                     intent(inout)   :: vel2_Fts_real(:,:,:)
        type(AD_D),                     intent(inout)   :: vel2_Fts_imag(:,:,:)
        type(AD_D),                     intent(inout)   :: vel3_Fts_real(:,:,:)
        type(AD_D),                     intent(inout)   :: vel3_Fts_imag(:,:,:)
        type(AD_D),                     intent(inout)   :: pressure_Fts_real(:,:,:)
        type(AD_D),                     intent(inout)   :: pressure_Fts_imag(:,:,:)

        type(AD_D), allocatable,    dimension(:)    ::  &
            density_bar, vel1_bar, vel2_bar, vel3_bar, pressure_bar
            
        type(AD_D)  :: density_bar_r, pressure_bar_r, c_bar_r

        integer(ik)             :: iradius, itheta, itime

        ! Get spatio-temporal average at radial stations
        density_bar  = density_Fts_real(:,1,1)
        vel1_bar     = vel1_Fts_real(:,1,1)
        vel2_bar     = vel2_Fts_real(:,1,1)
        vel3_bar     = vel3_Fts_real(:,1,1)
        pressure_bar = pressure_Fts_real(:,1,1)

        ! Convert characteristic Fourier modes back to primitive Fourier modes and store
        do iradius = 1,size(c1_hat_real,1)
            ! Get radius-local averages
            density_bar_r  = density_bar(iradius)
            pressure_bar_r = pressure_bar(iradius)
            c_bar_r        = sqrt(gam*pressure_bar_r/density_bar_r)
            do itheta = 1,size(c1_hat_real,2)
                do itime = 1,size(c1_hat_real,3)
                    density_Fts_real(iradius,itheta,itime) = -(ONE/(c_bar_r*c_bar_r))*c1_hat_real(iradius,itheta,itime) + (ONE/(TWO*c_bar_r*c_bar_r))*c4_hat_real(iradius,itheta,itime) + (ONE/(TWO*c_bar_r*c_bar_r))*c5_hat_real(iradius,itheta,itime)
                    vel1_Fts_real(iradius,itheta,itime) = (ONE/(density_bar_r*c_bar_r))*c2_hat_real(iradius,itheta,itime)
                    vel2_Fts_real(iradius,itheta,itime) = (ONE/(density_bar_r*c_bar_r))*c3_hat_real(iradius,itheta,itime)
                    vel3_Fts_real(iradius,itheta,itime) = (ONE/(TWO*density_bar_r*c_bar_r))*c4_hat_real(iradius,itheta,itime) - (ONE/(TWO*density_bar_r*c_bar_r))*c5_hat_real(iradius,itheta,itime)
                    pressure_Fts_real(iradius,itheta,itime) = HALF*c4_hat_real(iradius,itheta,itime) + HALF*c5_hat_real(iradius,itheta,itime)

                    density_Fts_imag(iradius,itheta,itime) = -(ONE/(c_bar_r*c_bar_r))*c1_hat_imag(iradius,itheta,itime) + (ONE/(TWO*c_bar_r*c_bar_r))*c4_hat_imag(iradius,itheta,itime) + (ONE/(TWO*c_bar_r*c_bar_r))*c5_hat_imag(iradius,itheta,itime)
                    vel1_Fts_imag(iradius,itheta,itime) = (ONE/(density_bar_r*c_bar_r))*c2_hat_imag(iradius,itheta,itime)
                    vel2_Fts_imag(iradius,itheta,itime) = (ONE/(density_bar_r*c_bar_r))*c3_hat_imag(iradius,itheta,itime)
                    vel3_Fts_imag(iradius,itheta,itime) = (ONE/(TWO*density_bar_r*c_bar_r))*c4_hat_imag(iradius,itheta,itime) - (ONE/(TWO*density_bar_r*c_bar_r))*c5_hat_imag(iradius,itheta,itime)
                    pressure_Fts_imag(iradius,itheta,itime) = HALF*c4_hat_imag(iradius,itheta,itime) + HALF*c5_hat_imag(iradius,itheta,itime)
                end do
            end do
        end do 

    end subroutine characteristics_to_primitive
    !*********************************************************************************











    !> Compute boundary average by averaging spatio-temporal average over
    !! radial stations.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/16/2018
    !!
    !--------------------------------------------------------------------------------
    subroutine compute_boundary_average(self,worker,bc_comm,density_bar,vel1_bar,vel2_bar,vel3_bar,pressure_bar,c_bar, &
                                                            density_avg,vel1_avg,vel2_avg,vel3_avg,pressure_avg,c_avg)
        class(giles_HB_base_t),  intent(inout)   :: self
        type(chidg_worker_t),                       intent(in)      :: worker
        type(mpi_comm),                             intent(in)      :: bc_comm
        type(AD_D),                                 intent(in)      :: density_bar(:)
        type(AD_D),                                 intent(in)      :: vel1_bar(:)
        type(AD_D),                                 intent(in)      :: vel2_bar(:)
        type(AD_D),                                 intent(in)      :: vel3_bar(:)
        type(AD_D),                                 intent(in)      :: pressure_bar(:)
        type(AD_D),                                 intent(in)      :: c_bar(:)
        type(AD_D),                                 intent(inout)   :: density_avg
        type(AD_D),                                 intent(inout)   :: vel1_avg
        type(AD_D),                                 intent(inout)   :: vel2_avg
        type(AD_D),                                 intent(inout)   :: vel3_avg
        type(AD_D),                                 intent(inout)   :: pressure_avg
        type(AD_D),                                 intent(inout)   :: c_avg

        real(rk)    :: area, dr
        integer(ik) :: irad


        density_avg  = ZERO*density_bar(1)
        vel1_avg     = ZERO*density_bar(1)
        vel2_avg     = ZERO*density_bar(1)
        vel3_avg     = ZERO*density_bar(1)
        pressure_avg = ZERO*density_bar(1)
        c_avg        = ZERO*density_bar(1)

        if (worker%coordinate_system() == 'Cartesian') then
            dr = self%r(2) - self%r(1)
            area = ZERO
            do irad = 1,size(density_bar)-1
                density_avg  = density_avg  + dr*(density_bar(irad+1) +density_bar(irad))/TWO
                vel1_avg     = vel1_avg     + dr*(vel1_bar(irad+1)    +vel1_bar(irad))/TWO
                vel2_avg     = vel2_avg     + dr*(vel2_bar(irad+1)    +vel2_bar(irad))/TWO
                vel3_avg     = vel3_avg     + dr*(vel3_bar(irad+1)    +vel3_bar(irad))/TWO
                pressure_avg = pressure_avg + dr*(pressure_bar(irad+1)+pressure_bar(irad))/TWO
                c_avg        = c_avg        + dr*(c_bar(irad+1)       +c_bar(irad))/TWO
                area = area + dr
            end do

        else if (worker%coordinate_system() == 'Cylindrical') then
            dr = self%r(2) - self%r(1)
            area = ZERO
            do irad = 1,size(density_bar)-1
                density_avg  = density_avg  + dr*(self%r(irad+1)*density_bar(irad+1)  + self%r(irad)*density_bar(irad))/TWO
                vel1_avg     = vel1_avg     + dr*(self%r(irad+1)*vel1_bar(irad+1)     + self%r(irad)*vel1_bar(irad))/TWO
                vel2_avg     = vel2_avg     + dr*(self%r(irad+1)*vel2_bar(irad+1)     + self%r(irad)*vel2_bar(irad))/TWO
                vel3_avg     = vel3_avg     + dr*(self%r(irad+1)*vel3_bar(irad+1)     + self%r(irad)*vel3_bar(irad))/TWO
                pressure_avg = pressure_avg + dr*(self%r(irad+1)*pressure_bar(irad+1) + self%r(irad)*pressure_bar(irad))/TWO
                c_avg        = c_avg        + dr*(self%r(irad+1)*c_bar(irad+1)        + self%r(irad)*c_bar(irad))/TWO
                area = area + dr*(self%r(irad+1)+self%r(irad))/TWO
            end do

        end if

        density_avg  = density_avg  / area
        vel1_avg     = vel1_avg     / area
        vel2_avg     = vel2_avg     / area
        vel3_avg     = vel3_avg     / area
        pressure_avg = pressure_avg / area
        c_avg        = c_avg        / area

    end subroutine compute_boundary_average
    !********************************************************************************





    !>  Determine rmin, rmax, and average theta.
    !!
    !!  Initialize radial stations and reference theta for Fourier transform.
    !!  
    !!  Defines:
    !!      self%r(:)
    !!      self%theta_ref
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/26/2018
    !!
    !--------------------------------------------------------------------------------
    subroutine analyze_bc_geometry(self,mesh,group_ID,bc_comm)
        class(giles_HB_base_t),  intent(inout)   :: self
        type(mesh_t),                               intent(in)      :: mesh
        integer(ik),                                intent(in)      :: group_ID
        type(mpi_comm),                             intent(in)      :: bc_comm

        integer(ik) :: bc_idomain_l, bc_ielement_l, patch_ID, face_ID, iface, ierr, inode
        real(rk)    :: face_rmin, face_rmax, local_rmin, local_rmax, global_rmin, global_rmax,  &
                       face_thetamin, face_thetamax, local_thetamin, local_thetamax, global_thetamin, global_thetamax
        real(rk)    :: ref_nodes(8,3), physical_nodes(8,3)

        ! Search for min/max radius on local processor
        local_rmin     =  HUGE(1._rk) ! any radius will be smaller than this, so it is guarunteed to be reset.
        local_rmax     = -HUGE(1._rk) ! any radius will be larger than this, so it is guarunteed to be reset.
        local_thetamin =  HUGE(1._rk) ! any theta will be smaller than this, so it is guarunteed to be reset.
        local_thetamax = -HUGE(1._rk) ! any theta will be larger than this, so it is guarunteed to be reset.
        do patch_ID = 1,mesh%bc_patch_group(group_ID)%npatches()
            do face_ID = 1,mesh%bc_patch_group(group_ID)%patch(patch_ID)%nfaces()

                iface = mesh%bc_patch_group(group_ID)%patch(patch_ID)%iface(face_ID)

                ! Pick points to evaluate coordinates
                if (iface == XI_MIN) then
                    ref_nodes(1,:) = [-ONE, -ONE, -ONE]
                    ref_nodes(2,:) = [-ONE, -ONE, ZERO]
                    ref_nodes(3,:) = [-ONE, -ONE,  ONE]
                    ref_nodes(4,:) = [-ONE, ZERO,  ONE]
                    ref_nodes(5,:) = [-ONE,  ONE,  ONE]
                    ref_nodes(6,:) = [-ONE,  ONE, ZERO]
                    ref_nodes(7,:) = [-ONE,  ONE, -ONE]
                    ref_nodes(8,:) = [-ONE, ZERO, -ONE]
                else if (iface == XI_MAX) then
                    ref_nodes(1,:) = [ONE, -ONE, -ONE]
                    ref_nodes(2,:) = [ONE, -ONE, ZERO]
                    ref_nodes(3,:) = [ONE, -ONE,  ONE]
                    ref_nodes(4,:) = [ONE, ZERO,  ONE]
                    ref_nodes(5,:) = [ONE,  ONE,  ONE]
                    ref_nodes(6,:) = [ONE,  ONE, ZERO]
                    ref_nodes(7,:) = [ONE,  ONE, -ONE]
                    ref_nodes(8,:) = [ONE, ZERO, -ONE]

                else if (iface == ETA_MIN) then
                    ref_nodes(1,:) = [ -ONE, -ONE, -ONE]
                    ref_nodes(2,:) = [ -ONE, -ONE, ZERO]
                    ref_nodes(3,:) = [ -ONE, -ONE,  ONE]
                    ref_nodes(4,:) = [ ZERO, -ONE,  ONE]
                    ref_nodes(5,:) = [  ONE, -ONE,  ONE]
                    ref_nodes(6,:) = [  ONE, -ONE, ZERO]
                    ref_nodes(7,:) = [  ONE, -ONE, -ONE]
                    ref_nodes(8,:) = [ ZERO, -ONE, -ONE]

                else if (iface == ETA_MAX) then
                    ref_nodes(1,:) = [ -ONE, ONE, -ONE]
                    ref_nodes(2,:) = [ -ONE, ONE, ZERO]
                    ref_nodes(3,:) = [ -ONE, ONE,  ONE]
                    ref_nodes(4,:) = [ ZERO, ONE,  ONE]
                    ref_nodes(5,:) = [  ONE, ONE,  ONE]
                    ref_nodes(6,:) = [  ONE, ONE, ZERO]
                    ref_nodes(7,:) = [  ONE, ONE, -ONE]
                    ref_nodes(8,:) = [ ZERO, ONE, -ONE]

                else if (iface == ZETA_MIN) then
                    ref_nodes(1,:) = [ -ONE, -ONE, -ONE]
                    ref_nodes(2,:) = [ -ONE, ZERO, -ONE]
                    ref_nodes(3,:) = [ -ONE,  ONE, -ONE]
                    ref_nodes(4,:) = [ ZERO,  ONE, -ONE]
                    ref_nodes(5,:) = [  ONE,  ONE, -ONE]
                    ref_nodes(6,:) = [  ONE, ZERO, -ONE]
                    ref_nodes(7,:) = [  ONE, -ONE, -ONE]
                    ref_nodes(8,:) = [ ZERO, -ONE, -ONE]

                else if (iface == ZETA_MAX) then
                    ref_nodes(1,:) = [ -ONE, -ONE, ONE]
                    ref_nodes(2,:) = [ -ONE, ZERO, ONE]
                    ref_nodes(3,:) = [ -ONE,  ONE, ONE]
                    ref_nodes(4,:) = [ ZERO,  ONE, ONE]
                    ref_nodes(5,:) = [  ONE,  ONE, ONE]
                    ref_nodes(6,:) = [  ONE, ZERO, ONE]
                    ref_nodes(7,:) = [  ONE, -ONE, ONE]
                    ref_nodes(8,:) = [ ZERO, -ONE, ONE]

                else
                    call chidg_signal(FATAL,"outlet_giles_quasi3d_unsteady_HB: analyze_bc_geometry, invalid face indec.")
                end if

                ! Evaluate physical coordinates on face edges
                bc_idomain_l  = mesh%bc_patch_group(group_ID)%patch(patch_ID)%idomain_l()
                bc_ielement_l = mesh%bc_patch_group(group_ID)%patch(patch_ID)%ielement_l(face_ID)
                do inode = 1,size(ref_nodes,1)
                    physical_nodes(inode,:) = mesh%domain(bc_idomain_l)%elems(bc_ielement_l)%physical_point(ref_nodes(inode,:),'Deformed')
                end do !inode

                ! Get face min/max radius
                face_rmin     = minval(physical_nodes(:,1))
                face_rmax     = maxval(physical_nodes(:,1))
                face_thetamin = minval(physical_nodes(:,2))
                face_thetamax = maxval(physical_nodes(:,2))

                ! Update processor-local value if new min/max values were found on the face
                if (face_rmin < local_rmin)         local_rmin     = face_rmin
                if (face_rmax > local_rmax)         local_rmax     = face_rmax
                if (face_thetamin < local_thetamin) local_thetamin = face_thetamin
                if (face_thetamax > local_thetamax) local_thetamax = face_thetamax

            end do !face_ID
        end do !patch_ID

        ! Reduce processor local values to determine boundary-global min/max values
        call MPI_AllReduce(local_rmin,    global_rmin,    1,MPI_REAL8,MPI_MIN,bc_comm,ierr)
        call MPI_AllReduce(local_rmax,    global_rmax,    1,MPI_REAL8,MPI_MAX,bc_comm,ierr)
        call MPI_AllReduce(local_thetamin,global_thetamin,1,MPI_REAL8,MPI_MIN,bc_comm,ierr)
        call MPI_AllReduce(local_thetamax,global_thetamax,1,MPI_REAL8,MPI_MAX,bc_comm,ierr)
        
        ! Create radial stations
        self%r = linspace(global_rmin,global_rmax,self%nr)

        ! Compute theta_ref
        self%theta_ref = (global_thetamin + global_thetamax)/TWO

    end subroutine analyze_bc_geometry
    !********************************************************************************




    !>  For each radial station, a theta-discretization exists upon which the 
    !!  discrete Fourier transform operation is computed. Since the theta-discretization
    !!  spans the entire boundary, a general interpolation procedure is used to fill
    !!  solution values that could come from many different elements on the boundary.
    !!
    !!  This procedure finds donor elements for each node in the theta-discretization
    !!  at each radial station. Additionally, the location of the physical coordinate
    !!  within the donor element's reference coordinate system is also found.
    !!
    !!  This procedure
    !!
    !!  Initializes:
    !!  ------------------------------
    !!  self%theta(:)
    !!  self%donor(:,:)
    !!  self%donor_node(:,:,:)
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/28/2018
    !!
    !----------------------------------------------------------------------------------
    subroutine initialize_fourier_discretization(self,mesh,group_ID,bc_comm)
        class(giles_HB_base_t),  intent(inout)   :: self
        type(mesh_t),                               intent(in)      :: mesh
        integer(ik),                                intent(in)      :: group_ID
        type(mpi_comm),                             intent(in)      :: bc_comm

        integer(ik)                 :: nmodes, ncoeff, nradius, ntheta, idomain_l, ielement_l, iface, &
                                       iradius, itheta, ierr, noverset
        real(rk)                    :: dtheta, dtheta_n, midpoint(3), try_offset(3), node(3), z
        real(rk),       allocatable :: pitch(:)
        character(:),   allocatable :: user_msg
        logical                     :: donor_found


        ! Determine z-location of some face on the boundary and assume 
        ! entire boundary is constant-z
        idomain_l  = mesh%bc_patch_group(group_ID)%patch(1)%idomain_l()
        ielement_l = mesh%bc_patch_group(group_ID)%patch(1)%ielement_l(1)
        iface      = mesh%bc_patch_group(group_ID)%patch(1)%iface(1)
        if (iface == XI_MIN) then
            midpoint = mesh%domain(idomain_l)%elems(ielement_l)%physical_point([-ONE,ZERO,ZERO],'Deformed')
        else if (iface == XI_MAX) then
            midpoint = mesh%domain(idomain_l)%elems(ielement_l)%physical_point([ ONE,ZERO,ZERO],'Deformed')
        else if (iface == ETA_MIN) then
            midpoint = mesh%domain(idomain_l)%elems(ielement_l)%physical_point([ZERO,-ONE,ZERO],'Deformed')
        else if (iface == ETA_MAX) then
            midpoint = mesh%domain(idomain_l)%elems(ielement_l)%physical_point([ZERO, ONE,ZERO],'Deformed')
        else if (iface == ZETA_MIN) then
            midpoint = mesh%domain(idomain_l)%elems(ielement_l)%physical_point([ZERO,ZERO,-ONE],'Deformed')
        else if (iface == ZETA_MAX) then
            midpoint = mesh%domain(idomain_l)%elems(ielement_l)%physical_point([ZERO,ZERO, ONE],'Deformed')
        end if
        z = midpoint(3)

        ! Define Fourier discretization
        nmodes  = self%nfourier_space
        ncoeff  = 1 + (nmodes-1)*2
        nradius = size(self%r)
        ntheta  = ncoeff
        
        ! Initialize theta discretization parameters
        pitch  = self%bcproperties%compute('Pitch',time=ZERO,coord=[point_t(ZERO,ZERO,ZERO)])
        dtheta = pitch(1)
        dtheta_n = dtheta/ntheta

        ! Construct theta discretization at each radius
        allocate(self%theta(size(self%r),ntheta), stat=ierr)
        if (ierr /= 0) call AllocationError
        do itheta = 1,ntheta
            self%theta(:,itheta) = self%theta_ref + (itheta-1)*dtheta_n
        end do

        ! Donor search offset, if needed
        try_offset = [ZERO, -pitch(1), ZERO]

        ! For each radial station, initialized donor for each node in theta grid
        allocate(self%donor(nradius,ntheta), self%donor_node(nradius,ntheta,3), stat=ierr)
        if (ierr /= 0) call AllocationError
        do iradius = 1,size(self%r)
            noverset = 0
            do itheta = 1,ntheta

                node = [self%r(iradius), self%theta(iradius,itheta), z]

                ! Try processor-LOCAL elements
                call find_gq_donor(mesh,                                    &
                                   node,                                    &
                                   [ZERO,ZERO,ZERO],                        &
                                   face_info_constructor(0,0,0,0,0),        &   ! we don't really have a receiver face
                                   self%donor(iradius,itheta),              &
                                   self%donor_node(iradius,itheta,1:3),     &
                                   donor_found)

                ! Try LOCAL elements with try_offset if still not found 
                if ( .not. donor_found ) then
                    call find_gq_donor(mesh,                                    &
                                       node,                                    &
                                       try_offset,                              &
                                       face_info_constructor(0,0,0,0,0),        &   ! we don't really have a receiver face
                                       self%donor(iradius,itheta),              &
                                       self%donor_node(iradius,itheta,1:3),     &
                                       donor_found)
                    if (donor_found) then
                        noverset=noverset+1
                        self%theta(iradius,itheta) = self%theta(iradius,itheta) - pitch(1)
                    end if
                end if

                ! Try PARALLEL_ELEMENTS if donor not found amongst local elements
                if (.not. donor_found) then
                    call find_gq_donor_parallel(mesh,                                   &
                                                node,                                   &
                                                [ZERO,ZERO,ZERO],                       &
                                                face_info_constructor(0,0,0,0,0),       &   ! we don't really have a receiver face
                                                self%donor(iradius,itheta),             &
                                                self%donor_node(iradius,itheta,1:3),    &
                                                donor_found)
                end if

                
                ! Try PARALLEL_ELEMENTS with try_offset if still not found 
                if ( .not. donor_found ) then
                    call find_gq_donor_parallel(mesh,                                   &
                                                node,                                   &
                                                try_offset,                             &
                                                face_info_constructor(0,0,0,0,0),       &   ! we don't really have a receiver face
                                                self%donor(iradius,itheta),             &
                                                self%donor_node(iradius,itheta,1:3),    &
                                                donor_found)
                    if (donor_found) then
                        noverset=noverset+1
                        self%theta(iradius,itheta) = self%theta(iradius,itheta) - pitch(1)
                    end if
                end if 


                ! Abort if we didn't find a donor
                user_msg = "bc_giles_HB_base%initialize_fourier_discretization: &
                            no donor element found for Fourier discretization node."
                if (.not. donor_found) call chidg_signal(FATAL,user_msg)

            end do !itheta

            ! Shift arrays so that we start with the theta_min point
            self%donor(iradius,:)        = cshift(self%donor(iradius,:), -noverset, dim=1)
            self%donor_node(iradius,:,:) = cshift(self%donor_node(iradius,:,:), -noverset, dim=1)
            self%theta(iradius,:)        = cshift(self%theta(iradius,:), -noverset, dim=1)

        end do !iradius

    end subroutine initialize_fourier_discretization
    !**************************************************************************************





    !>  Return the correct temporal frequency for the time index.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   5/28/2018
    !!
    !--------------------------------------------------------------------------------------
    function get_omega(worker,itime) result(omega)
        type(chidg_worker_t),   intent(in)  :: worker
        integer(ik),            intent(in)  :: itime

        real(rk) :: omega
        integer(ik) :: ntime

        ntime = worker%time_manager%ntime

        if (itime == 1) then
            omega = ZERO
        else if (itime <= ((ntime-1)/2 + 1)) then
            omega = worker%time_manager%freqs(itime-1)
        else
            omega = -worker%time_manager%freqs(ntime-itime+1)
        end if



!        if (itime == 1) then
!            omega = ZERO
!        else if (itime == 2) then
!            omega = worker%time_manager%freqs(1)
!        else if (itime == 3) then
!            omega = -worker%time_manager%freqs(1)
!        end if


    end function get_omega
    !**************************************************************************************


    !>  Return the correct spatial frequency for the space index.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   5/28/2018
    !!
    !--------------------------------------------------------------------------------------
    function get_lm(itheta,ntheta) result(lm)
        integer(ik),    intent(in)  :: itheta
        integer(ik),    intent(in)  :: ntheta

        real(rk) :: lm

        if (itheta <= ((ntheta-1)/2 + 1)) then
            lm = real(itheta-1,rk)  ! positive frequencies
        else
            lm = -real(ntheta-itheta+1,rk) ! negative frequencies
        end if

    end function get_lm
    !**************************************************************************************
    







end module bc_giles_HB_base
