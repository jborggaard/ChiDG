module fcn_mmd_ldiff
#include <messenger.h>
    use mod_kinds,      only: rk,ik
    use mod_constants,  only: ZERO, HALF, ONE, TWO, THREE, FIVE, EIGHT, PI
    use type_point,     only: point_t
    use type_function,  only: function_t
    implicit none
    private








    !>  @TODO NEEDS FIXED
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !-------------------------------------------------------------------------------
    type, extends(function_t), public :: mmd_ldiff_f

    contains

        procedure   :: init
        procedure   :: compute

    end type mmd_ldiff_f
    !********************************************************************************



contains




    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/2/2016
    !!
    !-------------------------------------------------------------------------
    subroutine init(self)
        class(mmd_ldiff_f),  intent(inout)  :: self

        !
        ! Set function name
        !
        !self%name = "isentropic vortex  ::  weeeee!"
        call self%set_name("Mesh Motion : Diffusion, constant viscosity")

        call self%add_option('dspl', ZERO)
        call self%add_option('frac', ZERO)
        call self%add_option('xstart', ZERO)
        call self%add_option('xlength', ZERO)
    end subroutine init
    !*************************************************************************







    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/2/2016
    !!
    !!
    !-----------------------------------------------------------------------------------------
    impure elemental function compute(self,time,coord) result(val)
        class(mmd_ldiff_f),     intent(inout)  :: self
        real(rk),                       intent(in)  :: time
        type(point_t),                  intent(in)  :: coord

        real(rk)                                    :: val

        real(rk)    :: x,   y,   z, &
                        dspl, frac, xstart, xlength

        dspl = self%get_option_value('dspl')
        frac = self%get_option_value('frac')
        xstart = self%get_option_value('xstart')
        xlength = self%get_option_value('xlength')
        x = (coord%c1_-xstart)/xlength
        y = coord%c2_
        z = coord%c3_

        val = dspl*(ONE-frac*x)

    end function compute
    !**********************************************************************************


end module fcn_mmd_ldiff
