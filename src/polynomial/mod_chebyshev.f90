module mod_chebyshev
    use mod_kinds,      only: rk,ik,rquad
    use mod_constants,  only: XI_DIR,ETA_DIR,ZETA_DIR, &
                              ZERO, ONE, TWO, THREE, FOUR, FIVE, EIGHTH, HALF
    use mod_ordering,   only: xi_order_2d, eta_order_2d, &
                              xi_order_3d, eta_order_3d, zeta_order_3d
    use ieee_arithmetic,    only: ieee_is_nan

    implicit none

contains


    !> Compute value of hierarchical Chebyshev polynomial expansion.
    !!
    !! based on mod_legendre.f90 by Nathan A. Wukie
    !!
    !!  @author Jeff Borggaard
    !!  @date   6/14/2022
    !!
    !!
    !!
    !-----------------------------------------------------------------------------------------
    function chebyshev_val(space_dim,currentmode,xpos,ypos,zpos) result(polyval)
        integer(ik),    intent(in)  :: space_dim
        integer(ik),    intent(in)  :: currentmode
        real(rk),       intent(in)  :: xpos
        real(rk),       intent(in)  :: ypos
        real(rk),       intent(in)  :: zpos

        real(rk)                    :: polyval

        select case (space_dim)
            case (1)    ! 1D
                polyval = chebyshev_val1D(currentmode,xpos)
            case (2)    ! 2D
                polyval = chebyshev_val2D(currentmode,xpos,ypos)
            case (3)    ! 3D
                polyval = chebyshev_val3D(currentmode,xpos,ypos,zpos)
            case default
                print*, "Error - chebyshev_val: valid space dimensions are (1,2,3)"
                stop
        end select

    end function chebyshev_val
    !*****************************************************************************************




    



    recursive function chebyshev_val1D(nterm,pos) result(polyval)
    !   Compute the value of the nterm
    !   chebyshev polynomial at the location pos
    !   between -1 and 1
    !   Edit list:  based on legendre code by Nathan A. Wukie - 2/11/2015
        integer(ik),    intent(in) :: nterm
        real(rk),       intent(in) :: pos

        real(rk)                   :: polyval, polyval_nm1, polyval_nm2

        select case (nterm)
            ! Start recursion terms
            case (1)
                polyval = ONE
            case (2)
                polyval = pos
            case (3 :)
                ! Recursive definition for norder >= 2
                polyval_nm1=chebyshev_val1D(nterm-1,pos)
                polyval_nm2=chebyshev_val1D(nterm-2,pos)
                polyval = TWO*pos*polyval_nm1 - polyval_nm2
        end select


    end function chebyshev_val1D
    !****************************************************************************************


    !>  Quad-precision version.
    !!
    !!  @author Jeff Borggaard, based on legendre version by Nathan A. Wukie
    !!  @date   6/14/2022
    !!
    !-----------------------------------------------------------------------------------------
    recursive function chebyshev_val1D_quad(nterm,pos) result(polyval)
    !   Compute the value of the nterm
    !   chebyshev polynomial at the location pos
    !   between -1 and 1
        integer(ik),    intent(in) :: nterm
        real(rquad),       intent(in) :: pos

        real(rquad)                   :: polyval, polyval_nm1, polyval_nm2

        select case (nterm)
            ! Start recursion terms
            case (1)
                polyval = 1._rquad
            case (2)
                polyval = pos
            case (3 :)
                ! Recursive definition for norder >= 2
                polyval_nm1=chebyshev_val1D_quad(nterm-1,pos)
                polyval_nm2=chebyshev_val1D_quad(nterm-2,pos)
                polyval = 2._rquad*pos*polyval_nm1 - polyval_nm2
        end select


    end function chebyshev_val1D_quad
    !****************************************************************************************









    !>  A set of 2D, Chebyshev polynomials is associated with the coordinates
    !!  in 'nodes'. This function computes the value of the Chebyshev polynomial
    !!  associated with the 'currentnode' at the coordinate '(xpos,ypos)'.
    !!  
    !!  @author Jeff Borggaard, based on legendre_val2D by Nathan A. Wukie
    !!  @date   6/14/2022
    !!
    !!
    !----------------------------------------------------------------------------------------
    function chebyshev_val2D(currentmode,xi,eta) result(polyval)
        integer(ik),    intent(in)  :: currentmode
        real(rk),       intent(in)  :: xi
        real(rk),       intent(in)  ::eta

        real(rk)    :: polyval, xi_polyval, eta_polyval
        integer(ik) :: xi_mode, eta_mode

        ! compute current x/y-node indices of 1D tensor product polynomials
        ! Example: for a 2D Chebyshev polynomial L(x,y)=L(x)*L(y), compute the
        ! x and y indices
        xi_mode  = xi_order_2d( currentmode)
        eta_mode = eta_order_2d(currentmode)

        xi_polyval  = chebyshev_val1D(xi_mode,  xi)
        eta_polyval = chebyshev_val1D(eta_mode,eta)

        polyval = xi_polyval*eta_polyval

    end function chebyshev_val2D
    !*****************************************************************************************










    !>  A set of 3D, Chebyshev polynomials is associated with the coordinates
    !!  in 'nodes'. This function computes the value of the Chebyshev polynomial
    !!  associated with the 'currentnode' at the coordinate '(xpos,ypos,zpos)'.
    !!  
    !!  @author Jeff Borggaard, based on the Legendre version by Nathan A. Wukie
    !!  @date   6/14/2022
    !!
    !!
    !-----------------------------------------------------------------------------------------
    function chebyshev_val3D(currentmode,xi,eta,zeta) result(polyval)
        integer(ik),    intent(in)  :: currentmode
        real(rk),       intent(in)  :: xi
        real(rk),       intent(in)  :: eta
        real(rk),       intent(in)  :: zeta

        real(rk)    :: polyval, xi_polyval, eta_polyval, zeta_polyval
        integer(ik) :: xi_mode, eta_mode, zeta_mode

        ! compute current x/y-node indices of 1D tensor product polynomials
        ! Example: for a 3D Chebyshev polynomial L(x,y,z)=L(x)*L(y)*L(z), 
        ! compute the x, y, and z indices
        xi_mode   = xi_order_3d(  currentmode)
        eta_mode  = eta_order_3d( currentmode)
        zeta_mode = zeta_order_3d(currentmode)

        xi_polyval   = chebyshev_val1D(xi_mode,  xi)
        eta_polyval  = chebyshev_val1D(eta_mode, eta)
        zeta_polyval = chebyshev_val1D(zeta_mode,zeta)

        polyval = xi_polyval*eta_polyval*zeta_polyval

    end function chebyshev_val3D
    !*****************************************************************************************










    !> Compute directional derivative of a Chebyshev polynomial.
    !!
    !!  @author Jeff Borggaard, based on the Legendre version by Nathan A. Wukie
    !!  @date   6/14/2022
    !!
    !!
    !!
    !-----------------------------------------------------------------------------------------
    function dchebyshev_val(space_dim,currentmode,xi,eta,zeta,dir) result(dpolyval)
        integer(ik),    intent(in)  :: space_dim
        integer(ik),    intent(in)  :: currentmode
        real(rk),       intent(in)  :: xi
        real(rk),       intent(in)  :: eta
        real(rk),       intent(in)  :: zeta
        integer(ik),    intent(in)  :: dir

        real(rk)    :: dpolyval

        select case (space_dim)
            case (1)    ! 1D
                dpolyval = dchebyshev_val1D(currentmode,xi)
            case (2)    ! 2D
                dpolyval = dchebyshev_val2D(currentmode,xi,eta,dir)
            case (3)    ! 3D
                dpolyval = dchebyshev_val3D(currentmode,xi,eta,zeta,dir)
            case default
                print*, "Error - dchebyshev_val: Valid space dimensions are (1,2,3)"
                stop
        end select

    end function dchebyshev_val
    !*****************************************************************************************







    !>  Compute the first derivative of the nterm Chebyshev polynomial at the location 
    !!  'pos' between -1 and 1.
    !!
    !!  @author Jeff Borggaard, based on the Legendre version by Nathan A. Wukie
    !!  @date   6/14/2022
    !!
    !-----------------------------------------------------------------------------------------
    recursive function dchebyshev_val1D(nterm,pos) result(dpolyval)
        integer(ik), intent(in)    :: nterm
        real(rk),    intent(in)    :: pos

        real(rk)                   :: dpolyval

        select case (nterm)
            ! Trivial evaluations
            case (1)
                dpolyval = ZERO
            case (2)
                dpolyval = ONE
            case (3 :)
                ! Recursive definition
                dpolyval = 2._rquad*chebyshev_val1D(nterm-1,pos) + &
                           2._rquad*pos*dchebyshev_val1D(nterm-1,pos) - &
                           dchebyshev_val1D(nterm-2,pos)

        end select

    end function dchebyshev_val1D
    !*****************************************************************************************


    !>  Compute the first derivative of the nterm Chebyshev polynomial at the location 
    !!  'pos' between -1 and 1.
    !!
    !!  Quad-precision version.
    !!
    !!  @author Jeff Borggaard, based on legendre version by Nathan A. Wukie
    !!  @date   6/14/2022
    !!
    !-----------------------------------------------------------------------------------------
    recursive function dchebyshev_val1D_quad(nterm,pos) result(dpolyval)
        integer(ik), intent(in)    :: nterm
        real(rquad),    intent(in)    :: pos

        real(rquad)                   :: dpolyval
        real(rk)                      :: pos_lp

        select case (nterm)
            ! Trivial evaluations
            case (1)
                dpolyval = 0._rquad
            case (2)
                dpolyval = 1._rquad
            case (3 :)
                ! Recursive definition
                pos_lp = pos
                dpolyval = 2._rquad*chebyshev_val1D(nterm-1,pos_lp) + &
                           2._rquad*pos*dchebyshev_val1D(nterm-1,pos_lp) - &
                           dchebyshev_val1D(nterm-2,pos_lp)

        end select

    end function dchebyshev_val1D_quad
    !*****************************************************************************************










    !>  A set of 1D-Chebyshev polynomials is associated with the coordinates
    !!  in 'nodes'. This function compute the derivative of the Chebyshev polynomial
    !!  associated with the 'currentnode' at the location 'pos'.
    !!
    !!  @author Jeff Borggaard, based on the Legendre version by Nathan A. Wukie
    !!  @date   6/14/2022
    !!
    !-----------------------------------------------------------------------------------------
    function dchebyshev_val2D(currentmode,xi,eta,dir) result(dpolyval)
        integer(ik),    intent(in)  :: currentmode
        real(rk),       intent(in)  :: xi
        real(rk),       intent(in)  :: eta
        integer(ik),    intent(in)  :: dir

        integer(ik) :: xi_mode, eta_mode
        real(rk)    :: dpolyval
        real(rk)    :: xi_val, eta_val, dxi_val, deta_val

        xi_mode  = xi_order_2d(currentmode)
        eta_mode = eta_order_2d(currentmode)

        select case (dir)
            case (XI_DIR)
                dxi_val = dchebyshev_val1D(xi_mode,xi)
                eta_val = chebyshev_val1D(eta_mode,eta)

                dpolyval = dxi_val*eta_val

            case (ETA_DIR)
                xi_val   = chebyshev_val1D(xi_mode,xi)
                deta_val = dchebyshev_val1D(eta_mode,eta)

                dpolyval = xi_val*deta_val

            case (ZETA_DIR)
                ! By definition of 2D polynomial, no derivative in ZETA dimension
                dpolyval = ZERO 

            case default
                print*, "valid derivative directions are - 'XI_DIR', 'ETA_DIR', 'ZETA_DIR'"
                stop
        end select

    end function dchebyshev_val2D
    !*****************************************************************************************






    !>  A set of 3D-Chebyshev polynomials is associated with the coordinates
    !!  in 'nodes'. This function computes the derivative of the Chebyshev polynomial
    !!  associated with the 'currentnode' at the location 'pos' my multiplying the
    !!  1D derivatives, since the 3D modes are constructed from a tensor product of 1D modes
    !!
    !!  @author Jeff Borggaard, based on the Legendre version by Nathan A. Wukie
    !!  @date   6/14/2022
    !!
    !----------------------------------------------------------------------------------------
    function dchebyshev_val3D(currentmode,xi,eta,zeta,dir) result(dpolyval)
        integer(ik),    intent(in)  :: currentmode
        real(rk),       intent(in)  :: xi
        real(rk),       intent(in)  :: eta
        real(rk),       intent(in)  :: zeta
        integer(ik),    intent(in)  :: dir

        integer(ik) :: xi_mode, eta_mode, zeta_mode
        real(rk)    :: dpolyval
        real(rk)    :: xi_val, eta_val, zeta_val, dxi_val, deta_val, dzeta_val



        xi_mode   = xi_order_3d(  currentmode)
        eta_mode  = eta_order_3d( currentmode)
        zeta_mode = zeta_order_3d(currentmode)

        select case (dir)
            case (XI_DIR)
                dxi_val   = dchebyshev_val1D(xi_mode,xi)
                eta_val   = chebyshev_val1D(eta_mode,eta)
                zeta_val  = chebyshev_val1D(zeta_mode,zeta)

                dpolyval  = dxi_val*eta_val*zeta_val
            case (ETA_DIR)
                xi_val    = chebyshev_val1D(xi_mode,xi)
                deta_val  = dchebyshev_val1D(eta_mode,eta)
                zeta_val  = chebyshev_val1D(zeta_mode,zeta)

                dpolyval  = xi_val*deta_val*zeta_val
            case (ZETA_DIR)
                xi_val    = chebyshev_val1D(xi_mode,xi)
                eta_val   = chebyshev_val1D(eta_mode,eta)
                dzeta_val = dchebyshev_val1D(zeta_mode,zeta)

                dpolyval  = xi_val*eta_val*dzeta_val
            case default
                print*, "valid derivative directions are - 'XI_DIR', 'ETA_DIR', 'ZETA_DIR'"
                stop
        end select

    end function dchebyshev_val3D
    !*****************************************************************************************










    !>  Second/mixed derivatives of modes in chebyshev tensor product basis.
    !!
    !!  @author Jeff Borggaard, based on the Legendre version by Nathan A. Wukie
    !!  @date   6/14/2022
    !!
    !-----------------------------------------------------------------------------------------
    function ddchebyshev_val(space_dim,currentmode,xi,eta,zeta,partial1,partial2) result(res)
        integer(ik),    intent(in)  :: space_dim
        integer(ik),    intent(in)  :: currentmode
        real(rk),       intent(in)  :: xi
        real(rk),       intent(in)  :: eta
        real(rk),       intent(in)  :: zeta
        integer(ik),    intent(in)  :: partial1
        integer(ik),    intent(in)  :: partial2

        real(rk)    :: res

        select case (space_dim)
            case (1)    ! 1D
                res = ddchebyshev_val1D(currentmode,xi)
            case (2)    ! 2D
                res = ddchebyshev_val2D(currentmode,xi,eta,partial1,partial2)
            case (3)    ! 3D
                res = ddchebyshev_val3D(currentmode,xi,eta,zeta,partial1,partial2)
            case default
                print*, "Error - ddchebyshev_val: Valid space dimensions are (1,2,3)"
                stop
        end select

    end function ddchebyshev_val
    !*****************************************************************************************














    !>  Compute the second derivative of the nterm Chebyshev polynomial at the location 
    !!  'pos' between -1 and 1.
    !!
    !!  The Chebyshev polynomials satisfy the recurrence relationship:
    !!        Tk = 2xTkm1 - Tkm2
    !!
    !!  Differentiating wrt x twice leads to
    !!      ddTk = 4dTkm1 + 2xddTkm1 - ddTkm2
    !!
    !!  @author Jeff Borggaard, based on the Legendre version by Nathan A. Wukie
    !!  @date   6/14/2022
    !!
    !-----------------------------------------------------------------------------------------
    recursive function ddchebyshev_val1d(nterm,pos) result(res)
        use mod_constants,  only: ZERO, TWO, FOUR
        integer(ik), intent(in)    :: nterm
        real(rk),    intent(in)    :: pos

        real(rk)    :: res

        select case(nterm)
            case(1)
                res = ZERO
            case(2)
                res = ZERO
            case(3 :)
                res = FOUR*dchebyshev_val1d(nterm-1,pos) + &
                      TWO*pos*ddchebyshev_val1d(nterm-1,pos) - &
                      ddchebyshev_val1d(nterm-2,pos)
        end select

    end function ddchebyshev_val1d
    !*****************************************************************************************




    !>  A set of 1D-Chebyshev polynomials is associated with the coordinates
    !!  in 'nodes'. This function computes the second derivative of the 2D Chebyshev modal tensor 
    !!  product associated with the 'currentnode' at the location 'pos'.
    !!
    !!  @author Jeff Borggaard, based on the Legendre version by Nathan A. Wukie
    !!  @date   6/14/2022
    !!
    !-----------------------------------------------------------------------------------------
    function ddchebyshev_val2d(currentmode,xi,eta,partial1,partial2) result(res)
        integer(ik),    intent(in)  :: currentmode
        real(rk),       intent(in)  :: xi
        real(rk),       intent(in)  :: eta
        integer(ik),    intent(in)  :: partial1
        integer(ik),    intent(in)  :: partial2

        integer(ik) :: xi_mode, eta_mode
        real(rk)    :: term1, term2, res

        !
        ! Check valid input for derivatives
        !
        if ((partial1 /= XI_DIR) .and. (partial1 /= ETA_DIR) .and. (partial1 /= ZETA_DIR) .or. &
            (partial2 /= XI_DIR) .and. (partial2 /= ETA_DIR) .and. (partial2 /= ZETA_DIR) ) then
            print*, "valid derivative directions are - 'XI_DIR', 'ETA_DIR', 'ZETA_DIR'"
            stop
        end if


        !
        ! Get indices of terms that construct tensor product mode
        !
        xi_mode  = xi_order_2d(currentmode)
        eta_mode = eta_order_2d(currentmode)


        !
        ! Pure second derivative
        !
        if (partial1 == partial2) then

            select case (partial1)
                case (XI_DIR)
                    term1 = ddchebyshev_val1D(xi_mode,xi)
                    term2 =   chebyshev_val1D(eta_mode,eta)

                case (ETA_DIR)
                    term1 =   chebyshev_val1D(xi_mode,xi)
                    term2 = ddchebyshev_val1D(eta_mode,eta)
            end select

        !
        ! Mixed derivative:  dd(phi)/dxideta = [d(phi)/dxi][d(phi)/deta]
        !
        else
            term1 = dchebyshev_val1D(xi_mode,xi)
            term2 = dchebyshev_val1D(eta_mode,eta)

        end if


        res = term1 * term2


        ! By definition of 2D polynomial, no derivative in ZETA dimension
        if ((partial1 == ZETA_DIR) .or. (partial2 == ZETA_DIR)) res = ZERO


    end function ddchebyshev_val2d
    !*****************************************************************************************






    !>  A set of 1D-chebyshev polynomials is associated with the coordinates
    !!  in 'nodes'. This function compute the derivative of the 3D chebyshev modal tensor 
    !!  product associated with the 'currentnode' at the location 'pos'.
    !!
    !!  @author Jeff Borggaard, based on the legendre version by Nathan A. Wukie
    !!  @date   6/14/2022
    !!
    !-----------------------------------------------------------------------------------------
    function ddchebyshev_val3d(currentmode,xi,eta,zeta,partial1,partial2) result(res)
        integer(ik),    intent(in)  :: currentmode
        real(rk),       intent(in)  :: xi
        real(rk),       intent(in)  :: eta
        real(rk),       intent(in)  :: zeta
        integer(ik),    intent(in)  :: partial1
        integer(ik),    intent(in)  :: partial2

        integer(ik) :: xi_mode, eta_mode, zeta_mode
        real(rk)    :: term1, term2, term3, res

        xi_mode   = xi_order_3d(currentmode)
        eta_mode  = eta_order_3d(currentmode)
        zeta_mode = zeta_order_3d(currentmode)


        !
        ! Check valid input for derivatives
        !
        if ((partial1 /= XI_DIR) .and. (partial1 /= ETA_DIR) .and. (partial1 /= ZETA_DIR) .or. &
            (partial2 /= XI_DIR) .and. (partial2 /= ETA_DIR) .and. (partial2 /= ZETA_DIR) ) then
            print*, "valid derivative directions are - 'XI_DIR', 'ETA_DIR', 'ZETA_DIR'"
            stop
        end if




        !
        ! Pure second derivative
        !
        if (partial1 == partial2) then

            select case (partial1)
                case (XI_DIR)
                    term1 = ddchebyshev_val1D(xi_mode,xi)
                    term2 =   chebyshev_val1D(eta_mode,eta)
                    term3 =   chebyshev_val1D(zeta_mode,zeta)

                case (ETA_DIR)
                    term1 =   chebyshev_val1D(xi_mode,xi)
                    term2 = ddchebyshev_val1D(eta_mode,eta)
                    term3 =   chebyshev_val1D(zeta_mode,zeta)

                case (ZETA_DIR)
                    term1 =   chebyshev_val1D(xi_mode,xi)
                    term2 =   chebyshev_val1D(eta_mode,eta)
                    term3 = ddchebyshev_val1D(zeta_mode,zeta)

            end select


        !
        ! Mixed derivative:  
        !   ex:  dd(phi)/dxideta   = [d(phi)/dxi] * [d(phi)/deta] * [phi(zeta)]
        !   ex:  dd(phi)/detadzeta = [phi(xi)]  *  [d(phi)/deta] * [d(phi)/dzeta]
        !
        else


            select case(partial1)
                case (XI_DIR)
                    term1 = dchebyshev_val1D(xi_mode,xi)
                case (ETA_DIR)
                    term1 = dchebyshev_val1D(eta_mode,eta)
                case (ZETA_DIR)
                    term1 = dchebyshev_val1D(zeta_mode,zeta)
            end select


            select case(partial2)
                case (XI_DIR)
                    term2 = dchebyshev_val1D(xi_mode,xi)
                case (ETA_DIR)
                    term2 = dchebyshev_val1D(eta_mode,eta)
                case (ZETA_DIR)
                    term2 = dchebyshev_val1D(zeta_mode,zeta)
            end select

            
            ! Determine third term by knowledge of first two terms in derivative
            select case(partial1 + partial2)
                case (XI_DIR + ETA_DIR)
                    term3 = chebyshev_val1D(zeta_mode,zeta)
                case (XI_DIR + ZETA_DIR)
                    term3 = chebyshev_val1D(eta_mode,eta)
                case (ETA_DIR + ZETA_DIR)
                    term3 = chebyshev_val1D(xi_mode,xi)
            end select



        end if


        ! Compute output
        res = term1 * term2 * term3

        if (ieee_is_nan(res)) print*, 'term ', currentmode, ' is nan'


    end function ddchebyshev_val3d
    !*****************************************************************************************










end module mod_chebyshev
