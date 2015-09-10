module mod_time_scheme
#include <messenger.h>
    use mod_kinds,          only: rk,ik
    use atype_time_scheme,  only: time_scheme_t
    use type_dict,          only: dict_t




    ! Import solverdata types
    use forward_euler,                  only: forward_euler_t
    use backward_euler,                 only: backward_euler_t
    use backward_euler_subiteration,    only: backward_euler_subiteration_t
    implicit none



    ! Instantiate solver types for sourcing
    type(forward_euler_t)               :: FORWARD_EULER_INSTANCE
    type(backward_euler_t)              :: BACKWARD_EULER_INSTANCE
    type(backward_euler_subiteration_t) :: BACKWARD_EULER_SUBITERATION_INSTANCE



    logical :: initialized = .false.



contains





    subroutine create_time_scheme(time_string,instance,options)
        character(*),                       intent(in)      :: time_string
        class(time_scheme_t), allocatable,  intent(inout)   :: instance
        type(dict_t), optional,             intent(inout)   :: options

        select case (trim(time_string))
            case ('ForwardEuler','forwardeuler','Forward_Euler','forward_euler','FE','fe')
                allocate(instance, source=FORWARD_EULER_INSTANCE)


            case ('BackwardEuler','backwardeuler','Backward_Euler','backward_euler','BE','be')
                allocate(instance, source=BACKWARD_EULER_INSTANCE)

            case ('BackwardEulerSub','backwardeulersub','Backward_Euler_Sub','backward_euler_sub','BES','bes')
                allocate(instance, source=BACKWARD_EULER_SUBITERATION_INSTANCE)


            case default
                call signal(FATAL,'create_time_scheme -- solver string not recognized')
        end select





        !
        ! Call options initialization if present
        !
        if (present(options)) then
            call instance%set(options)
        end if



        !
        ! Make sure the solver was allocated
        !
        if (.not. allocated(instance)) call signal(FATAL,"create_time_scheme: solver was not allocated. Check that the desired solver was registered and instantiated in the mod_time_scheme module")


    end subroutine








end module mod_time_scheme