@enum prima_message::UInt32 begin
    PRIMA_MSG_NONE = 0
    PRIMA_MSG_EXIT = 1
    PRIMA_MSG_RHO = 2
    PRIMA_MSG_FEVL = 3
end

@enum prima_rc::Int32 begin
    PRIMA_SMALL_TR_RADIUS = 0
    PRIMA_FTARGET_ACHIEVED = 1
    PRIMA_TRSUBP_FAILED = 2
    PRIMA_MAXFUN_REACHED = 3
    PRIMA_MAXTR_REACHED = 20
    PRIMA_NAN_INF_X = -1
    PRIMA_NAN_INF_F = -2
    PRIMA_NAN_INF_MODEL = -3
    PRIMA_NO_SPACE_BETWEEN_BOUNDS = 6
    PRIMA_DAMAGING_ROUNDING = 7
    PRIMA_ZERO_LINEAR_CONSTRAINT = 8
    PRIMA_INVALID_INPUT = 100
    PRIMA_ASSERTION_FAILS = 101
    PRIMA_VALIDATION_FAILS = 102
    PRIMA_MEMORY_ALLOCATION_FAILS = 103
end

function prima_get_rc_string(rc)
    @ccall libprimac.prima_get_rc_string(rc::Cint)::Cstring
end

# typedef void ( * prima_obj ) ( const double x [ ] , double * f )
const prima_obj = Ptr{Cvoid}

# typedef void ( * prima_objcon ) ( const double x [ ] , double * f , double constr [ ] )
const prima_objcon = Ptr{Cvoid}
