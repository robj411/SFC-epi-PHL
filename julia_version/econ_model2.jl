# ============================================================================
#  econ_model2.jl  —  OWNED BY THE SFC MODELLER
# ----------------------------------------------------------------------------
#  Model 2 = Model 1 + time-varying labour productivity λ.
#
#  Read this as a DIFF against econ_model.jl. Lines marked ## NEW or ## CHANGED
#  are the only differences. Everything else is identical — and note what is
#  ABSENT: there is no dot_cons_s / dot_cons_d / quotient-rule block. In the
#  odin Model 2 that block was ~15 lines of hand-derived calculus with the
#  fragile 1/(Y - λ·p1·available_lf) inversion. Here productivity tracks output
#  growth by writing D(Y) directly and letting MTK differentiate the algebraic
#  definition of Y during structural_simplify. That is the whole point.
# ============================================================================

using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

"""
    build_econ2(; name, alpha0, alpha1, alpha2, theta, G, cons0, emp0, lf,
                  lambda0, lambda_p0, lambda_p1, q1, q2, standalone=false)

SIM econ block with productivity `λ` promoted from a parameter to a state that
tracks output growth: dλ/dt = λ·(p0 + p1/Y · dY/dt).
"""
function build_econ2(; name,
                     alpha0, alpha1, alpha2, theta, G, cons0, emp0, lf,
                     lambda0, lambda_p0, lambda_p1, q1, q2,          ## CHANGED args
                     H0 = 0.0, standalone::Bool = false)

    @variables notifications(t) C_lf(t)

    @variables H_h(t)=H0
    @variables λv(t) =lambda0                                             ## NEW: productivity is now a state

    @variables scalar(t) available_lf(t) cons_d(t) cons_s(t) cons(t) Y(t) YD(t)
    @variables relative_consumption(t) relative_work(t)

    @parameters α0=alpha0 α1=alpha1 α2=alpha2 θ=theta Gp=G c0=cons0 e0=emp0 lfp=lf q1p=q1 q2p=q2
    @parameters p0=lambda_p0 p1=lambda_p1                           ## NEW productivity parameters
    # (λ is no longer a parameter — it is a state variable λv)

    α1t   = scalar * α1
    α2t   = scalar * α2
    denom = 1 - α1t * (1 - θ)

    eqs = [
        scalar       ~ q1p + (1 - q1p) / (1 + q2p * notifications),
        available_lf ~ lfp - C_lf / 1e6,

        cons_d ~ (α0 + α1t * Gp * (1 - θ) + α2t * H_h) / denom,
        cons_s ~ max(0, available_lf * λv - Gp),                    ## CHANGED: λv, the state
        # ifelse (not min) so MTK can differentiate through the switch
        # branch-wise when it forms D(Y) below:
        cons   ~ ifelse(cons_s < cons_d, cons_s, cons_d),          ## CHANGED

        Y      ~ cons + Gp,
        YD     ~ Y * (1 - θ),
        D(H_h) ~ YD - cons,

        # --- productivity dynamics, written the way the model reads --------
        # Y is algebraic; D(Y) is expanded symbolically by structural_simplify
        # into the implicit relation the odin code solved by hand. You never
        # write that inversion — the DAE solver resolves the implicit dλ.
        D(λv)  ~ λv * (p0 + p1 / Y * D(Y)),                         ## NEW

        relative_consumption ~ cons / c0,
        relative_work        ~ (cons + Gp) / λv / e0,              ## CHANGED: λv
    ]

    if standalone
        push!(eqs, notifications ~ 0)
        push!(eqs, C_lf ~ 0)
    end

    return ODESystem(eqs, t; name = name)

    # ---------------------------------------------------------------------
    # HONEST CAVEAT: MTK removes the hand-calculus, not the singularity. The
    # supply-determined regime is genuinely near-singular when p1·available_lf·λ
    # approaches Y (that is what the odin denominator was). The implicit solver
    # approaches it far more gracefully than dividing by it, and crucially YOU
    # are not the one writing and debugging the inverted expression — but if p1
    # is large you should still watch that regime.
    #
    # ROBUST FALLBACK if structural_simplify won't auto-expand D(Y): introduce
    # the derivative explicitly as its own variable,
    #     @variables dY(t)
    #     push!(eqs, dY ~ D(cons))            # MTK differentiates cons for you
    #     ... use dY in the D(λv) equation
    # This is still auto-derived — you are only naming the intermediate.
    # ---------------------------------------------------------------------
end
