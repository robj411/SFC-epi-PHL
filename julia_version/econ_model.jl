# ============================================================================
#  econ_model.jl  —  OWNED BY THE SFC MODELLER
# ----------------------------------------------------------------------------
#  Model 1: the SIM model of Godley & Lavoie (ch. 3), continuous time.
#  It knows nothing about the epidemic. It CONSUMES two inputs from whatever
#  epi model is attached:
#        notifications, C_lf
#  and PROVIDES two outputs back:
#        relative_consumption, relative_work
#
#  This is where the model grows. New behavioural equations, sectors, assets,
#  inventories, etc. are added HERE as extra state variables (D(x) ~ ...) and
#  algebraic identities (y ~ ...). The stocks are differential; the accounting
#  identities are algebraic — you never substitute them out or differentiate
#  them by hand. structural_simplify/mtkcompile does the index reduction and
#  hands the DAE to the solver.
#
#  Build with standalone=true to run it on its own: there is then no pandemic
#  signal (total_C = notifications = C_lf = 0), so it should sit at its
#  normal-times equilibrium — a useful "is the block well-posed?" check.
# ============================================================================

using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

"""
    build_econ(; name, alpha0, alpha1, alpha2, theta, G, cons0, emp0, lf,
                 lambda, q1, q2, standalone=false)

Return an ODESystem for the SIM econ block with the epi->econ behavioural
link. `lambda` is (constant) labour productivity in Model 1; promote it to a
state to obtain Model 2 (see note at the bottom).
"""
function build_econ(; name,
                     alpha0, alpha1, alpha2, theta, G, cons0, emp0, lf, lambda,
                     q1, q2, H0 = 0.0, standalone::Bool = false)

    # --- coupling inputs from the epi side --------------------------------
    @variables notifications(t) C_lf(t)

    # --- state: household wealth ------------------------------------------
    @variables H_h(t)=H0

    # --- algebraic intermediates (accounting identities + behaviour) ------
    @variables scalar(t) available_lf(t) cons_d(t) cons_s(t) cons(t) Y(t) YD(t)
    @variables relative_consumption(t) relative_work(t)

    # --- tunable econ parameters (defaults from the R setup) --------------
    @parameters α0=alpha0 α1=alpha1 α2=alpha2 θ=theta Gp=G c0=cons0 e0=emp0 lfp=lf λ=lambda q1p=q1 q2p=q2

    # behavioural response of propensity-to-consume to notifications
    α1t   = scalar * α1
    α2t   = scalar * α2
    denom = 1 - α1t * (1 - θ)

    eqs = [
        # epi -> econ behavioural link
        scalar       ~ q1p + (1 - q1p) / (1 + q2p * notifications),
        available_lf ~ lfp - C_lf / 1e6,

        # consumption: demand-determined vs supply-determined, take the min
        cons_d ~ (α0 + α1t * Gp * (1 - θ) + α2t * H_h) / denom,
        cons_s ~ max(0, available_lf * λ - Gp),
        cons   ~ min(cons_s, cons_d),          # <-- non-smooth; see note
        # to achieve soft min, choose large k:
        # cons_s ~ log1p(exp(k*(available_lf*λv - Gp))) / k
        # cons ~ -log(exp(-k*cons_s) + exp(-k*cons_d)) / k

        # accounting identities (SFC): output, disposable income, wealth flow
        Y      ~ cons + Gp,
        YD     ~ Y * (1 - θ),
        D(H_h) ~ YD - cons,                    # dH_h/dt = saving

        # econ -> epi outputs
        relative_consumption ~ cons / c0,
        relative_work        ~ (cons + Gp) / λ / e0,
    ]

    if standalone
        push!(eqs, notifications ~ 0)
        push!(eqs, C_lf ~ 0)
    end

    return ODESystem(eqs, t; name = name)

    # ---------------------------------------------------------------------
    # NOTE 1 (Model 2, no hand-calculus):
    #   Add   @variables λv(t)     (productivity as a STATE), replace the λ
    #   parameter with λv, and append
    #       D(λv) ~ λv * (lambda_p0 + lambda_p1 / Y * D(Y))
    #   You may write D(Y) directly: MTK differentiates the algebraic
    #   definition of Y symbolically during structural_simplify. There is NO
    #   need for the hand-derived dot_cons_s / quotient-rule branch that the
    #   odin version required — that whole block disappears.
    #
    # NOTE 2 (the min() kink):
    #   min(cons_s, cons_d) is non-smooth at the supply/demand crossing. For
    #   production, replace it with a ContinuousCallback root-find on
    #   (cons_s - cons_d) and let each regime be its own smooth equation,
    #   rather than integrating across the kink.
    # ---------------------------------------------------------------------
end
