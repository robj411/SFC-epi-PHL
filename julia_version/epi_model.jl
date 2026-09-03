# ============================================================================
#  epi_model.jl  —  OWNED BY THE EPIDEMIOLOGIST
# ----------------------------------------------------------------------------
#  A self-contained epidemic component. It knows nothing about the SFC model.
#  It CONSUMES two scalar inputs from whatever econ model is attached:
#        rel_cons  — consumption relative to normal times
#        rel_work  — employment relative to normal times
#  It PROVIDES three scalar outputs back to the econ model (the coupling
#  contract; see coupling in run.jl):
#        total_C        — total confirmed cases
#        notifications  — daily notification rate
#        C_lf           — confirmed cases among workers (labour-force loss)
#
#  Build with standalone=true to run it on its own: economic activity is
#  then pinned to normal-times levels (rel_cons = rel_work = 1), which is
#  exactly the epi-only counterfactual.
# ============================================================================

using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
#  NOTE (v10): `t_nounits`/`D_nounits` still exist; `ODESystem` is aliased to
#  `System` and `structural_simplify` to `mtkcompile`. Pin your versions.

"""
    build_epi(; name, NNs, Mww, Mcw, Mcc, Mcom, beta, TEtoI, TItoR,
                TItoC, TCtoR, prob_detected, standalone=false)

Return an ODESystem for the stratified S-E-Iu-Id-C-R model.
`NNs` is the stratum-population vector (length n); the four contact-component
matrices are n×n. Epidemic rate constants are exposed as tunable parameters.
"""

function build_epi(; name, NNs, Mww, Mcw, Mcc, Mcom,
                     beta, TEtoI, TItoR, TItoC, TCtoR, prob_detected,
                     S0, E0,
                     Iu0 = zeros(length(NNs)), Id0 = zeros(length(NNs)),
                     C0  = zeros(length(NNs)), R0  = zeros(length(NNs)),
                     standalone::Bool = false)
    n = length(NNs)
    @variables rel_cons(t) rel_work(t)
    @variables (S(t))[1:n]=S0 (E(t))[1:n]=E0 (Iu(t))[1:n]=Iu0 (Id(t))[1:n]=Id0 (C(t))[1:n]=C0 (R(t))[1:n]=R0

    # --- tunable epidemic parameters (defaults come from the RDS) ---------
    @parameters β=beta Tei=TEtoI Tir=TItoR Tic=TItoC Tcr=TCtoR pdet=prob_detected

    # Contact matrix reassembled every instant from its components, scaled by
    # current economic activity. Mww/Mcw/Mcc/Mcom are numeric n×n; rel_* are
    # symbolic scalars, so `contact` is an n×n array of expressions.
    contact = @. Mcom + rel_work^2 * Mww +
                 rel_work * rel_cons * Mcw + rel_cons^2 * Mcc

    Itot = [Id[i] + Iu[i] for i in 1:n]
    foi  = [β * sum(contact[i, j] * Itot[j] / NNs[j] for j in 1:n) for i in 1:n]

    latent_det   = [E[i]  * pdet       / Tei for i in 1:n]
    latent_undet = [E[i]  * (1 - pdet) / Tei for i in 1:n]
    undet_rec    = [Iu[i] / Tir              for i in 1:n]
    det_detected = [Id[i] / Tic              for i in 1:n]
    con_rec      = [C[i]  / Tcr              for i in 1:n]

    eqs = Equation[]
    for i in 1:n
        push!(eqs, D(S[i])  ~ -S[i] * foi[i])
        push!(eqs, D(E[i])  ~  S[i] * foi[i] - latent_det[i] - latent_undet[i])
        push!(eqs, D(Iu[i]) ~  latent_undet[i] - undet_rec[i])
        push!(eqs, D(Id[i]) ~  latent_det[i]   - det_detected[i])
        push!(eqs, D(C[i])  ~  det_detected[i] - con_rec[i])
        push!(eqs, D(R[i])  ~  undet_rec[i] + con_rec[i])
    end

    # --- outputs exposed to the econ side (the coupling contract) ---------
    @variables total_C(t) notifications(t) C_lf(t)
    push!(eqs, total_C       ~ sum(C[i] for i in 1:n))
    push!(eqs, notifications ~ sum(det_detected[i] for i in 1:n))
    # push!(eqs, notifications ~ total_C/Tcr)
    push!(eqs, C_lf          ~ C[1])            # workers are stratum 1

    if standalone
        # No econ attached → economic activity stays at normal-times level.
        push!(eqs, rel_cons ~ 1)
        push!(eqs, rel_work ~ 1)
    end

    return ODESystem(eqs, t; name = name)
end
