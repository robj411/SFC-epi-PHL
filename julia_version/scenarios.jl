# ============================================================================
#  scenarios.jl  —  reproduce the R script's 11-scenario metric table
# ----------------------------------------------------------------------------
#  Prints, per scenario, the same five numbers R's get_results_df() prints:
#     peak_infected(M) | peak_day | cum_incidence(%) | min_employment(%) | gdp_loss(%)
#  in the same "&"-separated layout, so you can diff against the R stdout.
#
#  Uses Model 2 (econ = model2 in R_script.R). Assumes epi_model.jl and
#  econ_model2.jl are already included and the data / econ-setup block from
#  run.jl has run (NNs, Mww, Mcw, Mcc, Mcom, beta, T*, pdet, alpha*, theta,
#  g0, cons0, emp0, lf, lambda0, h0, emprate all in scope).
# ============================================================================
include("get_started.jl")

# --- consumption_contact (η) only splits the community matrix; Mww, Mcw are
#     untouched, and the RDS was built at η = 0.4. Recover the full matrix once.
Mcom_full = Mcom .+ Mcc                     # 0.6·full + 0.4·full = full
contacts_at(cc) = (Mww, Mcw, cc .* Mcom_full, (1 - cc) .* Mcom_full)

days  = 0.0:1.0:365.0
trapz(y) = sum((y[1:end-1] .+ y[2:end]) ./ 2)          # unit day spacing

# --- scenario table (transcribed from R_script.R) ---------------------------
base = (q1 = 0.86273916, q2 = 0.02158386, γ1 = 0.5, cc = 0.4)
scenarios = [
    base,
    merge(base, (q1 = 0.73,)),
    merge(base, (q1 = 0.93,)),
    merge(base, (q2 = 0.00002,)),
    merge(base, (q2 = 0.2,)),
    merge(base, (γ1 = 0.25,)),
    merge(base, (γ1 = 0.75,)),
    merge(base, (cc = 0.2,)),
    merge(base, (cc = 0.6,)),
    (q1 = 0.93, q2 = 0.00002,   γ1 = 0.75, cc = 0.2),
    (q1 = 0.93, q2 = 0.0000002, γ1 = 0.75, cc = 0.2),
]

pop = sum(NNs)

# --- helpers to pull the epi aggregates from a solution ---------------------
strata = 1:length(NNs)
allinf(sol, epi) = [sum(sol(d, idxs = epi.C[i]) + sol(d, idxs = epi.Id[i]) +
                        sol(d, idxs = epi.Iu[i]) for i in strata) for d in days]
totR(sol, epi)   = [sum(sol(d, idxs = epi.R[i]) for i in strata) for d in days]

# ============================================================================
#  INDEPENDENT run (counterfactual): cc-invariant, computed once.
#  Y is flat at y0, so the GDP-loss denominator is y0 over the horizon.
# ============================================================================
Mww, Mcw_, Mcc0, Mcom0 = contacts_at(0.4)

@named epi_ind = build_epi(; NNs, Mww = Mww, Mcw = Mcw_, Mcc = Mcc0, Mcom = Mcom0,
                             beta, TEtoI, TItoR, TItoC, TCtoR, prob_detected = pdet,
                             S0 = NNs .- imported, E0 = imported, standalone = true)
epi_solo_s = structural_simplify(epi_ind)
prob_epi = ODEProblem(epi_solo_s, [], tspan)          # empty u0 map
sol_ind = solve(prob_epi, FBDF())
##
ai = allinf(sol_ind, epi_ind)
y0 = cons0 + g0                                   # flat counterfactual GDP

peak_inf_ind = maximum(ai) / 1e6
peak_day_ind = days[argmax(ai)]
cuminc_ind   = maximum(totR(sol_ind, epi_ind)) / pop * 100

println("model & peak_inf(M) & peak_day & cum_inc(%) & min_emp(%) & gdp_loss(%)")
println("Independent & $(round(peak_inf_ind, digits=1)) & $(Int(peak_day_ind)) & ",
        "$(round(Int, cuminc_ind)) & 100.0 & 0")

# ============================================================================
#  INTEGRATED runs, one per scenario
# ============================================================================
for (k, s) in enumerate(scenarios)
    local Mw, Mcw_, Mcc_s, Mcom_s = contacts_at(s.cc)

    local @named epic = build_epi(; NNs, Mww = Mw, Mcw = Mcw_, Mcc = Mcc_s, Mcom = Mcom_s,
                              beta, TEtoI, TItoR, TItoC, TCtoR, prob_detected = pdet,
                              S0 = NNs .- imported, E0 = imported, standalone = false)
    local @named econc = build_econ2(; alpha0, alpha1, alpha2, theta, G = g0, cons0, emp0,
                                 lf, lambda0, lambda_p0 = 0.0, lambda_p1 = s.γ1,
                                 q1 = s.q1, q2 = s.q2, H0 = h0, standalone = false)
    conns = [
        epic.rel_cons       ~ econc.relative_consumption,
        epic.rel_work       ~ econc.relative_work,
        econc.notifications ~ epic.notifications,
        econc.C_lf          ~ epic.C_lf,
    ]
    sys = structural_simplify(compose(ODESystem(conns, t; name = Symbol("cpl$k")), epic, econc))
    sol = solve(ODEProblem(sys, [], tspan;
                guesses = [econc.Y => y0, econc.cons => cons0, econc.scalar => 1.0]),
                FBDF())

    local ai   = allinf(sol, epic)
    Yint = [sol(d, idxs = econc.Y) for d in days]
    minemp = minimum(sol(d, idxs = econc.relative_work) for d in days) * emprate
    gdploss = (1 - trapz(Yint) / (y0 * (length(days) - 1))) * 100

    println("Integrated & $(round(maximum(ai)/1e6, digits=1)) & $(Int(days[argmax(ai)])) & ",
            "$(round(Int, maximum(totR(sol, epic))/pop*100)) & ",
            "$(round(minemp, digits=1)) & $(round(gdploss, digits=1))   % scenario $k ",
            "(q1=$(s.q1), q2=$(s.q2), γ1=$(s.γ1), η=$(s.cc))")
end
