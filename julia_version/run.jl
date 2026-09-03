# ============================================================================
#  run.jl  —  the only file that sees both sides
# ----------------------------------------------------------------------------
#  Reads data_file.Rds from the existing R pipeline, reproduces the Model-1
#  econ setup, then runs three things:
#     (1) epi alone      — the epidemic counterfactual (activity pinned to 1)
#     (2) econ alone      — the econ block with no pandemic (equilibrium check)
#     (3) coupled         — epi <-> econ, connected only through 5 equations
#
#  The point of this file: the ONLY place the two models meet is the
#  `connections` vector below. Neither modeller edits the other's file.
# ============================================================================

include("get_started.jl")

# ============================================================================
#  RUN (1): EPI ALONE  — economic activity pinned to normal times
# ============================================================================
# epi-alone
@named epi_solo = build_epi(; NNs, Mww, Mcw, Mcc, Mcom, beta, TEtoI, TItoR,
                              TItoC, TCtoR, prob_detected = pdet,
                              S0 = NNs .- imported, E0 = imported, standalone = true)
epi_solo_s = structural_simplify(epi_solo)
prob_epi = ODEProblem(epi_solo_s, [], tspan)          # empty u0 map

sol_epi = solve(prob_epi, FBDF())
@info "epi-alone done" peak_confirmed = maximum(sol_epi[epi_solo.total_C]) percent_infected = sum(maximum(sol_epi[epi_solo.R]))/sum(NNs)*100


# ============================================================================
#  RUN (2): ECON ALONE  — no pandemic signal, should hold at equilibrium
# ============================================================================
# econ-alone
@named econ_solo = build_econ(; alpha0, alpha1, alpha2, theta, G = g0, cons0,
                                emp0, lf, lambda, q1, q2, H0 = h0, standalone = true)
econ_solo_s = structural_simplify(econ_solo)
prob_econ = ODEProblem(econ_solo_s, [], tspan)

# prob_econ = ODEProblem(econ_solo_s, [econ_solo.H_h => h0], tspan, [])
sol_econ = solve(prob_econ, FBDF())
@info "econ-alone done" H_h_start = sol_econ[econ_solo.H_h][1] H_h_end = sol_econ[econ_solo.H_h][end]

# ============================================================================
#  RUN (3): COUPLED  — the two components connected by 4 equations
# ============================================================================

# --- coupled ---
@named epic  = build_epi(; NNs, Mww, Mcw, Mcc, Mcom, beta, TEtoI, TItoR,
                           TItoC, TCtoR, prob_detected = pdet,
                           S0 = NNs .- imported, E0 = imported,
                           standalone = false)

@named econc = build_econ(; alpha0, alpha1, alpha2, theta, G = g0, cons0,
                            emp0, lf, lambda, q1, q2, H0 = h0,
                            standalone = false)

connections = [
    epic.rel_cons        ~ econc.relative_consumption,
    epic.rel_work        ~ econc.relative_work,
    econc.notifications  ~ epic.notifications,
    econc.C_lf           ~ epic.C_lf,
]

@named coupled = compose(ODESystem(connections, t; name = :coupled), epic, econc)
coupled_s = structural_simplify(coupled)

prob_cpl = ODEProblem(coupled_s, [], tspan)
sol_cpl  = solve(prob_cpl, FBDF())
@info "coupled done" peak_confirmed = maximum(sol_cpl[epic.total_C]) H_h_end = sol_cpl[econc.H_h][end] percent_infected = sum(maximum(sol_cpl[epic.R]))/sum(NNs)*100


# ----------------------------------------------------------------------------
#  Quick comparison: pandemic impact on GDP path (coupled vs no epidemic)
# ----------------------------------------------------------------------------
# sol_cpl[econc.Y]  is the integrated GDP path;
# sol_econ[econ_solo.Y] is the counterfactual (flat) path.
using Plots
plot(sol_cpl.t,  sol_cpl[econc.Y],      label="GDP, integrated")
plot!(sol_econ.t, sol_econ[econ_solo.Y], label="GDP, counterfactual")
plot(sol_cpl.t,  sol_cpl[epic.total_C], label="Confirmed, integrated")

println("All three runs constructed. See the connections vector for the full")
println("epi<->econ interface — that is the only shared surface.")





## repeat for model 2

# ============================================================================
#  RUN (2): ECON ALONE  — no pandemic signal, should hold at equilibrium
# ============================================================================
# econ-alone
@named econ_solo = build_econ2(; alpha0, alpha1, alpha2, theta, G = g0, cons0, emp0, lf,
                  lambda0, lambda_p0, lambda_p1, q1, q2, H0 = h0, standalone = true)               
econ_solo_s = structural_simplify(econ_solo)
prob_econ = ODEProblem(econ_solo_s, [], tspan)

# prob_econ = ODEProblem(econ_solo_s, [econ_solo.H_h => h0], tspan, [])
sol_econ = solve(prob_econ, FBDF())
@info "econ-alone done" H_h_start = sol_econ[econ_solo.H_h][1] H_h_end = sol_econ[econ_solo.H_h][end]

# ============================================================================
#  RUN (3): COUPLED  — the two components connected by 4 equations
# ============================================================================

# --- coupled ---
@named epic  = build_epi(; NNs, Mww, Mcw, Mcc, Mcom, beta, TEtoI, TItoR,
                           TItoC, TCtoR, prob_detected = pdet,
                           S0 = NNs .- imported, E0 = imported,
                           standalone = false)

@named econc = build_econ2(; alpha0, alpha1, alpha2, theta, G = g0, cons0,
                            emp0, lf, lambda0, lambda_p0, lambda_p1, q1, q2, H0 = h0,
                            standalone = false)

connections = [
    epic.rel_cons        ~ econc.relative_consumption,
    epic.rel_work        ~ econc.relative_work,
    econc.notifications  ~ epic.notifications,
    econc.C_lf           ~ epic.C_lf,
]

@named coupled = compose(ODESystem(connections, t; name = :coupled), epic, econc)
coupled_s = structural_simplify(coupled)

prob_cpl = ODEProblem(coupled_s, [], tspan)
sol_cpl  = solve(prob_cpl, FBDF())
@info "coupled done" peak_confirmed = maximum(sol_cpl[epic.total_C]) H_h_end = sol_cpl[econc.H_h][end] percent_infected = sum(maximum(sol_cpl[epic.R]))/sum(NNs)*100


# ----------------------------------------------------------------------------
#  Quick comparison: pandemic impact on GDP path (coupled vs no epidemic)
# ----------------------------------------------------------------------------
##
# sol_cpl[econc.Y]  is the integrated GDP path;
# sol_econ[econ_solo.Y] is the counterfactual (flat) path.
using Plots
plot(sol_cpl.t,  sol_cpl[econc.Y],      label="GDP, integrated")
p = plot!(sol_econ.t, sol_econ[econ_solo.Y], label="GDP, counterfactual")
display(p)
##
strata = ["workers", "children", "non-workers", "retired"]   # NNs order
days = 0:365
I = reduce(hcat, [sol_cpl(days, idxs = epic.Id[i]).u .+ sol_cpl(days, idxs = epic.Iu[i]).u
                  for i in 1:4])
plot(days, I; label = permutedims(strata), xlabel = "day",
     ylabel = "infectious (Id + Iu)", lw = 2)
p = plot!(sol_epi.t, [sol_epi[epi_solo.Id[i]] .+ sol_epi[epi_solo.Iu[i]] for i in 1:4] |> x->reduce(hcat,x); ls=:dash, label=false)
display(p)

