using ModelingToolkit
using ModelingToolkit: t_nounits as t
using OrdinaryDiffEq        # solvers (FBDF / Rodas5 handle the DAE)

include("epi_model.jl")
include("econ_model.jl")
include("econ_model2.jl")

# ----------------------------------------------------------------------------
#  1. Load the country / epidemic data produced by the R pipeline
# ----------------------------------------------------------------------------
#  RData.load returns the R object as a nested Dict keyed by the list names.
#  The conmat S3 classes in full_cms/reduced_cms are read as plain matrices
#  with attributes attached; we only touch the plain-numeric parts below.
using CodecZlib, RData
raw = tempname()
open(raw, "w") do out
    s = GzipDecompressorStream(open("../data_file.Rds"))
    write(out, read(s)); close(s)
end
ld = RData.load(raw)

contacts = ld["contacts"]
Mww  = Matrix{Float64}(contacts["Mww"])
Mcw  = Matrix{Float64}(contacts["Mcw"])
Mcc  = Matrix{Float64}(contacts["Mcc"])
Mcom = Matrix{Float64}(contacts["Mcom"])
NNs  = Vector{Float64}(ld["NNs"])          # [workers, children, non-workers, retired]

ep    = ld["epidemic"]
beta  = ep["beta"];  #
TEtoI = ep["TEtoI"]; 
TItoR = ep["TItoR"]
TItoC = ep["TItoC"]; 
TCtoR = ep["TCtoR"]; 
pdet  = ep["prob_detected"]

# ----------------------------------------------------------------------------
#  2. Reproduce the Model-1 (SIM) econ setup  (ported from econ_models.R)
# ----------------------------------------------------------------------------
gdp = ld["gdp"]; 
tax = ld["tax"]; 
emprate = ld["employmentrate"]

y0    = gdp / 365                 # annual -> daily
tax0  = tax / 365
theta = tax0 / y0
yd0   = (1 - theta) * y0
cons0 = yd0
g0    = y0 - cons0
h0    = 19032.66 / 1e3            # sourced household wealth

alpha1 = 0.7695
alpha2 = 0.0105 * 4 / 365
alpha0 = cons0 - h0 * alpha2 - yd0 * alpha1

lf     = NNs[1] / 1e6            # labour force, millions
emp0   = lf * emprate / 100
lambda = y0 / emp0

# for econ model 2
lambda0 = lambda
lambda_p0 = 0
lambda_p1 = 0.5

q1 = 0.862752
q2 = 0.02158386

# ----------------------------------------------------------------------------
#  3. Initial epidemic conditions: 200 seed infections spread over strata
# ----------------------------------------------------------------------------
imported = 1 .* NNs ./ sum(NNs)
S0 = NNs .- imported
E0 = imported
Z4 = zeros(length(NNs))
tspan = (0.0, 365.0)
