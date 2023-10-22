using Revise
using ModelingToolkit
using DifferentialEquations
using Interpolations
using OrdinaryDiffEq: ReturnCode.Success
import ModelingToolkitFluidLibrary: CompressorMap, Compressor
using Plots


eta_table = [
    0.825 0.820 0.813 0.8115 0.81 0.805
    0.840 0.835 0.831 0.8301 0.827 0.82
    0.835 0.832 0.828 0.823 0.821 0.815
    0.825 0.820 0.8146 0.808 0.80 0.79
    0.795 0.790 0.783 0.777 0.772 0.765
]

phic_table = [
    68.534 73.0 76.64 79.14 82.14 85.864
    75 78.514 81.924 84.924 87.5 91.942
    78 81.788 85.5 88.2 90.5 94.5
    79.768 84.2 88.178 91 93.178 96
    80.454 85.454 89.668 92.668 94.668 97.314
]

PR_table = [
    12.99 14.39 15.85 16.85 18 20.1
    12.5 14 15.59 16.69 17.7 19.62
    11.9 13.52 15 16.2 17.175 18.9
    10.65 12.5 14.095 15.095 16.095 17.7
    9.35 10.75 12.275 13.275 14.275 16.0
]


beta_table = [1; 2; 3; 4; 5]
N_T_table = [90; 92; 95; 97; 100; 105]

interp_eta = interpolate((beta_table, N_T_table), eta_table, Gridded(Linear()))
interp_phic = interpolate((beta_table, N_T_table), phic_table, Gridded(Linear()))
interp_pr = interpolate((beta_table, N_T_table), PR_table, Gridded(Linear()))

Eta(beta, N_T) = interp_eta(beta, N_T)
Phic(beta, N_T) = interp_phic(beta, N_T)
PR(beta, N_T) = interp_pr(beta, N_T)

@register_symbolic Eta(beta, N_T)
@register_symbolic Phic(beta, N_T)
@register_symbolic PR(beta, N_T)

# график зависимости η от приведенного расхода при разных значения приведенной скорости и β
plot([[Phic(beta, N_T) for beta in beta_table] for N_T in N_T_table],
     [[Eta(beta, N_T) for beta in beta_table] for N_T in N_T_table],
     labels=reshape(["nᵣ=$N_T%" for N_T in N_T_table], 1, length(N_T_table)),
     xlabel = "phic",
     ylabel = "η")
plot!([[Phic(beta, N_T) for N_T in N_T_table] for beta in beta_table],
     [[Eta(beta, N_T) for N_T in N_T_table] for beta in beta_table],
     linestyle=:dash,
     labels=reshape(["β=$beta%" for beta in beta_table], 1, length(beta_table)),
     legend=:outerbottomright)

# график зависимости π от приведенного расхода при разных значения приведенной скорости и β
plot([[Phic(beta, N_T) for beta in beta_table] for N_T in N_T_table],
     [[PR(beta, N_T) for beta in beta_table] for N_T in N_T_table],
     labels=reshape(["nᵣ=$N_T%" for N_T in N_T_table], 1, length(N_T_table)),
     xlabel = "phic",
     ylabel = "π")
plot!([[Phic(beta, N_T) for N_T in N_T_table] for beta in beta_table],
     [[PR(beta, N_T) for N_T in N_T_table] for beta in beta_table],
     linestyle=:dash,
     labels=reshape(["β=$beta%" for beta in beta_table], 1, length(beta_table)),
     legend=:outerbottomright)


@parameters t

@named cm = CompressorMap(
    Phic=Phic,
    Eta=Eta,
    PR=PR,
    N_T_bounds=(90.0,105.0))

@parameters pr N_T

ns = compose(ODESystem([cm.pr ~ pr
    cm.N_T ~ N_T], t, [], [pr, N_T]; name=:connected), cm)


guess = [
    cm.beta => 2.0,
    cm.eta => 0.83,
    cm.phic => 85
    ]


ps = [
    pr => 17.0
    N_T => 100.
    ]


sys = structural_simplify(ns)

prob = ODEProblem(sys, ps, (0, 100.0))

sol = solve(prob)

sol[cm.beta][end]


plot(sol, idxs=[cm.beta])

