using Revise
using ModelingToolkit
using DifferentialEquations
using OrdinaryDiffEq: ReturnCode.Success
import ModelingToolkitFluidLibrary: Source, FlowNode, Media.Water
using Plots
using Test

# @testset "FlowNode demo" begin

@parameters t

@named source_a = Source(Medium=Water, P=2.5e5, T=30.0)
@named source_b = Source(Medium=Water, P=2.3e5, T=31.0)
@named flow_node = FlowNode(Medium=Water, Din=0.05, L=5, ke=0.00014)

subs = [source_a; source_b; flow_node]

flow_eqs = [
    connect(source_a.port, flow_node.port_a)
    connect(flow_node.port_b, source_b.port)
]

@named _flow_model = ODESystem(flow_eqs, t)

@named flow_model = compose(_flow_model, subs)

sys = structural_simplify(flow_model)

equations(sys)

u0 = [
]

prob = ODEProblem(sys, u0, (0, 10000.0))

sol = solve(prob, Rodas4())

@test sol.retcode == Success

# plot(sol, idxs=[node.h])

# end