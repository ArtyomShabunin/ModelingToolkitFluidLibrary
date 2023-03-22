using Revise

using ModelingToolkit
using DifferentialEquations

import ModelingToolkitFluidLibrary: Source, LinearMassFlow, Node, Media.Air

using Plots

@parameters t

@named source_a = Source(Medium=Air, P=2.5e5, T=300.0)
@named source_b = Source(Medium=Air, P=2.3e5, T=310.0)
@named lin_flow = LinearMassFlow()
@named node = Node()
@named lin_flow_2 = LinearMassFlow()

subs = [source_a; source_b; lin_flow; node; lin_flow_2]

flow_eqs = [
    connect(source_a.port, lin_flow.port_a)
    connect(lin_flow.port_b, node.port_a)
    connect(node.port_b, lin_flow_2.port_a)
    connect(lin_flow_2.port_b, source_b.port)
]

@named _flow_model = ODESystem(flow_eqs, t)

@named flow_model = compose(_flow_model, subs)

sys = structural_simplify(flow_model)

equations(sys)

u0 = [
    node.P => 2.31e5,
    node.h => 1300
]

prob = ODEProblem(sys, u0, (0, 10000.0))

sol = solve(prob, Rodas4())

plot(sol, idxs=[node.h])