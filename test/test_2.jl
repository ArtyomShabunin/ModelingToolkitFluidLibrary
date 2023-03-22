using Revise

using ModelingToolkit
using DifferentialEquations

import ModelingToolkitFluidLibrary: MassFlowSource_T, Source, LinearMassFlow, Node, Media.Air

using Plots


@parameters t

@named source_a = MassFlowSource_T(Medium=Air, T_in=293.0, m_flow_in=10)
@named source_b = Source(Medium=Air, P=2.3e5, T=310.0)
@named node = Node()
@named lin_flow = LinearMassFlow()

subs = [source_a; source_b; node; lin_flow]

flow_eqs = [
    connect(source_a.port, node.port_a)
    connect(node.port_b, lin_flow.port_a)
    connect(lin_flow.port_b, source_b.port)
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

sol = solve(prob, Rodas4(autodiff=false))

plot(sol, idxs=[node.P, node.h])