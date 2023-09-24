using Revise
using ModelingToolkit
using ModelingToolkitStandardLibrary.Mechanical.Rotational
import ModelingToolkitStandardLibrary.Blocks
using DifferentialEquations
using OrdinaryDiffEq: ReturnCode.Success
import ModelingToolkitFluidLibrary: Source, LinearMassFlow, Node, Media.Air
using Plots

@parameters t

@named torque = Torque(; use_support = false)
@named speed = Speed(; use_support = false)
@named torque_source = Blocks.Constant(; k=100)
@named speed_source = Blocks.Constant(; k=100)

connections = [
    connect(torque_source.output, torque.tau),
    connect(torque.flange, speed.flange),
    connect(speed.w_ref, speed_source.output)
]

@named model = ODESystem(connections, t,
        systems = [
            torque,
            speed,
            torque_source,
            speed_source
        ])

sys = structural_simplify(model)
prob = ODEProblem(sys, [], (0, 100.0))
sol = solve(prob, Rodas4())

plot(sol, idxs=[torque.flange.phi])