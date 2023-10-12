using Revise
using ModelingToolkit
using ModelingToolkitStandardLibrary.Mechanical.Rotational
import ModelingToolkitStandardLibrary.Blocks
using DifferentialEquations
using Interpolations
using OrdinaryDiffEq: ReturnCode.Success
import ModelingToolkitFluidLibrary: Source, LinearMassFlow, Node, Media.Air, TwoPorts
using Plots


EtaC = [
    0.825 0.820 0.813 0.8115 0.81 0.805
    0.840 0.835 0.831 0.8301 0.827 0.82
    0.835 0.832 0.828 0.823 0.821 0.815
    0.825 0.820 0.8146 0.808 0.80 0.79
    0.795 0.790 0.783 0.777 0.772 0.765
]

PhicC = [
    68.534 73.0 76.64 79.14 82.14 85.864
    75 78.514 81.924 84.924 87.5 91.942
    78 81.788 85.5 88.2 90.5 94.5
    79.768 84.2 88.178 91 93.178 96
    80.454 85.454 89.668 92.668 94.668 97.314
]

PRC = [
    12.99 14.39 15.85 16.85 18 20.1
    12.5 14 15.59 16.69 17.7 19.62
    11.9 13.52 15 16.2 17.175 18.9
    10.65 12.5 14.095 15.095 16.095 17.7
    9.35 10.75 12.275 13.275 14.275 16.0
]


β = [1; 2; 3; 4; 5]
N_T = [0; 92; 95; 97; 100; 105]

interp_eta = interpolate((β, N_T), EtaC, Gridded(Linear()))
interp_phi = interpolate((β, N_T), PhicC, Gridded(Linear()))
interp_pr = interpolate((β, N_T), PRC, Gridded(Linear()))

interp_eta_f(β, N_T) = interp_eta(β, N_T)
interp_phi_f(β, N_T) = interp_phi(β, N_T)
interp_pr_f(β, N_T) = interp_pr(β, N_T)

@register_symbolic interp_eta_f(β, N_T)
@register_symbolic interp_phi_f(β, N_T)
@register_symbolic interp_pr_f(β, N_T)

function Compressor(; name,
                    Medium, phi_f, eta_f, pr_f,
                    N_T_design)

    pars = @parameters begin
        N_T_design = N_T_design # Referred design velocity
    end

    vars = @variables begin
        phi(t)
        phic(t)
        eta(t)
        pr(t)
        β(t) = 2.0
        N_T(t) = 95
        hout_iso(t)
        omega(t) = 100
        pr_2(t)
    end

    @named two_ports = TwoPorts()
    @unpack dp, dh, m_flow, port_a, port_b  = two_ports

    # Mechanical TwoPorts
    @named flange_a = Flange()
    @named flange_b = Flange()

    @named gas_in = Medium.setState_ph()

    subs = [flange_a, flange_b, gas_in]

    D = Differential(t)

    eqns = Equation[]

    # w это m_flow

    push!(eqns, N_T ~ 100*omega/sqrt(gas_in.T/290)/N_T_design)

    push!(eqns, gas_in.P ~ port_a.P)
    push!(eqns, gas_in.h ~ instream(port_a.h_outflow))

    push!(eqns, hout_iso ~ Medium.isentropicEnthalpy(port_b.P, gas_in))
    push!(eqns, m_flow ~ phic)
    push!(eqns, dh ~ 1/eta*(hout_iso - gas_in.h))

    push!(eqns, phi ~ flange_a.phi)
    push!(eqns, phi ~ flange_b.phi)
    push!(eqns, 0 ~ flange_a.tau + flange_b.tau)

    push!(eqns, 1e1*D(β) ~ pr-pr_2)
    # push!(eqns, pr_2 ~ 1)
    push!(eqns, pr_2 ~ port_b.P / port_a.P)

    push!(eqns, phic ~ phi_f(β, N_T))
    push!(eqns, eta ~ eta_f(β, N_T))
    push!(eqns, pr ~ pr_f(β, N_T))

    push!(eqns, D(phi) ~ omega)

    extend(compose(ODESystem(eqns, t, vars, pars; name=name), subs), two_ports)
end

@parameters t

@named torque = Torque(; use_support=false)
@named speed = Speed(; use_support=false)
@named compressor = Compressor(; Medium=Air, phi_f=interp_phi_f, eta_f=interp_eta_f, pr_f=interp_pr_f, N_T_design=100)
@named torque_source = Blocks.Constant(; k=100)
@named speed_source = Blocks.Constant(; k=100)
@named source_a = Source(Medium=Air, P=2.5e5, T=30.0 + 273.15)
@named source_b = Source(Medium=Air, P=34e5, T=31.0 + 273.15)

connections = [
    connect(torque_source.output, torque.tau),
    connect(torque.flange, compressor.flange_a),
    connect(compressor.flange_b, speed.flange),
    connect(speed.w_ref, speed_source.output),
    connect(source_a.port, compressor.port_a),
    connect(compressor.port_b, source_b.port)
]

@named model = ODESystem(connections, t,
    systems=[
        torque,
        speed,
        torque_source,
        speed_source,
        compressor,
        source_a,
        source_b
    ])

sys = structural_simplify(model)
prob = ODAEProblem(sys, [], (0, 100.0))
sol = solve(prob, Tsit5())

sol[compressor.pr][end]

sol[t][29]

plot(sol, idxs=[compressor.β])