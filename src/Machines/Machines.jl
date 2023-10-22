"""
"""

function CompressorMap(; name, Phic, Eta, PR, N_T_bounds, beta_initial=3)
    
    vars = @variables begin
        beta(t) = beta_initial #, [description = "Number of beta line"]
        eta(t), [description = "Isentropic efficiency"]
        phic(t), [description = "Flow number"]
        N_Tb(t), [description = "Bounded referred speed"]
        N_T(t), [description = "Referred speed"]
        pr(t), [description = "Pressure ratio"]
    end
    
    eqs = [
        N_Tb ~ min(max(N_T, N_T_bounds[1]), N_T_bounds[2]),
        eta ~ Eta(beta, N_Tb),
        phic ~ Phic(beta, N_Tb),
        pr ~ PR(beta, N_Tb)
        ]

    ODESystem(eqs, t, vars, []; name=name)    
end
    
function Compressor(; name,
    Medium, cm,
    N_T_design)

    pars = @parameters begin
    N_T_design = N_T_design, [description = "Referred design velocity"]
    end

    vars = @variables begin
    phi(t)
    N_T(t)
    hout_iso(t)
    omega(t)
    end

    @named two_ports = TwoPorts()
    @unpack dp, dh, m_flow, port_a, port_b  = two_ports

    @named flange_a = Flange()
    @named flange_b = Flange()

    @named gas_in = Medium.setState_ph()

    subs = [flange_a, flange_b, gas_in, cm]

    D = Differential(t)

    eqns = Equation[]

    push!(eqns, cm.N_T ~ 100*omega/sqrt(gas_in.T/290)/N_T_design)    

    push!(eqns, gas_in.P ~ port_a.P)
    push!(eqns, gas_in.h ~ instream(port_a.h_outflow))

    push!(eqns, hout_iso ~ Medium.isentropicEnthalpy(port_b.P, gas_in))

    push!(eqns, m_flow ~ cm.phic / sqrt(gas_in.T/290)*1.013e5/(gas_in.P))

    push!(eqns, dh ~ 1/cm.eta*(hout_iso - gas_in.h))

    push!(eqns, phi ~ flange_a.phi)
    push!(eqns, phi ~ flange_b.phi)
    push!(eqns, 0 ~ flange_a.tau + flange_b.tau)

    push!(eqns, cm.pr ~ port_b.P / port_a.P)

    push!(eqns, D(phi) ~ omega)

    extend(compose(ODESystem(eqns, t, vars, pars; name=name), subs), two_ports)
end