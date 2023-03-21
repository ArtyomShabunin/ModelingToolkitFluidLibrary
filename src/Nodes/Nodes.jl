"""
"""

function Node(; name, Mh=1, Mp=1)
    pars = @parameters begin
        Mh = Mh
        Mp = Mp
    end

    vars = @variables P(t) = 1.0 Ha(t) = 1.0 Hb(t) = 1.0 h(t) = 1300

    @named port_a = FluidPort()
    @named port_b = FluidPort()

    subs = [port_a; port_b]

    D = Differential(t)

    eqns = Equation[]

    push!(eqns, P ~ port_a.P)
    push!(eqns, P ~ port_b.P)
    push!(eqns, Ha ~ IfElse.ifelse(port_a.m_flow > 0, port_a.m_flow * instream(port_a.h_outflow), port_a.m_flow * h))
    push!(eqns, Hb ~ IfElse.ifelse(port_b.m_flow > 0, port_b.m_flow * instream(port_b.h_outflow), port_b.m_flow * h))
    push!(eqns, Mh * D(h) ~ Ha + Hb)
    push!(eqns, Mp * D(P) ~ port_a.m_flow + port_b.m_flow)
    push!(eqns, port_a.h_outflow ~ h)
    push!(eqns, port_b.h_outflow ~ h)

    compose(ODESystem(eqns, t, vars, pars; name=name), subs)
end