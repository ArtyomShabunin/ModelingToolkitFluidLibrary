"""

"""

function Source(; name,
    P=1.013e5,
    T=293.0)
    pars = @parameters begin
        P = P
        T = T
    end

    @named port = FluidPort()

    subs = [port]

    eqns = Equation[]

    push!(eqns, port.P ~ P)
    push!(eqns, port.h_outflow ~ hpt(P, T))

    compose(ODESystem(eqns, t, [], pars; name=name), subs)
end

function MassFlowSource_T(; name,
    T_in=293.0,
    m_flow_in=-0.01)
    pars = @parameters begin
        T_in = T_in
        m_flow_in = m_flow_in
    end

    vars = @variables begin
        P(t)
    end

    @named port = FluidPort()

    subs = [port]

    eqns = Equation[]

    push!(eqns, port.P ~ P)
    push!(eqns, port.m_flow ~ -m_flow_in)
    push!(eqns, port.h_outflow ~ hpt(P, T_in))

    compose(ODESystem(eqns, t, vars, pars; name=name), subs)
end