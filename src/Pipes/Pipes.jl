"""
"""

function LinearMassFlow(; name)
    @named two_ports = TwoPorts()
    @unpack dp, dh, m_flow = two_ports

    eqns = Equation[]

    push!(eqns, m_flow ~ 1e-3 * dp)
    push!(eqns, dh ~ 0)

    extend(ODESystem(eqns, t, [], []; name=name), two_ports)
end