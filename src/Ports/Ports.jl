"""

"""
module Ports
using ModelingToolkit, CoolProp, IfElse

@parameters t
D = Differential(t)

hpt(P, T) = PropsSI("H", "P", P, "T", T, "Air")
@register hpt(P, T)

export FluidPort
@connector function FluidPort(; name, P=0.0, m_flow=0.0, h_outflow=0.0)
    vars = @variables h_outflow(t) = h_outflow [connect = Stream] m_flow(t) = m_flow [
        connect = Flow,
    ] P(t) = P
    ODESystem(Equation[], t, vars, []; name=name)
end

export TwoPorts
function TwoPorts(; name)
    vars = @variables dp(t) dh(t) m_flow(t)

    @named port_a = FluidPort()
    @named port_b = FluidPort()

    subs = [port_a; port_b]

    eqns = Equation[]

    push!(eqns, dp ~ port_a.P - port_b.P)
    push!(eqns, port_b.h_outflow ~ instream(port_a.h_outflow) + dh)
    push!(eqns, port_a.h_outflow ~ instream(port_b.h_outflow) + dh)
    push!(eqns, port_a.m_flow ~ port_b.m_flow)
    push!(eqns, m_flow ~ port_a.m_flow)

    compose(ODESystem(eqns, t, vars, []; name=name), subs)
end

end