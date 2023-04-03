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


function FlowNode(; name,
    Medium,
    Din,
    L,
    ke)

    pars = @parameters begin
        Din = Din
        L = L
        f = pi * Din^2 / 4
        ke = ke # абсолютная эквивалентная шероховатость (0.00014)
    end

    vars = @variables begin
        P(t)
        T(t)
        h(t)
        d(t)
        w(t)
        lambda_tr(t) # коэффициент трения
        Xi_flow(t) # коэффициент гидравлического сопротивления участка
        dp_fric(t) # потеря давления из-за сил трения
    end

    @named two_ports = TwoPorts()
    @unpack dp, dh, m_flow, port_a, port_b = two_ports

    @named state = Medium.setState_ph()

    subs = [two_ports, state]

    eqns = Equation[]

    push!(eqns, state.P ~ port_b.P)
    push!(eqns, state.h ~ instream(port_a.h_outflow))
    push!(eqns, w ~ m_flow / state.d / f)
    push!(eqns, lambda_tr ~ 1 / (1.14 + 2 * log10(Din / ke))^2)
    push!(eqns, Xi_flow ~ lambda_tr * L / Din)
    push!(eqns, dp_fric ~ abs(w) * w * Xi_flow * state.d / 2 / 9.81)
    push!(eqns, dp ~ dp_fric) # заменить на дифференциальное уравнение!!!???

    push!(eqns, dh ~ 0)

    extend(compose(ODESystem(eqns, t, vars, pars; name=name), state), two_ports)
end