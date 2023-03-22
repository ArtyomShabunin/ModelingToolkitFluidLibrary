"""
"""

module Air

using ModelingToolkit, CoolProp

@parameters t
D = Differential(t)

hpt(P, T) = PropsSI("H", "P", P, "T", T, "Air")
@register hpt(P, T)

export setState_pT
function setState_pT(; name)

    vars = @variables P(t) T(t) h(t)

    eqns = Equation[]

    push!(eqns, h ~ hpt(P, T))

    ODESystem(eqns, t, vars, []; name=name)
end

end