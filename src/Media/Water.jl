"""
"""

module Water

using ModelingToolkit, CoolProp

CoolProp

@parameters t
D = Differential(t)

hpt(P, T) = PropsSI("H", "P", P, "T", T, "Water")
@register hpt(P, T)

tph(P, H) = PropsSI("T", "P", P, "H", H, "Water")
@register tph(P, T)

dph(P, H) = PropsSI("D", "P", P, "H", H, "Water")
@register dph(P, H)

export setState_pT
function setState_pT(; name)

    vars = @variables P(t) T(t) h(t)

    eqns = Equation[]

    push!(eqns, h ~ hpt(P, T))

    ODESystem(eqns, t, vars, []; name=name)
end

export setState_ph
function setState_ph(; name)

    vars = @variables P(t) T(t) h(t) d(t)

    eqns = Equation[]

    push!(eqns, T ~ tph(P, h))
    push!(eqns, d ~ dph(P, h))

    ODESystem(eqns, t, vars, []; name=name)
end


end