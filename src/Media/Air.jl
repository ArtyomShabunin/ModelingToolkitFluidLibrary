"""
"""

module Air

using ModelingToolkit, CoolProp

@parameters t
D = Differential(t)

hpt(P, T) = PropsSI("H", "P", P, "T", T, "Air")
hps(P, S) = PropsSI("H", "P", P, "S", S, "Air")
tph(P, H) = PropsSI("T", "P", P, "H", H, "Air")
dpt(P, T) = PropsSI("D", "P", P, "T", T, "Air")
dph(P, H) = PropsSI("D", "P", P, "H", H, "Air")
spt(P, T) = PropsSI("S", "P", P, "T", T, "Air")
sph(P, H) = PropsSI("S", "P", P, "H", H, "Air")

@register hpt(P, T)
@register hps(P, S)
@register tph(P, H)
@register dpt(P, T)
@register dph(P, H)
@register spt(P, T)
@register sph(P, H)

export setState_pT
function setState_pT(; name)

    vars = @variables P(t) T(t) h(t)

    eqns = Equation[]

    push!(eqns, h ~ hpt(P, T))

    ODESystem(eqns, t, vars, []; name=name)
end

export setState_ph
function setState_ph(; name)

    vars = @variables P(t) T(t) h(t)

    eqns = Equation[]

    push!(eqns, T ~ tph(P, h))

    ODESystem(eqns, t, vars, [],; name=name)
end

export specificEntropy
function specificEntropy(state)
    sph(state.P, state.h)
end

export isentropicEnthalpy
function isentropicEnthalpy(p, state)
    hps(p, sph(state.P, state.h))
end

end