module ModelingToolkitFluidLibrary

using ModelingToolkit, CoolProp, IfElse

@parameters t
D = Differential(t)

hpt(P, T) = PropsSI("H", "P", P, "T", T, "Air")
@register hpt(P, T)

export FluidPort, TwoPorts
include("Ports/Ports.jl")

export Source, MassFlowSource_T
include("Sources/Sources.jl")

export LinearMassFlow
include("Pipes/Pipes.jl")

export Node
include("Nodes/Nodes.jl")

include("Media/Media.jl")

end # module
