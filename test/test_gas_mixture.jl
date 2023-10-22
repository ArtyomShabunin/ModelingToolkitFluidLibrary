using ModelingToolkit, CoolProp, IfElse


# @connector function FluidGasPort(; name, P=0.0, h_outflow=0.0)
#     vars = @variables h_outflow(t) = h_outflow [connect = Stream] m_flow(t) = mflow [connect = Flow] P(t) = P Xi_outflow[1:6] [connect = Stream]
#     ODESystem(Equation[], t, vars, []; name=name)
# end

PropsSI("H", "T", 300, "P", 101325,
        "Oxygen[0.23]&Argon[0.02]&Water[0.01]&CO2[0.04]&Nitrogen[0.7]")

PropsSI("H", "T", 300, "P", 101325, "Air")

HAPropsSI("H", "T", 300, "RH", 0.0, "P", 101325)