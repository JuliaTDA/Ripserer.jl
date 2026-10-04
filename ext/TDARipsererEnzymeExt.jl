module TDARipsererEnzymeExt
using TDARipserer
import Enzyme: EnzymeRules

# Freeze the combinatorial choice only. Endpoint arithmetic remains active.
EnzymeRules.inactive(::typeof(TDARipserer._persistence_pairing),args...) = true
end
