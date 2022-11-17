module Models4RapidFEM
using RapidFEM, Tensors, LargeDefs, Einsum
include("localBodyExternalForce.jl")
include("localAssemblyLinearElastic.jl")
include("localAssemblyHyperElastic.jl")
include("localAssemblyPoisson.jl")
include("localAssemblyUpdatedHyper.jl")
include("localAssemblyUpdatedHyperPlaneStrain2.jl")
end # module Models4RapidFEM
