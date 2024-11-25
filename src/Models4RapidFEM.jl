module Models4RapidFEM
using RapidFEM, Tensors, LargeDefs, Einsum
include("localBodyExternalForce.jl")
include("localAssemblyLinearElastic.jl")
include("localAssemblyHyperElastic.jl")
include("localAssemblyPoisson.jl")
include("localAssemblyUpdatedHyper.jl")
include("localAssemblyUpdatedHyperPlaneStrain2.jl")
include("localPressureProjectionV2.jl")
#From localBodyExternalForce
export v_S!, v_F!
#from localAssemblyLinearElastic
export ∇v_C_∇u!, gaussianStress, gaussianStrain, gaussTwiceLinStrainEnergy, createElasticTensor
export createPlaneStressElasticTensor
#from localAssemblyHyperElastic
export δE_Cᵀ_ΔE!, δE_S_ΔE!, δE_S!
#from localAssemblyPoisson
export ∇v_∇u!
#from localAssemblyUpdatedHyper
export δD_𝕔_ΔD!, δu_σ_Δu!, δD_σ!
#from localPressureProjectionV2
export getPpp_V2local
#from localAssemblyUpdatedHyperPlaneStrain2
export δD_𝕔_ΔD_PlaneStrain!, δu_σ_Δu_PlaneStrain!, δD_σ_PlaneStrain!
export gaussian_σ, gaussian_b, gaussian_DispGrad, gaussian_DefGrad
export StrainEnergy
end # module Models4RapidFEM
