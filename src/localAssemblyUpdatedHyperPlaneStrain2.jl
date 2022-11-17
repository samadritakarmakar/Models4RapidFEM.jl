"""Material stiffness matrix for plane strain problems, of updated lagrangian type."""
function δD_𝕔_ΔD_PlaneStrain!(𝕂::Array{Float64, 2}, passedData::T, problemDim::Int64, element::AbstractElement, elementNo::Int64, shapeFunction::Array{ShapeFunction}, coordArray::Array{Float64, 2}; kwargs4function...) where T
    noOfIpPoints = getNoOfElementIpPoints(shapeFunction)
    noOfNodes = getNoOfElementNodes(shapeFunction)
    (hyperModel, modelParams, ipState_∂uₙ_∂X, ipState_∂uₖ_∂X, lastSol_Δu) = passedData
    ∂Δu_∂xₙ_plnStrn = zeros(3, 3)
    ∂uₙ_∂X = zeros(3, 3)
    Δu_Nodes = getSolAtElement(lastSol_Δu, element, problemDim)
    currentCoordArray = getCurrentCoordArray(coordArray, Δu_Nodes)
    for ipNo::Int64 = 1:noOfIpPoints
        ∂x_∂ξ = get_∂x_∂ξ(currentCoordArray, shapeFunction, ipNo)
        ∂X_∂ξ = get_∂x_∂ξ(coordArray, shapeFunction, ipNo)
        dω = get_dΩ(element, ∂x_∂ξ, shapeFunction, ipNo)
        ∂ϕ_∂X = get_∂ϕ_∂x(element, ∂X_∂ξ, shapeFunction, ipNo)
        ∂ϕ_∂x = get_∂ϕ_∂x(element, ∂x_∂ξ, shapeFunction, ipNo)
        ∂Δu_∂xₙ = get_∂u_∂x(Δu_Nodes, ∂ϕ_∂X, Int64(length(Δu_Nodes) / size(∂ϕ_∂X, 1)))
        SaveFemIpState.getIpState!(∂uₙ_∂X, ipState_∂uₙ_∂X, elementNo, ipNo)
        fill!(∂Δu_∂xₙ_plnStrn, 0.0)
        ∂Δu_∂xₙ_plnStrn[1:2, 1:2] = ∂Δu_∂xₙ
        Δf = LargeDefs.getDeformationGradient(∂Δu_∂xₙ_plnStrn)
        Fₙ = LargeDefs.getDeformationGradient(∂uₙ_∂X)
        ∂uₖ_∂X = ∂Δu_∂xₙ_plnStrn * Fₙ + ∂uₙ_∂X
        updateIpStateDict!(∂uₖ_∂X, ipState_∂uₖ_∂X, elementNo, ipNo)
        𝕔3d = LargeDefs.spatialTangentTensor(hyperModel, Δf ⋅ Fₙ, modelParams)
        𝕔 = @view(𝕔3d[1:2, 1:2, 1:2, 1:2])
        𝕔_Complete = Tensor{4, 2, Float64}(((i, j, k, l)->begin
                        𝕔[i, j, k, l] + 𝕔[j, i, k, l] + 𝕔[i, j, l, k] + 𝕔[j, i, l, k]
                    end))
        for l = 1:size(𝕔_Complete, 4)
            for j = 1:size(∂ϕ_∂x, 2)
                for b = 1:noOfNodes
                    for k = 1:problemDim
                        for a = 1:noOfNodes
                            for i = 1:problemDim
                                𝕂[problemDim * (a - 1) + i, problemDim * (b - 1) + k] += 0.25 * ∂ϕ_∂x[a, j] * 𝕔_Complete[i, j, k, l] * ∂ϕ_∂x[b, l] * dω
                            end
                        end
                    end
                end
            end
        end
    end
end

"""Geometric stiffness matrix for plane strain problems, of updated lagrangian type."""
function δu_σ_Δu_PlaneStrain!(K::Array{Float64, 2}, passedData::T, problemDim::Int64, element::AbstractElement, elementNo::Int64, shapeFunction::Array{ShapeFunction}, coordArray::Array{Float64, 2}; kwargs4function...) where T
    noOfIpPoints = getNoOfElementIpPoints(shapeFunction)
    noOfNodes = getNoOfElementNodes(shapeFunction)
    (hyperModel, modelParams, ipState_∂uₙ_∂X, ipState_∂uₖ_∂X, lastSol_Δu) = passedData
    ∂Δu_∂xₙ_plnStrn = zeros(3, 3)
    ∂uₙ_∂X = zeros(3, 3)
    Δu_Nodes = getSolAtElement(lastSol_Δu, element, problemDim)
    currentCoordArray = getCurrentCoordArray(coordArray, Δu_Nodes)
    for ipNo::Int64 = 1:noOfIpPoints
        ∂x_∂ξ = get_∂x_∂ξ(currentCoordArray, shapeFunction, ipNo)
        ∂X_∂ξ = get_∂x_∂ξ(coordArray, shapeFunction, ipNo)
        dω = get_dΩ(element, ∂x_∂ξ, shapeFunction, ipNo)
        ∂ϕ_∂X = get_∂ϕ_∂x(element, ∂X_∂ξ, shapeFunction, ipNo)
        ∂ϕ_∂x = get_∂ϕ_∂x(element, ∂x_∂ξ, shapeFunction, ipNo)
        ∂Δu_∂xₙ = get_∂u_∂x(Δu_Nodes, ∂ϕ_∂X, Int64(length(Δu_Nodes) / size(∂ϕ_∂X, 1)))
        SaveFemIpState.getIpState!(∂uₙ_∂X, ipState_∂uₙ_∂X, elementNo, ipNo)
        fill!(∂Δu_∂xₙ_plnStrn, 0.0)
        ∂Δu_∂xₙ_plnStrn[1:2, 1:2] = ∂Δu_∂xₙ
        Δf = LargeDefs.getDeformationGradient(∂Δu_∂xₙ_plnStrn)
        Fₙ = LargeDefs.getDeformationGradient(∂uₙ_∂X)
        σ3d = LargeDefs.cauchyStress(hyperModel, Δf ⋅ Fₙ, modelParams)
        σ = @view(σ3d[1:2, 1:2])
        for k = 1:size(∂ϕ_∂x, 2)
            for b = 1:noOfNodes
                for l = 1:size(∂ϕ_∂x, 2)
                    for a = 1:noOfNodes
                        for j = 1:problemDim
                            K[problemDim * (a - 1) + j, problemDim * (b - 1) + j] += ∂ϕ_∂x[a, l] * ∂ϕ_∂x[b, k] * σ[l, k] * dω
                        end
                    end
                end
            end
        end
    end
end

"""Internal force vector for plane strain problems, of updated lagrangian type."""
function δD_σ_PlaneStrain!(f::Vector, passedData::T, problemDim::Int64, element::AbstractElement, elementNo::Int64, shapeFunction::Array{ShapeFunction}, coordArray::Array{Float64, 2}; kwargs4function...) where T
    noOfIpPoints = getNoOfElementIpPoints(shapeFunction)
    noOfNodes = getNoOfElementNodes(shapeFunction)
    (hyperModel, modelParams, ipState_∂uₙ_∂X, ipState_∂uₖ_∂X, lastSol_Δu) = passedData
    ∂Δu_∂xₙ_plnStrn = zeros(3, 3)
    ∂uₙ_∂X = zeros(3, 3)
    Δu_Nodes = getSolAtElement(lastSol_Δu, element, problemDim)
    currentCoordArray = getCurrentCoordArray(coordArray, Δu_Nodes)
    for ipNo::Int64 = 1:noOfIpPoints
        ∂x_∂ξ = get_∂x_∂ξ(currentCoordArray, shapeFunction, ipNo)
        ∂X_∂ξ = get_∂x_∂ξ(coordArray, shapeFunction, ipNo)
        dω = get_dΩ(element, ∂x_∂ξ, shapeFunction, ipNo)
        ∂ϕ_∂X = get_∂ϕ_∂x(element, ∂X_∂ξ, shapeFunction, ipNo)
        ∂ϕ_∂x = get_∂ϕ_∂x(element, ∂x_∂ξ, shapeFunction, ipNo)
        ∂Δu_∂xₙ = get_∂u_∂x(Δu_Nodes, ∂ϕ_∂X, Int64(length(Δu_Nodes) / size(∂ϕ_∂X, 1)))
        SaveFemIpState.getIpState!(∂uₙ_∂X, ipState_∂uₙ_∂X, elementNo, ipNo)
        fill!(∂Δu_∂xₙ_plnStrn, 0.0)
        ∂Δu_∂xₙ_plnStrn[1:2, 1:2] = ∂Δu_∂xₙ
        Δf = LargeDefs.getDeformationGradient(∂Δu_∂xₙ_plnStrn)
        Fₙ = LargeDefs.getDeformationGradient(∂uₙ_∂X)
        σ3d = LargeDefs.cauchyStress(hyperModel, Δf ⋅ Fₙ, modelParams)
        σ = @view(σ3d[1:2, 1:2])
        for j = 1:size(∂ϕ_∂x, 2)
            for a = 1:noOfNodes
                for i = 1:problemDim
                    f[problemDim * (a - 1) + i] += 0.5 * ∂ϕ_∂x[a, j] * (σ[i, j] + σ[i, j]) * dω
                end
            end
        end
    end
end

"""Cauchy stress at gauss point."""
function gaussian_σ(passedData::T, solAtNodes::Array{Float64, 1}, problemDim::Int64, element::AbstractElement, elementNo::Int64, shapeFunction::Array{ShapeFunction}, coordArray::Array{Float64, 2}; kwargs4function...) where T
    begin
        noOfIpPoints = getNoOfElementIpPoints(shapeFunction)
        noOfNodes = getNoOfElementNodes(shapeFunction)
        (ipState, hyperModel, modelParams) = passedData
        σ_g = zeros(noOfIpPoints, 3, 3)
        for ipNo::Int64 = 1:noOfIpPoints
            ∂uₖ_∂X = SaveFemIpState.getIpState(ipState, elementNo, ipNo)
            F = LargeDefs.getDeformationGradient(∂uₖ_∂X)
            σ = LargeDefs.cauchyStress(hyperModel, F, modelParams)
            for j = 1:problemDim
                for i = 1:problemDim
                    σ_g[ipNo, i, j] = σ[i, j]
                end
            end
        end
    end
    return σ_g
end

"""Left Cauchy Tensor at gauss point."""
function gaussian_b(passedData::T, solAtNodes::Array{Float64, 1}, problemDim::Int64, element::AbstractElement, elementNo::Int64, shapeFunction::Array{ShapeFunction}, coordArray::Array{Float64, 2}; kwargs4function...) where T
    begin
        noOfIpPoints = getNoOfElementIpPoints(shapeFunction)
        noOfNodes = getNoOfElementNodes(shapeFunction)
        (ipState, hyperModel, modelParams) = passedData
        b_g = zeros(noOfIpPoints, 3, 3)
        for ipNo::Int64 = 1:noOfIpPoints
            ∂uₖ_∂X = SaveFemIpState.getIpState(ipState, elementNo, ipNo)
            F = LargeDefs.getDeformationGradient(∂uₖ_∂X)
            b = getLeftCauchyTensor(F)
            for j = 1:problemDim
                for i = 1:problemDim
                    b_g[ipNo, i, j] = b[i, j]
                end
            end
        end
    end
    return b_g
end

"""Displacement gradient at gauss point."""
function gaussian_DispGrad(passedData::T, solAtNodes::Array{Float64, 1}, problemDim::Int64, element::AbstractElement, elementNo::Int64, shapeFunction::Array{ShapeFunction}, coordArray::Array{Float64, 2}; kwargs4function...) where T
    begin
        noOfIpPoints = getNoOfElementIpPoints(shapeFunction)
        noOfNodes = getNoOfElementNodes(shapeFunction)
        (ipState, hyperModel, modelParams) = passedData
        ∂u_∂X_g = zeros(noOfIpPoints, 3, 3)
        for ipNo::Int64 = 1:noOfIpPoints
            ∂uₖ_∂X = SaveFemIpState.getIpState(ipState, elementNo, ipNo)
            for j = 1:problemDim
                for i = 1:problemDim
                    ∂u_∂X_g[ipNo, i, j] = ∂uₖ_∂X[i, j]
                end
            end
        end
    end
    return ∂u_∂X_g
end

"""Deformation gradient at gauss point."""
function gaussian_DefGrad(passedData::T, solAtNodes::Array{Float64, 1}, problemDim::Int64, element::AbstractElement, elementNo::Int64, shapeFunction::Array{ShapeFunction}, coordArray::Array{Float64, 2}; kwargs4function...) where T
    begin
        noOfIpPoints = getNoOfElementIpPoints(shapeFunction)
        noOfNodes = getNoOfElementNodes(shapeFunction)
        (ipState, hyperModel, modelParams) = passedData
        F_g = zeros(noOfIpPoints, 3, 3)
        for ipNo::Int64 = 1:noOfIpPoints
            ∂uₖ_∂X = SaveFemIpState.getIpState(ipState, elementNo, ipNo)
            F = LargeDefs.getDeformationGradient(∂uₖ_∂X)
            for j = 1:problemDim
                for i = 1:problemDim
                    F_g[ipNo, i, j] = F[i, j]
                end
            end
        end
    end
    return F_g
end

"""Strain Energy at gauss point."""
function StrainEnergy(StrainE::Array{Float64, 1}, passedData::T, problemDim::Int64, element::AbstractElement, elementNo::Int64, shapeFunction::Array{ShapeFunction}, coordArray::Array{Float64, 2}; kwargs4function...) where T
    noOfIpPoints = getNoOfElementIpPoints(shapeFunction)
    noOfNodes = getNoOfElementNodes(shapeFunction)
    (ipState, hyperModel, modelParams) = passedData
    for ipNo::Int64 = 1:noOfIpPoints
        ∂X_∂ξ = get_∂x_∂ξ(coordArray, shapeFunction, ipNo)
        dΩ = get_dΩ(element, ∂X_∂ξ, shapeFunction, ipNo)
        ∂uₙ_∂X = SaveFemIpState.getIpState(ipState, elementNo, ipNo)
        F = LargeDefs.getDeformationGradient(∂uₙ_∂X)
        E = LargeDefs.getGreenLagrangeStrain(F)
        Ψ_local = hyperModel.strainEnergyDensity(E, modelParams)
        Jacobian = LargeDefs.getJacobianDeformationGradient(F)
        dωₙ = dΩ
        for i = 1:problemDim
            StrainE[i] += if i == 1
                    Ψ_local * ((1.0 / Jacobian) * dωₙ)
                else
                    0.0
                end
        end
    end
end
