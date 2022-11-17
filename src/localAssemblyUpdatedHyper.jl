"""Material Stiffness Matrix for hyper elastic problems of Updated Lagrangian type"""
function δD_𝕔_ΔD!(𝕂::Array{Float64, 2}, passedData::T, problemDim::Int64, element::AbstractElement, elementNo::Int64, shapeFunction::Array{ShapeFunction}, coordArray::Array{Float64, 2}; kwargs4function...) where T
    noOfIpPoints = getNoOfElementIpPoints(shapeFunction)
    noOfNodes = getNoOfElementNodes(shapeFunction)
    (hyperModel, modelParams, ipState_∂uₙ_∂X, ipState_∂uₖ_∂X, lastSol_Δu) = passedData
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
        Δf = LargeDefs.getDeformationGradient(∂Δu_∂xₙ)
        Fₙ = LargeDefs.getDeformationGradient(∂uₙ_∂X)
        ∂uₖ_∂X = ∂Δu_∂xₙ * Fₙ + ∂uₙ_∂X
        updateIpStateDict!(∂uₖ_∂X, ipState_∂uₖ_∂X, elementNo, ipNo)
        𝕔 = LargeDefs.spatialTangentTensor(hyperModel, Δf ⋅ Fₙ, modelParams)
        𝕔_Complete = Tensor{4, 3, Float64}(((i, j, k, l)->begin
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

"""Geometric Stiffness Matrix for hyper elastic problems of Updated Lagrangian type"""
function δu_σ_Δu!(K::Array{Float64, 2}, passedData::T, problemDim::Int64, element::AbstractElement, elementNo::Int64, shapeFunction::Array{ShapeFunction}, coordArray::Array{Float64, 2}; kwargs4function...) where T
    noOfIpPoints = getNoOfElementIpPoints(shapeFunction)
    noOfNodes = getNoOfElementNodes(shapeFunction)
    (hyperModel, modelParams, ipState_∂uₙ_∂X, ipState_∂uₖ_∂X, lastSol_Δu) = passedData
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
        Δf = LargeDefs.getDeformationGradient(∂Δu_∂xₙ)
        Fₙ = LargeDefs.getDeformationGradient(∂uₙ_∂X)
        σ = LargeDefs.cauchyStress(hyperModel, Δf ⋅ Fₙ, modelParams)
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

"""Internal Force Vector for hyper elastic problems of Updated Lagrangian type"""
function δD_σ!(f::Vector, passedData::T, problemDim::Int64, element::AbstractElement, elementNo::Int64, shapeFunction::Array{ShapeFunction}, coordArray::Array{Float64, 2}; kwargs4function...) where T
    noOfIpPoints = getNoOfElementIpPoints(shapeFunction)
    noOfNodes = getNoOfElementNodes(shapeFunction)
    (hyperModel, modelParams, ipState_∂uₙ_∂X, ipState_∂uₖ_∂X, lastSol_Δu) = passedData
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
        Δf = LargeDefs.getDeformationGradient(∂Δu_∂xₙ)
        Fₙ = LargeDefs.getDeformationGradient(∂uₙ_∂X)
        σ = LargeDefs.cauchyStress(hyperModel, Δf ⋅ Fₙ, modelParams)
        for j = 1:size(∂ϕ_∂x, 2)
            for a = 1:noOfNodes
                for i = 1:problemDim
                    f[problemDim * (a - 1) + i] += 0.5 * ∂ϕ_∂x[a, j] * (σ[i, j] + σ[i, j]) * dω
                end
            end
        end
    end
end
