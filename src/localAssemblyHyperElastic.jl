"""Material Stiffness tensor for hyper elastic problems of total lagrangian strain type."""
function δE_Cᵀ_ΔE!(𝕂::Array{Float64, 2}, passedData::T, problemDim::Int64, element::AbstractElement, elementNo::Int64, shapeFunction::Array{ShapeFunction}, coordArray::Array{Float64, 2}; kwargs4function...) where T
    noOfIpPoints = getNoOfElementIpPoints(shapeFunction)
    noOfNodes = getNoOfElementNodes(shapeFunction)
    (model, modelParams, lastSol_u) = passedData
    ∂ϕ_∂X_F_ℂ = zeros(3, noOfNodes, 3, 3)
    ∂ϕ_∂X_F = zeros(3, noOfNodes, 3, 3)
    u_Nodes = getSolAtElement(lastSol_u, element, problemDim)
    for ipNo::Int64 = 1:noOfIpPoints
        ∂X_∂ξ = get_∂x_∂ξ(coordArray, shapeFunction, ipNo)
        dΩ = get_dΩ(element, ∂X_∂ξ, shapeFunction, ipNo)
        ∂ϕ_∂X = get_∂ϕ_∂x(element, ∂X_∂ξ, shapeFunction, ipNo)
        ∂u_∂X = get_∂u_∂x(u_Nodes, ∂ϕ_∂X, Int64(length(u_Nodes) / size(∂ϕ_∂X, 1)))
        F = LargeDefs.getDeformationGradient(∂u_∂X)
        E = LargeDefs.getGreenLagrangeStrain(F)
        S = model.secondPiolaStress(E, modelParams)
        ℂ = model.materialTangentTensor(E, modelParams) * dΩ
        @einsum ∂ϕ_∂X_F[i, a, I, J] = ∂ϕ_∂X[a, I] * F[i, J]
        @einsum ∂ϕ_∂X_F_ℂ[i, a, K, L] = ∂ϕ_∂X_F[i, a, I, J] * ℂ[I, J, K, L]
        for K = 1:size(∂ϕ_∂X_F_ℂ, 3)
            for L = 1:size(∂ϕ_∂X_F_ℂ, 4)
                for b = 1:noOfNodes
                    for j = 1:problemDim
                        for a = 1:noOfNodes
                            for i = 1:problemDim
                                𝕂[problemDim * (a - 1) + i, problemDim * (b - 1) + j] += 0.5 * ∂ϕ_∂X_F_ℂ[i, a, K, L] * (∂ϕ_∂X_F[j, b, L, K] + ∂ϕ_∂X_F[j, b, K, L])
                            end
                        end
                    end
                end
            end
        end
    end
end

"""Geometric Stiffness tensor for hyper elastic problems of total lagrangian strain type."""
function δE_S_ΔE!(𝕂::Array{Float64, 2}, passedData::T, problemDim::Int64, element::AbstractElement, elementNo::Int64, shapeFunction::Array{ShapeFunction}, coordArray::Array{Float64, 2}; kwargs4function...) where T
    noOfIpPoints = getNoOfElementIpPoints(shapeFunction)
    noOfNodes = getNoOfElementNodes(shapeFunction)
    (model, modelParams, lastSol_u) = passedData
    ∂ϕ_∂X_S = zeros(noOfNodes, 3, 3)
    u_Nodes = getSolAtElement(lastSol_u, element, problemDim)
    for ipNo::Int64 = 1:noOfIpPoints
        ∂X_∂ξ = get_∂x_∂ξ(coordArray, shapeFunction, ipNo)
        dΩ = get_dΩ(element, ∂X_∂ξ, shapeFunction, ipNo)
        ∂ϕ_∂X = get_∂ϕ_∂x(element, ∂X_∂ξ, shapeFunction, ipNo)
        ∂ϕ_∂X = get_∂ϕ_∂x(element, ∂X_∂ξ, shapeFunction, ipNo)
        ∂u_∂X = get_∂u_∂x(u_Nodes, ∂ϕ_∂X, Int64(length(u_Nodes) / size(∂ϕ_∂X, 1)))
        F = LargeDefs.getDeformationGradient(∂u_∂X)
        E = LargeDefs.getGreenLagrangeStrain(F)
        S = model.secondPiolaStress(E, modelParams) * dΩ
        @einsum ∂ϕ_∂X_S[a, I, J] = ∂ϕ_∂X[a, I] * S[I, J]
        for J = 1:size(∂ϕ_∂X_S, 3)
            for b = 1:noOfNodes
                for I = 1:size(∂ϕ_∂X_S, 2)
                    for a = 1:noOfNodes
                        for i = 1:problemDim
                            𝕂[problemDim * (a - 1) + i, problemDim * (b - 1) + i] += ∂ϕ_∂X_S[a, I, J] * ∂ϕ_∂X[b, J]
                        end
                    end
                end
            end
        end
    end
end

"""Internal Force Vector for hyper elastic problems of total lagrangian strain type."""
function δE_S!(B::Vector, passedData::T, problemDim::Int64, element::AbstractElement, elementNo::Int64, shapeFunction::Array{ShapeFunction}, coordArray::Array{Float64, 2}; kwargs4function...) where T
    noOfIpPoints = getNoOfElementIpPoints(shapeFunction)
    noOfNodes = getNoOfElementNodes(shapeFunction)
    (model, modelParams, lastSol_u) = passedData
    F_S = zeros(3, 3)
    u_Nodes = getSolAtElement(lastSol_u, element, problemDim)
    for ipNo::Int64 = 1:noOfIpPoints
        ∂X_∂ξ = get_∂x_∂ξ(coordArray, shapeFunction, ipNo)
        dΩ = get_dΩ(element, ∂X_∂ξ, shapeFunction, ipNo)
        ∂ϕ_∂X = get_∂ϕ_∂x(element, ∂X_∂ξ, shapeFunction, ipNo)
        ∂u_∂X = get_∂u_∂x(u_Nodes, ∂ϕ_∂X, Int64(length(u_Nodes) / size(∂ϕ_∂X, 1)))
        F = LargeDefs.getDeformationGradient(∂u_∂X)
        E = LargeDefs.getGreenLagrangeStrain(F)
        S = model.secondPiolaStress(E, modelParams) * dΩ
        @einsum F_S[i, I] = F[i, J] * S[I, J]
        for I = 1:size(∂ϕ_∂X, 2)
            for a = 1:noOfNodes
                for i = 1:problemDim
                    B[problemDim * (a - 1) + i] += ∂ϕ_∂X[a, I] * F_S[i, I]
                end
            end
        end
    end
end