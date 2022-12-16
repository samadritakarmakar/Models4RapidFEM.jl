"""Function to create the stiffness matrix in Linear Elastic Problems"""
function ∇v_C_∇u!(K::Array{Float64, 2}, passedData::T, problemDim::Int64, element::AbstractElement, elementNo::Int64, shapeFunction::Array{ShapeFunction}, coordArray::Array{Float64, 2}; kwargs4function...) where T
    noOfIpPoints = getNoOfElementIpPoints(shapeFunction)
    noOfNodes = getNoOfElementNodes(shapeFunction)
    C = passedData
    for ipNo::Int64 = 1:noOfIpPoints
        ∂X_∂ξ = get_∂x_∂ξ(coordArray, shapeFunction, ipNo)
        dΩ = get_dΩ(element, ∂X_∂ξ, shapeFunction, ipNo)
        ∂ϕ_∂X = get_∂ϕ_∂x(element, ∂X_∂ξ, shapeFunction, ipNo)
        for l = 1:size(∂ϕ_∂X, 2)
            for j = 1:size(∂ϕ_∂X, 2)
                for b = 1:noOfNodes
                    for k = 1:problemDim
                        for a = 1:noOfNodes
                            for i = 1:problemDim
                                K[problemDim * (a - 1) + i, problemDim * (b - 1) + k] += 0.25 * ∂ϕ_∂X[a, j] * C[i, j, k, l] * ∂ϕ_∂X[b, l] * dΩ
                                K[problemDim * (a - 1) + j, problemDim * (b - 1) + l] += 0.25 * ∂ϕ_∂X[a, i] * C[i, j, k, l] * ∂ϕ_∂X[b, k] * dΩ
                                K[problemDim * (a - 1) + j, problemDim * (b - 1) + k] += 0.25 * ∂ϕ_∂X[a, i] * C[i, j, k, l] * ∂ϕ_∂X[b, l] * dΩ
                                K[problemDim * (a - 1) + i, problemDim * (b - 1) + l] += 0.25 * ∂ϕ_∂X[a, j] * C[i, j, k, l] * ∂ϕ_∂X[b, k] * dΩ
                            end
                        end
                    end
                end
            end
        end
    end
end

"""Function to stress at gauss points in Linear Elastic Problems."""
function gaussianStress(passedData::T, solAtNodes::Array{Float64, 1}, problemDim::Int64, element::AbstractElement, elementNo::Int64, shapeFunction::Array{ShapeFunction}, coordArray::Array{Float64, 2}; kwargs4function...) where T
    begin
        noOfIpPoints = getNoOfElementIpPoints(shapeFunction)
        #noOfNodes = getNoOfElementNodes(shapeFunction)
        C = passedData
        σ = zeros(noOfIpPoints, 3, 3)
        u_Nodes = solAtNodes
        for ipNo::Int64 = 1:noOfIpPoints
            ∂X_∂ξ = get_∂x_∂ξ(coordArray, shapeFunction, ipNo)
            ∂ϕ_∂X = get_∂ϕ_∂x(element, ∂X_∂ξ, shapeFunction, ipNo)
            ∂u_∂X = get_∂u_∂x(u_Nodes, ∂ϕ_∂X, Int64(length(u_Nodes) / size(∂ϕ_∂X, 1)))
            for l = 1:size(∂ϕ_∂X, 2)
                for k = 1:size(∂ϕ_∂X, 2)
                    for j = 1:problemDim
                        for i = 1:problemDim
                            σ[ipNo, i, j] += C[i, j, k, l] * 0.5 * (∂u_∂X[k, l] + ∂u_∂X[l, k])
                        end
                    end
                end
            end
        end
    end
    return σ
end

"""Function to strain at gauss points in Linear Elastic Problems."""
function gaussianStrain(passedData::T, solAtNodes::Array{Float64, 1}, problemDim::Int64, element::AbstractElement, elementNo::Int64, shapeFunction::Array{ShapeFunction}, coordArray::Array{Float64, 2}; kwargs4function...) where T
    begin
        noOfIpPoints = getNoOfElementIpPoints(shapeFunction)
        #noOfNodes = getNoOfElementNodes(shapeFunction)
        #(C, lastSol_u) = passedData
        ϵ = zeros(noOfIpPoints, 3, 3)
        u_Nodes = solAtNodes
        for ipNo::Int64 = 1:noOfIpPoints
            ∂X_∂ξ = get_∂x_∂ξ(coordArray, shapeFunction, ipNo)
            ∂ϕ_∂X = get_∂ϕ_∂x(element, ∂X_∂ξ, shapeFunction, ipNo)
            ∂u_∂X = get_∂u_∂x(u_Nodes, ∂ϕ_∂X, Int64(length(u_Nodes) / size(∂ϕ_∂X, 1)))
            for l = 1:problemDim
                for k = 1:problemDim
                    ϵ[ipNo, k, l] += 0.5 * (∂u_∂X[k, l] + ∂u_∂X[l, k])
                end
            end
        end
    end
    return ϵ
end

"""Function to create Elastic Tensor for Linear Elastic Isotropic Materials"""
function createElasticTensor(E::Float64, ν::Float64)
    λ = (ν * E) / ((1 + ν) * (1 - 2ν))
    μ = E / (2 * (1 + ν))
    C = λ * (one(SymmetricTensor{2, 3, Float64}) ⊗ one(SymmetricTensor{2, 3, Float64}))
    C += 2 * μ * one(SymmetricTensor{4, 3, Float64})
    return C
end
