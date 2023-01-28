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

"""Function to find stress at gauss points in Linear Elastic Problems."""
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

"""Function to find strain at gauss points in Linear Elastic Problems."""
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


"""Function to find twice strain energy at gauss points in Linear Elastic Problems."""
function gaussTwiceLinStrainEnergy(passedData::T, solAtNodes::Array{Float64, 1}, problemDim::Int64, 
    element::AbstractElement, elementNo::Int64, shapeFunction::Array{ShapeFunction}, 
    coordArray::Array{Float64, 2}; kwargs4function...) where T

    noOfIpPoints = getNoOfElementIpPoints(shapeFunction)
    #noOfNodes = getNoOfElementNodes(shapeFunction)
    C = passedData
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
        @einsum σ[i,j] := C[i,j,k,l] * ϵ[k,l]
    return dot(ϵ, σ)
end

"""Function to create Elastic Tensor for Linear Elastic Isotropic Materials"""
function createElasticTensor(E::Float64, ν::Float64)
    λ = (ν * E) / ((1 + ν) * (1 - 2ν))
    μ = E / (2 * (1 + ν))
    C = λ * (one(SymmetricTensor{2, 3, Float64}) ⊗ one(SymmetricTensor{2, 3, Float64}))
    C += 2 * μ * one(SymmetricTensor{4, 3, Float64})
    return C
end

"""Function to create Plane Stress 2D Elastic Tensor for Linear Elastic Isotropic Materials"""
function createPlaneStressElasticTensor(E::Float64, ν::Float64)
    m_n(m::Int,n::Int) = 10*m+n
    c = E/(1.0-ν^2)
    C(i::Int,j::Int,k::Int,l::Int) = begin
        ij = m_n(i,j)
        kl = m_n(k,l)
        if i == j && k == l
            if ij == kl
                return c
            else
                return c*ν
            end
        elseif i != j && k != l && ij == kl
            return c*(1.0-ν)/2.0
        end
        return 0.0
    end
    return SymmetricTensor{4,2, Float64}(C)
end

