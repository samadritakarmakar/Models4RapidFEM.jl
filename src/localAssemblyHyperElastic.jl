"""Material Stiffness tensor for hyper elastic problems of total lagrangian strain type."""
function Ξ΄E_Cα΅_ΞE!(π::Array{Float64, 2}, passedData::T, problemDim::Int64, element::AbstractElement, elementNo::Int64, shapeFunction::Array{ShapeFunction}, coordArray::Array{Float64, 2}; kwargs4function...) where T
    noOfIpPoints = getNoOfElementIpPoints(shapeFunction)
    noOfNodes = getNoOfElementNodes(shapeFunction)
    (model, modelParams, lastSol_u) = passedData
    βΟ_βX_F_β = zeros(3, noOfNodes, 3, 3)
    βΟ_βX_F = zeros(3, noOfNodes, 3, 3)
    u_Nodes = getSolAtElement(lastSol_u, element, problemDim)
    for ipNo::Int64 = 1:noOfIpPoints
        βX_βΞΎ = get_βx_βΞΎ(coordArray, shapeFunction, ipNo)
        dΞ© = get_dΞ©(element, βX_βΞΎ, shapeFunction, ipNo)
        βΟ_βX = get_βΟ_βx(element, βX_βΞΎ, shapeFunction, ipNo)
        βu_βX = get_βu_βx(u_Nodes, βΟ_βX, Int64(length(u_Nodes) / size(βΟ_βX, 1)))
        F = LargeDefs.getDeformationGradient(βu_βX)
        E = LargeDefs.getGreenLagrangeStrain(F)
        S = model.secondPiolaStress(E, modelParams)
        β = model.materialTangentTensor(E, modelParams) * dΞ©
        @einsum βΟ_βX_F[i, a, I, J] = βΟ_βX[a, I] * F[i, J]
        @einsum βΟ_βX_F_β[i, a, K, L] = βΟ_βX_F[i, a, I, J] * β[I, J, K, L]
        for K = 1:size(βΟ_βX_F_β, 3)
            for L = 1:size(βΟ_βX_F_β, 4)
                for b = 1:noOfNodes
                    for j = 1:problemDim
                        for a = 1:noOfNodes
                            for i = 1:problemDim
                                π[problemDim * (a - 1) + i, problemDim * (b - 1) + j] += 0.5 * βΟ_βX_F_β[i, a, K, L] * (βΟ_βX_F[j, b, L, K] + βΟ_βX_F[j, b, K, L])
                            end
                        end
                    end
                end
            end
        end
    end
end

"""Geometric Stiffness tensor for hyper elastic problems of total lagrangian strain type."""
function Ξ΄E_S_ΞE!(π::Array{Float64, 2}, passedData::T, problemDim::Int64, element::AbstractElement, elementNo::Int64, shapeFunction::Array{ShapeFunction}, coordArray::Array{Float64, 2}; kwargs4function...) where T
    noOfIpPoints = getNoOfElementIpPoints(shapeFunction)
    noOfNodes = getNoOfElementNodes(shapeFunction)
    (model, modelParams, lastSol_u) = passedData
    βΟ_βX_S = zeros(noOfNodes, 3, 3)
    u_Nodes = getSolAtElement(lastSol_u, element, problemDim)
    for ipNo::Int64 = 1:noOfIpPoints
        βX_βΞΎ = get_βx_βΞΎ(coordArray, shapeFunction, ipNo)
        dΞ© = get_dΞ©(element, βX_βΞΎ, shapeFunction, ipNo)
        βΟ_βX = get_βΟ_βx(element, βX_βΞΎ, shapeFunction, ipNo)
        βΟ_βX = get_βΟ_βx(element, βX_βΞΎ, shapeFunction, ipNo)
        βu_βX = get_βu_βx(u_Nodes, βΟ_βX, Int64(length(u_Nodes) / size(βΟ_βX, 1)))
        F = LargeDefs.getDeformationGradient(βu_βX)
        E = LargeDefs.getGreenLagrangeStrain(F)
        S = model.secondPiolaStress(E, modelParams) * dΞ©
        @einsum βΟ_βX_S[a, I, J] = βΟ_βX[a, I] * S[I, J]
        for J = 1:size(βΟ_βX_S, 3)
            for b = 1:noOfNodes
                for I = 1:size(βΟ_βX_S, 2)
                    for a = 1:noOfNodes
                        for i = 1:problemDim
                            π[problemDim * (a - 1) + i, problemDim * (b - 1) + i] += βΟ_βX_S[a, I, J] * βΟ_βX[b, J]
                        end
                    end
                end
            end
        end
    end
end

"""Internal Force Vector for hyper elastic problems of total lagrangian strain type."""
function Ξ΄E_S!(B::Vector, passedData::T, problemDim::Int64, element::AbstractElement, elementNo::Int64, shapeFunction::Array{ShapeFunction}, coordArray::Array{Float64, 2}; kwargs4function...) where T
    noOfIpPoints = getNoOfElementIpPoints(shapeFunction)
    noOfNodes = getNoOfElementNodes(shapeFunction)
    (model, modelParams, lastSol_u) = passedData
    F_S = zeros(3, 3)
    u_Nodes = getSolAtElement(lastSol_u, element, problemDim)
    for ipNo::Int64 = 1:noOfIpPoints
        βX_βΞΎ = get_βx_βΞΎ(coordArray, shapeFunction, ipNo)
        dΞ© = get_dΞ©(element, βX_βΞΎ, shapeFunction, ipNo)
        βΟ_βX = get_βΟ_βx(element, βX_βΞΎ, shapeFunction, ipNo)
        βu_βX = get_βu_βx(u_Nodes, βΟ_βX, Int64(length(u_Nodes) / size(βΟ_βX, 1)))
        F = LargeDefs.getDeformationGradient(βu_βX)
        E = LargeDefs.getGreenLagrangeStrain(F)
        S = model.secondPiolaStress(E, modelParams) * dΞ©
        @einsum F_S[i, I] = F[i, J] * S[I, J]
        for I = 1:size(βΟ_βX, 2)
            for a = 1:noOfNodes
                for i = 1:problemDim
                    B[problemDim * (a - 1) + i] += βΟ_βX[a, I] * F_S[i, I]
                end
            end
        end
    end
end