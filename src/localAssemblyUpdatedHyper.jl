"""Material Stiffness Matrix for hyper elastic problems of Updated Lagrangian type"""
function Ξ΄D_π_ΞD!(π::Array{Float64, 2}, passedData::T, problemDim::Int64, element::AbstractElement, elementNo::Int64, shapeFunction::Array{ShapeFunction}, coordArray::Array{Float64, 2}; kwargs4function...) where T
    noOfIpPoints = getNoOfElementIpPoints(shapeFunction)
    noOfNodes = getNoOfElementNodes(shapeFunction)
    (hyperModel, modelParams, ipState_βuβ_βX, ipState_βuβ_βX, lastSol_Ξu) = passedData
    βuβ_βX = zeros(3, 3)
    Ξu_Nodes = getSolAtElement(lastSol_Ξu, element, problemDim)
    currentCoordArray = getCurrentCoordArray(coordArray, Ξu_Nodes)
    for ipNo::Int64 = 1:noOfIpPoints
        βx_βΞΎ = get_βx_βΞΎ(currentCoordArray, shapeFunction, ipNo)
        βX_βΞΎ = get_βx_βΞΎ(coordArray, shapeFunction, ipNo)
        dΟ = get_dΞ©(element, βx_βΞΎ, shapeFunction, ipNo)
        βΟ_βX = get_βΟ_βx(element, βX_βΞΎ, shapeFunction, ipNo)
        βΟ_βx = get_βΟ_βx(element, βx_βΞΎ, shapeFunction, ipNo)
        βΞu_βxβ = get_βu_βx(Ξu_Nodes, βΟ_βX, Int64(length(Ξu_Nodes) / size(βΟ_βX, 1)))
        SaveFemIpState.getIpState!(βuβ_βX, ipState_βuβ_βX, elementNo, ipNo)
        Ξf = LargeDefs.getDeformationGradient(βΞu_βxβ)
        Fβ = LargeDefs.getDeformationGradient(βuβ_βX)
        βuβ_βX = βΞu_βxβ * Fβ + βuβ_βX
        updateIpStateDict!(βuβ_βX, ipState_βuβ_βX, elementNo, ipNo)
        π = LargeDefs.spatialTangentTensor(hyperModel, Ξf β Fβ, modelParams)
        π_Complete = Tensor{4, 3, Float64}(((i, j, k, l)->begin
                        π[i, j, k, l] + π[j, i, k, l] + π[i, j, l, k] + π[j, i, l, k]
                    end))
        for l = 1:size(π_Complete, 4)
            for j = 1:size(βΟ_βx, 2)
                for b = 1:noOfNodes
                    for k = 1:problemDim
                        for a = 1:noOfNodes
                            for i = 1:problemDim
                                π[problemDim * (a - 1) + i, problemDim * (b - 1) + k] += 0.25 * βΟ_βx[a, j] * π_Complete[i, j, k, l] * βΟ_βx[b, l] * dΟ
                            end
                        end
                    end
                end
            end
        end
    end
end

"""Geometric Stiffness Matrix for hyper elastic problems of Updated Lagrangian type"""
function Ξ΄u_Ο_Ξu!(K::Array{Float64, 2}, passedData::T, problemDim::Int64, element::AbstractElement, elementNo::Int64, shapeFunction::Array{ShapeFunction}, coordArray::Array{Float64, 2}; kwargs4function...) where T
    noOfIpPoints = getNoOfElementIpPoints(shapeFunction)
    noOfNodes = getNoOfElementNodes(shapeFunction)
    (hyperModel, modelParams, ipState_βuβ_βX, ipState_βuβ_βX, lastSol_Ξu) = passedData
    βuβ_βX = zeros(3, 3)
    Ξu_Nodes = getSolAtElement(lastSol_Ξu, element, problemDim)
    currentCoordArray = getCurrentCoordArray(coordArray, Ξu_Nodes)
    for ipNo::Int64 = 1:noOfIpPoints
        βx_βΞΎ = get_βx_βΞΎ(currentCoordArray, shapeFunction, ipNo)
        βX_βΞΎ = get_βx_βΞΎ(coordArray, shapeFunction, ipNo)
        dΟ = get_dΞ©(element, βx_βΞΎ, shapeFunction, ipNo)
        βΟ_βX = get_βΟ_βx(element, βX_βΞΎ, shapeFunction, ipNo)
        βΟ_βx = get_βΟ_βx(element, βx_βΞΎ, shapeFunction, ipNo)
        βΞu_βxβ = get_βu_βx(Ξu_Nodes, βΟ_βX, Int64(length(Ξu_Nodes) / size(βΟ_βX, 1)))
        SaveFemIpState.getIpState!(βuβ_βX, ipState_βuβ_βX, elementNo, ipNo)
        Ξf = LargeDefs.getDeformationGradient(βΞu_βxβ)
        Fβ = LargeDefs.getDeformationGradient(βuβ_βX)
        Ο = LargeDefs.cauchyStress(hyperModel, Ξf β Fβ, modelParams)
        for k = 1:size(βΟ_βx, 2)
            for b = 1:noOfNodes
                for l = 1:size(βΟ_βx, 2)
                    for a = 1:noOfNodes
                        for j = 1:problemDim
                            K[problemDim * (a - 1) + j, problemDim * (b - 1) + j] += βΟ_βx[a, l] * βΟ_βx[b, k] * Ο[l, k] * dΟ
                        end
                    end
                end
            end
        end
    end
end

"""Internal Force Vector for hyper elastic problems of Updated Lagrangian type"""
function Ξ΄D_Ο!(f::Vector, passedData::T, problemDim::Int64, element::AbstractElement, elementNo::Int64, shapeFunction::Array{ShapeFunction}, coordArray::Array{Float64, 2}; kwargs4function...) where T
    noOfIpPoints = getNoOfElementIpPoints(shapeFunction)
    noOfNodes = getNoOfElementNodes(shapeFunction)
    (hyperModel, modelParams, ipState_βuβ_βX, ipState_βuβ_βX, lastSol_Ξu) = passedData
    βuβ_βX = zeros(3, 3)
    Ξu_Nodes = getSolAtElement(lastSol_Ξu, element, problemDim)
    currentCoordArray = getCurrentCoordArray(coordArray, Ξu_Nodes)
    for ipNo::Int64 = 1:noOfIpPoints
        βx_βΞΎ = get_βx_βΞΎ(currentCoordArray, shapeFunction, ipNo)
        βX_βΞΎ = get_βx_βΞΎ(coordArray, shapeFunction, ipNo)
        dΟ = get_dΞ©(element, βx_βΞΎ, shapeFunction, ipNo)
        βΟ_βX = get_βΟ_βx(element, βX_βΞΎ, shapeFunction, ipNo)
        βΟ_βx = get_βΟ_βx(element, βx_βΞΎ, shapeFunction, ipNo)
        βΞu_βxβ = get_βu_βx(Ξu_Nodes, βΟ_βX, Int64(length(Ξu_Nodes) / size(βΟ_βX, 1)))
        SaveFemIpState.getIpState!(βuβ_βX, ipState_βuβ_βX, elementNo, ipNo)
        Ξf = LargeDefs.getDeformationGradient(βΞu_βxβ)
        Fβ = LargeDefs.getDeformationGradient(βuβ_βX)
        Ο = LargeDefs.cauchyStress(hyperModel, Ξf β Fβ, modelParams)
        for j = 1:size(βΟ_βx, 2)
            for a = 1:noOfNodes
                for i = 1:problemDim
                    f[problemDim * (a - 1) + i] += 0.5 * βΟ_βx[a, j] * (Ο[i, j] + Ο[i, j]) * dΟ
                end
            end
        end
    end
end
