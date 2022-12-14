"""Material stiffness matrix for plane strain problems, of updated lagrangian type."""
function Ξ΄D_π_ΞD_PlaneStrain!(π::Array{Float64, 2}, passedData::T, problemDim::Int64, element::AbstractElement, elementNo::Int64, shapeFunction::Array{ShapeFunction}, coordArray::Array{Float64, 2}; kwargs4function...) where T
    noOfIpPoints = getNoOfElementIpPoints(shapeFunction)
    noOfNodes = getNoOfElementNodes(shapeFunction)
    (hyperModel, modelParams, ipState_βuβ_βX, ipState_βuβ_βX, lastSol_Ξu) = passedData
    βΞu_βxβ_plnStrn = zeros(3, 3)
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
        fill!(βΞu_βxβ_plnStrn, 0.0)
        βΞu_βxβ_plnStrn[1:2, 1:2] = βΞu_βxβ
        Ξf = LargeDefs.getDeformationGradient(βΞu_βxβ_plnStrn)
        Fβ = LargeDefs.getDeformationGradient(βuβ_βX)
        βuβ_βX = βΞu_βxβ_plnStrn * Fβ + βuβ_βX
        updateIpStateDict!(βuβ_βX, ipState_βuβ_βX, elementNo, ipNo)
        π3d = LargeDefs.spatialTangentTensor(hyperModel, Ξf β Fβ, modelParams)
        π = @view(π3d[1:2, 1:2, 1:2, 1:2])
        π_Complete = Tensor{4, 2, Float64}(((i, j, k, l)->begin
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

"""Geometric stiffness matrix for plane strain problems, of updated lagrangian type."""
function Ξ΄u_Ο_Ξu_PlaneStrain!(K::Array{Float64, 2}, passedData::T, problemDim::Int64, element::AbstractElement, elementNo::Int64, shapeFunction::Array{ShapeFunction}, coordArray::Array{Float64, 2}; kwargs4function...) where T
    noOfIpPoints = getNoOfElementIpPoints(shapeFunction)
    noOfNodes = getNoOfElementNodes(shapeFunction)
    (hyperModel, modelParams, ipState_βuβ_βX, ipState_βuβ_βX, lastSol_Ξu) = passedData
    βΞu_βxβ_plnStrn = zeros(3, 3)
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
        fill!(βΞu_βxβ_plnStrn, 0.0)
        βΞu_βxβ_plnStrn[1:2, 1:2] = βΞu_βxβ
        Ξf = LargeDefs.getDeformationGradient(βΞu_βxβ_plnStrn)
        Fβ = LargeDefs.getDeformationGradient(βuβ_βX)
        Ο3d = LargeDefs.cauchyStress(hyperModel, Ξf β Fβ, modelParams)
        Ο = @view(Ο3d[1:2, 1:2])
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

"""Internal force vector for plane strain problems, of updated lagrangian type."""
function Ξ΄D_Ο_PlaneStrain!(f::Vector, passedData::T, problemDim::Int64, element::AbstractElement, elementNo::Int64, shapeFunction::Array{ShapeFunction}, coordArray::Array{Float64, 2}; kwargs4function...) where T
    noOfIpPoints = getNoOfElementIpPoints(shapeFunction)
    noOfNodes = getNoOfElementNodes(shapeFunction)
    (hyperModel, modelParams, ipState_βuβ_βX, ipState_βuβ_βX, lastSol_Ξu) = passedData
    βΞu_βxβ_plnStrn = zeros(3, 3)
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
        fill!(βΞu_βxβ_plnStrn, 0.0)
        βΞu_βxβ_plnStrn[1:2, 1:2] = βΞu_βxβ
        Ξf = LargeDefs.getDeformationGradient(βΞu_βxβ_plnStrn)
        Fβ = LargeDefs.getDeformationGradient(βuβ_βX)
        Ο3d = LargeDefs.cauchyStress(hyperModel, Ξf β Fβ, modelParams)
        Ο = @view(Ο3d[1:2, 1:2])
        for j = 1:size(βΟ_βx, 2)
            for a = 1:noOfNodes
                for i = 1:problemDim
                    f[problemDim * (a - 1) + i] += 0.5 * βΟ_βx[a, j] * (Ο[i, j] + Ο[i, j]) * dΟ
                end
            end
        end
    end
end

"""Cauchy stress at gauss point."""
function gaussian_Ο(passedData::T, solAtNodes::Array{Float64, 1}, problemDim::Int64, element::AbstractElement, elementNo::Int64, shapeFunction::Array{ShapeFunction}, coordArray::Array{Float64, 2}; kwargs4function...) where T
    begin
        noOfIpPoints = getNoOfElementIpPoints(shapeFunction)
        noOfNodes = getNoOfElementNodes(shapeFunction)
        (ipState, hyperModel, modelParams) = passedData
        Ο_g = zeros(noOfIpPoints, 3, 3)
        for ipNo::Int64 = 1:noOfIpPoints
            βuβ_βX = SaveFemIpState.getIpState(ipState, elementNo, ipNo)
            F = LargeDefs.getDeformationGradient(βuβ_βX)
            Ο = LargeDefs.cauchyStress(hyperModel, F, modelParams)
            for j = 1:problemDim
                for i = 1:problemDim
                    Ο_g[ipNo, i, j] = Ο[i, j]
                end
            end
        end
    end
    return Ο_g
end

"""Left Cauchy Tensor at gauss point."""
function gaussian_b(passedData::T, solAtNodes::Array{Float64, 1}, problemDim::Int64, element::AbstractElement, elementNo::Int64, shapeFunction::Array{ShapeFunction}, coordArray::Array{Float64, 2}; kwargs4function...) where T
    begin
        noOfIpPoints = getNoOfElementIpPoints(shapeFunction)
        noOfNodes = getNoOfElementNodes(shapeFunction)
        (ipState, hyperModel, modelParams) = passedData
        b_g = zeros(noOfIpPoints, 3, 3)
        for ipNo::Int64 = 1:noOfIpPoints
            βuβ_βX = SaveFemIpState.getIpState(ipState, elementNo, ipNo)
            F = LargeDefs.getDeformationGradient(βuβ_βX)
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
        βu_βX_g = zeros(noOfIpPoints, 3, 3)
        for ipNo::Int64 = 1:noOfIpPoints
            βuβ_βX = SaveFemIpState.getIpState(ipState, elementNo, ipNo)
            for j = 1:problemDim
                for i = 1:problemDim
                    βu_βX_g[ipNo, i, j] = βuβ_βX[i, j]
                end
            end
        end
    end
    return βu_βX_g
end

"""Deformation gradient at gauss point."""
function gaussian_DefGrad(passedData::T, solAtNodes::Array{Float64, 1}, problemDim::Int64, element::AbstractElement, elementNo::Int64, shapeFunction::Array{ShapeFunction}, coordArray::Array{Float64, 2}; kwargs4function...) where T
    begin
        noOfIpPoints = getNoOfElementIpPoints(shapeFunction)
        noOfNodes = getNoOfElementNodes(shapeFunction)
        (ipState, hyperModel, modelParams) = passedData
        F_g = zeros(noOfIpPoints, 3, 3)
        for ipNo::Int64 = 1:noOfIpPoints
            βuβ_βX = SaveFemIpState.getIpState(ipState, elementNo, ipNo)
            F = LargeDefs.getDeformationGradient(βuβ_βX)
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
        βX_βΞΎ = get_βx_βΞΎ(coordArray, shapeFunction, ipNo)
        dΞ© = get_dΞ©(element, βX_βΞΎ, shapeFunction, ipNo)
        βuβ_βX = SaveFemIpState.getIpState(ipState, elementNo, ipNo)
        F = LargeDefs.getDeformationGradient(βuβ_βX)
        E = LargeDefs.getGreenLagrangeStrain(F)
        Ξ¨_local = hyperModel.strainEnergyDensity(E, modelParams)
        Jacobian = LargeDefs.getJacobianDeformationGradient(F)
        dΟβ = dΞ©
        for i = 1:problemDim
            StrainE[i] += if i == 1
                    Ξ¨_local * ((1.0 / Jacobian) * dΟβ)
                else
                    0.0
                end
        end
    end
end
