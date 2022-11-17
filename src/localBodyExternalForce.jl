"""External Body Force Vector for hyper elastic problems of total lagrangian strain type."""
function v_S!(B::Vector, passedData::T, problemDim::Int64, element::AbstractElement, elementNo::Int64, shapeFunction::Array{ShapeFunction}, coordArray::Array{Float64, 2}; kwargs4function...) where T
    noOfIpPoints = getNoOfElementIpPoints(shapeFunction)
    noOfNodes = getNoOfElementNodes(shapeFunction)
    sourceFunc = passedData
    for ipNo::Int64 = 1:noOfIpPoints
        ∂X_∂ξ = get_∂x_∂ξ(coordArray, shapeFunction, ipNo)
        dΩ = get_dΩ(element, ∂X_∂ξ, shapeFunction, ipNo)
        ϕ = get_ϕ(shapeFunction, ipNo)
        X = getInterpolated_x(coordArray, ϕ)
        s = sourceFunc(X)
        for a = 1:noOfNodes
            for i = 1:problemDim
                B[problemDim * (a - 1) + i] += ϕ[a] * s[i] * dΩ
            end
        end
    end
end

"""External Force Vector for hyper elastic problems of total lagrangian strain type."""
function v_F!(B::Vector, passedData::T, problemDim::Int64, element::AbstractElement, elementNo::Int64, shapeFunction::Array{ShapeFunction}, coordArray::Array{Float64, 2}; kwargs4function...) where T
    noOfIpPoints = getNoOfElementIpPoints(shapeFunction)
    noOfNodes = getNoOfElementNodes(shapeFunction)
    neumFunc = passedData
    for ipNo::Int64 = 1:noOfIpPoints
        ∂X_∂ξ = get_∂x_∂ξ(coordArray, shapeFunction, ipNo)
        dS = get_dS(element, ∂X_∂ξ, shapeFunction, ipNo)
        ϕ = get_ϕ(shapeFunction, ipNo)
        X = getInterpolated_x(coordArray, ϕ)
        f = neumFunc(X)
        for a = 1:noOfNodes
            for i = 1:problemDim
                B[problemDim * (a - 1) + i] += ϕ[a] * f[i] * dS
            end
        end
    end
end