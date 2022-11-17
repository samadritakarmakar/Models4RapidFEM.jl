"""Stiffness Matrix for Poission Problem."""
function ∇v_∇u!(K::Array{Float64, 2}, passedData::T, problemDim::Int64, element::AbstractElement, elementNo::Int64, shapeFunction::Array{ShapeFunction}, coordArray::Array{Float64, 2}; kwargs4function...) where T
    noOfIpPoints = getNoOfElementIpPoints(shapeFunction)
    noOfNodes = getNoOfElementNodes(shapeFunction)
    for ipNo::Int64 = 1:noOfIpPoints
        ∂X_∂ξ = get_∂x_∂ξ(coordArray, shapeFunction, ipNo)
        dΩ = get_dΩ(element, ∂X_∂ξ, shapeFunction, ipNo)
        ∂ϕ_∂X = get_∂ϕ_∂x(element, ∂X_∂ξ, shapeFunction, ipNo)
        for b = 1:noOfNodes
            for j = 1:size(∂ϕ_∂X, 2)
                for a = 1:noOfNodes
                    for i = 1:problemDim
                        K[problemDim * (a - 1) + i, problemDim * (b - 1) + i] += ∂ϕ_∂X[a, j] * ∂ϕ_∂X[b, j] * dΩ
                    end
                end
            end
        end
    end
end


