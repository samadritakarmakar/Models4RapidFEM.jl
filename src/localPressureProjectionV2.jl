"""This function calculates the second stabilization matrix of Pressure projection method for the given element.
To be called locally when calculating matrix V 
(see: "The Finite Element Method: Its Basis and Fundamentals by Olek C Zienkiewicz, Robert L Taylor, JZ Zhu")
WARNING: This function has only been tested for triangles and tets. It may not work for other elements.

    V2::Matrix{Float64} = getPpp_V2local(element, shapeFunction, coordArray)

A per element constant (generally, α/μ, where α = 1 or 2 and μ is shear modulus) needs to be multiplied with V2 to get the final V2 matrix.
First proposed by Dohrmann and Bochev in 2004. 
"""

function getPpp_V2local(element::AbstractElement, shapeFunction::Array{ShapeFunction}, coordArray::Array{Float64,2})
    noOfIpPoints::Int64 = length(shapeFunction)
    noOfNodes = getNoOfElementNodes(shapeFunction)
    h = RapidFEM.getPoly(element, coordArray, element.order-1)
    H = zeros(size(h,2), size(h,2))
    G = zeros(size(h,2), noOfNodes)
    vol = 0.0
    for ipNo ∈ 1:noOfIpPoints
        ϕ = get_ϕ(shapeFunction, ipNo)
        ∂x_∂ξ = get_∂x_∂ξ(coordArray, shapeFunction, ipNo)
        dΩ = get_dΩ(element, ∂x_∂ξ, shapeFunction, ipNo)
        vol += dΩ
        x = getInterpolated_x(coordArray, ϕ)
        h = RapidFEM.getPoly(element, x, element.order-1)
        H_temp = h'*h
        if H_temp isa Float64
            H[1] += H_temp*dΩ
        else
            H += H_temp*dΩ
        end
        #for a ∈ 1:size(h,2)
        for a ∈ eachindex(h)
            for b ∈ 1:noOfNodes
                G[a, b] += h[a]*ϕ[b]*dΩ
            end
        end
    end
    return G'*(H\G)
end