function ud_grade(input_map, output_nside; threshold=abs(1e-6UNSEEN), pess=false, power=0)
    map_out = udgrade(input_map, output_nside; threshold=threshold, pess=pess)
    map_out.pixels .*= ((input_map.resolution.nside//output_nside)^-power)
    return map_out
end 

function convert(::Type{T}, v::AbstractVector{W}) where {T<:Matrix, W<:HealpixMap}
    return convert(T, stack(v))
end

function Matrix(v::AbstractVector{W}) where {W<:HealpixMap}
    return convert(Matrix{Float64}, v)
end

function isapprox_in(element, list; kwargs...)
    # Check if the element is approximately equal to any element in the list.
    any(isapprox(element, el, kwargs...) for el in list)
end

function approx_geq(x, y; kwargs...)
    # Check if x is approximately greater than or equal to y.
    isapprox(x, y, kwargs...) | (x > y)
end

function approx_leq(x, y; kwargs...)
    # Check if x is approximately less than or equal to y.
    isapprox(x, y, kwargs...) | (x < y)
end