function ud_grade(input_map, output_nside; threshold=abs(1e-6UNSEEN), pess=false, power=0)
    map_out = udgrade(input_map, output_nside; threshold=threshold, pess=pess)
    map_out.pixels .*= ((input_map.resolution.nside//output_nside)^-power)
    return map_out
end 