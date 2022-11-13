export get_EEPs, get_secondary_EEP, interpolate_grid_quantity

function get_EEPs(track)
    return get_EEPs_internal(track.c_h1, track.c_he4, size(track,1))
end

function get_EEPs_internal(central_h1, central_he4, nrows)
    initial_central_h1 = central_h1[1]
    ZAMS_EEP = 0
    IAMS_EEP = 0
    TAMS_EEP = 0
    RGBTip_EEP = 0
    ZACHeB_EEP = 0
    TACHeB_EEP = 0

    EEPs = [0,0,0,0,0,0]

    # ZAMS_EEP
    for i in 1:nrows
        if central_h1[i] < central_h1[1] - 0.0015
            ZAMS_EEP = i
            break
        end
    end
    if ZAMS_EEP == 0
        return EEPs
    end
    EEPs[1] = ZAMS_EEP
    # IAMS_EEP
    for i in ZAMS_EEP+1:nrows
        if central_h1[i] < 0.3
            IAMS_EEP = i
            break
        end
    end
    if IAMS_EEP == 0
        return EEPs
    end
    EEPs[2] = IAMS_EEP
    # TAMS_EEP
    for i in IAMS_EEP+1:nrows
        if central_h1[i] < 1e-12
            TAMS_EEP = i
            break
        end
    end
    if TAMS_EEP == 0
        return EEPs
    end
    EEPs[3] = TAMS_EEP
    # RGBTip_EEP
    for i in TAMS_EEP+1:nrows
        if central_he4[i] < central_he4[TAMS_EEP] - 0.01
            RGBTip_EEP = i
            break
        end
    end
    if RGBTip_EEP == 0
        return EEPs
    end
    EEPs[4] = RGBTip_EEP
    # ZACHeB_EEP
    for i in RGBTip_EEP+1:nrows
        if central_he4[i] < central_he4[RGBTip_EEP] - 0.03
            ZACHeB_EEP = i
            break
        end
    end
    if ZACHeB_EEP == 0
        return EEPs
    end
    EEPs[5] = ZACHeB_EEP
    # TACHeB_EEP
    for i in ZACHeB_EEP+1:nrows
        if central_he4[i] < 1e-4
            TACHeB_EEP = i
            break
        end
    end
    EEPs[6] = TACHeB_EEP
    return EEPs
end

function get_secondary_EEP(x, track, EEPs, interpolated_quantity)
    start_index = 0
    end_index = 0
    #rescale x to go from zero to one in the relevant phase
    if x>=0 && x<1 # ZAMS to IAMS
        start_index = EEPs[1]
        end_index = EEPs[2]
    elseif x>=1 && x<2 # IAMS to TAMS
        start_index = EEPs[2]
        end_index = EEPs[3]
        x = x-1
    elseif x>=2 && x<3 # TAMS to RGBTip
        start_index = EEPs[3]
        end_index = EEPs[4]
        x = x-2
    elseif x>=3 && x<4 # RGBTip to ZACHeB
        start_index = EEPs[4]
        end_index = EEPs[5]
        x = x-3
    elseif x>=4 && x<=5 # ZACHeB to TACHeB
        start_index = EEPs[5]
        end_index = EEPs[6]
        x = x-4
    end

    if start_index == 0 || end_index == 0
        return NaN
    end

    desired_distance = get_desired_distance(start_index, end_index, track.distance, x)

    (distances, vals) = get_values_around_desired_distance(start_index, end_index, 
                                        desired_distance, track.distance, getproperty(track, interpolated_quantity))

    return vals[1] + (desired_distance-distances[1])/(distances[2]-distances[1])*(vals[2]-vals[1]) # interpolate in x
end

function get_desired_distance(start_index, end_index, distance, x)
    return distance[start_index] + x*(distance[end_index] - distance[start_index])
end

function get_values_around_desired_distance(start_index, end_index, desired_distance, distance, required_quantity)
    lower_i = 0
    for i in start_index:(end_index-1)
        if desired_distance >= distance[i] && desired_distance <= distance[i+1]
            lower_i = i
            break
        end
    end
    return ((distance[lower_i], distance[lower_i+1]),(required_quantity[lower_i], required_quantity[lower_i+1]))
end

function interpolate_grid_quantity(grid, grid_parameters, interpolated_quantity, x)
    lower_i = zeros(Int, length(grid.input_names))
    grid_values = zeros(length(grid.input_names),2)
    # get values of parameters at edges
    for j in 1:length(grid.input_names)
        for i in 1:(length(grid.input_values[j])-1)
            if grid.input_values[j][i] <= grid_parameters[j] && grid.input_values[j][i+1] >= grid_parameters[j]
                grid_values[j,1] = grid.input_values[j][i]
                grid_values[j,2] = grid.input_values[j][i+1]
                lower_i[j] = i
                break
            end
        end
    end
    for i in lower_i
        if i == 0
            return NaN
        end
    end

    # get value of quantity to interpolate at vertices at a given x
    yvalues = zeros(typeof(x),[2 for i in 1:length(grid.input_names)]...)
    df_index = zeros(Int, length(grid.input_names))
    for index in Base.product([1:2 for i in 1:length(grid.input_names)]...)
        df_index .= lower_i.+index.-1
        if !isassigned(grid.dfs,df_index...)
            return NaN
        end
        yvalues[index...] = get_secondary_EEP(x,grid.dfs[df_index...], grid.EEPs[df_index...,:], interpolated_quantity)
    end

    for i in 1:length(grid.input_names)
        for index in Base.product([1:2 for i in 1:(length(grid.input_names)-i)]...)
            yval1 = yvalues[[1 for j in 1:(i-1)]...,1,index...]
            yval2 = yvalues[[1 for j in 1:(i-1)]...,2,index...]
            xval1 = grid_values[i,1]
            xval2 = grid_values[i,2]
            xtarget = grid_parameters[i]
            yvalues[[1 for j in 1:(i-1)]...,1,index...] = yval1 +(yval2-yval1)*(xtarget-xval1)/(xval2-xval1)
        end
    end

    return yvalues[[1 for i in 1:length(grid.input_names)]...]

end