export get_EEPs, get_secondary_EEP, interpolate_grid_quantity

function get_EEPs(track, Xc_TAMS)
    return get_EEPs_internal(track.c_h1, track.c_he4, size(track,1), Xc_TAMS)
end

function get_EEPs_internal(central_h1, central_he4, nrows, Xc_TAMS)
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
        if central_h1[i] < Xc_TAMS
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

function get_secondary_EEP(x::T, track, interpolated_quantity)::T where{T}
    get_secondary_EEP_internal(x, track.x, getproperty(track, interpolated_quantity))
end
function get_secondary_EEP_internal(x::T, trackx, interp_values)::T where{T}
    #dirty check to not overflow at edge
    lower_i = min(searchsorted(trackx,x).stop, length(trackx)-1)
    return interp_values[lower_i] + (x-trackx[lower_i])/(trackx[lower_i+1]-trackx[lower_i])*(interp_values[lower_i+1]-interp_values[lower_i]) # interpolate in x
end

function interpolate_grid_quantity(grid, grid_parameters, interpolated_quantity, x::T) where{T}
    return interpolate_grid_quantity_internal(grid.dfs, grid.input_values, grid_parameters, interpolated_quantity, x)
end

function df_index_element(i, lower_i, cartesian_index)
    return cartesian_index[i]+lower_i[i]-1
end
function interpolator_index1_element(j,i)
    if j<= i
        return 1
    end
    return 2
end
function interpolator_index2_element(j,i,index1)
    if j< i
        return 1
    elseif j==i
        return 2
    end
    return index1[j]
end

function interpolate_grid_quantity_internal(dfs, input_values, grid_parameters, interpolated_quantity, x::T) where{T}
    dimensions = ndims(dfs)
    lower_i = zeros(Int, dimensions)
    grid_values = zeros(dimensions,2)
    # get values of parameters at edges
    for j in 1:dimensions
        for i in 1:(length(input_values[j])-1)
            if input_values[j][i] <= grid_parameters[j] && input_values[j][i+1] >= grid_parameters[j]
                grid_values[j,1] = input_values[j][i]
                grid_values[j,2] = input_values[j][i+1]
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
    yvalues::Array{T,dimensions} = zeros(T,ntuple(_->2, dimensions)...)
    for index in CartesianIndices(yvalues)
        df_index = ntuple(i->df_index_element(i,lower_i,index),dimensions)
        try
            df = dfs[CartesianIndex(df_index)]
            yvalues[index] = get_secondary_EEP(x, df, interpolated_quantity)
        catch
            return NaN
        end
    end

    for i in 1:dimensions
        xval1 = grid_values[i,1]
        xval2 = grid_values[i,2]
        xtarget = grid_parameters[i]
        for index1 in CartesianIndices(ntuple(j->interpolator_index1_element(j,i),dimensions))
            index2 = CartesianIndex(ntuple(j->interpolator_index2_element(j,i, index1),dimensions))
            yval1 = yvalues[index1]
            yval2 = yvalues[index2]
            yvalues[index1] = yval1 +(yval2-yval1)*(xtarget-xval1)/(xval2-xval1)
        end
    end

    return yvalues[1]
end