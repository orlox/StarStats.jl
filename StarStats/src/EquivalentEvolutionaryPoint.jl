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

function get_secondary_EEP(x, track, interpolated_quantity)
    get_secondary_EEP_internal(x, track.x, getproperty(track, interpolated_quantity))
end
function get_secondary_EEP_internal(x, trackx, interp_values)
    #dirty check to not overflow at edge
    lower_i = min(searchsorted(trackx,x).stop, length(trackx)-1)
    return interp_values[lower_i] + (x-trackx[lower_i])/(trackx[lower_i+1]-trackx[lower_i])*(interp_values[lower_i+1]-interp_values[lower_i]) # interpolate in x
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
    for index in collect(Base.product([1:2 for i in 1:length(grid.input_names)]...))
        df_index .= lower_i.+index.-1
        if !isassigned(grid.dfs,df_index...)
            return NaN
        end
        yvalues[index...] = get_secondary_EEP(x,grid.dfs[df_index...], interpolated_quantity)
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