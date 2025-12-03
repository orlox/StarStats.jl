export get_EEPs, get_secondary_EEP, compute_distance_and_EEPs!

using DataFrames

function get_EEPs(track, Xc_TAMS)
    return get_EEPs_internal(track.c_h1, track.c_he4, size(track,1), Xc_TAMS)
end

function names_of_EEPs_single_star!(df)
    names = []
    push!(names,"ZAMS_EEP")
    push!(names,"IAMS_EEP")
    push!(names, "TAMS_EEP")
    push!(names, "RGBTip_EEP")
    push!(names,"ZACHeB_EEP" )
    push!(names, "TACHeB_EEP ")
    #println()
    names_new = Symbol.(names)
    return names_new
end
function get_EEPs_internal(central_h1, central_he4, nrows, Xc_TAMS)
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

"""
   compute_distance_and_EEPs

Function 

# Fields: 
-   
"""
function compute_distance_and_EEPs!(df::DataFrame; Xc_TAMS=1e-10)

    EEPs = get_EEPs(df, Xc_TAMS)

    distance = zeros(size(df,1))
    delta_log_Teff = 0
    delta_log_L = 0

    for i in 2:size(df,1)
        delta_log_Teff = df.logTeff[i]-df.logTeff[i-1]
        delta_log_L = df.logL[i]-df.logL[i-1]

        distance[i] = distance[i-1] + sqrt(delta_log_Teff^2 + delta_log_L^2)
    end
    df[!,:distance] = distance

    x = zeros(size(df,1))
    for j in 1:(length(EEPs)-1)
        start_EEP = EEPs[j]
        end_EEP = EEPs[j+1]
        if start_EEP==0 || end_EEP==0
            break
        end
        start_distance = df.distance[start_EEP]
        end_distance = df.distance[end_EEP]
        for i in start_EEP:end_EEP
            x[i] = (df.distance[i]-start_distance)/(end_distance-start_distance) + (j-1)
        end
    end
    for i in size(df,1):-1:1
        if x[i] > 0
            x[i:end] .= x[i]
            break
        end
    end
    df[!,:x] = x

    dtdx = zeros(size(df,1))
    for i in 2:size(df,1)
        if (df.x[i]-df.x[i-1])<=0
            continue
        end
        dtdx[i] = (df.star_age[i]-df.star_age[i-1])/(df.x[i]-df.x[i-1])
    end
    df[!,:dtdx] = dtdx

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

function check_coords_of_simplex(simplex, models)
    barycenter = simplex.barycenter
    coords = barycentric_coords(barycenter, simplex)
    number_of_eeps = []
    all_good = 1
    for i in 1:length(coords)          
        if coords[i] == 0
            continue
        else
            append!(number_of_eeps,length(models[simplex.point_indeces[i]].EEPs))
        end
    end

    if Statistics.mean(number_of_eeps) != number_of_eeps[1]
        all_good = 0
    else
        all_good = 1
    end

    return all_good
end

function check_names_of_EEPs(simplex, models)
    barycenter = simplex.barycenter
    coords = barycentric_coords(barycenter, simplex)
    all_good = 1
    if coords[1] != 0
        model1 = models[simplex.point_indeces[1]].EEPs_type
    end
    
    for i in 2:length(coords)
        if coords[i] != 0
            model_test = models[simplex.point_indeces[i]].EEPs_type
            for k in 1:length(model_test)
                if model1[k] == model_test[k]
                    continue
                else
                    println("EEPs names are different, intorpolations is not recommended")
                    println("model1 $(model1[k]) and model $(i) $(model_test[k])")
                    all_good = 0
                    break
                end
            end
        end
    end
    return all_good
end