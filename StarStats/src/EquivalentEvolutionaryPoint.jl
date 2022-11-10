export get_EEPs, get_secondary_EEP

function get_EEPs(track)
    initial_central_h1 = track.c_h1[1]
    ZAMS_EEP = 0
    IAMS_EEP = 0
    TAMS_EEP = 0
    RGBTip_EEP = 0
    ZACHeB_EEP = 0
    TACHeB_EEP = 0
    # ZAMS_EEP
    for (i, central_h1) in enumerate(track.c_h1)
        if central_h1 < initial_central_h1 - 0.0015
            ZAMS_EEP = i
            break
        end
    end
    # IAMS_EEP
    for (i, central_h1) in enumerate(track.c_h1[ZAMS_EEP+1:end])
        if central_h1 < 0.3
            IAMS_EEP = i + ZAMS_EEP
            break
        end
    end
    # TAMS_EEP
    for (i, central_h1) in enumerate(track.c_h1[IAMS_EEP+1:end])
        if central_h1 < 1e-12
            TAMS_EEP = i + IAMS_EEP
            break
        end
    end
    # RGBTip_EEP
    TAMS_center_he4 = track.c_he4[TAMS_EEP]
    for (i, central_he4) in enumerate(track.c_he4[TAMS_EEP+1:end])
        if central_he4 < TAMS_center_he4 - 0.01
            RGBTip_EEP = i + TAMS_EEP
            break
        end
    end
    # ZACHeB_EEP
    RGBTip_center_he4 = track.c_he4[RGBTip_EEP]
    for (i, central_he4) in enumerate(track.c_he4[RGBTip_EEP+1:end])
        if central_he4 < RGBTip_center_he4 - 0.03
            ZACHeB_EEP = i + RGBTip_EEP
            break
        end
    end
    # TACHeB_EEP
    for (i, central_he4) in enumerate(track.c_he4[ZACHeB_EEP+1:end])
        if central_he4 < 1e-4
            TACHeB_EEP = i + ZACHeB_EEP
            break
        end
    end
    return (ZAMS_EEP, IAMS_EEP, TAMS_EEP, RGBTip_EEP, ZACHeB_EEP, TACHeB_EEP)
end

function get_secondary_EEP(x, track, interpolated_quantity)
    EEPs = get_EEPs(track)
    start_index = 0
    end_index = 0
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
    total_distance = 0
    for i in start_index:(end_index-1)
        delta_log_Teff = track.Teff[i+1]-track.Teff[i]
        delta_log_L = track.logL[i+1]-track.logL[i]
        total_distance += sqrt(delta_log_Teff^2 + delta_log_L^2)
    end
    desired_distance = x*total_distance

    distance1 = 0
    distance2 = 0
    lower_i = 0
    for i in start_index:(end_index-1)
        delta_log_Teff = track.Teff[i+1]-track.Teff[i]
        delta_log_L = track.logL[i+1]-track.logL[i]
        distance1 = distance2
        distance2 += sqrt(delta_log_Teff^2 + delta_log_L^2)
        if desired_distance >= distance1 && desired_distance <= distance2
            lower_i = i
            break
        end
    end

    val1 = getproperty(track,interpolated_quantity)[lower_i]
    val2 = getproperty(track,interpolated_quantity)[lower_i+1]
    return val1 + (desired_distance-distance1)/(distance2-distance1)*(val2-val1) # interpolate in x
end