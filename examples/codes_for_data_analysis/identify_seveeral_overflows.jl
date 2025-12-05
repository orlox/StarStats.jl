using Pkg
using DataFrames
using CSV
using Plots
using LinearAlgebra
using Interpolations
using Base.Filesystem
using CairoMakie
using Statistics


function identify_overflow(h_donor)
    overflow = []
    arr_index = []
    arr_stage = []
    stage = 0
    a = 0
    for i in 1:length(h_donor.log_Teff)-1
        
        #if (h_donor.rl_relative_overflow_1[i] > -1e-4).&& (h_donor.rl_relative_overflow_1[i+1] > -1e-4)
        if (h_donor.lg_mtransfer_rate[i]>-10) .&& (h_donor.lg_mtransfer_rate[i+1]>-10)
              #overflow ends before end of track
            append!(overflow,i)
            if (h_donor.center_h1[i] > 1e-12)
                stage = 0 #case A (MS)
            elseif (h_donor.center_h1[i] < 1e-12) .&& (h_donor.he_core_mass[i] > 0)
                stage = 1 #case B (NO MS but with He)
            else
                stage = 2 #case C (NO MS NO He)
            end
        elseif (h_donor.lg_mtransfer_rate[i] < -11)
            continue
        else
            if length(overflow) == 0
                continue
            #println(string(Integer(length(arr_index)/2+1))*" OVERFLOW ENDED")
            else
                append!(arr_index, overflow[1]) #index начала overflow
                append!(arr_index, overflow[end]) #index конец overflow
                append!(arr_stage, stage) #stage начлаа overflow
                append!(arr_stage, stage) #stage конец overflow
                overflow = []
            end
            a = a+1
        end

    end
    if (a == 0) #this is for outflow once started and never ended
        append!(arr_index, overflow[1]) #index начала overflow
        append!(arr_index, overflow[end]) #index конец overflow
        append!(arr_stage, stage) #stage начлаа overflow
        append!(arr_stage, stage) #stage конец overflow
    end

    df = DataFrame(
        index_overflow = arr_index,
        stage = arr_stage)
    return df
end

function find_number_overflows(df)
    casea = 0
    caseb =0
    casec = 0
    for i in 1:Integer(length(df.stage)/2)
        if (df.stage[2*i] == 0)
            casea = casea+1
        elseif (df.stage[2*i] == 1)
            caseb = caseb+1
        else
            casec = casec+1
        end
    end
    return casea,caseb,casec
end 

