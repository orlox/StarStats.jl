using StarStats
using CairoMakie
using Random
using CSV
using DataFrames

include("../../../condor_outputs/find_peaks.jl")
include("../../../condor_outputs/identify_seveeral_overflows.jl")
include("../../../condor_outputs/merge_files.jl")
include("../../../condor_outputs/make_EEPs.jl")

param_names = [:q, :logP, :logM1]

function get_subfolders(path::String)
    items = readdir(path)
    subfolders = String[]
    for item in items
        item_path = joinpath(path, item)
        if isdir(item_path)
            push!(subfolders, item)
        end
    end
    return subfolders
end


name_list3d = get_subfolders("../../condor_outputs/3d_data/results")

##

path_consructor_3d = x-> "../../condor_outputs/3d_data/results/$(x[1])_$(x[2])_$(x[3])/LOGS1/history.data"


inputs = Matrix{String}(undef,3, length(name_list3d))
input_values = Matrix{Float64}(undef, 3, length(name_list3d))

for i in 1:length(name_list3d)
    name = name_list3d[i]
    mass_ratio_str, period_str, mass_pr_str = split(name,"_")
    mass_ratio = parse(Float64, mass_ratio_str)
    period = parse(Float64, period_str)
    mass =  parse(Float64, mass_pr_str)
    inputs[1,i] = mass_ratio_str
    inputs[2,i] = period_str
    inputs[3,i] =  mass_pr_str
    input_values[1,i] = mass_ratio
    input_values[2,i] = period
    input_values[3,i] = mass
end

function binary_dataframe_loader(path)

    h_donor =  DataFrame(CSV.File(path, delim=" ", ignorerepeated=true, skipto=7, header=6))
 
    for name in names(h_donor)
        if name == "log_Teff"
       
            h_donor[!, :logTeff] = copy(h_donor.log_Teff)  # создаем колонку logTeff
        end
        if name == "log_L"
            h_donor[!, :logL] = copy(h_donor.log_L)  # создаем колонку logL
        end
    end

    return h_donor

end

function compute_distance_and_EEPs_binaries!(df::DataFrame)

    ########## loading EEPS
    h_donor = df
    df1 = DataFrame(
        time =  h_donor.star_age ,
        signal =h_donor.thermal_eq_difference
    )
    smoothed_data = smooth_data(df1.signal, 20)
    index_max, index_min= find_peacs_indexes(df1.time, log10.(smoothed_data))
    left, right = find_edges_of_the_peak(index_max, index_min, log10.(smoothed_data), df1.time)
    df_overflows = identify_overflow(h_donor)

    EEPs, EEPs_names = construct_EEPs_vector(left, right, df_overflows, h_donor)
    EEPs_int = convert(Vector{Int}, EEPs)
    ########## Calculationg distance
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
        start_EEP = EEPs_int[j]
        end_EEP = EEPs_int[j+1]
        
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

    return EEPs_int
end

function EEPs_symbols!(df)
    h_donor = df
    df1 = DataFrame(
        time =  h_donor.star_age ,
        signal =h_donor.thermal_eq_difference
    )
    smoothed_data = smooth_data(df1.signal, 20)
    index_max, index_min= find_peacs_indexes(df1.time, log10.(smoothed_data))
    left, right = find_edges_of_the_peak(index_max, index_min, log10.(smoothed_data), df1.time)
    df_overflows = identify_overflow(h_donor)
    EEPs, EEPs_names = construct_EEPs_vector(left, right, df_overflows, h_donor)
    EEPs_names = Symbol.(EEPs_names)
    return EEPs_names
end
##

model_set = StellarModelSet(inputs, [:q, :logP, :logM1], path_consructor_3d, binary_dataframe_loader, compute_distance_and_EEPs_binaries!,EEPs_symbols!);
