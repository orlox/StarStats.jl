using StarStats
using CairoMakie
using Random
using CSV
using DataFrames

include("codes_for_binary_data_analysis/find_peaks.jl")
include("codes_for_binary_data_analysis/identify_seveeral_overflows.jl")
include("codes_for_binary_data_analysis/merge_files.jl")
include("codes_for_binary_data_analysis/make_EEPs.jl")

param_names = [:q, :logP]

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

name_list = get_subfolders("../../condor_outputs/100models/results")

path_constructor = x-> "../../condor_outputs/100models/results/$(x[1])_$(x[2])/LOGS1/history.data"


inputs = Matrix{String}(undef,2, length(name_list))
input_values = Matrix{Float64}(undef, 2, length(name_list))

for i in 1:length(name_list)
    name = name_list[i]
    mass_ratio_str, period_str = split(name,"_")
    mass_ratio = parse(Float64, mass_ratio_str)
    period = parse(Float64, period_str)
    inputs[1,i] = mass_ratio_str
    inputs[2,i] = period_str
    input_values[1,i] = mass_ratio
    input_values[2,i] = period
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
    EEPs_names = Symbol.(EEPs_names)
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

    return EEPs_int, EEPs_names         #this line added
end

##

model_set = StellarModelSet(inputs, [:q, :logP], path_constructor, binary_dataframe_loader, compute_distance_and_EEPs_binaries!);

##
x_min = 0.0
x_max = maximum([maximum(model.df.x[.!isnan.(model.df.x)]) for model in model_set.models])

#x = ["0000.87500","0000.47500"]
x = ["0000.86250","0000.46250"]

p = "../../condor_outputs/100models/0000.86250_0000.46250/LOGS1/history.data"
donor = DataFrame(CSV.File(p, delim=" ", ignorerepeated=true, skipto=7, header=6))
logTdata = donor.log_Teff
logLdata = donor.log_L

#q = 0.875
#logP = 0.4875
q = 0.86250
logP = 0.46250
coords, indeces, simplex_id = StarStats.interpolation_info([q,logP],model_set.simplex_interpolant)
marker = model_set.check_possibility_of_interpolation[simplex_id]
eeps_number =length(model_set.models[indeces[1]].EEPs_type)
xvals = LinRange(x_min, eeps_number-1, 1000)

##
if marker == 1
    logTeff = interpolate_grid_quantity.(Ref(model_set),Ref([q, logP]),:logTeff, xvals)
    logL = interpolate_grid_quantity.(Ref(model_set),Ref([q, logP]),:logL, xvals)

    f = Figure()
    gl = GridLayout(f[1, 1])
    ax = Axis(gl[1, 1],xlabel = "log T_eff", ylabel = "log L")
    lines!(ax, logTeff, logL, linewidth=5, alpha=0.5, color = "blue", label="interpolated")
    lines!(ax, logTdata,logLdata,linewidth=3,color = "orange", label="original")
    legend = Legend(gl[1, 2], ax)
    ax.xreversed[] = true
    #save("HR_test.png", f)
    f
else
    println("ERROR: number of EEPs in intorpolation curves are different")
end


##
#=
function visualize_simplex_interpolant(ax, si::StarStats.SimplexInterpolant{N,P,LU,E,V}) where {N,P,LU,E,V}
    if size(si.points)[1] != 2
        return
    end
    for simplex in si.simplexes
        lines!(ax, [simplex.points[1,1], simplex.points[1,2], simplex.points[1,3], simplex.points[1,1]],
                   [simplex.points[2,1], simplex.points[2,2], simplex.points[2,3], simplex.points[2,1]],
                   linewidth=3)
    end
end

fig = Figure()
ax = Axis(fig[1,1], xlabel=L"\log_{10}\left(M_\mathrm{i}/M_\odot\right)", ylabel=L"v_\mathrm{rot,i}\;\mathrm{km\;s^{-1}}")
visualize_simplex_interpolant(ax, model_set.simplex_interpolant)
save("testing_simpl.png", fig)
fig
=#