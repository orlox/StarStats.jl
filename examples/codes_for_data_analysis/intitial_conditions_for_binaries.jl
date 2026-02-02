using StarStats
using CairoMakie
using Random
using CSV
using DataFrames
using CairoMakie
import CairoMakie: scatter!

include("find_peaks.jl")
include("identify_seveeral_overflows.jl")
include("merge_files.jl")
include("make_EEPs.jl")
include("verify_result.jl")
include("Statistic.jl")

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
metric = StarStats.EuclideanMetric()
##
model_set = StellarModelSet(inputs, [:q, :logP], path_constructor, binary_dataframe_loader, compute_distance_and_EEPs_binaries!, metric);

##
suggest_ref_failed = StarStats.refine_model_set_failed_interpolation(model_set)
number_of_params = length(model_set.input_names)
#verify_simplexes(model_set)
##
df,evol_colors = array_of_colors(model_set)
##
suggest_extra_ref = suggested_refinment_after_IQR(model_set)
##

function visualize_simplex_interpolant(ax, model_set,  si::StarStats.SimplexInterpolant{N,P,LU,E,V}, suggest_ref_failed,suggest_extra_ref,df ) where {N,P,LU,E,V}
    if size(si.points)[1] != 2
        return
    end
    for i in eachindex(model_set.simplex_interpolant.simplexes)
        simplex = model_set.simplex_interpolant.simplexes[i]
        
        if model_set.can_interpolate_simplex[i] == 1
            lines!(ax, [simplex.points[1,1], simplex.points[1,2], simplex.points[1,3], simplex.points[1,1]],
                [simplex.points[2,1], simplex.points[2,2], simplex.points[2,3], simplex.points[2,1]],
                linewidth=3, color = "green")
        else
            lines!(ax, [simplex.points[1,1], simplex.points[1,2], simplex.points[1,3], simplex.points[1,1]],
                [simplex.points[2,1], simplex.points[2,2], simplex.points[2,3], simplex.points[2,1]],
                linewidth=3, color = "black")
        end
    end

    #This is for bad simplexes
    name_1 = model_set.input_names[1]
    name_2 = model_set.input_names[2]
    for i in 1:length(suggest_ref_failed[name_1])
        scatter!(ax,suggest_ref_failed[name_1][i], suggest_ref_failed[name_2][i], color = "blue", markersize=7)
    end
    #This is for extra refinment
   
    for m in 1:length(suggest_extra_ref[name_1])
        scatter!(ax,suggest_extra_ref[name_1][m], suggest_extra_ref[name_2][m], color = "pink", markersize=7)
    end

    for i in 1:length(df.simplex_id)
        for k in 1:length(df.submodel_id[1])
            coord1 = parse(Float64, model_set.models[df.submodel_id[i][k]].input_params[1])
            coord2 = parse(Float64, model_set.models[df.submodel_id[i][k]].input_params[2])
            color_model = df.color_of_model[i][k]
            scatter!(coord1,coord2,color = color_model)
        end
    end
    
end

##
fig = Figure()
ax = Axis(fig[1,1], xlabel="q", ylabel="logP")
visualize_simplex_interpolant(ax, model_set, model_set.simplex_interpolant, suggest_ref_failed,suggest_extra_ref,df )
labels = ["different_evol","different_evol","different_evol", "bad interpolation","extra_refinment"]
my_colors = evol_colors
##
append!(my_colors, ["pink"])
append!(my_colors, ["blue"])

leg = Legend(fig[1, 2],
    [PolyElement(color=color) for color in my_colors],
    labels,
    valign=:center,  # вертикальное выравнивание по центру
    patchsize=(30, 20),
    labelsize=12
)
#save("testing_simpl.png", fig)
fig

##
#=
#HERE IS JUST A BARPLOT TO HAVE SOME SENCE OF DATASET

matrix_total_for_statistic, matrix_total_for_statistic_point_to_point = statistics_for_dataset(model_set)
upper_outlier, upper_outlier_point_to_point  = precompute_data_for_extra_refinment(model_set)
fig = Figure()
ax = Axis(fig[1, 1], title="Hist of differences", xlabel="EEPs", ylabel="differences")    
barplot!(ax, 2, maximum(matrix_total_for_statistic), color=:blue, strokecolor=:black, strokewidth=1, alpha = 0.5)
barplot!(ax, 4, maximum(matrix_total_for_statistic_point_to_point), color=:red, strokecolor=:black, strokewidth=1, alpha = 0.5)
scatter!(ax, 2, median(matrix_total_for_statistic), color=:red, strokecolor=:black, strokewidth=1, alpha = 0.5)
scatter!(ax, 4, median(matrix_total_for_statistic_point_to_point), color=:red, strokecolor=:black, strokewidth=1, alpha = 0.5)
hlines!(ax, upper_outlier, color=:blue,linewidth=1)
hlines!(ax, upper_outlier_point_to_point, color=:red,linewidth=1)


my_colors  = [ "blue", "red"]
labels = ["difference of length","distance EEPs"]
leg = Legend(fig[1, 2],
    [PolyElement(color=color) for color in my_colors],
    labels,
    valign=:center,  # вертикальное выравнивание по центру
    patchsize=(30, 20),
    labelsize=12
)

fig
save("2methods.png", fig) 
##
#########JUST A FEW TESTS ON EXTRA refinment
upper_outlier, upper_outlier_point_to_point  = precompute_data_for_extra_refinment(model_set)

array_of_models_for_ref_1 = []
array_of_models_for_ref_2 = []
array_of_differences_1 = []
array_of_differences_2 = []
for i in eachindex(model_set.simplex_interpolant.simplexes)
    if model_set.can_interpolate_simplex[i] == 1
        simplex = model_set.simplex_interpolant.simplexes[i]
        models = model_set.models
        df = StarStats.length_between_EEPs(simplex, models, model_set.metric)
        df_new = filter_data(df,upper_outlier,upper_outlier_point_to_point)
        for m in 1:length(df_new.extra_refinment_flag)
            if df_new.extra_refinment_flag[m] == 1
                model_1 = df_new.model_pair[m][1]
                model_2 = df_new.model_pair[m][2]
                append!(array_of_models_for_ref_1 , model_1)
                append!(array_of_models_for_ref_2 , model_2)
                append!(array_of_differences_1, df_new.difference[m])
                append!(array_of_differences_2, df_new.distance_between_eeps[m])

            end
        end
    end
end

index_max_1 = argmax(array_of_differences_1)
index_max_2 = argmax(array_of_differences_2)
println(index_max_1)
println(index_max_2)

logT1 = model_set.models[array_of_models_for_ref_1[5]].df.log_Teff
logL1 = model_set.models[array_of_models_for_ref_1[5]].df.logL
logT2 = model_set.models[array_of_models_for_ref_2[5]].df.log_Teff
logL2 = model_set.models[array_of_models_for_ref_2[5]].df.logL

fig = Figure()
ax = Axis(fig[1,1], xlabel="logT", ylabel="logL")

lines!(logT1, logL1)
lines!(logT2, logL2)

ax.xreversed = true
fig

##
#EXTRA TEST WITH BROTT MODELS

# code below can be used after the data has been downloaded

model_set_brott = StarStats.load_brott_data("brott_models", metric, :lmc);
##

number_of_params_brott, id_s_brott, matrix_of_params_brott = StarStats.suggest_new_simulations(model_set_brott, "output.txt")

#If you want ot verify the bad simplexes
verify_simplexes(model_set_brott)

df_brott,evol_colors_brott = array_of_colors(model_set_brott)

matrix_of_params_for_ref_brott, id_s_for_ref_brott  = suggest_refinment(model_set_brott)

fig = Figure()
ax = Axis(fig[1,1], xlabel="q", ylabel="logP")
visualize_simplex_interpolant(ax, model_set_brott, model_set_brott.simplex_interpolant, matrix_of_params_brott, matrix_of_params_for_ref_brott)
labels = ["different_evol","different_evol","different_evol", "extra refinment", "bad interpolation"]
my_colors = evol_colors_brott

append!(my_colors, ["pink"])
append!(my_colors, ["blue"])

#=
leg = Legend(fig[1, 2],
    [PolyElement(color=color) for color in my_colors],
    labels,
    valign=:center,  # вертикальное выравнивание по центру
    patchsize=(30, 20),
    labelsize=12
)
save("testing_simpl.png", fig)
=#
fig
=#