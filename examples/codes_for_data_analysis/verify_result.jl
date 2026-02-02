
using StarStats
using CairoMakie
using Random
using CSV
using DataFrames
using CairoMakie
import CairoMakie: scatter!


function verify_simplexes(model_set::StellarModelSet)

    println("THIS FUNCTION IS A VERIFICATION OF BAD SIMPLEXES, BELOW IT WILL STATE ALL INFROMATION ABOUT SIMPLEXES")
    for i in eachindex(model_set.simplex_interpolant.simplexes)
        simplex = model_set.simplex_interpolant.simplexes[i]
        if !model_set.can_interpolate_simplex[i]
            println("In simpllex number $(simplex.id) there are models used: ")
            for i in 1:length(simplex.point_indeces)
                i_th_model = simplex.point_indeces[i]
                i_th_params = [simplex.points[1,i],simplex.points[2,i]]
                println("$(i_th_model) with parameters $(i_th_params)")
                i_th_EEPS = model_set.models[i_th_model].EEP_names
                number_of_EEPS = length(model_set.models[i_th_model].EEP_names)
                println("And EEPS are $(i_th_EEPS), number of EEPS is $(number_of_EEPS)")
            end
            println("------------------------------------------------------")
        end
    

    end

end

#function used to make vector of all eeps used for coloring array
function make_huge_vector_of_eeps(model_set::StellarModelSet)

    huge_vector_of_eeps = []
    for i in eachindex(model_set.simplex_interpolant.simplexes)
        simplex = model_set.simplex_interpolant.simplexes[i]
        for k in 1:length(simplex.point_indeces)
            k_th_model = simplex.point_indeces[k]
            k_th_EEPS = model_set.models[k_th_model].EEP_names
            push!(huge_vector_of_eeps, k_th_EEPS)
        end
    end
    markers = unique(huge_vector_of_eeps) # identifies how many unique evolution behaviers there are
    
    return markers

end

#this fuction give diffrent colors for different evolutions
function array_of_colors(model_set::StellarModelSet)
    color_list = ["cyan","red","yellow","orange","brown","gray","white","magenta",
    "violet","indigo","gold","silver","maroon","olive","navy","teal","coral","lavender"]
    df = DataFrame(
    simplex_id = [],
    submodel_id = [],
    point_name = [],
    marker_of_model = [],
    color_of_model = [])

    arr_id = []
    arr_submodel = []
    arr_eeps = []
    arr_marker = []
    arr_color = []
    
    markers = make_huge_vector_of_eeps(model_set)

    for i in eachindex(model_set.simplex_interpolant.simplexes)
        simplex = model_set.simplex_interpolant.simplexes[i]
        append!(arr_id, simplex.id)
        for k in 1:length(simplex.point_indeces)
                k_th_model = simplex.point_indeces[k]
                k_th_EEPS = model_set.models[k_th_model].EEP_names
                append!(arr_submodel,k_th_model)
                push!(arr_eeps,k_th_EEPS )
                for m in 1:length(markers)
                    if k_th_EEPS == markers[m]
                        append!(arr_marker, m) # this marker of model set 0 from the start will then be changed
                        push!(arr_color,color_list[m] )
                    end
                end
        end
        push!(df, (arr_id, arr_submodel, arr_eeps, arr_marker, arr_color))
        arr_id = []
        arr_submodel = []
        arr_eeps = []
        arr_marker = []
        arr_color = []
    end
    
    evol_colors =  color_list[1:length(markers)]
    return df, evol_colors

end