using Base.Threads, LinearAlgebra

export StellarModelSet, interpolate_grid_quantity

mutable struct StellarModelSet{N,P,LU,E,V}
    models::Vector{SimulationData}
    inputs::Matrix{String}
    input_names::Vector{Symbol}
    input_values::Matrix{Float64}
    simplex_interpolant::SimplexInterpolant{N,P,LU,E,V}
    can_interpolate_simplex::Vector{Bool}
end
#Extra field to StellarModelSet a vector of booleans

function StellarModelSet(inputs, input_names, path_constructor, dataframe_loader, EEP_and_distance_calculator!; input_values = nothing) 
    if isnothing(input_values)
        input_values = parse.(Float64, inputs)
    else
        if size(input_values) != size(inputs)
            throw(DomainError(input_values, "Length of input_names and input_values must match"))
        end
    end
    models = Vector{SimulationData}(undef,size(inputs)[2])
    # load up data
    @threads for i in 1:(size(inputs)[2])
        strings = Vector{String}(undef, length(input_names))
        for j in eachindex(input_names)
            strings[j] = inputs[j, i]
        end
        models[i] = SimulationData(strings, input_names, path_constructor, dataframe_loader, EEP_and_distance_calculator!)
    end
    simplex_interpolant = SimplexInterpolant(input_values)
    all_good = 1
    check_interpolation = zeros(Bool,length(simplex_interpolant.simplexes))
   
    for i in eachindex(simplex_interpolant.simplexes)
        simplex = simplex_interpolant.simplexes[i]
        all_good = check_names_of_EEPs(simplex, models)
        check_interpolation[i] =  all_good
    end

    return StellarModelSet(models,inputs,input_names, input_values, simplex_interpolant, check_interpolation)
end

function interpolate_grid_quantity(grid::StellarModelSet{N,P,LU,E,V}, grid_parameters, interpolated_quantity, x::T) where{N,P,LU,E,V,T}
    models = grid.models
    coords, indeces, simplex_id  = interpolation_info(grid_parameters,grid.simplex_interpolant)

    if maximum(coords)==0
        return NaN
    end

    # get value of quantity to interpolate at vertices at a given x
    yvalues::Vector{T} = zeros(T,length(grid.input_names)+1) # if data dimensions are N, interpolating simplex has N+1 points (eg triangle in 2D)
    for i in eachindex(indeces)
        model = models[indeces[i]]
        try
            yvalues[i] = get_secondary_EEP(x, model.df, interpolated_quantity)
        catch
            return NaN
        end
    end

    return dot(coords, yvalues)
end

"""
    coords_of_bad_symplex(all_good, simplex, models)
If interpolation is not recommended in the current simplex it returnes a file with params for the simulations that needs to be run!
"""

function suggest_new_simulations(model_set::StellarModelSet, file_name::String)
    f = open(file_name, "a") 
    models = model_set.models
    for i in eachindex(model_set.simplex_interpolant.simplexes)
        simplex = model_set.simplex_interpolant.simplexes[i]

        if !model_set.can_interpolate_simplex[i]  
            sorted_indexes = sort(simplex.point_indeces) #sorts the array with indexes of models
            dist_max = sorted_indexes[2] - sorted_indexes[1]
            dict = [sorted_indexes[1], sorted_indexes[2]]
            #calculate the largest edge
            for i in 1:length(sorted_indexes)
                for k in i+1:length(sorted_indexes)
                    dist_calc = sorted_indexes[k] - sorted_indexes[i]
                    if (dist_calc > dist_max)
                        dist_max = sorted_indexes[k] - sorted_indexes[i]
                        dict = [sorted_indexes[k],sorted_indexes[i]]
                    end
                end
            end
                #The idea is to devide this edge into half and provide new parameters to calculated models  
            for s in 1:length(models[dict[1]].input_params)
                mean_param = (parse(Float64,models[dict[1]].input_params[s])+parse(Float64,models[dict[2]].input_params[s]))/2
                mean_param = round(mean_param, digits=5) 
                write(f, string(mean_param), " ")
            end
            write(f, "\n")          
        end
    end
    close(f)
end

"""
    check_files_existing(simplex_interpolant , models)
checks the names of parameters which would be passed to run simmulations. Writes their names as strings to output file.
"""
#=
function check_files_existing(simplex_interpolant , models)

    simplex1 = simplex_interpolant.simplexes[1]
    names = string.(models[simplex1.point_indeces[1]].input_names)
    if isfile("output.txt")
        rm("output.txt")  
    end

    open("output.txt", "a") do f
        for i in 1:length(names)
            write(f,names[i], " ")
        end
        write(f, "\n")  
    end

end
=#