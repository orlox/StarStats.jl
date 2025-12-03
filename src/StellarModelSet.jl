using Base.Threads, LinearAlgebra

export StellarModelSet, interpolate_grid_quantity

mutable struct StellarModelSet{N,P,LU,E,V}
    models::Vector{SimulationData}
    inputs::Matrix{String}
    input_names::Vector{Symbol}
    input_values::Matrix{Float64}
    simplex_interpolant::SimplexInterpolant{N,P,LU,E,V}
    check_possibility_of_interpolation::Vector{Bool}
end
#Extra field to StellarModelSet a vector of booleans

function StellarModelSet(inputs, input_names, path_constructor, dataframe_loader, EEP_and_distance_calculator!, EEPs_symbols!; input_values = nothing)
    if isnothing(input_values)
        input_values = parse.(Float64, inputs)
    else
        if size(input_values) != size(inputs)
            throw(DomainError(input_values, "Length of input_names and input_values must match"))
        end
    end
    models = Vector{SimulationData}(undef,size(inputs)[2])
    # load up data
    #@threads for i in 1:(size(inputs)[2])
    for i in 1:(size(inputs)[2])
        strings = Vector{String}(undef, length(input_names))
        for j in eachindex(input_names)
            strings[j] = inputs[j, i]
        end
        models[i] = SimulationData(strings, input_names, path_constructor, dataframe_loader, EEP_and_distance_calculator!, EEPs_symbols!)
    end
    simplex_interpolant = SimplexInterpolant(input_values)
    all_good = 1
    check_interpolation = Bool[]
    #if file already exists -> delete it
    if isfile("output.txt")
        rm("output.txt")  
    end
    for simplex in simplex_interpolant.simplexes
        all_good = check_coords_of_simplex(simplex, models)
        all_good = chaeck_names_of_EEPs(simplex, models)
        coords_of_bad_symplex(all_good, simplex, models)
        push!(check_interpolation , all_good)
       # println(all_good, " ",simplex.id)
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

function coords_of_bad_symplex(all_good, simplex, models)
    if all_good == 0
        println(simplex.id, simplex.point_indeces)
        sorted_indexes = sort(simplex.point_indeces) #sort the array with indexes of models
        dist_max = sorted_indexes[2] - sorted_indexes[1]
        dict = [sorted_indexes[1], sorted_indexes[2]]
        for i in 1:length(sorted_indexes)
            for k in i+1:length(sorted_indexes)
                dist_calc = sorted_indexes[k] - sorted_indexes[i]
                if (dist_calc > dist_max)
                    dist_max = sorted_indexes[k] - sorted_indexes[i]
                    dict = [sorted_indexes[k],sorted_indexes[i]]
                end
            end
        end
        println(dist_max, dict)
        #println(models[dict[1]].input_params, models[dict[2]].input_params)
       
        open("output.txt", "a") do f
            for s in 1:length(models[dict[1]].input_params)
                mean_param = (parse(Float64,models[dict[1]].input_params[s])+parse(Float64,models[dict[2]].input_params[s]))/2
                mean_param = round(mean_param, digits=5) 
                write(f, string(mean_param), " ")
            end
            write(f, "\n")      
        end  
    end


end