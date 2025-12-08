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

    if isfile(file_name)
        rm(file_name)  
    end
    f = open(file_name, "a") 
    models = model_set.models
    names_1 = models[1].input_names
    for i in 1:length(names_1)
            write(f,names_1[i], " ")
    end
    write(f, "\n")  

    params_array = []
    id_s = []
    for i in eachindex(model_set.simplex_interpolant.simplexes)
        simplex = model_set.simplex_interpolant.simplexes[i]

        if !model_set.can_interpolate_simplex[i]  
            append!(id_s, simplex.id)
            sum_sqrt = 0
            dict = [simplex.point_indeces[1], simplex.point_indeces[2]]
            for i in 1:length(models[1].input_params)
                first = simplex.point_indeces[1]
                second = simplex.point_indeces[2]
                sum_sqrt  = sum_sqrt +sqrt((parse(Float64,models[first].input_params[i])-parse(Float64,models[second].input_params[i]))^2)
            end
            println(simplex.id," ", sum_sqrt, " ", 0000)
            println(simplex.id, " ", dict)
            for k in 1: length(simplex.point_indeces)
                for s in k+1: length(simplex.point_indeces)
                    sum_sqr_calc = 0
                    for m in 1:length(models[1].input_params)
                        k_th_model = simplex.point_indeces[k]
                        s_th_model = simplex.point_indeces[s]
                        sum_sqr_calc = sum_sqr_calc + sqrt((parse(Float64,models[k_th_model].input_params[m])-parse(Float64,models[s_th_model].input_params[m]))^2)
                    end
                    println(simplex.id," ", sum_sqr_calc, " ", 1111)
                    if sum_sqr_calc > sum_sqrt
                        dict = [simplex.point_indeces[k], simplex.point_indeces[s]]
                        sum_sqrt = sum_sqr_calc
                    end
                    println(simplex.id, " ", dict)
                end
            end
            #The idea is to devide this edge into half and provide new parameters to calculated models  
            for s in 1:length(models[dict[1]].input_params)
                mean_param = (parse(Float64,models[dict[1]].input_params[s])+parse(Float64,models[dict[2]].input_params[s]))/2
                mean_param = round(mean_param, digits=5) 
                append!(params_array, mean_param)
                write(f, string(mean_param), " ")
            end
            write(f, "\n")          
        end
    end
    close(f)
    number_of_params = length(models[1].input_params)
    return number_of_params, id_s, params_array
end

 """   
    params_for_simulations(number_of_params, params_array)
    Independently of dimensions cunstracts a matrix with the paramenters for simulations to be run. This can be used for plotting 
 """

function params_for_simulations(number_of_params, params_array)
    length_of_data = Integer(length(params_array)/number_of_params)
    matrix_of_params = Matrix{Float64}(undef,number_of_params,length_of_data)

    sch = 1
    for i in 1:Integer(length(params_array)/number_of_params)
        for k in 1:number_of_params
            matrix_of_params[k,i] = params_array[sch+(k-1)]
        end
        println(sch)
        sch = sch+number_of_params
    end

    return matrix_of_params
end





