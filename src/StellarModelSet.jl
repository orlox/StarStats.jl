using Base.Threads, LinearAlgebra

export StellarModelSet, interpolate_grid_quantity

mutable struct StellarModelSet{N,P,LU,E,V,METRIC}
    models::Vector{SimulationData}
    inputs::Matrix{String}
    input_names::Vector{Symbol}
    input_values::Matrix{Float64}
    simplex_interpolant::SimplexInterpolant{N,P,LU,E,V}
    can_interpolate_simplex::Vector{Bool}
    metric::METRIC
end
#Extra field to StellarModelSet a vector of booleans

function StellarModelSet(inputs, input_names, path_constructor, dataframe_loader, EEP_and_distance_calculator!, metric; input_values = nothing) 
  
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
        distance_along_curve!(models[i], metric)
    end
    simplex_interpolant = SimplexInterpolant(input_values)
    all_good = 1
    check_interpolation = zeros(Bool,length(simplex_interpolant.simplexes))
   
    for i in eachindex(simplex_interpolant.simplexes)
        simplex = simplex_interpolant.simplexes[i]
        all_good = check_names_of_EEPs(simplex, models)
        check_interpolation[i] =  all_good
    end

    return StellarModelSet(models,inputs,input_names, input_values, simplex_interpolant, check_interpolation, metric)
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
#refine_model_set_working_interpolation
function refine_model_set_failed_interpolation(model_set::StellarModelSet)#, file_name::String)
    ##delete file if one already exists in the directory
    #if isfile(file_name)
    #    rm(file_name)  
    #end
    #f = open(file_name, "a") 

    ##create modelset 
    #models = model_set.models
    #names = model_set.input_names
    ##write parameter names in the file
    #for i in 1:length(names)
    #        write(f,names[i], " ")
    #end
    #write(f, "\n")

    suggested_sims = Dict()
    for input_param in model_set.input_params
        suggested_sims[input_param] = zeros(Float64, 0)
    end
    suggested_sims[:distance] = zeros(Float64, 0)
    suggested_sims[:simplex_id] = zeros(Int, 0)

    #create arrays to return
    for i in eachindex(model_set.simplex_interpolant.simplexes)
        simplex = model_set.simplex_interpolant.simplexes[i]
        if !model_set.can_interpolate_simplex[i]  
            append!(suggested_sims[simplex_id], simplex.id)
            model_ids, dist = calculate_longest_edge(models, simplex) #TODO
        
            #The idea is to devide this edge into half and provide new parameters to calculated models  
            for s in 1:length(model_set.input_params)
                mean_param = (parse(Float64,models[models_ids[1]].input_params[s])+parse(Float64,models[model_ids[2]].input_params[s]))/2
                append!(suggested_sims[model_set.input_params[s]], mean_param)  
            end
            append!(suggested_sims[:distance], dist)
        end
    end

    new_params = [[suggested_sims[name][i] for name in model_set.input_names] for i in 1:length(suggested_sims[:simplex_id])]
    unique_filter = unique(i -> new_params[i], 1:length(new_params))

    for input_param in model_set.input_params
        suggested_sims[input_param] = suggested_sims[input_param][unique_filter]
    end
    suggested_sims[:distance] = suggested_sims[:distance][unique_filter]
    suggested_sims[:simplex_id] = suggested_sims[:simplex_id][unique_filter]

    return suggested_sims
    
    #number_of_params = length(models[1].input_params)
    #matrix_of_params = params_for_simulations(number_of_params, params_array) #make a matrix with data instead of arrays
    #matrix_of_params_uniq,id_s_uniq = delete_same_models_2D(matrix_of_params,id_s) #delete same elements
    #number_of_params_uniq = length(id_s_uniq)

    #for i in 1:length(matrix_of_params_uniq[1,:])
    #    write(f, string(matrix_of_params_uniq[1,i]), " ",string(matrix_of_params_uniq[2,i]))
    #    write(f, "\n")          
    #end
    #close(f)
    #return number_of_params_uniq, id_s_uniq, matrix_of_params_uniq
end

 """   
    params_for_simulations(number_of_params, params_array)
    Independently of dimensions cunstracts a matrix with the paramenters for simulations to be run. This can be used for plotting 
 """

function calculate_longest_edge(models, simplex)
    max_distance = -1
    # TODO: add squared values, take sqrt after
    #calculate other edges and compare 
    for k in 1: length(simplex.point_indeces)
        for s in k+1: length(simplex.point_indeces)
            sum_sqr_calc = 0

            for m in 1:length(models[1].input_params)
                k_th_model = simplex.point_indeces[k]
                s_th_model = simplex.point_indeces[s]
                sum_sqr_calc = sum_sqr_calc + sqrt((parse(Float64,models[k_th_model].input_params[m])-parse(Float64,models[s_th_model].input_params[m]))^2)
            end
            if max_distance == -1 && sum_sqr_calc > max_distance
                model_ids = [simplex.point_indeces[k], simplex.point_indeces[s]]
                max_distance = sum_sqr_calc
            end
        end
    end
    return model_ids, max_distance
end

#function used for statistics and refinment
#returns df with the largest difference between distance between 2eeps on 2 tracs for 2 models, 
#as well returns laregs difference between 2 eeps on 2 curves
#and the number of 2 model for which these distances are calculated
function length_between_EEPs(simplex, models, metric)
    
    number_of_comparision = Integer(length(simplex.point_indeces)*(length(simplex.point_indeces)-1)/2)
    length_of_differences = length(models[simplex.point_indeces[1]].EEPs) -1
    array_of_difference = zeros(number_of_comparision) # here we keep difference of distances between 2 eeps for each curve between 2 curves 
    array_of_difference_point_to_point = zeros(number_of_comparision) #here we keep differences of distances between 2 points on 2 curves
    temp_array =  zeros(length_of_differences)
    temp_array2 =  zeros(length_of_differences+1) #+1 because we dont calculate difference between eeps on one curve, we calculate difference between curves
    temp_array_models =  Vector{Any}(undef, number_of_comparision)

    number = 0
    for k in 1: length(simplex.point_indeces)
        for s in k+1: length(simplex.point_indeces)
            number = number +1
            k_th_model = simplex.point_indeces[k]
            s_th_model = simplex.point_indeces[s]
            for m in 1:length(models[k_th_model].EEPs) -1
                EEP_k_index = models[k_th_model].EEPs[m]
                EEP_k1_index = models[k_th_model].EEPs[m+1]
                EEP_s_index = models[s_th_model].EEPs[m]
                EEP_s1_index = models[s_th_model].EEPs[m+1]
               
                #special check for brott  models
                if (EEP_k_index != 0) && (EEP_k1_index !=0) && (EEP_s_index != 0) && (EEP_s1_index != 0)
                    # Here we calculate difference between distances in 2 curves 
                    # curve1(EEP2-EEP1) - curve2(EEP2-EEP1)
                    k_th_model_dist = models[k_th_model].df.distance[EEP_k1_index] - models[k_th_model].df.distance[EEP_k_index]
                    s_th_model_dist = models[s_th_model].df.distance[EEP_s1_index] - models[s_th_model].df.distance[EEP_s_index]
                    temp_array[m] = abs(k_th_model_dist - s_th_model_dist ) # abs because we dont know which has longer dist between eeps
                    # Here we calcaulte distance between EEPs in 2 curves
                    #p1(Tk,Lk) (kth model), p2(Ts, Ls) (sth model)

                    temp_array2[m] = dist(models[k_th_model], models[s_th_model], EEP_k_index, EEP_s_index, metric)
                end

            end
            
            #As we have length(models[k_th_model].EEPs) -1 we need to calculate distance between 2 end points
            EEP_k_end = models[k_th_model].EEPs[end]
            EEP_s_end = models[s_th_model].EEPs[end]
            
            #Broat loader check
            if (EEP_k_end != 0) && (EEP_s_end != 0)
                temp_array2[end] = dist(models[k_th_model], models[s_th_model], EEP_k_end, EEP_s_end, metric)
            end


            temp_array_models[number] = [k_th_model, s_th_model]
            array_of_difference[number] = maximum(temp_array)
            array_of_difference_point_to_point[number] = maximum(temp_array2)
            temp_array = zeros(length_of_differences)

        end
    end
    df = DataFrame(
      difference = array_of_difference,   # biggest difference 
      distance_between_eeps = array_of_difference_point_to_point, # biggest distance between eeps
      model_pair = temp_array_models      #models numbers
    )
    return df
end


