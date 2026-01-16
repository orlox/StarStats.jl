using Base.Threads, LinearAlgebra

export StellarModelSet, interpolate_grid_quantity

mutable struct StellarModelSet{N,P,LU,E,V}
    models::Vector{SimulationData}
    inputs::Matrix{String}
    input_names::Vector{Symbol}
    input_values::Matrix{Float64}
    simplex_interpolant::SimplexInterpolant{N,P,LU,E,V}
    can_interpolate_simplex::Vector{Bool}
    metric
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
        make_distance_along_curve(models[i], metric)
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

function suggest_new_simulations(model_set::StellarModelSet, file_name::String)
    #delete file if one already exists in the directory
    if isfile(file_name)
        rm(file_name)  
    end
    f = open(file_name, "a") 

    #create modelset 
    models = model_set.models
    names_1 = models[1].input_names
    #write parameter names in the file
    for i in 1:length(names_1)
            write(f,names_1[i], " ")
    end
    write(f, "\n")  

    #create arrays to return
    params_array = []
    id_s = []
    for i in eachindex(model_set.simplex_interpolant.simplexes)
        simplex = model_set.simplex_interpolant.simplexes[i]
        if !model_set.can_interpolate_simplex[i]  
            append!(id_s, simplex.id)
            dict  = calculate_longest_edge(models, simplex)
        
            #The idea is to devide this edge into half and provide new parameters to calculated models  
            for s in 1:length(models[dict[1]].input_params)
                mean_param = (parse(Float64,models[dict[1]].input_params[s])+parse(Float64,models[dict[2]].input_params[s]))/2
                mean_param = round(mean_param, digits=5) 
                append!(params_array, mean_param)  
            end
            
        end
    end
    
    number_of_params = length(models[1].input_params)
    matrix_of_params = params_for_simulations(number_of_params, params_array) #make a matrix with data instead of arrays
    matrix_of_params_uniq,id_s_uniq = delete_same_models_2D(matrix_of_params,id_s) #delete same elements
    number_of_params_uniq = length(id_s_uniq)

    for i in 1:length(matrix_of_params_uniq[1,:])
        write(f, string(matrix_of_params_uniq[1,i]), " ",string(matrix_of_params_uniq[2,i]))
        write(f, "\n")          
    end
    close(f)
    return number_of_params_uniq, id_s_uniq, matrix_of_params_uniq
end


 """   
    params_for_simulations(number_of_params, params_array)
    Independently of dimensions cunstracts a matrix with the paramenters for simulations to be run. This can be used for plotting 
 """

function calculate_longest_edge(models, simplex)
    sum_sqrt = 0
    dict = [simplex.point_indeces[1], simplex.point_indeces[2]]
    #calculate the 1st edge of simplex
    for i in 1:length(models[1].input_params)
        first = simplex.point_indeces[1]
        second = simplex.point_indeces[2]
        sum_sqrt  = sum_sqrt +sqrt((parse(Float64,models[first].input_params[i])-parse(Float64,models[second].input_params[i]))^2)
    end
    #calculate other edges and compare 
    for k in 1: length(simplex.point_indeces)
        for s in k+1: length(simplex.point_indeces)
            sum_sqr_calc = 0

            for m in 1:length(models[1].input_params)
                k_th_model = simplex.point_indeces[k]
                s_th_model = simplex.point_indeces[s]
                sum_sqr_calc = sum_sqr_calc + sqrt((parse(Float64,models[k_th_model].input_params[m])-parse(Float64,models[s_th_model].input_params[m]))^2)
            end
            if sum_sqr_calc > sum_sqrt
                dict = [simplex.point_indeces[k], simplex.point_indeces[s]]
                sum_sqrt = sum_sqr_calc
            end
        end
    end
    return dict
end
#Function used for suggesting new simulations
#IS IT MULIDIMENSIONAL
function params_for_simulations(number_of_params, params_array)
    length_of_data = Integer(length(params_array)/number_of_params)
    matrix_of_params = Matrix{Float64}(undef,number_of_params,length_of_data)

    sch = 1
    for i in 1:Integer(length(params_array)/number_of_params)
        for k in 1:number_of_params
            matrix_of_params[k,i] = params_array[sch+(k-1)]
        end
        sch = sch+number_of_params
    end

    return matrix_of_params
end
#function used for suggesting new simulations
function delete_same_models_2D(matrix_of_params,id_s)
    arr1 = []
    arr2 = []
    arr3 = []
   
    for i in 1:length(matrix_of_params[1,:])
        a = 0
        for k in i+1:length(matrix_of_params[1,:])
            if (matrix_of_params[1,i] == matrix_of_params[1,k] ) .&& (matrix_of_params[2,i] == matrix_of_params[2,k] )
                a = 1 # marker of coindsiding paors
            end
        end
        if a == 0
            append!(arr1, matrix_of_params[1,i])
            append!(arr2, matrix_of_params[2,i])
            append!(arr3, id_s[i])
        end
    end
    number_of_params = 2
    length_of_data = length(arr1)
    matrix_of_params_new = Matrix{Float64}(undef,number_of_params,length_of_data)
    matrix_of_params_new[1,:] = arr1
    matrix_of_params_new[2,:] = arr2
    id_s_new = arr3
    return matrix_of_params_new,id_s_new

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
                    point1 = (models[k_th_model].df.logTeff[EEP_k_index], models[k_th_model].df.logL[EEP_k_index])
                    point2 = (models[s_th_model].df.logTeff[EEP_s_index], models[s_th_model].df.logL[EEP_s_index])
                    pp = PointPair2(point1, point2, metric)
                    
                    temp_array2[m] = distance_point_pair(pp)    
                end

            end
            
            #As we have length(models[k_th_model].EEPs) -1 we need to calculate distance between 2 end points
            EEP_k_end = models[k_th_model].EEPs[end]
            EEP_s_end = models[s_th_model].EEPs[end]
            
            #Broat loader check
            if (EEP_k_end != 0) && (EEP_s_end != 0)

                point1 = (models[k_th_model].df.logTeff[EEP_k_end], models[k_th_model].df.logL[EEP_k_end])
                point2 = (models[s_th_model].df.logTeff[EEP_s_end], models[s_th_model].df.logL[EEP_s_end])
                pp = PointPair2(point1, point2, metric)
                temp_array2[end] = distance_point_pair(pp) 
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


