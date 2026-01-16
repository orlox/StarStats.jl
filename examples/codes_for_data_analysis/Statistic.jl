using Statistics


# this function returnes a matrix with values of diff of distances along the curves with same evolution
function statistics_for_dataset(model_set)

indices = findall(x -> x == 1, model_set.can_interpolate_simplex)
number_of_simplexes = length(indices) # we take only simplexes where we can interpolate
#we construc dataframe with comparisions for 1st simplex and take the len so we know number of comparisions
number_of_comparision = length(StarStats.length_between_EEPs(model_set.simplex_interpolant.simplexes[1], models, model_set.metric).difference)
matrix_total_for_statistic = zeros(number_of_simplexes*number_of_comparision )  # Here we store difference between distacenses of 2 eeps on one curve between 2 curves
matrix_total_for_statistic_point_to_point = zeros(number_of_simplexes*number_of_comparision )  # Here we store distance between EEPs in 2 curves
number = 0
for i in eachindex(model_set.simplex_interpolant.simplexes)
    if model_set.can_interpolate_simplex[i] == 1 #this checks that all all models have same number of eeps in simplex
        models = model_set.models
        simplex = model_set.simplex_interpolant.simplexes[i]
        df = StarStats.length_between_EEPs(simplex, models, model_set.metric)
        for m in 1:length(df.difference)
            number = number+1
            matrix_total_for_statistic[number] = df.difference[m]
            matrix_total_for_statistic_point_to_point[number] = df.distance_between_eeps[m]
        end
    end
end
return  matrix_total_for_statistic, matrix_total_for_statistic_point_to_point

end

#function provides which number of eeps can be in subset
function find_uniq_number_of_eeps(model_set)
    markers = make_huge_vector_of_eeps(model_set)

    array_of_lengths = zeros(length(markers))
    for i in 1:length(markers)
        array_of_lengths[i] = length(markers[i])
    end
    return array_of_lengths
end

function IQR(data)

    Q1 = quantile(data, 0.25)
    Q3 = quantile(data, 0.75)
    iqr = Q3 - Q1

    upper = Q3 + 1.5 * iqr #k = 1.5
    return upper
end

#data return dataframe with flags that suggest extra refinmnet

function filter_data(df,upper_outlier,upper_outlier_point_to_point)
    number_of_comparisions = length(df.difference)
    extra_refinment_flag = zeros(number_of_comparisions)
    for m in 1:length(df.difference)
        if (df.difference[m] > upper_outlier) || (df.distance_between_eeps[m] > upper_outlier_point_to_point) 
            extra_refinment_flag[m] = 1
        else
            extra_refinment_flag[m] = 0
        end
    end
    #HERE WE ACTUALLY CAN ADD IT TO THE SIMPLEX
    df[!, :extra_refinment_flag] = extra_refinment_flag
    return df
end

#precompute everything for extra refinment
function precompute_data_for_extra_refinment(model_set)

    matrix_total_for_statistic, matrix_total_for_statistic_point_to_point = statistics_for_dataset(model_set)
    upper_outlier = IQR(matrix_total_for_statistic)
    upper_outlier_point_to_point = IQR(matrix_total_for_statistic_point_to_point)

    return  upper_outlier, upper_outlier_point_to_point 
end

function suggest_refinment(model_set)

    upper_outlier, upper_outlier_point_to_point  =  precompute_data_for_extra_refinment(model_set)
    id_s = []
    params_array = []
    for i in eachindex(model_set.simplex_interpolant.simplexes)
        if model_set.can_interpolate_simplex[i] == 1
            models = model_set.models
            simplex = model_set.simplex_interpolant.simplexes[i]
            df = StarStats.length_between_EEPs(simplex, models, model_set.metric)
            df_new = filter_data(df,upper_outlier,upper_outlier_point_to_point)
            for m in 1:length(df_new.extra_refinment_flag)
                if df_new.extra_refinment_flag[m] == 1 #need extra refinment
                    append!(id_s, simplex.id)
                    for s in 1:length(models[df_new.model_pair[1][1]].input_params)

                        model1 = models[df_new.model_pair[m][1]] #mth model which need ref and 1 is 1st model in pair
                        model2 = models[df_new.model_pair[m][2]] #mth model which need ref and 2 is 2nd model in pair
                        mean_param = (parse(Float64,model1.input_params[s])+parse(Float64,model2.input_params[s]))/2
                        mean_param = round(mean_param, digits=5) 
                        append!(params_array, mean_param)  
                        #params array is array with mean params to show where refinment will be needed 
                    end
                end
            end
        end
    end
    number_of_params = length(models[1].input_params)
    matrix_of_params = StarStats.params_for_simulations(number_of_params, params_array)

    matrix_of_params_uniq,id_s_uniq = StarStats.delete_same_models_2D(matrix_of_params,id_s)   
    return   matrix_of_params_uniq, id_s_uniq        
end
