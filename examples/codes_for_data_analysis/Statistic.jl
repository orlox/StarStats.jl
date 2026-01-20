using Statistics

function IQR(data)

    Q1 = quantile(data, 0.25)
    Q3 = quantile(data, 0.75)
    iqr = Q3 - Q1

    upper = Q3 + 1.5 * iqr #k = 1.5
    return upper
end

#data return dataframe with flags that suggest extra refinmnet


function filter_data_new!(model_set)
    
    #extra_refinment_flag = zeros(length(refine_model_set_bad_resolution[:difference]))
    refine_dict = StarStats.refine_model_set_bad_resolution(model_set)
    upper_outlier = IQR(refine_dict[:difference])
    upper_outlier_point_to_point = IQR(refine_dict[:distance_between_eeps])
    refine_dict[:extra_refinment_flag] = zeros(Int,0)
    for m in 1:length(refine_dict[:difference])
        if (refine_dict[:difference][m] > upper_outlier) || (refine_dict[:distance_between_eeps][m] > upper_outlier_point_to_point) 
            append!(refine_dict[:extra_refinment_flag],1)
        else
            append!(refine_dict[:extra_refinment_flag],0)
        end
    end
    #HERE WE ACTUALLY CAN ADD IT TO THE SIMPLEX
    
    return refine_dict
end

function suggested_refinment_after_IQR(model_set)

    refine_dict = filter_data_new!(model_set)
    params_array = Dict()
    for input_param in model_set.input_names
        params_array[input_param] = zeros(Float64, 0)
    end
    models = model_set.models
    for m in 1:length(refine_dict[:extra_refinment_flag])
        println(m)
        if refine_dict[:extra_refinment_flag][m] == 1 #need extra refinment
            for s in 1:length(model_set.input_names)
                model1 = models[refine_dict[:models][m][1]] #mth model which need ref and 1 is 1st model in pair
                model2 = models[refine_dict[:models][m][2]] #mth model which need ref and 2 is 2nd model in pair
                mean_param = (parse(Float64,model1.input_params[s])+parse(Float64,model2.input_params[s]))/2
                mean_param = round(mean_param, digits=5) 
                append!(params_array[model_set.input_names[s]], mean_param)  
                #params array is array with mean params to show where refinment will be needed 
            end
        end
    end
   

    return  params_array    
end
