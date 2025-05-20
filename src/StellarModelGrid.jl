using Base.Threads

export StellarModelGrid, interpolate_grid_quantity

"""
    ModelDataGrid

Struct containing the definition of a grid including the input variables defining it and their values 

# Fields: 
-   dfs: Multidimensional array of DataFrames containing individual evolutionary track information from a grid of models.
        This is populated with the load_grid function.    
- inputs: Vector of String vectors containing the parth of the path from a track. (TODO give an example)
- input_names: Vector of Symbol containing the names of each parameter of the grid.  
- input_values: Result of parsing the inputs into floats.   
- EEPs: Array of integers representing the Equivelant Evolutionary Points (see Dotter, 2016)
- Xc_TAMS: Float denoting the limit  in central hydrogen at which one defines the TAMS
"""
mutable struct StellarModelGrid
    models::Array{SimulationData}
    inputs::Vector{Vector{String}}
    input_names::Vector{Symbol}
    input_values::Vector{Vector{Float64}}
    function StellarModelGrid(inputs, input_names, path_constructor, dataframe_loader, EEP_and_distance_calculator!; input_values = nothing)
        dimensions = [length(input) for input in inputs]
        if isnothing(input_values)
            input_values = [parse.(Float64, input) for input in inputs]
        else
            if length(input_values) != length(input_names)
                throw(DomainError(input_values, "Length of input_names and input_values must match"))
            end
        end
        models = Array{SimulationData}(undef,dimensions...)
        # load up data
        @threads for index in collect(Base.product([1:length(input) for input in inputs]...))
            strings = Vector{String}(undef, length(index))
            for (i,j) in enumerate(index)
                strings[i] = inputs[i][j]
            end

            models[index...] = SimulationData(strings, input_names, path_constructor, dataframe_loader, EEP_and_distance_calculator!)
        end
        new(models,inputs,input_names, input_values)
    end
end

function interpolate_grid_quantity(grid, grid_parameters, interpolated_quantity, x::T) where{T}
    return interpolate_grid_quantity_internal(grid.models, grid.input_values, grid_parameters, interpolated_quantity, x)
end

function model_index_element(i, lower_i, cartesian_index)
    return cartesian_index[i]+lower_i[i]-1
end
function interpolator_index1_element(j,i)
    if j<= i
        return 1
    end
    return 2
end
function interpolator_index2_element(j,i,index1)
    if j< i
        return 1
    elseif j==i
        return 2
    end
    return index1[j]
end

function interpolate_grid_quantity_internal(models, input_values, grid_parameters, interpolated_quantity, x::T) where{T}
    dimensions = ndims(models)
    lower_i = zeros(Int, dimensions)
    grid_values = zeros(dimensions,2)
    # get values of parameters at edges
    for j in 1:dimensions
        for i in 1:(length(input_values[j])-1)
            if input_values[j][i] <= grid_parameters[j] && input_values[j][i+1] >= grid_parameters[j]
                grid_values[j,1] = input_values[j][i]
                grid_values[j,2] = input_values[j][i+1]
                lower_i[j] = i
                break
            end
        end
    end
    for i in lower_i
        if i == 0
            return NaN
        end
    end

    # get value of quantity to interpolate at vertices at a given x
    yvalues::Array{T,dimensions} = zeros(T,ntuple(_->2, dimensions)...)
    for index in CartesianIndices(yvalues)
        model_index = ntuple(i->model_index_element(i,lower_i,index),dimensions)
        try
            model = models[CartesianIndex(model_index)]
            yvalues[index] = get_secondary_EEP(x, model.df, interpolated_quantity)
        catch
            return NaN
        end
    end

    for i in 1:dimensions
        xval1 = grid_values[i,1]
        xval2 = grid_values[i,2]
        xtarget = grid_parameters[i]
        for index1 in CartesianIndices(ntuple(j->interpolator_index1_element(j,i),dimensions))
            index2 = CartesianIndex(ntuple(j->interpolator_index2_element(j,i, index1),dimensions))
            yval1 = yvalues[index1]
            yval2 = yvalues[index2]
            yvalues[index1] = yval1 +(yval2-yval1)*(xtarget-xval1)/(xval2-xval1)
        end
    end

    return yvalues[1]
end