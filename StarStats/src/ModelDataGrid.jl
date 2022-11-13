using CSV, DataFrames, CodecZlib

export ModelDataGrid, load_grid, compute_distances_and_EEPs

"""
    ModelDataGrid

Struct containing the definition of a grid including the input variables defining it and their values 

# Fields: 
-   dfs: Multidimensional array of DataFrames containing individual evolutionary track information from a grid of models.
        This is populated with the load_grid function.    
- inputs: Vector of String vectors containing the parth of the path from a track. (TODO give an example)
- input_names: Vector of Symbol containing the names of each parameter of the grid.  
- input_values: Result of parsing the inputs into floats.   
"""
mutable struct ModelDataGrid
    dfs::Array{DataFrame}
    inputs::Vector{Vector{String}}
    input_names::Vector{Symbol}
    input_values::Vector{Vector{Float64}}
    EEPs::Array{Int64}
    function ModelDataGrid(inputs, input_names)
        dimensions = [length(input) for input in inputs]
        input_values = [parse.(Float64, input) for input in inputs  ]
        dfs = Array{DataFrame}(undef,dimensions...)
        EEPs = zeros(dimensions...,6) # We consider six EEPs right now
        new(dfs,inputs,input_names, input_values, EEPs)
    end
end

function load_grid(grid::ModelDataGrid, path_constructor)
    for index in Base.product([1:length(input) for input in grid.inputs]...)
        strings = Vector{String}(undef, length(index))
        for (i,j) in enumerate(index)
            strings[i] = grid.inputs[i][j]
        end
        path = path_constructor(strings)
        if !isfile(path)
            continue
        end
        file = GzipDecompressorStream(open(path))
        grid.dfs[index...] = CSV.read(file, DataFrame, delim=" ", ignorerepeated=true)
    end
end

function compute_distances_and_EEPs(grid::ModelDataGrid)
    for index in Base.product([1:length(input) for input in grid.inputs]...)
        if !isassigned(grid.dfs,index...)
            continue
        end

        df = grid.dfs[index...]
        grid.EEPs[index...,:] = get_EEPs(df)

        distance = zeros(size(df,1))
        delta_log_Teff = 0
        delta_log_L = 0

        for i in 2:size(df,1)
            delta_log_Teff = df.Teff[i]-df.Teff[i-1]
            delta_log_L = df.logL[i]-df.logL[i-1]

            distance[i] = distance[i-1] + sqrt(delta_log_Teff^2 + delta_log_L^2)
        end

        df[!,:distance] = distance
    end
end
