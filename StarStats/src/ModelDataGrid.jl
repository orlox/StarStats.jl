using CSV, DataFrames, Base.Threads

export ModelDataGrid, load_grid, compute_distances_and_EEPs, gz_dataframe_loader_with_Teff_and_star_age_fix

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
mutable struct ModelDataGrid
    dfs::Array{DataFrame}
    inputs::Vector{Vector{String}}
    input_names::Vector{Symbol}
    input_values::Vector{Vector{Float64}}
    EEPs::Array{Int64}
    Xc_TAMS::Float64
    function ModelDataGrid(inputs, input_names; input_values = nothing, Xc_TAMS = 1e-12)
        dimensions = [length(input) for input in inputs]
        if isnothing(input_values)
            input_values = [parse.(Float64, input) for input in inputs]
        else
            if length(input_values) != length(input_names)
                throw(DomainError(input_values, "Length of input_names and input_values must match"))
            end
        end
        dfs = Array{DataFrame}(undef,dimensions...)
        EEPs = zeros(dimensions...,6) # We consider six EEPs right now
        new(dfs,inputs,input_names, input_values, EEPs, Xc_TAMS)
    end
end

"""
    load_grid

Function that does bla

# Fields: 
-  grid:
-  path_constructor
-  dataframe_loader:
"""

function load_grid(grid::ModelDataGrid, path_constructor, dataframe_loader)
    @threads for index in collect(Base.product([1:length(input) for input in grid.inputs]...))
        strings = Vector{String}(undef, length(index))
        for (i,j) in enumerate(index)
            strings[i] = grid.inputs[i][j]
        end
        path = path_constructor(strings)
        if !isfile(path)
            continue
        end
        grid.dfs[index...] = dataframe_loader(path)
    end
end

"""
    mesa_dataframe_loader

Function 

# Fields: 
-   
"""

function mesa_dataframe_loader(path)
    df = DataFrame(CSV.File(path, delim=" ", ignorerepeated=true, skipto=7, header=6))
    df[!,:c_h1] = copy(df.center_h1)
    df[!,:c_he4] = copy(df.center_he4)
    df[!,:logTeff] = copy(df.log_Teff)
    df[!,:logL] = copy(df.log_L)
    return df
end

"""
   compute_distances_and_EEPs

Function 

# Fields: 
-   
"""
function compute_distances_and_EEPs(grid::ModelDataGrid)
    for index in Base.product([1:length(input) for input in grid.inputs]...)
        if !isassigned(grid.dfs,index...)
            continue
        end

        df = grid.dfs[index...]
        grid.EEPs[index...,:] = get_EEPs(df, grid.Xc_TAMS)

        distance = zeros(size(df,1))
        delta_log_Teff = 0
        delta_log_L = 0

        for i in 2:size(df,1)
            delta_log_Teff = df.logTeff[i]-df.logTeff[i-1]
            delta_log_L = df.logL[i]-df.logL[i-1]

            distance[i] = distance[i-1] + sqrt(delta_log_Teff^2 + delta_log_L^2)
        end
        df[!,:distance] = distance

        x = zeros(size(df,1))
        for j in 1:(length(grid.EEPs[index...,:])-1)
            start_EEP = grid.EEPs[index...,j]
            end_EEP = grid.EEPs[index...,j+1]
            if start_EEP==0 || end_EEP==0
                break
            end
            start_distance = df.distance[start_EEP]
            end_distance = df.distance[end_EEP]
            for i in start_EEP:end_EEP
                x[i] = (df.distance[i]-start_distance)/(end_distance-start_distance) + (j-1)
            end
        end
        for i in size(df,1):-1:1
            if x[i] > 0
                x[i:end] .= x[i]
                break
            end
        end
        df[!,:x] = x

        dtdx = zeros(size(df,1))
        for i in 2:size(df,1)
            if (df.x[i]-df.x[i-1])<=0
                continue
            end
            dtdx[i] = (df.star_age[i]-df.star_age[i-1])/(df.x[i]-df.x[i-1])
        end
        df[!,:dtdx] = dtdx
    end
end
