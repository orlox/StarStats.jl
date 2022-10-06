using CSV, DataFrames, CodecZlib

export ModelDataGrid, load_grid

mutable struct ModelDataGrid
    dfs::Array{DataFrame}
    inputs::Vector{Vector{String}}
    input_names::Vector{Symbol}
    input_values::Vector{Vector{Float64}}
    function ModelDataGrid(inputs, input_names)
        dimensions = [length(input) for input in inputs]
        input_values = [parse.(Float64, input) for input in inputs  ]
        dfs = Array{DataFrame}(undef,dimensions...)
        new(dfs,inputs,input_names, input_values)
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
