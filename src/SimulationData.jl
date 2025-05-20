using DataFrames, CSV

export mesa_dataframe_loader

mutable struct SimulationData{A,B}
    df::DataFrame
    EEPs::Vector{Int}
    input_params::A
    input_names::B
end

function SimulationData(input_strings, input_names, path_constructor, dataframe_loader, EEP_and_distance_calculator!)
    path = path_constructor(input_strings)
    if !(isfile(path) || isdir(path))
        error("Constructed path $(path) does not exist")
    end
    df = dataframe_loader(path)
    EEPs = EEP_and_distance_calculator!(df)
    return SimulationData(df, EEPs, input_strings, input_names)
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