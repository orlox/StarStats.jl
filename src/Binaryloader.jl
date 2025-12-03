using Pkg
using DataFrames
using CSV
using Base.Filesystem

struct container_mesa_data
    mesa_data_donor:: Vector{DataFrame}
    mesa_data_accretor:: Vector{DataFrame}
    mesa_data_q:: Vector{Float64}   
    mesa_data_logP:: Vector{Float64}   
    mesa_folder_name:: Vector{String} #through this field we find the needed model
end

function get_subfolders(path::String)
    items = readdir(path)
    subfolders = String[]
    for item in items
        item_path = joinpath(path, item)
        if isdir(item_path)
            push!(subfolders, item)
        end
    end
    return subfolders
end

function construct_path()
    file_folder = "mesa_data/results/"
    name_list = get_subfolders(file_folder)
    path = file_folder.*name_list
    return path, name_list
end

function construct_path_to_model_n(n)
    file_folder = "mesa_data/results/"
    name_list = get_subfolders(file_folder)
    path = file_folder.*name_list[n]
    return path,name_list[n]
end


function download_binary_data(path, name_list)
    
    donor_tables = DataFrame[]
    accretor_tables = DataFrame[]
    q = []
    logP = []
    folder_name = []
    for number in 1:length(name_list)
        #mesa data read dfs to stuct
        path_donor = path[number]*"/LOGS1/history.data"
        path_accretor = path[number]*"/LOGS2/history.data"
        h_donor = DataFrame(CSV.File(path_donor, delim=" ", ignorerepeated=true, skipto=7, header=6))
        h_accretor = DataFrame(CSV.File(path_accretor, delim=" ", ignorerepeated=true, skipto=7, header=6))
        #from file_names read q and logP
        q_val = parse(Float64,split(name_list[number],"_")[1])
        logP_val = parse(Float64,split(name_list[number],"_")[2])

        push!(q,q_val)
        push!(logP,logP_val)
        push!(donor_tables, h_donor)
        push!(accretor_tables,h_accretor)
        push!(folder_name,path[number])
    end
    return container_mesa_data(donor_tables, accretor_tables, q, logP,folder_name)
end

path, name_list = construct_path()
container = download_binary_data(path, name_list)

struct simulation_data
    df::DataFrame
    EEPs::Vector{Int}
    q_and_logP::Vector{Float64}
    input_names::Vector{Symbol}
end

function construct_simulation_data(path,name,donor_or_accretor)

    if (donor_or_accretor == "donor")
        path_donor = path*"/LOGS1/history.data"
        h = DataFrame(CSV.File(path_donor, delim=" ", ignorerepeated=true, skipto=7, header=6))
    else
        path_accretor = path*"/LOGS2/history.data"
        h = DataFrame(CSV.File(path_accretor, delim=" ", ignorerepeated=true, skipto=7, header=6))
    end
    
    #from file_names read q and logP
    q = parse(Float64,split(name,"_")[1])
    logP = parse(Float64,split(name,"_")[2])
    EEPs = []
    q_and_logP =[q,logP]
    input_names = [:q, :logP]
    return simulation_data(h,EEPs,q_and_logP, input_names)
   
    
end

path_to_model, name = construct_path_to_model_n(1)
s1 = construct_simulation_data(path_to_model, name, "donor")

#construct_simulation_data(container,"mesa_data/results/0000.50000_0000.00000" )
