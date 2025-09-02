using FTPClient, DataFrames, CSV

#type can be :all, :smc, :lmc or :mw
function download_brott_data(folder_name; type=:all)
    ftp = FTP("ftp://cdsarc.cds.unistra.fr/ftp/J/A+A/530/A115/evol/")
    files = readdir(ftp)

    if !isdir(folder_name)
        mkdir(folder_name)
    end
    for file in files
        if type==:smc
            if !occursin("smc", file)
                continue
            end
        elseif type==:lmc
            if !occursin("lmc", file)
                continue
            end
        elseif type==:mw
            if !occursin("mw", file)
                continue
            end
        elseif type!=:all
            error("type must be either :all, :smc, :lmc or :mw, got $(type) instead")
        end
        download(ftp, file, "$(folder_name)/$(file)")
    end
end

function brott_path_constructor(strings::Vector{String}, folder_name, type)
    # this tells me how to construct a path to the data from the simulation input
    return "$(folder_name)/f$(strings[1])-$(strings[2]).$type.track.dat"
end

function brott_dataframe_loader(path)
    df = DataFrame(CSV.File(path, delim="|", header=[
        :star_age, :mass, :Teff, :logL, :radius, :log_mdot, :logg, :vsurf, :Prot, :vcrit, :Ge,
        :ab_H, :ab_He, :ab_Li, :ab_Be, :ab_B, :ab_C, :ab_N, :ab_O, :ab_F, :ab_Ne, :ab_Na, :ab_Mg, :ab_Al, :ab_Si, :ab_Fe,
        :s_h1, :s_he3, :s_he4, :s_li7, :s_be9, :s_b10, :s_b11, :s_c12, :s_c13, :s_n14, :s_n15, :s_o16, :s_o17, :s_o18,
        :s_f19, :s_ne20, :s_ne21, :s_ne22, :s_na23, :s_mg24, :s_mg25, :s_mg26, :s_al26, :s_al27, :s_si28, :s_si29, :s_si30, :s_fe56,
        :c_h1, :c_he3, :c_he4, :c_li7, :c_be9, :c_b10, :c_b11, :c_c12, :c_c13, :c_n14, :c_n15, :c_o16, :c_o17, :c_o18,
        :c_f19, :c_ne20, :c_ne21, :c_ne22, :c_na23, :c_mg24, :c_mg25, :c_mg26, :c_al26, :c_al27, :c_si28, :c_si29, :c_si30, :c_fe56,
        ]))
    df[!,:logTeff] = log10.(df.Teff)
    return df
end

function load_brott_data(folder_name, type; Xc_TAMS=1e-6)
    # linear interpolation works better with logM
    input_names = [:logM, :vrot]
    file_names = readdir(folder_name)
    # remove everything that does not match the desired set
    if type==:smc
        filter = "smc"
    elseif type==:lmc
        filter = "lmc"
    elseif type==:mw
        filter = "gal"
    else
        error("type must be either :smc, :lmc or :mw, got $(type) instead")
    end
    subset_file_names = []
    for file in file_names
        if occursin(r"f*"*"$(filter).track.dat", file)
            push!(subset_file_names, file)
        end
    end
    inputs = Matrix{String}(undef,2, length(subset_file_names))
    input_values = zeros(2, length(subset_file_names))

    for (i,file) in enumerate(subset_file_names)
        mass_str, rest = split(file,"-")
        mass_str = mass_str[2:end]
        logM = log10(parse(Float64, mass_str))
        vrot_str = split(rest,".")[1]
        vrot = parse(Float64, vrot_str)
        inputs[1,i] = mass_str
        inputs[2,i] = vrot_str
        input_values[1,i] = logM
        input_values[2,i] = vrot
    end

    path_constructor = x->brott_path_constructor(x, folder_name, filter)
    my_EEP_and_distance_calculator! = (df)->compute_distance_and_EEPs!(df; Xc_TAMS=Xc_TAMS)
    return StellarModelSet(inputs,input_names, path_constructor, brott_dataframe_loader, my_EEP_and_distance_calculator!; input_values=input_values);
end