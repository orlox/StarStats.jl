using Pkg
using DataFrames
using CSV
using Plots
using LinearAlgebra
using Interpolations
using Base.Filesystem
using CairoMakie
using Statistics

function construct_EEPs_vector(left, right, df_overflows, h_donor)
    index_in_tables = []
    meaning_of_index_in_tables = []
    new_table = DataFrame(
        index = index_in_tables,
        meaning_of_index = meaning_of_index_in_tables
    )
    index_in_tables = []
    meaning_of_index_in_tables = []

    for i in 1:length(left)
        append!(index_in_tables, left[i])
        push!(meaning_of_index_in_tables, "TE_brocken")
        append!(index_in_tables, right[i])
        push!(meaning_of_index_in_tables, "TE_regained")
    end

    for k in 1:Integer(length(df_overflows.index_overflow)/2)
        append!(index_in_tables, df_overflows.index_overflow[2*k-1])
        push!(meaning_of_index_in_tables, "Start_overflow")
        append!(index_in_tables, df_overflows.index_overflow[2*k])
        push!(meaning_of_index_in_tables, "End_overflow")
    end

   # append!(index_in_tables,1 )
   # append!(meaning_of_index_in_tables, "ZAMS")
    append!(index_in_tables,length(h_donor.log_Teff))
    push!(meaning_of_index_in_tables, "End_simulation")

    
    #########Conditions below are taken from paper Dotter 2012
    for s in 1:length(h_donor.log_Teff)
        if  h_donor.log_LH[s] / h_donor.log_L[s] > 0.99
            append!(index_in_tables,s)
            push!(meaning_of_index_in_tables, "ZAMS")
            break
        end
        
    end
    
    end_of_H_burning = 0
    for s in 1:length(h_donor.log_Teff)
        if(h_donor.center_h1[s] < 1e-4)
            
            append!(index_in_tables,s)
            push!(meaning_of_index_in_tables, "TAMS")
            end_of_H_burning = s
            break
        end
        
    end

    RGBTip_EEP = 0
    TAMS = index_in_tables[length(index_in_tables)]

    for i in TAMS+1:length(h_donor.log_Teff)
        if h_donor.center_he4[i] < h_donor.center_he4[TAMS] - 0.01
            RGBTip_EEP = i
            break
        end
    end
    for i in RGBTip_EEP+1:length(h_donor.log_Teff)
        if h_donor.center_he4[i] < h_donor.center_he4[RGBTip_EEP] - 0.03
            ZACHeB_EEP = i
            append!(index_in_tables,ZACHeB_EEP)
            push!(meaning_of_index_in_tables, "ZACHeB")
            break
        end
    end
    for s in 1:length(h_donor.log_Teff)
        if (h_donor.center_he4[s] < 1e-4)
            
            append!(index_in_tables,s)
            push!(meaning_of_index_in_tables, "TACHeB")
            
           # end_of_He_burning  = s
            break
        end
    end
    indexes, meaning_of_indexes = sort_vector(index_in_tables, meaning_of_index_in_tables)
    return modify_vector(indexes, meaning_of_indexes )

end


function sort_vector(index_in_tables, meaning_of_index_in_tables)

    for i in 1: length(index_in_tables)
        for k in i+1:length(index_in_tables)
            if (index_in_tables[i] > index_in_tables[k])
                a = index_in_tables[i]
                b = index_in_tables[k]
                index_in_tables[i] = b
                index_in_tables[k] = a
                c = meaning_of_index_in_tables[i]
                d = meaning_of_index_in_tables[k]
                meaning_of_index_in_tables[i] =d
                meaning_of_index_in_tables[k] = c
            end
        end
    end
    
    if index_in_tables[1] == index_in_tables[2] 
        if (meaning_of_index_in_tables[1] == "ZAMS")
            deleteat!(index_in_tables, 2)
            deleteat!(meaning_of_index_in_tables,2)
        else
            deleteat!(index_in_tables, 1)
            deleteat!(meaning_of_index_in_tables,1)
        end
    end


    if index_in_tables[length(index_in_tables) - 1] == index_in_tables[length(index_in_tables)]
        if(meaning_of_index_in_tables[length( meaning_of_index_in_tables)] == "End_simulation")
            deleteat!(index_in_tables, length(index_in_tables)-1)
            deleteat!(meaning_of_index_in_tables,length( meaning_of_index_in_tables)-1)
        else
            deleteat!(index_in_tables, length(index_in_tables))
            deleteat!(meaning_of_index_in_tables,length( meaning_of_index_in_tables))
        end
    end

    return index_in_tables, meaning_of_index_in_tables

end

#This function modifies vector if there is a situation where the code places 2 same points (Te brocken and regained (sometimes))
function modify_vector(index_in_tables, meaning_of_index_in_tables)
    index_to_remove = []
    for i in 1:length(index_in_tables)-1
        if index_in_tables[i] == index_in_tables[i+1]
            if (meaning_of_index_in_tables[i] == "TE_brocken") .&& (meaning_of_index_in_tables[i+1] == "TE_regained")
                append!(index_to_remove, i)
            end
        end
    end
    
    for i in 1:length(index_to_remove)
        deleteat!(index_in_tables, index_to_remove[i])
        deleteat!(meaning_of_index_in_tables, index_to_remove[i])
    end
    return index_in_tables, meaning_of_index_in_tables
        
end

