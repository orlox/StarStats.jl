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



    #=
    end_of_He_burning = 0
    for s in 1:length(h_donor.log_Teff)
        if (h_donor.center_he4[s] < 1e-4)
            append!(index_in_tables,s)
            push!(meaning_of_index_in_tables, "TACHeB")
            end_of_He_burning  = s
            break
        end
    end

    #this is a complicated way of finding end of RGBtip

    min_temp = minimum(h_donor.log_cntr_T[end_of_H_burning:end_of_He_burning])
    index = findfirst(x -> x == min_temp,h_donor.log_cntr_T[end_of_H_burning:end_of_He_burning])
    index = index+end_of_H_burning

  
    
    for s in 1:length(h_donor.log_cntr_T)
        if h_donor.center_he4[s] > h_donor.center_he4[index] -0.03
            global index_new = s
            break
        end
    end

    
    new_min_T = minimum(h_donor.log_cntr_T[index:length(h_donor.log_cntr_T)])
    index = findfirst(x -> x == new_min_T,h_donor.log_cntr_T[index_new:length(h_donor.log_cntr_T)])
    index = index+index_new

    append!(index_in_tables,index)
    push!(meaning_of_index_in_tables, "ZACHeB")
  
    =#

    return sort_vector(index_in_tables, meaning_of_index_in_tables)

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
        deleteat!(index_in_tables, 1)
        deleteat!(meaning_of_index_in_tables,1)
    end


    if index_in_tables[length(index_in_tables) - 1] == index_in_tables[length(index_in_tables)]
        deleteat!(index_in_tables, length(index_in_tables)-1)
        deleteat!(meaning_of_index_in_tables,length( meaning_of_index_in_tables)-1)
    end

    return index_in_tables, meaning_of_index_in_tables

end

