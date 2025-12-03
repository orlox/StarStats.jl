using Pkg
using DataFrames
using CSV
using Plots
using LinearAlgebra
using Interpolations
using Base.Filesystem
using CairoMakie
using Statistics

# if there is TE => 0, if there is NO TE 1
function merge_tables(left_index_peak, right_index_peak, data_overflows)
    data_peaks =  DataFrame(
        left = left_index_peak,
        right = right_index_peak)
    outflow = [] #0 if No outflow
    TE_left = [] #1 peak, 0 no Peak
    TE_right = []
    index_merged_left = [] #all indexes
    index_merged_right = []
    for i in 1:length(data_peaks.left )
        for k in 1:Integer(length(data_overflows.index_overflow)/2)
            #случай когда выход и приход обратно в ТЕ равновесие вне outflow, 2k-1 левая граница
            if k != 1
                if (data_peaks.left[i]<data_overflows.index_overflow[2*k-1]) && 
                    (data_peaks.right[i]<data_overflows.index_overflow[2*k-1]) &&
                     (data_peaks.left[i] > data_overflows.index_overflow[2*k-2])
                     #левый левее левого и правый левее левого => нет пересечения
                    append!(index_merged_left, data_peaks.left[i])

                    append!(index_merged_right,data_peaks.right[i])

                    append!(outflow, 0)
                    append!(TE_left,1)
                    append!(TE_right,0)
                end
            else
                if (data_peaks.left[i]<data_overflows.index_overflow[2*k-1]) && 
                    (data_peaks.right[i]<data_overflows.index_overflow[2*k-1]) 
                     #левый левее левого и правый левее левого => нет пересечения
                    append!(index_merged_left, data_peaks.left[i])

                    append!(index_merged_right,data_peaks.right[i])

                    append!(outflow, 0)
                    append!(TE_left,1)
                    append!(TE_right,0)
                end
            end

            if i != 1
                    if (data_overflows.index_overflow[2*k-1] < data_peaks.left[i]) &&
                        (data_overflows.index_overflow[2*k] < data_peaks.right[i]) &&
                            (data_overflows.index_overflow[2*k-1] > data_peaks.right[i-1])
                            append!(index_merged_left, data_overflows.index_overflow[2*k-1])
                            append!(index_merged_right, data_overflows.index_overflow[2*k])
                            append!(outfflow, 1)
                            append!(TE_left,0)
                            append!(TE_right,0)
                    #левый правее последнего правого
                    elseif (data_overflows.index_overflow[2*k-1] > data_peaks.right[end])
                        append!(index_merged_left, data_overflows.index_overflow[2*k-1])
                        append!(index_merged_right, data_overflows.index_overflow[2*k])
                        append!(outfflow, 1)
                        append!(TE_left,0)
                        append!(TE_right,0)
                        break
                    end
            else
                    if (data_overflows.index_overflow[2*k-1] < data_peaks.left[i]) &&
                        (data_overflows.index_overflow[2*k] < data_peaks.right[i])
                            append!(index_merged_left, data_overflows.index_overflow[2*k-1])
                            append!(index_merged_right, data_overflows.index_overflow[2*k])
                            append!(outfflow, 1)
                            append!(TE_left,0)
                            append!(TE_right,0)
                    end
               
            end

            if (data_peaks.left[i] >data_overflows.index_overflow[2*k-1]) &&
                (data_peaks.right[i] < data_overflows.index_overflow[2*k]) #левый правее левого, правый левее правого (между outflow)
                #######от левого края outflow к левому неравновесию
                append!(index_merged_left,data_overflows.index_overflow[2*k-1])
                append!(index_merged_right, data_peaks.left[i])
                append!(outflow, 1)
                append!(TE_left,0)
                append!(TE_right,1)

                #######от левого края неравновесия к правому неравновесию
                append!(index_merged_left, data_peaks.left[i])
                append!(index_merged_right, data_peaks.right[i])
                append!(outflow,1)
                append!(TE_left,1)
                append!(TE_right,0)

                #######от правого неравновесия к правому outflow
                append!(index_merged_left, data_peaks.right[i])
                append!(index_merged_right, data_overflows.index_overflow[2*k])
                append!(outflow, 1)
                append!(TE_left,0)
                append!(TE_right,0)
            
            # левый край левее левого outflow, правый правее левого, но правй правее правого outglow
            elseif  (data_peaks.left[i] < data_overflows.index_overflow[2*k-1]) &&
                (data_peaks.right[i] > data_overflows.index_overflow[2*k-1]) &&
                (data_peaks.right[i] < data_overflows.index_overflow[2*k]) 

                #потеря ТЕ к началу overflow
                append!(index_merged_left, data_peaks.left[i] )
                append!(index_merged_right, data_overflows.index_overflow[2*k-1])
                append!( outflow, 0)
                append!(TE_left,1)
                append!(TE_right,1)

                #начало overflow и все еще нет ТЕ
                append!(index_merged_left, data_overflows.index_overflow[2*k-1] )
                append!(index_merged_right, data_peaks.right[i] )
                append!(outflow, 1)
                append!(TE_left,1)
                append!(TE_right,0)

                # начало ТЕ и все еще есть overflow
                append!(index_merged_left,data_peaks.right[i])
                append!(index_merged_right, data_overflows.index_overflow[2*k])
                append!(outflow, 1)
                append!(TE_left,0)
                append!(TE_right,0)
            
            #левый правее левого overflow, но левее правого, правый правее правого overflow
            elseif (data_peaks.left[i] > data_overflows.index_overflow[2*k-1]) &&
                (data_peaks.left[i] < data_overflows.index_overflow[2*k]) &&
                (data_peaks.right[i] > data_overflows.index_overflow[2*k]) 

                append!(index_merged_left,data_overflows.index_overflow[2*k-1])
                append!(index_merged_right,data_peaks.left[i])
                append!(outflow, 1)
                append!(TE_left,0)
                append!(TE_right,1)

                append!(index_merged_left,data_peaks.left[i])
                append!(index_merged_right,data_overflows.index_overflow[2*k]) 
                append!(outflow, 1)
                append!(TE_left,1)
                append!(TE_right,1)

                append!(index_merged_left,data_overflows.index_overflow[2*k])
                append!(index_merged_right,data_peaks.right[i])
                append!(outflow, 0)
                append!(TE_left,1)
                append!(TE_right,0)
            
            elseif (data_peaks.left[i] > data_overflows.index_overflow[end])
                append!(index_merged_left,data_peaks.left[i])
                append!(index_merged_right,data_peaks.right[i])
                append!(outflow, 0)
                append!(TE_left,1)
                append!(TE_right,0)
                break
            
            #outflow внутри неравновесия лежит
            elseif (data_overflows.index_overflow[2*k-1] > data_peaks.left[i]) &&
                    (data_overflows.index_overflow[2*k] < data_peaks.right[i])
                    append!(index_merged_left,data_peaks.left[i])
                    append!(index_merged_right,  data_overflows.index_overflow[2*k-1])
                    append!(outflow, 0)
                    append!(TE_left,1)
                    append!(TE_right,1)

                    append!(index_merged_left,  data_overflows.index_overflow[2*k-1] )
                    append!(index_merged_right,  data_overflows.index_overflow[2*k])
                    append!(outflow,1)
                    append!(TE_left,1)
                    append!(TE_right,1)

                    append!(index_merged_left,data_overflows.index_overflow[2*k])
                    append!(index_merged_right,data_peaks.right[i] )
                    append!(outflow, 0)
                    append!(TE_left,1)
                    append!(TE_right,0)
                
            end

        end
    end

    
    data_merged = DataFrame(
        index_left = index_merged_left,
        index_right =index_merged_right,
        overflow = outflow,
        th_eq_left = TE_left,
        th_eq_right = TE_right)
    
    return data_merged
end