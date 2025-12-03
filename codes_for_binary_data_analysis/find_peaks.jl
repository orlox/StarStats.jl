using Pkg
using DataFrames
using CSV
using Plots
using LinearAlgebra
using Interpolations
using Base.Filesystem
using CairoMakie
using Statistics


function identify_pics_caseA(h_donor)
    dTh_eq = diff(log10.(h_donor.thermal_eq_difference))
    dt = diff(h_donor.star_age)
    dTh_dt = dTh_eq  ./ dt
    max_dTh_dt = maximum(dTh_dt )

    fig = Figure()
    ax = Axis(fig[1, 1])
    lines!(ax,h_donor.star_age, log10.(h_donor.thermal_eq_difference))
    CairoMakie.ylims!(ax,-5,5)  
    hlines!(ax, Statistics.median(log10.(h_donor.thermal_eq_difference)))

    ########IF OUTFLOW ENDS BEFORE THE PEAK IN THERMAL_NONEQ => ONLY OUTFLOW WILL BE MARKED
    arr_0 = []
    for m in 1:length(h_donor.rl_relative_overflow_1)
         if (h_donor.rl_relative_overflow_1[m] >0) .&&
            (h_donor.center_h1[m] > 1e-12)
            append!(arr_0,m)
         end
    end

    arr_1 = []

    for i in 1:length(dTh_dt)-1
        #if (h_donor.rl_relative_overflow_1[i] >0) .&& (h_donor.center_h1[i] > 1e-12)
            
        #    CairoMakie.scatter!(ax,h_donor.star_age[i], log10.(h_donor.thermal_eq_difference[i]))
             #println(abs(dTh_dt[i] > abs(max_dTh_dt)/100000), " ",h_donor.rl_relative_overflow_1[i] >0, " ", h_donor.center_h1[i] > 1e-12)
             #println( h_donor.center_h1[i])
        #end

        if (abs(dTh_dt[i]) > (abs(max_dTh_dt)/10000)) .&&
            (h_donor.rl_relative_overflow_1[i] >0) .&&
            (h_donor.center_h1[i] > 1e-12)
          
            CairoMakie.scatter!(ax,h_donor.star_age[i], log10.(h_donor.thermal_eq_difference[i]))
            append!(arr_1,i)
        end
    end
   
    if length(arr_1) == 0
        return  arr_0[1], arr_0[end]
    else
        return arr_1[1], arr_1[end]
    end
    fig 
    #return arr[1], arr[end]
    
end



####SMOOTHS THE DATA
function smooth_data(data::Vector{Float64}, window_size::Int)
    """Сглаживает данные скользящим средним"""
    n = length(data)
    smoothed = similar(data)
    
    for i in 1:n
        start_idx = max(1, i - window_size ÷ 2)
        end_idx = min(n, i + window_size ÷ 2)
        smoothed[i] = mean(data[start_idx:end_idx])
    end
    
    return smoothed
end

################# FINDS THE INDEX OF THE MAXIMUM and 2 INDEXES OF mins AROUND
function find_peacs_indexes(time::Vector{Float64}, th_eq_smoothed::Vector{Float64})

    length_data = length(time)
    mediana = Statistics.median(th_eq_smoothed)
    length_of_window =  Integer(round(length_data /100))
    number_of_windows = 100
    index_of_maximum = []
    for i in 1:Integer(number_of_windows)
        a = i*length_of_window
        b = (i+1)*length_of_window
        if (b > length_data)
            b = length_data
        else

            subset = th_eq_smoothed[a : b]
            max = maximum(subset)

            if (max == subset[end]) || (max < mediana)
                continue
            else
                index = findfirst(x -> x == max,th_eq_smoothed )
                append!(index_of_maximum, index)
            end
        end
    end

    min_in_data = [] # WE FIND WHERE THE DATA INTERCEPTS WITH LINE  mediana/2 
    for i in 1:length(th_eq_smoothed)-1
        if  (th_eq_smoothed[i] - mediana/2) >0  && (th_eq_smoothed[i+1] - mediana/2) <0
            append!(min_in_data, i)
        elseif (th_eq_smoothed[i] - mediana/2) <0  && (th_eq_smoothed[i+1] - mediana/2) >0
            append!(min_in_data, i)
        else
            continue
        end
    end

    max_in_data = [] # AND WE FIND MAXIMUM between the mins in data
    if (iseven(length(min_in_data)) == true )
        len = Integer(length(min_in_data)/2)
    else
        len = Integer((length(min_in_data)-1)/2)+1
    end
    
    for i in 1:len
        if (iseven(length(min_in_data)) == true )
            a = min_in_data[2*i-1]
            b = min_in_data[2*i]
            subset = th_eq_smoothed[a : b]
            max = maximum(subset)
            index = findfirst(x -> x == max,th_eq_smoothed )
            append!(max_in_data,index)
        else
            if 2*i-1 != length(min_in_data)
                a = min_in_data[2*i-1]
                b = min_in_data[2*i]
                subset = th_eq_smoothed[a : b]
                max = maximum(subset)
                index = findfirst(x -> x == max,th_eq_smoothed )
                append!(max_in_data,index)
            else
                max = th_eq_smoothed[end]
                index = findfirst(x -> x == max,th_eq_smoothed )
                append!(max_in_data,index)
                break
            end
        end

    end
    return max_in_data, min_in_data
end

#this function finds the interception with the plot and from these 2 points we will surch for the derivative
# START_Index is maximums in data
function find_interception_with_data_left_and_right(hight,start_index,th_eq_smoothed)
    
    for i in 1:length(th_eq_smoothed[1:start_index])
        if start_index - (i+1) == 0
            println("You are on the left edge there is no beggining of the peak")
            global index_left = 1
            break
        elseif (th_eq_smoothed[start_index - i] > hight ) && (th_eq_smoothed[start_index - (i+1)] < hight)
            global index_left = start_index - (i+1)
            break
        end
    end

    for i in start_index:length(th_eq_smoothed)
        if i - length(th_eq_smoothed) == 0
            println("You are on the right edge there is no end of the peak")
            global index_right = length(th_eq_smoothed)
            break
        elseif (th_eq_smoothed[i] > hight ) && (th_eq_smoothed[i+1] < hight)
            global index_right = i+1
            break
        end
    end


    return  index_left, index_right
end       


function find_edges_of_the_peak(max_in_data, min_in_data, th_eq_smoothed, time)

    edges_left = []
    edges_right = []
    mins = []
    minimum_in_data = minimum(th_eq_smoothed)
    ######IN THIS BLOCK WE FIND EDGES ON 1/3 OF THE HIGHT OF THE PEAK
    ######WE DO THIS BECAUSE FROM THA POINT THE DERIVATIVE CAN BE USED AS A CRITERIA

    if (iseven(length(min_in_data)) == false)
        mins = min_in_data[1:end-1] 
        for i in 1:Integer(length(mins)/2)
            hight = minimum_in_data + (th_eq_smoothed[max_in_data[i]] - minimum_in_data)/ 3 #hight of pick / 10
            current_peak_index = max_in_data[i]
            left, right = find_interception_with_data_left_and_right(hight,current_peak_index,th_eq_smoothed)
            append!(edges_left,left)
            append!(edges_right, right)
        end

        hight = minimum_in_data + (th_eq_smoothed[max_in_data[end]] - minimum_in_data)/ 3 #hight of pick / 10
        current_peak_index = max_in_data[end]
        left, right = find_interception_with_data_left_and_right(hight,current_peak_index,th_eq_smoothed)
        append!(edges_left,left)
        append!(edges_right, right)

    else
        mins = min_in_data
        for i in 1:Integer(length(mins)/2)
            
            hight = minimum_in_data + (th_eq_smoothed[max_in_data[i]] - minimum_in_data)/ 3 #hight of pick / 10
            
            current_peak_index = max_in_data[i]
            left, right = find_interception_with_data_left_and_right(hight,current_peak_index,th_eq_smoothed)
            
            append!(edges_left,left)
            append!(edges_right, right)
        end

        
    end
    #############NEWXT WE CALCULATE WITH THE EDGES THE DERIVATIVE CRITERIA
    result_left = []
    result_right = []
    right_index_stop = 1
    for i in 1:length(edges_left)
        lef_index_peak, right_index_peak = make_derivative_criteria(edges_left[i], edges_right[i], th_eq_smoothed, time,right_index_stop)
        append!(result_left, left_index_peak)
        append!(result_right, right_index_peak)
        right_index_stop = right_index_peak
    end

    return result_left, result_right 
end

function make_derivative_criteria(left_index, right_index, th_eq_smoothed, time, right_index_stop)
    dTh_eq = diff(th_eq_smoothed)
    dt = diff(time)
    dTh_dt = dTh_eq  ./ dt 

    if left_index == 1
        println("This is the beggining of the peak")
        global left_index_peak = 1
    else
        for i in 1:left_index #we go to the left from the left point
            if (abs(dTh_dt[left_index - i]) <= abs(dTh_dt[left_index]/50))
                if left_index - i < right_index_stop
                    global  left_index_peak= right_index_stop
                    break
                else
                    global left_index_peak = left_index - i
                    break
                end
            end
        end   
    end

    if right_index == length(th_eq_smoothed)
        println("This is the end of the peak")
        global right_index_peak = length(th_eq_smoothed)
    else
        for i in right_index: length(th_eq_smoothed) # now we go to the right from the right point
            if (abs(dTh_dt[i]) <= abs(dTh_dt[right_index]/50))
                global right_index_peak = i
                break
            end
        end
    end

    return  left_index_peak, right_index_peak
end

