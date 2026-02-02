function particle_in_box(x_0, y_0, vx_0, vy_0; max_steps=5_000, max_dt=0.1, min_dt =1e-4,friction_coeff=0.5, stop_velocity=0.1)
    x = zeros(max_steps)
    y = zeros(max_steps)
    vx = zeros(max_steps)
    vy = zeros(max_steps)
    time = zeros(max_steps)
    bounces = zeros(Int, max_steps)

    x[1] = x_0
    y[1] = y_0
    vx[1] = vx_0
    vy[1] = vy_0

    stop_index = max_steps
    for i in 2:max_steps
        dt = max_dt

        #find time it would take to bounce of edges
        if vx[i-1]>0
            dt_bouncex = (1-x[i-1])/vx[i-1]
        else
            dt_bouncex = -x[i-1]/vx[i-1]
        end
        if vy[i-1]>0
            dt_bouncey = (1-y[i-1])/vy[i-1]
        else
            dt_bouncey = -y[i-1]/vy[i-1]
        end
        #find time it will take to stop
        #simple friction proportional to v^2
        velocity = sqrt(vx[i-1]^2 + vy[i-1]^2)
        # acceleration is -\vec{v}/|v|*(|v|^2*friction_coeff)
        friction = velocity^2*friction_coeff
        dt_stop = (velocity-stop_velocity)/friction
        dt_stop = max(dt_stop, min_dt)

        bouncex = false
        bouncey = false
        stop = false
        if dt_bouncex <= dt_bouncey && dt_bouncex <= dt_stop && dt_bouncex <= dt
            bouncex = true
            dt = dt_bouncex
        elseif dt_bouncey <= dt_bouncex && dt_bouncey <= dt_stop && dt_bouncey <= dt
            bouncey = true
            dt = dt_bouncey
        elseif dt_stop <= dt_bouncex && dt_stop <= dt_bouncey && dt_stop <= dt
            stop = true
            dt = dt_stop
        end

        time[i] = time[i-1] + dt

        if bouncex || bouncey
            bounces[i] = bounces[i-1] + 1
        else
            bounces[i] = bounces[i-1]
        end

        # update position
        x[i] = x[i-1] + dt*vx[i-1]
        y[i] = y[i-1] + dt*vy[i-1]

        # update velocity
        friction_x = -vx[i-1]*velocity*friction_coeff
        friction_y = -vy[i-1]*velocity*friction_coeff
        if !stop
            vx[i] = vx[i-1] + friction_x*dt
            vy[i] = vy[i-1] + friction_y*dt
        else
            vx[i] = stop_velocity*sign(vx[i])
            vy[i] = stop_velocity*sign(vy[i])
            stop_index=i
            break
        end

        if bouncex
            vx[i] = -vx[i]
        elseif bouncey
            vy[i] = -vy[i]
        end
    end

    #truncate vectors to stop point if needed
    if stop_index < max_steps
        x = x[1:stop_index]
        y = y[1:stop_index]
        vx = vx[1:stop_index]
        vy = vy[1:stop_index]
        time = time[1:stop_index]
        bounces = bounces[1:stop_index]
    end

    return x,y,vx,vy,time, bounces
end

##

using StarStats, DataFrames

function ball_path_constructor(strings::Vector{String})
    try
        mkdir("sims")
    catch e
    end
    mkdir("sims/$(strings[1])_$(strings[2])")
    return "sims/$(strings[1])_$(strings[2])"
end

function ball_dataframe_loader(path)
    #extract parameters from path and run simulation
    vals = split(path[6:end],'_') # remove 'sims' from name
    vx_0 = parse(Float64, vals[1])
    vy_0 = parse(Float64,vals[2])
    x_0 = 0.5
    y_0 = 0.5
    x,y,vx,vy,time, bounces = particle_in_box(x_0, y_0, vx_0, vy_0)
    return DataFrame(xpos=x,ypos=y,vx=vx,vy=vy,time=time,bounces=bounces) # use xpos since x is reserved for the interpolating variable
end

# EEP loader, each bounce is an EEP in addition to start and stopping points
function get_EEPs_ball(track)
    return get_EEPs_ball_internal(track.bounces)
end
function get_EEPs_ball_internal(bounces)
    if bounces[end] == bounces[end-1]
        EEPs = zeros(Int,maximum(bounces)+2) # we need to have two more EEPs for start and finish
    else
        #unlikely case where simulation finishes at the same time there is a bounce
        EEPs = zeros(Int,maximum(bounces)+1)
    end
    names = Vector{Symbol}(undef, length(EEPs))
    EEPs[1] = 1
    names[1] = :start
    num_EEP = 2
    for i in 2:length(bounces)
        if bounces[i] > bounces[i-1]
            EEPs[num_EEP] = i
            names[num_EEP] = :bounce
            num_EEP += 1
        end
    end
    if bounces[end] == bounces[end-1]
        EEPs[end] = length(bounces)
        names[end] = :stop
    end
    return EEPs, names
end
function compute_distance_and_EEPs_ball!(df::DataFrame)

    EEPs, EEP_names = get_EEPs_ball(df)

    distance = zeros(size(df,1))
    delta_x = 0
    delta_y = 0

    for i in 2:size(df,1)
        delta_x = df.xpos[i]-df.xpos[i-1]
        delta_y = df.ypos[i]-df.ypos[i-1]

        distance[i] = distance[i-1] + sqrt(delta_x^2 + delta_y^2)
    end
    df[!,:distance] = distance

    x = zeros(size(df,1))
    for j in 1:(length(EEPs)-1)
        start_EEP = EEPs[j]
        end_EEP = EEPs[j+1]
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
        dtdx[i] = (df.time[i]-df.time[i-1])/(df.x[i]-df.x[i-1])
    end
    df[!,:dtdx] = dtdx

    return EEPs, EEP_names
end



struct BallMetric <: StarStats.AbstractMetric
end

function StarStats.dist(sim1::StarStats.SimulationData, sim2::StarStats.SimulationData, index1, index2, metric::BallMetric)
    x1 = sim1.df.xpos[index1]
    x2 = sim2.df.xpos[index2]
    y1 = sim1.df.ypos[index1]
    y2 = sim2.df.ypos[index2]

    return sqrt((x1-x2)^2 + (y1-y2)^2)
end
##

function remove_duplicate_pairs_simple(mat; eps=1e-10)
    unique_cols = []
    
    for i in 1:size(mat, 2)
        col = mat[:, i]
        is_duplicate = false
        
        for existing in unique_cols
            if all(abs.(col - existing) .< eps)
                is_duplicate = true
                break
            end
        end
        
        if !is_duplicate
            push!(unique_cols, col)
        end
    end
    
    return hcat(unique_cols...)
end

function make_inputs(input_values,suggest_ref_failed_ball,suggest_extra_ref_ball)
    n_failed = length(suggest_ref_failed_ball[:vx_i])
    n_extra = length(suggest_extra_ref_ball[:vx_i])
    nsims = n_failed + n_extra

    new_data = zeros(2, nsims)

    for i in 1:n_failed
        new_data[1, i] = suggest_ref_failed_ball[:vx_i][i]
        new_data[2, i] = suggest_ref_failed_ball[:vy_i][i]
    end

    for k in 1:n_extra
        col_index = n_failed + k
        new_data[1, col_index] = suggest_extra_ref_ball[:vx_i][k]
        new_data[2, col_index] = suggest_extra_ref_ball[:vy_i][k]
    end
    input_values_total = hcat(input_values, new_data)
    input_values_new = remove_duplicate_pairs_simple(input_values_total; eps=1e-10)
    return input_values_new
end

function make_random_inputs(input_values_old,input_values_new)
    
    number_to_append =  size(input_values_new,2) - size(input_values_old,2)
    new_data = rand(2,number_to_append)
    input_values_total = hcat(input_values_old, new_data)
    input_values_new = remove_duplicate_pairs_simple(input_values_total; eps=1e-10)
    return input_values_new
end

#makes new  modelset with new data 
function make_new_modelset(input_values)
    nsims = length(input_values[1,:])
    inputs = Matrix{String}(undef,2,nsims)
    input_names = [:vx_i, :vy_i]
    for i in 1:nsims
        inputs[1,i] = string(input_values[1,i])
        inputs[2,i] = string(input_values[2,i])
    end
    model_set_new = StellarModelSet(inputs,input_names, ball_path_constructor, ball_dataframe_loader, compute_distance_and_EEPs_ball!, BallMetric());
    rm("sims", recursive=true)

    return model_set_new
end


#makes new data and returns dataframes with points for refinment and colors for the plot
function make_new_data(model_set_new)
    
    suggest_ref_failed_new= StarStats.refine_model_set_failed_interpolation(model_set_new)
 
    df_new,evol_colors_new= array_of_colors(model_set_new)
   
    suggest_extra_ref_new = suggested_refinment_after_IQR(model_set_new)

    return suggest_ref_failed_new, df_new,evol_colors_new, suggest_extra_ref_new

end

#probably will work only in 2d
function calculate_simplex_area(model_set)
    good_area = 0
    bad_area = 0
    sides_of_simplex = zeros(Float64,0)
    for i in eachindex(model_set.simplex_interpolant.simplexes)
        simplex = model_set.simplex_interpolant.simplexes[i]
        for k in 1:length(simplex.point_indeces)
            for m in k+1:length(simplex.point_indeces)
                k_th_model = simplex.point_indeces[k]
                m_th_model = simplex.point_indeces[m]
                dist = 0
                for l in 1: length(model_set.models[k_th_model].input_params)
                    dist = dist+ (parse(Float64,model_set.models[k_th_model].input_params[l])-parse(Float64,model_set.models[m_th_model].input_params[l]))^2
                end
                append!(sides_of_simplex, sqrt(dist))
            end
        end
        #we assume simplex has 3 sides?
        p = sum(sides_of_simplex)/2
        area = p*(p-sides_of_simplex[1])*(p-sides_of_simplex[2])*(p-sides_of_simplex[3])  
        if area > 0
            area = sqrt(area)
               
            if model_set.can_interpolate_simplex[i] == 1
                good_area = good_area+area
            else
                bad_area = bad_area+area
            end
        end
        sides_of_simplex = zeros(Float64,0)  

    end
    return good_area, bad_area
end



##
good_area_arr = zeros(Float64,0)
bad_area_arr = zeros(Float64,0)
good_area_arr_evol = zeros(Float64,0)
bad_area_arr_evol = zeros(Float64,0)
good_area_arr_random = zeros(Float64,0)
bad_area_arr_random = zeros(Float64,0)
number_of_input_values = zeros(0)
number_of_input_values_evol = zeros(0)
n_sims_max = 5
max_res_level = 3 #+1 for 0 element
for sim_number in 1:n_sims_max
    nsims = 100000
    # sample vx and vy, keep initial x and y fixed at 0.5
    input_values = rand(2,nsims)
    input_values = hcat(input_values, [0.00001, 0.00001])
    input_values = hcat(input_values, [0.00001, 0.99999])
    input_values = hcat(input_values, [0.99999, 0.00001])
    input_values = hcat(input_values, [0.99999, 0.99999])   
    
    inputs = Matrix{String}(undef,2,nsims+4)
    input_names = [:vx_i, :vy_i]
    for i in eachindex(input_values)
        inputs[i] = string(input_values[i])
    end
    model_set_ball = StellarModelSet(inputs,input_names, ball_path_constructor, ball_dataframe_loader, compute_distance_and_EEPs_ball!, BallMetric());
    rm("sims", recursive=true)
    suggest_ref_failed_ball = StarStats.refine_model_set_failed_interpolation(model_set_ball)
    suggest_ref_failed_diff_evol_ball = StarStats.refine_model_set_failed_interpolation_diff_evol(model_set_ball)
    df_ball,evol_colors_ball = array_of_colors(model_set_ball)
    suggest_extra_ref_ball = suggested_refinment_after_IQR(model_set_ball)

    #fig = Figure()
    #ax = Axis(fig[1,1], xlabel="v1", ylabel="v2")
    #visualize_simplex_interpolant(ax, model_set_ball, model_set_ball.simplex_interpolant, suggest_ref_failed_diff_evol_ball,suggest_extra_ref_ball,df_ball )
    #save("sim_pics/testing_simpl_evol_level0_sim_number_1.png", fig)
    #fig

    suggest_ref_failed_old = suggest_ref_failed_ball
    suggest_extra_ref_old = suggest_extra_ref_ball
    input_values_old = input_values
    input_values_old_random = input_values

    suggest_ref_failed_old_evol = suggest_ref_failed_diff_evol_ball
    suggest_extra_ref_old_evol = suggest_extra_ref_ball
    input_values_old_evol = input_values
    
    
    good_area, bad_area = calculate_simplex_area(model_set_ball)
    append!(good_area_arr, good_area)
    append!(bad_area_arr, bad_area)

    append!(good_area_arr_evol, good_area)
    append!(bad_area_arr_evol, bad_area)
    #to random we as well put same area for start
    append!(good_area_arr_random, good_area)
    append!(bad_area_arr_random, bad_area)
    append!(number_of_input_values,length(input_values_old[1,:]))
    append!(number_of_input_values_evol,length(input_values_old_evol[1,:]))

    for i in 1:max_res_level

        input_values_new = make_inputs(input_values_old,suggest_ref_failed_old,suggest_extra_ref_old)
        model_set_new = make_new_modelset(input_values_new)
        good_area, bad_area = calculate_simplex_area(model_set_new)
        
        input_values_new_evol = make_inputs(input_values_old_evol,suggest_ref_failed_old_evol,suggest_extra_ref_old_evol)
        model_set_new_evol = make_new_modelset(input_values_new_evol)
        good_area_evol, bad_area_evol = calculate_simplex_area(model_set_new_evol)

        input_values_new_random = make_random_inputs(input_values_old_random, input_values_new) #never change input_values_old_random so the points will always be random
        model_set_random = make_new_modelset(input_values_new_random)
        good_area_random, bad_area_random = calculate_simplex_area(model_set_random)

        append!(number_of_input_values,length(input_values_new[1,:]))
        append!(good_area_arr, good_area)
        append!(bad_area_arr, bad_area)

        append!(number_of_input_values_evol,length(input_values_new_evol[1,:]))
        append!(good_area_arr_evol, good_area_evol)
        append!(bad_area_arr_evol, bad_area_evol)

        append!(good_area_arr_random, good_area_random)
        append!(bad_area_arr_random, bad_area_random)

        suggest_ref_failed_new, df_new,evol_colors_new, suggest_extra_ref_new = make_new_data(model_set_new)
        suggest_ref_failed_new_evol, df_new_evol,evol_colors_new_evol, suggest_extra_ref_new_evol = make_new_data(model_set_new_evol)
        #fig = Figure()
        #ax = Axis(fig[1,1], xlabel="v1", ylabel="v2")
        #visualize_simplex_interpolant(ax, model_set_new_evol, model_set_new_evol.simplex_interpolant, suggest_ref_failed_new_evol,suggest_extra_ref_new_evol,df_new_evol )
        #name = "sim_pics/testing_simpl_evol_level($i)_simnumber_$sim_number.png"
        #save(name, fig)
        #fig 
        #fig = Figure()
        #ax = Axis(fig[1,1], xlabel="v1", ylabel="v2")
        #visualize_simplex_interpolant(ax, model_set_new, model_set_new.simplex_interpolant, suggest_ref_failed_new,suggest_extra_ref_new,df_new )
        #name = "sim_pics/testing_simpl_level($i)_simnumber_$sim_number.png"
        #save(name, fig)
        #fig
        suggest_ref_failed_old = suggest_ref_failed_new
        suggest_extra_ref_old = suggest_extra_ref_new
        input_values_old = input_values_new
        
        suggest_ref_failed_old_evol = suggest_ref_failed_new_evol
        suggest_extra_ref_old_evol = suggest_extra_ref_new_evol
        input_values_old_evol = input_values_new_evol

        input_values_old_random = input_values_new
        
    end
   
end 
##
# MAKE LOG DATA
number_of_input_values = log10.(number_of_input_values)
number_of_input_values_evol = log10.(number_of_input_values_evol)
bad_area_arr = log10.(bad_area_arr)
bad_area_arr_random = log10.(bad_area_arr_random)
bad_area_arr_evol = log10.(bad_area_arr_evol)  
## 
fig = Figure()
    
    ax = Axis(fig[1,1], xlabel="Nsims", ylabel="area")
    for i in 1:n_sims_max
        #lines!(good_area_arr[(max_res_level+1)*i-(max_res_level+1)+1:(max_res_level+1)*i], number_of_input_values[(max_res_level+1)*i-(max_res_level+1)+1:(max_res_level+1)*i], color="green")
        lines!( number_of_input_values[(max_res_level+1)*i-(max_res_level+1)+1:(max_res_level+1)*i], bad_area_arr[(max_res_level+1)*i-(max_res_level+1)+1:(max_res_level+1)*i],  color="red")
        #lines!(good_area_arr_random[(max_res_level+1)*i-(max_res_level+1)+1:(max_res_level+1)*i], number_of_input_values[(max_res_level+1)*i-(max_res_level+1)+1:(max_res_level+1)*i],  color="blue")
        lines!( number_of_input_values[(max_res_level+1)*i-(max_res_level+1)+1:(max_res_level+1)*i], bad_area_arr_random[(max_res_level+1)*i-(max_res_level+1)+1:(max_res_level+1)*i], color="magenta")
        lines!( number_of_input_values_evol[(max_res_level+1)*i-(max_res_level+1)+1:(max_res_level+1)*i], bad_area_arr_evol[(max_res_level+1)*i-(max_res_level+1)+1:(max_res_level+1)*i], color="green")
    end

    labels = ["random dots", "middle of longest edge", "diff evol"]
    my_colors = ["magenta", "red", "green"]


    leg = Legend(fig[1, 2],
        [PolyElement(color=color) for color in my_colors],
        labels,
        valign=:center,  # вертикальное выравнивание по центру
        patchsize=(30, 20),
        labelsize=12
    )

fig
name = "sim_pics/number_of_runs_$(n_sims_max)_numer_ref_$(max_res_level+1)_new.png"
save(name, fig)
##

using BenchmarkTools


#@profview for i in 1:1000
#    suggested_refinment_after_IQR(model_set_ball)
#end
@code_warntype filter_data_new!(model_set_ball)

#@benchmark filter_data_new!($model_set_ball)
