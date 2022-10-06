
using QuadGK, Interpolations

export ModelDataGridLikelihood, marginalized_likelihood, credible_interval

mutable struct ModelDataGridLikelihood
    grid::ModelDataGrid
    likelihood::Array{Float64}
    observable_names::Vector{Symbol}
    observable_values::Vector{Float64}
    observable_errors::Vector{Float64}
    
    function ModelDataGridLikelihood(grid, observable_names, observable_values, observable_errors)
        likelihood = zeros(size(grid.dfs)...)
      
        for index in Base.product([1:length(input) for input in grid.inputs]...)
            if !isassigned(grid.dfs,index...)
                continue
            end
            for x in eachrow(grid.dfs[index...])
                step_likelihood = x.dt
                         
                for (i,name) in enumerate(observable_names)
                    step_likelihood *= exp(-((getproperty(x,name)-observable_values[i])/
                                                                  observable_errors[i])^2)
                end    
            likelihood[index...] += step_likelihood                                       
            end
        end
        new(grid, likelihood, observable_names, observable_values, observable_errors)
    end
end


function marginalized_likelihood(grid_likelihood::ModelDataGridLikelihood, non_marginalized_names::Vector{Symbol})
    sizes = []
    for name in non_marginalized_names
       index = findfirst(isequal(name),grid_likelihood.grid.input_names)  
       push!(sizes, length(grid_likelihood.grid.input_values[index]))
       
    end
    
    marginalized_likelihood = zeros(sizes...)
    dimensions_to_marginalize = []
    for (i, name) in enumerate(grid_likelihood.grid.input_names)
        if name ∉ non_marginalized_names
            push!(dimensions_to_marginalize, i)
        end     
    end 
    summed_likelihood = sum(grid_likelihood.likelihood,dims=dimensions_to_marginalize)
    summed_likelihood = dropdims(summed_likelihood, dims = (findall(size(summed_likelihood) .== 1)...,))
    permutation_array = []
    for name in non_marginalized_names
        index = 0
        for input_name in grid_likelihood.grid.input_names
            if input_name ∈ non_marginalized_names
                index += 1 
                if name == input_name
                     push!(permutation_array, index)
                     break 
                end    
            end     
        end 
    end
    summed_likelihood = permutedims(summed_likelihood, permutation_array)

    return summed_likelihood
end 

function credible_interval(grid_likelihood::ModelDataGridLikelihood, observable_name, 
        fraction, num_points)
    index = 0
    for (i, name) in enumerate(grid_likelihood.grid.input_names)
        if name == observable_name
            index = i
        end
    end
    minval = minimum(grid_likelihood.grid.input_values[index])
    maxval = maximum(grid_likelihood.grid.input_values[index])
    
    ml= marginalized_likelihood(grid_likelihood,[observable_name])
    itp = interpolate((grid_likelihood.grid.input_values[index],), ml, Gridded(Linear()))
    
    integral = quadgk(x -> itp(x), minval, maxval, rtol=1e-8)[1]
    
    xvalues = collect(LinRange(minval,maxval,num_points))
    yvalues = itp(xvalues)
    max_valy = findmax(yvalues)
    maxx = xvalues[max_valy[2]]
    maxy = max_valy[1]
    
    function itp2(x, itp, bound)
        value = itp(x)
        if value>bound
            return value
        else
            return 0
        end
    end
    minbound = 0
    maxbound = maxy
    newbound = 0
    for i in 1:50
        newbound = 0.5*(minbound+maxbound)
        integral2 = quadgk(x -> itp2(x,itp,newbound), minval, maxval, rtol=1e-8)[1]
        newfraction = integral2/integral
        if newfraction>fraction
            minbound = newbound
        else
            maxbound = newbound
        end
    end
    filtered_xvalues = xvalues[yvalues.>newbound]
    return ([maximum(filtered_xvalues), maxx, minimum(filtered_xvalues)], newbound)
end 
