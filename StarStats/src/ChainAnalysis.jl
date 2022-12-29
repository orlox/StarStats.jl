using Makie
using StatsBase
using CairoMakie

export create_corner_plot, create_2D_density, create_1D_density

function create_corner_plot(chain_values,names,label_names, chain_weights, fractions, figure)
    ga = figure[1, 1] = GridLayout()
   
    num_col = length(chain_values)-1
    for i in 1:num_col
        for  j in i+1:num_col+1
            axis = Axis(ga[j,i],xtickalign=1,xtickcolor = :white,ytickalign=1,ytickcolor = :white, 
                   xlabel=label_names[i], ylabel=label_names[j] )
            
            create_2D_density(chain_values[names[i]], chain_values[names[j]], chain_weights,fractions, axis)
            if i>1
                hideydecorations!(axis, ticks=false, minorticks=false)
            end
            if j!=num_col+1
                hidexdecorations!(axis,ticks=false, minorticks=false)
            end         
        end  
   
    end 
    for i in 1:num_col+1
        axis = Axis(ga[i,i], xgridvisible = false, ygridvisible = false,xtickalign=1,
               xlabel=label_names[i])
        create_1D_density(chain_values[names[i]], chain_weights,fractions,axis)
        hideydecorations!(axis)
        if i !=num_col+1
            hidexdecorations!(axis,ticks=false, minorticks=false)
        end    
    end     
    rowgap!(ga,10)
    colgap!(ga,10)

    return figure 
end

function create_2D_density(values1,values2, chain_weights,fractions,axis)

    h = fit(Histogram,(values1,values2),weights(chain_weights),nbins=100)
    x= (h.edges[2][2:end].+h.edges[2][1:end-1])./2
    y =(h.edges[1][2:end].+h.edges[1][1:end-1])./2
    heatmap!(axis,y,x, h.weights)
    bounds = get_bounds_for_fractions(h,fractions)
   
    contour!(axis,y,x, h.weights, levels=bounds,color=:black,linewidth=2)

    return axis
end  

function get_bounds_for_fractions(h,fractions)
    integral = sum(h.weights)
    #smartely decide when the bisection is ended 
    bounds =zeros(length(fractions))
    for (j,fraction) in enumerate(fractions)
        minbound = 0
        maxbound = maximum(h.weights)
        newbound = 0
        for i in 1:15
            newbound = 0.5*(minbound+maxbound)
            integral2 = sum(h.weights[h.weights.>newbound])
            newfraction = integral2/integral
            if newfraction>fraction
                minbound = newbound
            else
                maxbound = newbound
            end

        end
        bounds[j] = newbound
        
    end 
    return bounds
end

function create_1D_density(values, chain_weigths,fractions,axis)
    h = fit(Histogram,(values),weights(chain_weights),nbins=100)
    x =(h.edges[1][2:end].+h.edges[1][1:end-1])./2
    lines!(axis,x, h.weights)

    return axis
   
end