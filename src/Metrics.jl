# abstract metric
abstract type AbstractMetric end

struct EuclideanMetric <: AbstractMetric
end

function dist(sim1::SimulationData, sim2::SimulationData, index1, index2, metric::EuclideanMetric)
    logTeff1 = sim1.df.logTeff[index1]
    logTeff2 = sim2.df.logTeff[index2]
    logL1 = sim1.df.logL[index1]
    logL2 = sim2.df.logL[index2]

    return sqrt((logTeff1-logTeff2)^2 + (logL1-logL2)^2)
end

function distance_along_curve!(model::SimulationData, metric::METRIC) where{METRIC<:AbstractMetric}
    distance = zeros(size(model.df,1))

    for i in 2:size(model.df,1)
        distance[i] = distance[i-1] + dist(model, model, i, i-1, metric)
    end
    model.df[!,:distance] = distance

    x = zeros(size(model.df,1))
    for j in 1:(length(model.EEPs)-1)
        start_EEP = model.EEPs[j]
        end_EEP = model.EEPs[j+1]
        if start_EEP==0 || end_EEP==0
            break
        end
        start_distance = model.df.distance[start_EEP]
        end_distance = model.df.distance[end_EEP]
        for i in start_EEP:end_EEP
            x[i] = (model.df.distance[i]-start_distance)/(end_distance-start_distance) + (j-1)
        end
    end
    for i in size(model.df,1):-1:1
        if x[i] > 0
            x[i:end] .= x[i]
            break
        end
    end
    model.df[!,:x] = x   
end
