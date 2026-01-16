
struct EuclideanMetric 
end

function dist(point1, point2, metric::EuclideanMetric)
    return sqrt((point1[1]-point2[1])^2 + (point1[2]-point2[2])^2)
end

struct PointPair2{TMETRIC}
    point1::Tuple{Float64,Float64}
    point2::Tuple{Float64,Float64}
    metric::TMETRIC
end

function distance_point_pair(pp::PointPair2)
    return dist(pp.point1,pp.point2,pp.metric)
end


##
function make_distance_along_curve(model, metric)
    
    distance = zeros(size(model.df,1))

    #HERE WE CALCULAE METRICS ARBITRARY 
    for i in 2:size(model.df,1)
        ##########ВОЗМОЖНО БАГ В РАСПОЛОЖЕНИИИТОЧЕК!!!! ПРОВЕРИТЬ!!!!
        point1 = ( model.df.logTeff[i], model.df.logL[i])
        point2 = ( model.df.logTeff[i-1],  model.df.logL[i-1])
        pp = PointPair2(point1, point2, metric)
        delta_logT_L = distance_point_pair(pp)

        distance[i] = distance[i-1] + delta_logT_L

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

