using StarStats
using CairoMakie
using Random
using CSV
using DataFrames

# simulate a cannon shot in a 2d plane, starting at an initial height y_i, an initial velocity v_i and an angle theta (in radians)
# theta=0 means a horizontal shot, theta=pi/2 is a vertical one
function simulate_cannon(h_i,v_i,theta_i; sim_time=20, steps = 1000)
    t = collect(LinRange(0,sim_time,steps)) # collect turns a range into a vector
    dt = sim_time/(steps-1)

    x = zeros(steps)
    y = zeros(steps)
    vx = zeros(steps)
    vy = zeros(steps)

    # set initial conditions
    x[1] = 0
    y[1] = h_i
    vx[1] = cos(theta_i)*v_i
    vy[1] = sin(theta_i)*v_i

    g = 10 # gravity in m/s

    #integrate motion with the Euler method
    for i in 2:steps
        x[i] = x[i-1] + dt*vx[i-1]
        y[i] = y[i-1] + dt*vy[i-1]
        vx[i] = vx[i-1]
        vy[i] = vy[i-1] - g*dt
    end
    return t,x,y,vx,vy
end

##
#test function
t,x,y,vx,vy = simulate_cannon(100,100,pi/5)
f = Figure()
ax = Axis(f[1,1])
lines!(ax,x,y)
f

##
# construct a grid with random initial conditions
Random.seed!(10) # fix random seed to get the same result next
# we have 3 initial parameters (h_i, v_i, theta_i), the following
# array will have in each column the value of each of these for one
# simulation and we'll include nsim simulations
nsim = 1000
points = rand(3,nsim)
# first parameter will be h_i, we'll make it go from 0 to 100 meters
points[1,:] = points[1,:].*100
# next we'll have v_i which will go from 0 to 100 m/s
points[2,:] = points[2,:].*100
# and then we have theta_i, which will go from zero to pi/2
points[3,:] = points[3,:].*(pi/2)

if isdir("cannon_simulations")
    rm("cannon_simulations", recursive=true)
end
mkdir("cannon_simulations")
#run sims and store data
for i in 1:nsim
    h_i = points[1,i]
    v_i = points[2,i]
    theta_i = points[3,i]
    t,x,y,vx,vy = simulate_cannon(h_i,v_i,theta_i)
    #create and save DataFrame
    df = DataFrame(t=t, xpos=x, ypos=y, vx=vx, vy=vy) # I call the column xpos and ypos because "x" is reserved for the interpolating variable in StarStats
    filename = "$(h_i)_$(v_i)_$(theta_i).txt"
    CSV.write("cannon_simulations/$filename", df)
end

##
# Now we want to read the data
file_names = readdir("cannon_simulations")
inputs = Matrix{String}(undef,3, length(file_names))
for (i,file_name) in enumerate(file_names)
    params = split(file_name[1:end-4],"_") # the end-4 is to remove .txt
    inputs[:,i] .= params
end
# create the basics needed to load a model set
path_constructor = x->"cannon_simulations/$(x[1])_$(x[2])_$(x[3]).txt"

function cannon_dataframe_loader(path)
    return DataFrame(CSV.File(path))
end

function cannon_distance_and_EEPs!(df::DataFrame)
    # we are just gonna put an EEP at the start and one at the end
    y = df[!,:ypos]
    x = df[!,:xpos]
    npoints = length(y)
    EEPs = [1, npoints]

    # as a metric, I'm gonna use distance traveled
    distance = zeros(npoints)
    for i in 2:npoints
        dx = x[i] - x[i-1]
        dy = y[i] - y[i-1]
        distance[i] = distance[i-1] + sqrt(dx^2 + dy^2)
    end
    df[!,:distance] = distance # add this column to the DataFrame

    # we now create the interpolating variable. This will go from zero to one
    x = distance./distance[end]
    df[!,:x] = x # add this column to the DataFrame

    return EEPs
end
# create the model set
model_set = StellarModelSet(inputs,[:h_i, :v_i, :theta_i], path_constructor, cannon_dataframe_loader, cannon_distance_and_EEPs!);

##
# test interpolation
xvals = LinRange(0,1,1000)
h_i = 50
v_i = 80
theta_i = pi/4
# get interpolated values for x and y
x_interp = interpolate_grid_quantity.(Ref(model_set),Ref([h_i, v_i, theta_i]),:xpos, xvals)
y_interp = interpolate_grid_quantity.(Ref(model_set),Ref([h_i, v_i, theta_i]),:ypos, xvals)
#plot them together with a simulation using those exact inputs
f = Figure()
ax = Axis(f[1,1])
lines!(ax, x_interp, y_interp, linewidth=5)
t,x,y,vx,vy = simulate_cannon(h_i,v_i,theta_i)
lines!(ax,x,y,linewidth=5,linestyle=:dash)
f
