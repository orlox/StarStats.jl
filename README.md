# StarStats.jl

This package is designed to read grids of stellar evolution models from different evolutionary codes and perform interpolation and Bayesian inference against observed systems.

The following code  block shows an example of loading a grid with three variable parameters

* masses: sampled logarithmically $\log M/M_{\odot}=$0.9-2.1 in steps of 0.025
* r

```julia
using StarStats
using Printf

function path_constructor(strings::Vector{String})
    DATA_FOLDER = ENV["STARSTATS_TEST_DATA_FOLDER"]
    return DATA_FOLDER*"/LMC/LMC_$(strings[1])_$(strings[2])_$(strings[3]).track.gz"
end
masses = [@sprintf("%.3f", x) for x in range(0.9,2.1,step=0.025)]
rotation = [@sprintf("%.2f", x) for x in range(0.0,0.9,step=0.1)]
overshoot = [@sprintf("%.2f
", x) for x in range(0.5,4.5,step=0.5)]

star_grid = ModelDataGrid([rotation,masses,overshoot],[:rotation,:logM,:overshoot])
load_grid(star_grid,path_constructor,gz_dataframe_loader_with_Teff_and_star_age_fix); 
compute_distances_and_EEPs(grid)
```

The code
