# this only needs to be done once to get the data
using StarStats
StarStats.download_brott_data("brott_models"; type=:lmc)

##
# code below can be used after the data has been downloaded
using StarStats
##
model_set = StarStats.load_brott_data("brott_models", :lmc);

##
# show some HR diagrams from interpolation
using LaTeXStrings, CairoMakie, MathTeXEngine
# pretty theme
basic_theme = Theme(fonts=(regular=texfont(:text), bold=texfont(:bold),
                           italic=texfont(:italic), bold_italic=texfont(:bolditalic)),
                    fontsize=30, size=(1000, 750), linewidth=7,
                    Axis=(xlabelsize=40, ylabelsize=40, titlesize=40, xgridvisible=false, ygridvisible=false,
                          spinewidth=2.5, xminorticksvisible=true, yminorticksvisible=true, xtickalign=1, ytickalign=1,
                          xminortickalign=1, yminortickalign=1, xticksize=14, xtickwidth=2.5, yticksize=14,
                          ytickwidth=2.5, xminorticksize=7, xminortickwidth=2.5, yminorticksize=7, yminortickwidth=2.5,
                          xticklabelsize=35, yticklabelsize=35, xticksmirrored=true, yticksmirrored=true),
                    Legend=(patchsize=(70, 10), framevisible=false, patchlabelgap=20, rowgap=10))
set_theme!(basic_theme)


xvals = LinRange(0.0,2.0, 100)
f = Figure()
ax = Axis(f[1,1], ylabel=L"\log_{10}(L/L_\odot)", xlabel=L"\log_{10}(T_\mathrm{eff}/[\mathrm{K}])", xreversed=true)
for logM in LinRange(log10(5.0),log10(50.0),10)
    vrot = 400
    logTeff = interpolate_grid_quantity.(Ref(model_set),Ref([logM, vrot]),:logTeff, xvals)
    logL = interpolate_grid_quantity.(Ref(model_set),Ref([logM, vrot]),:logL, xvals)
    lines!(ax, logTeff, logL)
end
save("HR_test.png", f)
f

##
# visualize triangulation
function visualize_simplex_interpolant(ax, si::StarStats.SimplexInterpolant{N,P,LU,E,V}) where {N,P,LU,E,V}
    if size(si.points)[1] != 2
        return
    end
    for simplex in si.simplexes
        lines!(ax, [simplex.points[1,1], simplex.points[1,2], simplex.points[1,3], simplex.points[1,1]],
                   [simplex.points[2,1], simplex.points[2,2], simplex.points[2,3], simplex.points[2,1]],
                   linewidth=3)
    end
end
f = Figure()
ax = Axis(f[1,1], xlabel=L"\log_{10}\left(M_\mathrm{i}/M_\odot\right)", ylabel=L"v_\mathrm{rot,i}\;\mathrm{km\;s^{-1}}")
visualize_simplex_interpolant(ax, model_set.simplex_interpolant)
f
