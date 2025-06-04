## beware docstrings are not fully consistent
using StaticArrays
using LinearAlgebra
using Quickhull
using StatsBase

"""
    struct InterpSimplex
Structure to describe a simplex together with precomputed properties used to determine barycentric coordinates. In ``N`` dimensions, given ``N+1`` vertices on a simplex ``\\vec{x}_i`` and a point ``\\vec{a}`` for which we want to determine the barycentric coordinates ``\\alpha_i``, it should be satisfied that
``\\vec{a} = \\sum_{i=1}^{N+1}\\alpha_i \\vec{x}_i``
together with the condition
``\\sum_{i=1}^{N+1}\alpha_i = 1.``
This is equivalent to solving a linear system ``A\vec{\alpha}=\vec{b}`` where
``A=[(\\vec{x}_1-\\vec{x}_N)\\;\\dots\\;(\\vec{x}_{N-1}-\\vec{x}_N)],\\qquad \\vec{b}=\\vec{a}-\\vec{x}_N``
and ``\\vec{\\alpha}`` contains th barycentric coordinates from ``1`` to ``(N-1)``. To efficiently determine barycentric coordinates we use `StaticArrays` with precomputed LU factorizations of the ``A`` matrix. We also store the vertex ``\\vec{x}_N`` separately.
- `points`: matrix containing the points with dimensions ndims x npoints
- `A_lu`: factorized matrix in the linear system to determine the barycentric coordinate coefficients
- `end_point`: Coordinates of ``\\vec{x}_N``, used to set up the value of ``\\vec{b}`` in the linear system when a point ``\\vec{a}`` is provided.
"""
struct InterpSimplex{N,P,LU,E}
    points::P
    A_lu::LU
    end_point::E
    max_span::E
    min_span::E
    barycenter::E
    point_indeces::Vector{Int}
end

"""
    function InterpSimplex(points)
Construct a simplex based on the given `points`. `points` must be a matrix of dimensions ndims x npoints. 
"""
function InterpSimplex(points::Matrix, point_indeces::Vector{Int})
    ndims = size(points,1)
    npoints = size(points,2)
    if ndims != npoints -1
        throw(DimensionMismatch(["Input matrix of points must be of size N x N+1"]))
    end

    A = zeros(ndims,ndims)
    for i in 1:ndims
        A[:,i] .= points[:,i] .- points[:,end]
    end

    sv_points = SizedMatrix{ndims, ndims+1}(points)
    sv_A = SizedMatrix{ndims,ndims}(A)
    A_lu = lu(sv_A)
    end_point = SVector{ndims}(points[:,npoints])

    max_span = points[:,1]
    min_span = points[:,1]
    for i in 2:ndims+1 # iterate over remaining simplex points
        for j in 1:ndims # iterate over dimensions
            if points[j,i] > max_span[j]
                max_span[j] = points[j,i]
            end
            if points[j,i] < min_span[j]
                min_span[j] = points[j,i]
            end
        end
    end
    max_span = SVector{ndims}(max_span)
    min_span = SVector{ndims}(min_span)

    barycenter = dropdims(sum(points,dims=2);dims=2)./npoints
    barycenter = SVector{ndims}(barycenter)

    InterpSimplex{ndims, typeof(sv_points), typeof(A_lu), typeof(end_point)}(sv_points, A_lu, end_point, max_span, min_span, barycenter, point_indeces)
end

"""
    function barycentric_coords!(point, s::InterpSimplex{N,P,LU,E}, coords, bcache) where {N,P,LU,E}
Non-allocating calculation of barycenter coordinates. Results are stored in place in `coords`.
"""
function barycentric_coords!(point, s::InterpSimplex{N,P,LU,E}, coords, bcache) where {N,P,LU,E}
    bcache .= point .- s.end_point
    res = s.A_lu\bcache
    for i in 1:N
        coords[i] = res[i]
    end
    coords[N+1] = 1-sum(res)
    return coords
end

"""
    function barycentric_coords(point, s::InterpSimplex{N,P,LU,E}) where {N,P,LU,E}
Wrapper for `barycentric_coords!` that is non-mutating.
"""
function barycentric_coords(point, s::InterpSimplex{N,P,LU,E}) where {N,P,LU,E}
    coords = zeros(N+1)
    bcache = zeros(N)
    barycentric_coords!(point, s, coords, bcache)
    return coords
end

function plot_simplex_2d(ax,s,color)
    p1 = s.points[:,1]
    p2 = s.points[:,2]
    p3 = s.points[:,3]
    lines!(ax, [p1[1],p2[1],p3[1],p1[1]], [p1[2],p2[2],p3[2],p1[2]], color=color)
end

struct NullSimplexTree end
struct SimplexTree{V}
    simplexes::V
    tree_lower::Union{SimplexTree{V},NullSimplexTree}
    tree_higher::Union{SimplexTree{V},NullSimplexTree}
    max_span::Vector{Float64}
    min_span::Vector{Float64}
    cut::Float64
    cut_dim::Int
    depth::Int
    finish::Bool
end

function SimplexTree(simplexes,depth,maxdepth,ndims)
    max_span = zeros(ndims)
    min_span = zeros(ndims)
    max_span .= simplexes[1].max_span
    min_span .= simplexes[1].min_span
    for simplex in simplexes
        for idim in 1:ndims
            if (max_span[idim] < simplex.max_span[idim])
                max_span[idim] = simplex.max_span[idim]
            end
            if (min_span[idim] < simplex.min_span[idim])
                min_span[idim] = simplex.min_span[idim]
            end
        end
    end
    if depth==maxdepth
        tree_lower = NullSimplexTree()
        tree_higher = NullSimplexTree()
        cut = NaN
        cut_dim = 0
        finish = true
    else
        barycenters = zeros(length(simplexes))
        median_barys = zeros(ndims)
        indeces_lower = zeros(Bool,ndims,length(simplexes))
        indeces_higher = zeros(Bool,ndims,length(simplexes))
        for idim in 1:ndims # do cuts in all dimensions, check for best after
            for j in eachindex(simplexes)
                barycenters[j] = simplexes[j].barycenter[idim]
            end
            median_bary = median(barycenters)
            median_barys[idim] = median_bary
            for j in eachindex(simplexes)
                if simplexes[j].max_span[idim] > median_barys[idim]
                    indeces_higher[idim, j] = 1
                end
                if simplexes[j].min_span[idim] < median_barys[idim]
                    indeces_lower[idim, j] = 1
                end
            end
        end
        best_cut = length(simplexes)
        i_best_cut = 0
        for idim in 1:ndims #find best cut as the one that produces the smaller biggest set
            num_higher = sum(indeces_higher[idim,:])
            num_lower = sum(indeces_lower[idim,:])
            biggest_set = max(num_higher,num_lower)
            if biggest_set <= best_cut
                best_cut = biggest_set
                i_best_cut = idim
            end
        end
        indeces_higher = indeces_higher[i_best_cut,:]
        indeces_lower = indeces_lower[i_best_cut,:]
        num_higher = sum(indeces_higher)
        num_lower = sum(indeces_lower)
        #check if cut meaningfully reduces candidates
        if max(num_higher,num_lower)/length(simplexes) > 0.8
            tree_lower = NullSimplexTree()
            tree_higher = NullSimplexTree()
            cut = NaN
            cut_dim = 0
            finish = true
        else
            tree_lower = SimplexTree(simplexes[indeces_lower], depth+1, maxdepth,ndims)
            tree_higher = SimplexTree(simplexes[indeces_higher], depth+1, maxdepth,ndims)
            cut = median_barys[i_best_cut]
            cut_dim = i_best_cut
            finish = false
        end
    end
    return SimplexTree(simplexes, tree_lower,tree_higher,max_span,min_span,cut,cut_dim,depth,finish)
end

function find_subset(point,tree)
    if tree.finish
        return tree.simplexes
    else
        if point[tree.cut_dim] > tree.cut
            return find_subset(point, tree.tree_higher)
        else
            return find_subset(point, tree.tree_lower)
        end
    end
end

"""
    struct SimplexInterpolant{N,P,LU,E}
Structure containing information on a Delaunay triangulation to perform interpolation.
- `simplexes`: Instances of `Simplex` for each of the simplexes in the triangulation
- `points`: Points used for the triangulation. `points` must be a matrix of dimensions ndims x npoints.
- `points_indeces`: For each simplex, contains the indeces of the points that compose it. Dimensions are `ndims+1` x `nsimplexes`.
"""
struct SimplexInterpolant{N,P,LU,E,V}
    simplexes::Vector{InterpSimplex{N,P,LU,E}}
    points::Matrix
    simplex_tree::SimplexTree{V}
end

"""
    function SimplexInterpolant(points)
Constructor for a `SimplexInterpolant` based on a set of `points`. `points` must be a matrix of dimensions ndims x npoints. 
"""
function SimplexInterpolant(points; maxdepth=12)
    tri = delaunay(points)
    facets_eval = facets(tri)
    ndims = size(points,1)
    simplexes = []
    for (i, facet) in enumerate(facets_eval)
        indeces = facet.data
        simplex_points = zeros(ndims, ndims+1)
        for j in 1:(ndims+1)
            simplex_points[:,j] .= points[:,indeces[j]]
        end
        indeces_vec = zeros(Int, ndims+1)
        indeces_vec .= [indeces...]
        push!(simplexes, InterpSimplex(simplex_points, indeces_vec)) #turn indeces from tuple to vector
    end
    simplexes_typed::Vector{typeof(simplexes[1])} = [simplex for simplex in simplexes]

    simplex_tree = SimplexTree(simplexes_typed,1,maxdepth,ndims)

    SimplexInterpolant(simplexes_typed,points, simplex_tree)
end

"""
    function interpolation_info(point, si::SimplexInterpolant{N,P,LU,E}) where {N,P,LU,E}
Iterate through all simplexes in the `SimplexInterpolant` `si` to find the simplex containing the given `point`. Returns the Barycentric coordinates of the point for the containing simplex and the indeces of the points that form the vertices of the simplex. If no simplex is found that contains the point zero values are returned for the indeces and coefficients.
"""
function interpolation_info(point::Vector{T}, si::SimplexInterpolant{N,P,LU,E,V}) where {T,N,P,LU,E,V}
    coords = zeros(T,N+1)
    bcache = zeros(T,N)
    simplexes = find_subset(point,si.simplex_tree)
    for i in eachindex(simplexes)
        simplex = simplexes[i]
        barycentric_coords!(point, simplex, coords, bcache)
        if minimum(coords) >= 0
            return (coords, simplex.point_indeces)
        end
    end
    coords .= 0
    point_indeces = zeros(N+1)
    return (coords, point_indeces) # will be zero if nothing was found
end

"""
    function compute_simplex_interpolation(info, values)
Given the results from the function `interpolation_info` and the `values` of a function evaluated in the points used to create the `SimplexInterpolant`, return the interpolated value. Returns `NaN` if no simplex contains the point.
"""
function compute_simplex_interpolation(info, values)
    val = 0
    coeffs = info[1]
    indeces = info[2]
    if indeces[1]==0
        return NaN # in case nothing is found return NaN
    end
    for i in eachindex(coeffs)
        val = val + coeffs[i]*values[indeces[i]]
    end
    return val
end