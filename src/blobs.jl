using Unitful, UnitfulAstro, DimensionfulAngles
using Cosmology
using Statistics: mean, median
using CairoMakie
using Clustering
using LinearAlgebra: norm

"""Find blobs in a 2D matrix"""
function find_blobs_2D(data::AbstractMatrix{T}, threshold::T, minimum_neighbours, point_source_threshold) where {T}

    # init matrix with ghost cells at edges
    blobby = zeros(size(data, 1) + 2, size(data, 2) + 2)
    blobby[2:end-1, 2:end-1] = data

    # remove point sources
    # iterate over matrix, except for the ghost elements
    @views for i in 2:size(blobby, 2)-1
        for j in 2:size(blobby, 1)-1
            # skip over removed or missing pixels
            if blobby[j, i] == 0
                continue
            elseif ismissing(blobby[j, i])
                continue
            end
            # take a 3x3 slice centred on the target pixel
            chunk = blobby[j-1:j+1, i-1:i+1]

            # remove pixels if their fractional contribution to the chunk is greater than point_source_threshold
            # this is designed to eliminate point sources by detecting high and localised unlikelihood
            # substructure should be more consistent across multiple pixels
            if (blobby[j, i] / sum(chunk)) > point_source_threshold
                blobby[j, i] = 0
            end
        end
    end
    d = count(data .> 0)
    c = count(blobby .> 0)
    @info "$(d - c)/$d pixels removed by point source threshold" point_source_threshold

    # apply threshold
    blobby[blobby.<threshold] .= 0 # use zero for removed pixels
    c1 = count(blobby .> 0)
    @info "$(c - c1)/$c pixels removed by threshold" threshold

    # maybe point source detection should be a seperate, earlier pass.
    @info "Removing pixels based on minimum neighbours" minimum_neighbours
    diff = 1
    # loop until the data stops changing
    @views while diff > 0
        # tracking changes
        c2 = count(blobby .> 0)

        # iterate over matrix, except for the ghost elements
        for i in 2:size(blobby, 2)-1
            for j in 2:size(blobby, 1)-1
                # skip over removed or missing pixels
                if blobby[j, i] == 0
                    continue
                elseif ismissing(blobby[j, i])
                    continue
                end

                # take a 3x3 slice centred on the target pixel
                chunk = blobby[j-1:j+1, i-1:i+1]

                # check that sufficent neighbouring pixels exist
                # we could check for nonzero-ness instead of > threshold but this opens up possibility for changing theshholds
                blobby[j, i] = count(chunk .> threshold) > minimum_neighbours ? blobby[j, i] : 0
            end
        end
        c3 = count(blobby .> 0)
        diff = c2 - c3
        @info "$(c2 - c3)/$c2 pixels removed by neighbour check"
    end

    return blobby[2:end-1, 2:end-1]
end

"""Get points from the boundry of the ellipse"""
function getellipsepoints(e::Ellipse)
    t = range(0, 2 * pi, length=100)
    ellipse_x_r = @. e.r1 / 2 * cos(t)
    ellipse_y_r = @. e.r2 * sin(t)
    R = [cos(e.θ) sin(e.θ); -sin(e.θ) cos(e.θ)]
    r_ellipse = [ellipse_x_r ellipse_y_r] * R
    x = @. e.x + r_ellipse[:, 1]
    y = @. e.y + r_ellipse[:, 2]
    (x, y)
end

"""Fit an ellipse to a set of pixels"""
function fit_ellipse(pixels::Matrix)

    if size(pixels, 1) < 1
        @warn "Insufficent pixels" pixels
    end

    p1 = [0, 0]
    p2 = [0, 0]
    dist = 0.0

    for i in axes(pixels, 1)
        for j in axes(pixels, 1)
            if j > i
                if norm(pixels[i, :] - pixels[j, :]) > dist
                    p1 = pixels[i, :]
                    p2 = pixels[j, :]
                    dist = norm(pixels[i, :] - pixels[j, :])
                end
            end
        end
    end

    majoraxis = dist
    centre = (p1 + p2) / 2
    endpoint = p1 - centre

    @assert norm(endpoint) == majoraxis / 2

    θ = atan(endpoint[2], endpoint[1])

    @debug "Major axis search complete" majoraxis centre[1] centre[2] endpoint[1] endpoint[2] θ pixels

    if majoraxis < 1
        @warn "Improbably small ellipse" majoraxis centre endpoint
    end

    e = Ellipse(
        centre[1],
        centre[2],
        majoraxis,
        2.0,
        θ
    )

    while all(test_point.(Ref(e), eachrow(pixels))) == false
        e = Ellipse(
            centre[1],
            centre[2],
            majoraxis,
            e.r2 + 2.0,
            θ
        )
    end

    if (e.r2 > e.r1)
        # @warn "Major and minor axis inverted" e.r1 e.r2 e.θ
        # mostly affects the very small ones
        e = Ellipse(
            centre[1],
            centre[2],
            e.r2,
            e.r1,
            θ
        )
    end

    return e #, ([p1[1] p1[2]; p2[1] p2[2]])
end

"""Cluster pixels from blobs"""
function cluster(data::Matrix, dbscan_radius::Real=3)
    cluster_input = Matrix{Float64}(undef, (2, count(data .> 0)))

    k = 1
    for j in axes(data, 2)
        for i in axes(data, 1)
            if data[i, j] > 0
                cluster_input[:, k] = [i, j]
                k += 1
            end
        end
    end

    if size(cluster_input)[2] == 0
        return [], []
    end

    @info "Clustering points"
    clusters::Clustering.DbscanResult = dbscan(cluster_input, dbscan_radius, min_cluster_size=3, min_neighbors=2)
    cluster_mat = hcat(cluster_input[1, :], cluster_input[2, :], clusters.assignments)

    @info "Fitting ellipses to clusters"
    trimmed = cluster_mat[findall(>(1), cluster_mat[:, 3]), :]
    ellipses = [(fit_ellipse(trimmed[trimmed[:, 3].==i, 1:2]), i) for i in unique(trimmed[:, 3])]

    return trimmed, ellipses
end

function plot_blobs(
    data::Matrix,
    raw::Matrix,
    obs::Matrix,
    centre::NTuple{2,<:Real},
    # r200::Real,
    dbscan_radius::Real,
    colormap=:binary
)
    # generate clusters
    clustered, ellipses = cluster(data, dbscan_radius)

    @info "Generating plots"
    f = Figure(size=(3400, 1600))

    # heatmap of blobs
    ax1 = Axis(f[1, 1], aspect=1, width=600, height=600, title="Blob Finder")
    h1 = heatmap!(data, colormap=colormap, colorscale=Makie.pseudolog10)
    scatter!(centre..., marker=:cross, color=:white)
    # poly!(Circle(Point2f(centre...), r200), color=:transparent, strokecolor=:white, strokewidth=1, overdraw=true, transparency=true)
    Colorbar(f[1, 2], h1)

    # annotate with ellipses around clusters
    for tup in ellipses
        e, _ = tup
        lines!(getellipsepoints(e)..., color=:green)
    end

    # original likelihood data
    ax2 = Axis(f[1, 3], aspect=1, width=600, height=600, title="Absolute Value of Log-likelihood")
    h2 = heatmap!(raw, colormap=colormap, colorscale=Makie.pseudolog10)
    scatter!(centre..., marker=:cross, color=:white)
    # poly!(Circle(Point2f(centre...), r200), color=:transparent, strokecolor=:white, strokewidth=1, overdraw=true, transparency=true)
    Colorbar(f[1, 4], h2)

    # original obs
    ax3 = Axis(f[1, 5], aspect=1, width=600, height=600, title="Original Observation")
    h3 = heatmap!(obs, colormap=colormap, colorscale=Makie.pseudolog10)
    scatter!(centre..., marker=:cross, color=:white)
    # poly!(Circle(Point2f(centre...), r200), color=:transparent, strokecolor=:white, strokewidth=1, overdraw=true, transparency=true)
    Colorbar(f[1, 6], h3)

    linkaxes!(ax1, ax2, ax3)

    xlims!(0, size(data, 2) + 1)
    ylims!(0, size(data, 1) + 1)

    # plot clusters in 3d
    if length(clustered) > 0
        ax4 = Axis3(f[2, 1], width=600, height=600, title="Clustered ($(trunc(Int64, maximum(clustered[:,3]))) clusters)")
        scatter!(clustered[:, 1], clustered[:, 2], clustered[:, 3], color=clustered[:, 3])
        for tup in ellipses
            e, i = tup
            el = getellipsepoints(e)
            z = repeat([i], size(el[1], 1))
            lines!(el..., z, color=:orange)
        end
    end

    # histogram of data dist.
    ax5 = Axis(f[2, 3:6], yscale=log10, xlabel="log-likelihood", ylabel="Count in bin",)
    bins = range(extrema(data)..., length=30)
    hist!(vec(raw), bins=bins, label="Abs. Log-likelihood")
    hist!([d for d in data if d > 0], bins=bins, label="Blobs")
    Legend(f[2, 6], ax5, valign=:top)


    return f
end


function run_blob_analysis(
    observation::Array{<:Number,3},
    log_likelihood::Array{<:Union{Missing,<:Real},3},
    ij::NTuple{2,<:Real},
    # r200::Unitful.Length,
    minimum_neighbours=3,
    point_source_threshold=0.5,
)
    replace!(log_likelihood, missing => 0)
    ll_2d = abs.(dropdims(sum(log_likelihood, dims=1), dims=1))
    obs_2d = dropdims(sum(observation, dims=1), dims=1)

    b = find_blobs_2D(ll_2d, mean(ll_2d), minimum_neighbours, point_source_threshold)
    # p = plot_blobs(b, ll_2d, obs_2d, ij, r200, 4, :plasma)
    p = plot_blobs(b, ll_2d, obs_2d, ij, 4, :plasma)
    return p
end