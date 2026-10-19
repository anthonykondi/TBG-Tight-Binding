using StaticArrays
using LinearAlgebra
import PyPlot as plt
using BenchmarkTools
using ChunkSplitters
using Distributed
using SparseArrays


# note that the graphene indices swap roles on the second layer (should only need to worry about this in geometry construction)
# need to double check that the thresh_connections are being determined corretly

# ===============================================================================
# Geometry
# ===============================================================================

# NOTE: maybe I could also calculate the position vector at runtime and store it here to save speed for memory
"""Label the sites in tbg with 5 dimentional multi-indices (position calculated by constructor)"""
struct SiteIdx         
    moire::SVector{2, Int}
    graphene::SVector{2, Int}
    sublattice::Int               # 1 or 2 (for b1 or b2 vectors)
    layer::Int                    # 1 or 2
    spin::Int                     # +1 or -1
    pos::SVector{3, Float64}      
end


"""Using this for the thresh_connections values (relative moire-graphene-index)"""
struct RelativeIdx
    rel_moire::SVector{2, Int}
    graphene::SVector{2, Int}
    layer::Int
end


"""Geometry and energy constants stored in a common place"""
struct TBGConstants
    d::Float64
    h::Float64 
    mn::SVector{2, Int}
    graphene_cells::Dict{Int, Vector{SVector{2, Int}}}
    θ::Float64
    t0::Float64 
    t1::Float64 
    q_sigma::Float64 
    q_pi::Float64
    r_c::Float64 
    l_c::Float64
    a1_down::SVector{3, Float64}
    a1_up::SVector{3, Float64} 
    a2_down::SVector{3, Float64} 
    a2_up::SVector{3, Float64} 
    b1_down::SVector{3, Float64} 
    b1_up::SVector{3, Float64} 
    b2_down::SVector{3, Float64} 
    b2_up::SVector{3, Float64} 
    L1::SVector{3, Float64}
    L2::SVector{3, Float64}
    r_thresh::Float64 
    n_thresh::Int                 
    thresh_connections::Dict{RelativeIdx, Vector{RelativeIdx}}
end


"""Get the site indices of graphene cells in a moire cell, for layers 1 and 2"""
function get_graphene_cells(mn::SVector{2, Int})::Dict{Int, Vector{SVector{2,Int}}}
    m, n = mn
    R = n + m 
    N = m^2 + m*n + n^2
    graphene_in_moire = Dict{Int, Vector{SVector{2,Int}}}()

    # layer 1
    layer_1_graphene = []
    for i in -R:R, j in 0:2*R
        s1 = (n + m) * i + m * j
        s2 = -m * i + n * j
        if 0 <= s1 < N && 0 <= s2 < N
            push!(layer_1_graphene, [i, j])
        end
    end
    graphene_in_moire[1] = layer_1_graphene

    # layer 2
    layer_2_graphene = []
    for i in -R:R, j in 0:2*R
        s1 = (m + n) * i + n * j
        s2 = -n * i + m * j
        if 0 <= s1 < N && 0 <= s2 < N
            push!(layer_2_graphene, [i, j])
        end
    end
    graphene_in_moire[2] = layer_2_graphene

    return graphene_in_moire
end


"""Get the twist angle based on the commensurate lattice numbers (m,n)"""
function get_twist_angle(mn::SVector{2, Int})::Float64
    m, n = mn
    θ = acos((n^2 + 4*n*m + m^2) / (2 * (n^2 + n*m + m^2)))
    return θ
end


function get_n_thresh(r_thresh::Float64, L1::SVector{3, Float64}, L2::SVector{3, Float64})::Int
    n_thresh = ceil(r_thresh / min(norm(L1), norm(L2)))
    return n_thresh
end


"""
Converts from site indices to real space position vector using tbg geometry.
Meant to be only used by the TBGConstants as it is in the process 
of building the TBGConstants object.
"""
function get_position_tmp(site_idx::RelativeIdx, L1::SVector{3, Float64}, L2::SVector{3, Float64}, 
                          a1::SVector{3, Float64}, a2::SVector{3, Float64}, h::Float64)::SVector{3, Float64}
    # building the position vector
    position_vector = site_idx.rel_moire[1] * L1 + site_idx.rel_moire[2] * L2 
    position_vector += site_idx.graphene[1] * a1 + site_idx.graphene[2] * a2

    if site_idx.layer == 2
        position_vector += [0, 0, h]
    end

    return position_vector
end


"""
Get the connection dictionary for each site in a moire cell (sites connected 
    within the distance threshold). Note that they include themselves in the count.
"""
function get_thresh_connections(r_thresh::Float64, n_thresh::Int,
                                graphene_cells::Dict{Int, Vector{SVector{2, Int}}},
                                L1::SVector{3, Float64}, L2::SVector{3, Float64}, 
                                a1_down::SVector{3, Float64}, a2_down::SVector{3, Float64}, 
                                a1_up::SVector{3, Float64}, a2_up::SVector{3, Float64},
                                h::Float64)::Dict{RelativeIdx, Vector{RelativeIdx}}
    thresh_connections = Dict()
    
    # first looping through the origin graphene cells (keys) 
    for layer_orig in 1:2
        a1_o, a2_o = (layer_orig == 1) ? (a1_down, a2_down) : (a1_up, a2_up)

        for g_orig in graphene_cells[layer_orig]
            orig_site = RelativeIdx((0, 0), g_orig, layer_orig)
            connected_sites = Vector()
            
            # then looping though all the cells in within the moire cells in the threshold (values)
            for rel_moire_1 in (-2 * n_thresh):(2 * n_thresh), rel_moire_2 in (-2 * n_thresh):(2 * n_thresh)
                rel_moire = (rel_moire_1, rel_moire_2)

                for layer_compare in 1:2
                    a1_c, a2_c = (layer_compare == 1) ? (a1_down, a2_down) : (a1_up, a2_up)

                    for g_compare in graphene_cells[layer_compare]
                        compare_site = RelativeIdx(rel_moire, g_compare, layer_compare)
                        pos_orig = get_position_tmp(orig_site, L1, L2, a1_o, a2_o, h)
                        pos_compare = get_position_tmp(compare_site, L1, L2, a1_c, a2_c, h)
                        r_separation = norm(pos_orig - pos_compare)

                        # if the comparison site falls within threshold, add it to the connected sites vector
                        if r_separation < r_thresh 
                            push!(connected_sites, compare_site)
                        end
                    end 
                end
            end
            
            thresh_connections[orig_site] = connected_sites
        end
    end

    return thresh_connections
end


function rot_matrix_xy(θ::Float64)::Matrix{Float64}
    R_θ = [cos(θ) -sin(θ) 0;
           sin(θ) cos(θ) 0;
           0 0 1]
    return R_θ
end


"""constructor for the TBGConstants object"""
function TBGConstants(m::Int, n::Int)::TBGConstants
    0 < m < n || throw(ArgumentError("require 0 < m < n"))
    println("Beginning construction of TBGConstants")

    d = 1.42
    h = 3.34
    mn = SVector{2}(m, n)
    graphene_cells = get_graphene_cells(mn)
    θ = get_twist_angle(mn)
    t0 = 2.7
    t1 = 0.48
    q_sigma = 7.41
    q_pi = 3.15
    r_c = 5.0
    l_c = 0.265
    a1_down = SVector{3}((d*sqrt(3)) * [1, 0, 0])
    a1_up = SVector{3}(rot_matrix_xy(-θ) * a1_down)
    a2_down = SVector{3}((d*sqrt(3)) * [1/2, sqrt(3)/2, 0])
    a2_up = SVector{3}(rot_matrix_xy(-θ) * a2_down)
    b1_down = SVector{3}(d * [sqrt(3)/2, 1/2, 0])
    b1_up = SVector{3}(rot_matrix_xy(-θ) * b1_down)
    b2_down = SVector{3}(d * [sqrt(3)/2, -1/2, 0])
    b2_up = SVector{3}(rot_matrix_xy(-θ) * b2_down)
    L1 = SVector{3}(n * a1_down + m * a2_down)
    L2 = SVector{3}(-m * a1_down + (n + m) * a2_down)
    r_thresh = 7.0                                  # need to think carefully about this value
    n_thresh = get_n_thresh(r_thresh, L1, L2)       # number of moire cells within proximity (in L1/L2 directions, for a total of 4n^2 + 4n + 1)
    thresh_connections = get_thresh_connections(r_thresh, n_thresh, graphene_cells, L1, L2, 
                                                a1_down, a2_down, a1_up, a2_up, h)

    println("Finished construction of TBGConstants")
    return TBGConstants(d, h, mn, graphene_cells, θ, t0, t1, q_sigma, q_pi, r_c, 
                        l_c, a1_down,a1_up, a2_down, a2_up, b1_down, b1_up, b2_down, 
                        b2_up, L1, L2, r_thresh, n_thresh, thresh_connections)    # goodness gracious                
end


"""Converts from site indices to real space position vector using tbg geometry"""
function get_position(moire::SVector{2, Int}, graphene::SVector{2, Int}, sublattice::Int,
                      layer::Int, tbg_csts::TBGConstants)::SVector{3, Float64}
    h = tbg_csts.h
    L1 = tbg_csts.L1
    L2 = tbg_csts.L2

    if layer == 1
        a1 = tbg_csts.a1_up
        a2 = tbg_csts.a2_up
        b1 = tbg_csts.b1_up
        b2 = tbg_csts.b2_up
    else
        a1 = tbg_csts.a1_down
        a2 = tbg_csts.a2_down
        b1 = tbg_csts.b1_down
        b2 = tbg_csts.b2_down
    end

    # building the position vector
    position_vector = moire[1] * L1 + moire[2] * L2 
    position_vector += graphene[1] * a1 + graphene[2] * a2
    if sublattice == 1
        position_vector += b1
    else
        position_vector += b2
    end

    if layer == 2
        position_vector += [0, 0, h]
    end

    return SVector{3}(position_vector)
end


# NOTE: I should also add some kind of wrapping function (mod lattice) to keep graphene idxs within the allowed range
"""Constructor for the SiteIdx object"""
function SiteIdx(moire::SVector{2, Int}, graphene::SVector{2, Int}, sublattice::Int,
                 layer::Int, spin::Int, tbg_csts::TBGConstants)::SiteIdx
    pos = get_position(moire, graphene, sublattice, layer, tbg_csts)
    return SiteIdx(moire, graphene, sublattice, layer, spin, pos)
end


# ===============================================================================
# Physics
# ===============================================================================



# ===============================================================================
# Testing/Benchmarking Functions
# ===============================================================================

"""
benchmarking tools and reminders:
- use @benchmark and display(_) to get the statistics of a function
- use @btime to get time and memory statistics (it prints automatically)
- use @time for time of a single run
"""


"""Plot the moire superlattice unit cell"""
function plot_cells(m::Int, n::Int)
    tbg_csts = TBGConstants(m, n)

    lattice_sites = []
    for layer in 1:2
        for moire_x in 1:5, moire_y in 1:5
            moire = SVector{2}(moire_x, moire_y)
            for graphene in tbg_csts.graphene_cells[layer]
                for sublattice in 1:2
                    site_idx = SiteIdx(moire, graphene, sublattice, layer, 1, tbg_csts)
                    site_pos = site_idx.pos
                    push!(lattice_sites, site_pos)
                end
            end
        end
    end

    xs = [v[1] for v in lattice_sites]
    ys = [v[2] for v in lattice_sites]

    fig, ax = plt.subplots()
    ax.scatter(xs, ys, marker=".", alpha=0.5)

    ax.set_aspect("equal")
    ax.set_title("θ = $(tbg_csts.θ * (180 / pi))")

    plt.show()
end


function plot_thresh_sites(orig_site::SiteIdx, tbg_csts::TBGConstants)
    fig, ax = plt.subplots()

    # =================================================================
    # plot the entire lattice over top as well

    lattice_sites = []
    for layer in 1:2
        for moire_x in -1:1, moire_y in -1:1
            moire = SVector{2}(moire_x, moire_y)
            for graphene in tbg_csts.graphene_cells[layer]
                for sublattice in 1:2
                    site_idx = SiteIdx(moire, graphene, sublattice, layer, 1, tbg_csts)
                    site_pos = site_idx.pos
                    push!(lattice_sites, site_pos)
                end
            end
        end
    end

    xs_2 = [v[1] for v in lattice_sites]
    ys_2 = [v[2] for v in lattice_sites]

    ax.scatter(xs_2, ys_2, marker=".", alpha=0.1)

    # =================================================================

    # convert orig_site from SiteIdx to RelativeIdx object
    orig_site_r = RelativeIdx(orig_site.moire, orig_site.graphene, orig_site.layer)
    connected_sites = tbg_csts.thresh_connections[orig_site_r]
    
    connected_sites_pos = []
    for other_site_r in connected_sites
        for sublattice in 1:2
            other_site = SiteIdx(other_site_r.rel_moire, other_site_r.graphene, sublattice, 
                                 other_site_r.layer, 1, tbg_csts)
            push!(connected_sites_pos, other_site.pos)
        end
    end

    xs = [v[1] for v in connected_sites_pos]
    ys = [v[2] for v in connected_sites_pos]

    ax.scatter(xs, ys, marker=".", alpha=0.5)

    ax.scatter(orig_site.pos[1], orig_site.pos[2], marker="o")

    ax.set_aspect("equal")
    ax.set_title("θ = $(tbg_csts.θ * (180 / pi)), N_thresh_sites = $(size(connected_sites_pos)[1])")

    plt.show()
end


# ===============================================================================
# User Area
# ===============================================================================


function main()
    m, n = (8, 9)
    tbg_csts = TBGConstants(m, n)
    orig_site = SiteIdx(SVector(0, 0), SVector(4, 6), 1, 1, 1, tbg_csts)
    @time plot_thresh_sites(orig_site, tbg_csts)
end


if abspath(PROGRAM_FILE) == @__FILE__
    main()
end


"""
Questions:

- currently, im going to write the code such that it runs through all moire
cells and treats the whole DNxDN matrix, so it can only work for small moire 
superlattice cells (ex: (m,n)=(2,3)). To treat larger moire cells we would need
to decouple the cells and treat them individually, reducing to N DxD matrices 
Note: for our case of (8, 9) --> D ~ 1e3 and for decent resolution in 

- For the threshold distance, I could either set up one distance for the whole 
program, or set up distances for hopping, hartree, and fock separately so that 
each has its own length scale. Then for each threshold distance, I would generate
a dictionary with key: graphene site, and value: (relative moire site, graphene site),
so that I can use sparce matrices and just fill non-zero terms for sites connected to 
the site we are looking at (also note that I won't also account for sublattice as the 
sepration between A and B is small wrt the threshold length scale, so this simplifies it)
"""