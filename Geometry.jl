using StaticArrays

# use export keyword for sending out structs and functions

# ======================================================================================
# Parameters 
# ======================================================================================

const LATTICE_CONSTANTS = Dict("d" => 1, "h" => 0.3, "t0" => 1, "t1" => 1, "mn" => (17,18))    # set these based on the paper's settings


struct SiteIdx
    moire::SVector{2, Int}
    graphene::SVector{2, Int}     # (i, j) cannot exceed (the moire unit cell) -> need to add error raising
    sublattice::Int               # 1 or 2 (for b1 or b2 vectors)
    layer::Int                    # 1 or 2
    spin::Int                     # +1 or -1
end


# ======================================================================================
# Functions
# ======================================================================================


"""Get the allowed graphene indices based on the moire unit cell size, returns a dictionary"""
function get_graphene_in_moire_idxs()::Dict{Int, Vector{SVector{2,Int}}} 
    m, n = LATTICE_CONSTANTS["mn"]
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


function get_twist_angle()::Float64 
    m, n = LATTICE_CONSTANTS["mn"]
    θ = acos((n^2 + 4*n*m + m^2) / (2 * (n^2 + n*m + m^2)))
    return θ
end


function rot_matrix_xy(θ::Float64)::Matrix{Float64}
    R_θ = [cos(θ) -sin(θ) 0;
           sin(θ) cos(θ) 0;
           0 0 1]
    return R_θ
end


function get_position(site_idx::SiteIdx)::SVector{3, Float64}
    d = LATTICE_CONSTANTS["d"]
    h = LATTICE_CONSTANTS["h"]
    m, n = LATTICE_CONSTANTS["mn"]
    θ = get_twist_angle()

    # bottom layer vectors
    a1 = (d*sqrt(3)) * [1, 0, 0]
    a2 = (d*sqrt(3)) * [1/2, sqrt(3)/2, 0]
    b1 = d * [sqrt(3)/2, 1/2, 0]
    b2 = d * [sqrt(3)/2, -1/2, 0]
    L1 = n * a1 + m * a2
    L2 = -m * a1 + (n + m) * a2

    # switch to top layer vectors if on second layer
    if site_idx.layer == 2
        a1 = rot_matrix_xy(-θ) * a1
        a2 = rot_matrix_xy(-θ) * a2
        b1 = rot_matrix_xy(-θ) * b1
        b2 = rot_matrix_xy(-θ) * b2
    end

    # building the position vector
    position_vector = site_idx.moire[1] * L1 + site_idx.moire[2] * L2 
    position_vector += site_idx.graphene[1] * a1 + site_idx.graphene[2] * a2
    if site_idx.sublattice == 1
        position_vector += b1
    else
        position_vector += b2
    end

    if site_idx.layer == 2
        position_vector += [0, 0, h]
    end

    return position_vector
end


# want to test this by plotting 

# calculated at runtime
const GRAPHENE_IDXS = get_graphene_in_moire_idxs()

# lattice_sites = []
# for layer in 1:2
#     for moire_x in 1:1, moire_y in 1:1
#         for graphene in GRAPHENE_IDXS[layer]
#             for sublattice in 1:2
#                 site_idx = SiteIdx([moire_x, moire_y], graphene, sublattice, layer, 1)
#                 site_pos = get_position(site_idx)
#                 push!(lattice_sites, site_pos)
#             end
#         end
#     end
# end

# lattice_sites_x = [v[1] for v in lattice_sites]
# lattice_sites_y = [v[2] for v in lattice_sites]
# lattice_sites_z = [v[3] for v in lattice_sites]

# plt.scatter3D(lattice_sites_x, lattice_sites_y, lattice_sites_z, marker=".")
# plt.show()