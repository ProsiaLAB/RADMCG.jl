#=
==========================================================================
                    Adaptive Mesh Refinement Module

 This is a simple adaptive mesh refinement module. Nothing fancy, but
 functional. It does not support ghost cells in each branch, nor does it do
 zone decomposition for MPI parallel computing. For such complexity we
 refer to e.g. the PARAMESH or CHOMBO library. The objective of the present
 AMR module is to provide a simple mesh routine which allows to refine
 anywhere at will by dividing a cell into 2x2x2 (for 3-D) subcells.  In the
 end this module delivers a base grid in which each grid cell can be either
 a true cell (leaf) or a tree of sub-cells. All the branches are linked via
 pointers (a pointer to the parent and pointers to the children) and
 neighbors are also linked via pointers.  Moreover one can let the module
 produce a 1-D array of pointers to all cells (leafs), so that referencing
 to cells can be done in a do-loop rather than requiring the necessity of
 recursion. A similar thing can be done to all branches (i.e. leafs AND
 nodes). 

 The objective of this module is to provide a meshing for e.g.  radiative
 transfer or for diffusion-type equations. It is not quite efficient,
 however, for the explicit solving of hydrodynamics equations, since
 the module is not adapted for parallelization. For that we refer to
 the PARAMESH library or CHOMBO library.

--------------------------------------------------------------------------
 BASE GRID AND TREE STRUCTURED REFINEMENTS

 The base grid (forest) with possibly hierarchically structured 
 refinements (tree). The end-members of the branches are called leafs.
 The leafs are the "true" cells! Each symbol (o or *) is a branch 
 (the amr_branch type). A branch can be a node (consisting of children)
 or a leaf (end-member = AMR grid cell).

  Base grid (level=0)        o----------*----------o----------o----------o
                                      /||\
                                    / |  | \
  First level refinement           *  o  o  o
                                 /||\\   
                               / |  | \
  Second level refinement     o  o  o  o

 Two types of branches:
   o = leaf
   * = node

 NOTE: Compared to previous versions of amr_module, the base level is now
       0 instead of 1. This is in better comparison to the levels of the
       layers, and is in many respects more consistent.

--------------------------------------------------------------------------
 HOW THE REFEINEMENT WORKS IN 2-D (EXAMPLE)

  |-----------|-----------|-----------------------|
  |           |           |                       |
  |           |           |                       |
  |           |           |                       |
  |-----|-----|-----------|                       |
  |     |     |           |                       |
  |-----|-----|           |                       |
  |     |     |           |                       |
  |-----|-----|-----------|-----------------------|
  |                       |                       |
  |                       |                       |
  |                       |                       |
  |                       |                       |
  |                       |                       |
  |                       |                       |
  |                       |                       |
  |-----------------------|-----------------------|

--------------------------------------------------------------------------
 Author: C.P. Dullemond
 Date:   07.09.07
==========================================================================
=#

include("Constants.jl")
include("Utils.jl")



mutable struct AMRBranch
    level::Int
    id::Int
    branch_index::Int
    leaf_index::Int
    leaf::Bool
    neighbor::Array{AMRBranch,2}
    child::Array{AMRBranch,3}
    parent::Union{AMRBranch,Nothing}
    parent_slot::Array{Int,1}
    ixyzf::Array{Int,1}

    # Constructor
    function AMRBranch(
        level::Int,
        id::Int,
        branch_index::Int=0,
        leaf_index::Int=0,
        leaf::Bool=false,
        parent::Union{AMRBranch,Nothing}=nothing,
        parent_slot::Array{Int,1}=[0, 0, 0],
        ixyzf::Array{Int,1}=[0, 0, 0],
    )
        neighbor = Array{AMRBranch,2}(undef, 0, 0)
        child = Array{AMRBranch,3}(undef, 0, 0, 0)

        new(level, id, branch_index, leaf_index, leaf, neighbor, child, parent, parent_slot, ixyzf)
    end
end


mutable struct AMRBranchLink
    link::Union{AMRBranch,Nothing}
end


function amr_initialize(
    include_x::Bool,
    include_y::Bool,
    include_z::Bool,
    nx::Int, ny::Int, nz::Int,
    xi::Array{Float64},
    yi::Array{Float64},
    zi::Array{Float64},
    level_max::Int,
    fill_base_grid::Int,
    nr_leaves_max::Int,
    xcyclic::Bool,
    ycyclic::Bool,
    zcyclic::Bool,
    always_amr_tree::Bool,
    nr_branches_max::Int=0,
)
    amr_leaf_index_free_index = 1
    amr_branch_index_free_index = 1

    if level_max > 0
        amr_tree_present = true
    else
        amr_tree_present = false
        if always_amr_tree
            amr_tree_present = true
        end
    end

    if nx <= 0
        throw(ArgumentError("nx must be >= 1"))
    end

    if ny <= 0
        throw(ArgumentError("ny must be >= 1"))
    end

    if nz <= 0
        throw(ArgumentError("nz must be >= 1"))
    end

    amr_dim = 0
    amr_grid_nx = nx
    amr_grid_ny = ny
    amr_grid_nz = nz

    if include_x
        amr_xdim = 1
        amr_dim += 1
    else
        if nx != 1
            throw(ArgumentError("nx must be 1 if include_x is false"))
        end
        amr_xdim = 0
    end

    if include_y
        amr_ydim = 1
        amr_dim += 1
    else
        if ny != 1
            throw(ArgumentError("ny must be 1 if include_y is false"))
        end
        amr_ydim = 0
    end

    if include_z
        amr_zdim = 1
        amr_dim += 1
    else
        if nz != 1
            throw(ArgumentError("nz must be 1 if include_z is false"))
        end
        amr_zdim = 0
    end

    if amr_dim == 0
        throw(ArgumentError("Zero-dimensional AMR error"))
    end

    amr_nrchildref = (1 + amr_xdim) * (1 + amr_ydim) * (1 + amr_zdim)
    amr_xyzdim = [amr_xdim, amr_ydim, amr_zdim]

    amr_cyclic_xyz = [false, false, false]

    if xcyclic
        amr_cyclic_xyz[1] = true
    end

    if ycyclic
        amr_cyclic_xyz[2] = true
    end

    if zcyclic
        amr_cyclic_xyz[3] = true
    end

    amr_leafcount = 0
    amr_branchcount = 0
    amr_nrleaves = 0
    amr_nrbranches = 0
    amr_level_max = level_max
    amr_last_branch_id = 0

    if amr_tree_present
        if (nr_branches_max != 0)
            amr_use_index = true
            amr_nrbranches_max = nr_branches_max
            amr_index_to_branch = Array{AMRBranchLink}(undef, nr_branches_max)
            amr_index_to_leaf = Array{AMRBranchLink}(undef, nr_branches_max)
            amr_branch_index_holes = zeros(Int, nr_branches_max)
            amr_leaf_index_holes = zeros(Int, nr_branches_max)
            amr_branch_index_nrholes = 0
            amr_leaf_index_nrholes = 0
            if nr_leaves_max != 0
                if nr_leaves_max > nr_branches_max
                    throw(ArgumentError("Cannot have more leaves than branches"))
                end
                amr_nrleaves_max = nr_leaves_max
            else
                amr_nrleaves_max = nr_branches_max
            end
        else
            if nr_leaves_max != 0
                throw(ArgumentError("Cannot have leaves without branches"))
            end
            amr_use_index = false
            amr_nrbranches_max = 0
        end

    else
        amr_use_index = true
        amr_level_max = 0
        amr_nrleaves_max = amr_grid_nx * amr_grid_ny * amr_grid_nz
        amr_nrbranches_max = amr_grid_nx * amr_grid_ny * amr_grid_nz
    end

    if amr_tree_present
        amr_grid_branch = Array{AMRBranch,3}(undef, amr_grid_nx, amr_grid_ny, amr_grid_nz)
    end

    nxyzmax = max(amr_grid_nx, amr_grid_ny, amr_grid_nz)
    amr_grid_xi = zeros(Float64, nxyzmax + 1, 3)

    for ix = 1:(nx+1)
        amr_grid_xi[ix, 1] = xi[ix]
    end

    for ix = 2:(nx+1)
        if xi[ix] <= xi[ix-1]
            throw(ArgumentError("xi must be monotonically increasing"))
        end
    end

    for iy = 1:(ny+1)
        amr_grid_xi[iy, 2] = yi[iy]
    end

    for iy = 2:(ny+1)
        if yi[iy] <= yi[iy-1]
            throw(ArgumentError("yi must be monotonically increasing"))
        end
    end

    for iz = 1:(nz+1)
        amr_grid_xi[iz, 3] = zi[iz]
    end

    for iz = 2:(nz+1)
        if zi[iz] <= zi[iz-1]
            throw(ArgumentError("zi must be monotonically increasing"))
        end
    end

    nxmax = nx * (2^level_max)
    nymax = ny * (2^level_max)
    nzmax = nz * (2^level_max)

    amr_nxyzfmax = max(nxmax, nymax, nzmax)

    amr_finegrid_xi = Array{Float64}(undef, amr_nxyzfmax + 1, 3, level_max + 1)
    amr_finegrid_xc = Array{Float64}(undef, amr_nxyzfmax, 3, level_max + 1)

    amr_finegrid_xi[1:nx+1, 1, 1] = amr_grid_xi[1:nx+1, 1]
    amr_finegrid_xi[1:ny+1, 2, 1] = amr_grid_xi[1:ny+1, 2]
    amr_finegrid_xi[1:nz+1, 3, 1] = amr_grid_xi[1:nz+1, 3]

    for ix = 1:nx
        amr_finegrid_xc[ix, 1, 1] = 0.5 * (
            amr_finegrid_xi[ix, 1, 1] + amr_finegrid_xi[ix+1, 1, 1]
        )
    end

    for iy = 1:ny
        amr_finegrid_xc[iy, 2, 1] = 0.5 * (
            amr_finegrid_xi[iy, 2, 1] + amr_finegrid_xi[iy+1, 2, 1]
        )
    end

    for iz = 1:nz
        amr_finegrid_xc[iz, 3, 1] = 0.5 * (
            amr_finegrid_xi[iz, 3, 1] + amr_finegrid_xi[iz+1, 3, 1]
        )
    end

    nnx = nx
    nny = ny
    nnz = nz

    if amr_tree_present
        for ilevel = 1:level_max
            if include_x
                for ix = 1:(nnx+1)
                    amr_finegrid_xi[2*(ix-1)+1, 1, ilevel+1] = amr_finegrid_xi[ix, 1, ilevel]
                end
                for ix = 1:nnx
                    dx = (amr_finegrid_xi[ix+1, 1, ilevel] - amr_finegrid_xi[ix, 1, ilevel]) / 2
                end
                nnx = 2 * nnx
            else
                if nnx != 1
                    throw(ArgumentError("nnx must be 1 if include_x is false"))
                end
                amr_finegrid_xi[1, 1, ilevel+1] = amr_finegrid_xi[1, 1, ilevel]
                amr_finegrid_xi[2, 1, ilevel+1] = amr_finegrid_xi[2, 1, ilevel]
            end
            if include_y
                for iy = 1:(nny+1)
                    amr_finegrid_xi[2*(iy-1)+1, 2, ilevel+1] = amr_finegrid_xi[iy, 2, ilevel]
                end
                for iy = 1:nny
                    dy = (amr_finegrid_xi[iy+1, 2, ilevel] - amr_finegrid_xi[iy, 2, ilevel]) / 2
                end
                nny = 2 * nny
            else
                if nny != 1
                    throw(ArgumentError("nny must be 1 if include_y is false"))
                end
                amr_finegrid_xi[1, 2, ilevel+1] = amr_finegrid_xi[1, 2, ilevel]
                amr_finegrid_xi[2, 2, ilevel+1] = amr_finegrid_xi[2, 2, ilevel]
            end
            if include_z
                for iz = 1:(nnz+1)
                    amr_finegrid_xi[2*(iz-1)+1, 3, ilevel+1] = amr_finegrid_xi[iz, 3, ilevel]
                end
                for iz = 1:nnz
                    dz = (amr_finegrid_xi[iz+1, 3, ilevel] - amr_finegrid_xi[iz, 3, ilevel]) / 2
                end
                nnz = 2 * nnz
            else
                if nnz != 1
                    throw(ArgumentError("nnz must be 1 if include_z is false"))
                end
                amr_finegrid_xi[1, 3, ilevel+1] = amr_finegrid_xi[1, 3, ilevel]
                amr_finegrid_xi[2, 3, ilevel+1] = amr_finegrid_xi[2, 3, ilevel]
            end
            for ix = 1:nnx
                amr_finegrid_xc[ix, 1, ilevel+1] = 0.5 * (
                    amr_finegrid_xi[ix, 1, ilevel+1] + amr_finegrid_xi[ix+1, 1, ilevel+1]
                )
            end
            for iy = 1:nny
                amr_finegrid_xc[iy, 2, ilevel+1] = 0.5 * (
                    amr_finegrid_xi[iy, 2, ilevel+1] + amr_finegrid_xi[iy+1, 2, ilevel+1]
                )
            end
            for iz = 1:nnz
                amr_finegrid_xc[iz, 3, ilevel+1] = 0.5 * (
                    amr_finegrid_xi[iz, 3, ilevel+1] + amr_finegrid_xi[iz+1, 3, ilevel+1]
                )
            end
        end
    end

    a = AMRBranch(0, 0)
    b = Nothing
    if amr_tree_present
        slot = Array{Int,1}(undef, 3)
        if fill_base_grid == 1
            for iz = 1:amr_grid_nz
                for iy = 1:amr_grid_ny
                    for ix = 1:amr_grid_nx
                        slot[1] = ix
                        slot[2] = iy
                        slot[3] = iz
                        a = amr_branch_construct(
                            a,
                            Nothing,
                            slot,
                            amr_nrleaves,
                            amr_nrbranches,
                            amr_last_branch_id,
                            amr_xyzdim,
                            amr_branch_index_nrholes,
                            amr_leaf_index_free_index,
                            amr_branch_index_free_index,
                        )
                    end
                end
            end
        end
    end


end


function amr_branch_construct(
    a::Union{AMRBranch,Nothing},
    b::Union{AMRBranch,Nothing},
    slot::Array{Int,1},
    amr_nrleaves::Int,
    amr_nrbranches::Int,
    amr_last_branch_id::Int,
    amr_xyzdim::Array{Int,1},
    amr_branch_index_nrholes::Int,
    amr_leaf_index_free_index::Int,
    amr_branch_index_free_index::Int,
)
    amr_nrbranches = amr_nrbranches + 1
    amr_last_branch_id = amr_last_branch_id + 1

    a.id = amr_last_branch_id
    a.leaf = true

    amr_nrleaves = amr_nrleaves + 1

    if amr_use_index
        index = amr_assign_branch_index(a, amr_branch_index_nrholes, amr_branch_index_free_index)
        a.branch_index = index
        index = amr_assign_leaf_index(a, amr_leaf_index_nrholes, amr_leaf_index_free_index)
        a.leaf_index = index
    end

    a.child = Nothing

    if b isa AMRBranch
        a.parent = b
        b.child[slot[1], slot[2], slot[3]] = a
        a.level = b.level + 1
        for idir = 1:3
            if amr_xyzdim[idir] == 1
                a.ixyzf[idir] = 2 * (b.ixyzf[idir] - 1) + slot[idir]
            else
                a.ixyzf[idir] = 1
            end
        end
    else
        a.parent = Nothing
        a.level = 0
        if (slot[1] < 1 || slot[1] > amr_grid_nx) || (slot[2] < 1 || slot[2] > amr_grid_ny) || (slot[3] < 1 || slot[3] > amr_grid_nz)
            throw(ArgumentError("Invalid slot"))
        end

        if amr_grid_branch[slot[1], slot[2], slot[3]] isa AMRBranch
            throw(ArgumentError("Branch already exists"))
        end

        amr_grid_branch[slot[1], slot[2], slot[3]] = a

        for idir = 1:3
            a.ixyzf[idir] = slot[idir]
        end
    end

    for i = 1:3
        a.parent_slot[i] = slot[i]
    end

    amr_find_and_link_all_neighbors_branch()

    return a
end

function amr_assign_branch_index(
    a::Union{AMRBranch,Nothing},
    amr_branch_index_nrholes::Int,
    amr_branch_index_free_index::Int,
)
    if !amr_use_index
        throw(ArgumentError("Cannot assign branch index if not using index"))
    end

    if amr_branch_index_nrholes
        if amr_branch_index_free_index > amr_nrbranches_max
            throw(ArgumentError("Cannot allocate more holes due to lack of memory"))
        end
        index = amr_branch_index_free_index
        amr_branch_index_free_index = amr_branch_index_nrholes + 1
    else
        index = amr_branch_index_holes[amr_branch_index_nrholes]
        amr_branch_index_nrholes = amr_branch_index_nrholes - 1
    end

    if amr_index_to_branch[index] isa AMRBranch
        throw(ArgumentError("Branch index already in use"))
    end

    amr_index_to_branch[index] = a

end


function amr_assign_leaf_index(
    a::Union{AMRBranch,Nothing},
    amr_leaf_index_nrholes::Int,
    amr_leaf_index_free_index::Int,
)
    if !amr_use_index
        throw(ArgumentError("Cannot assign leaf index if not using index"))
    end

    if amr_leaf_index_nrholes
        if amr_leaf_index_free_index > amr_nrleaves_max
            throw(ArgumentError("Cannot allocate more holes due to lack of memory"))
        end
        index = amr_leaf_index_free_index
        amr_leaf_index_free_index = amr_leaf_index_nrholes + 1
    else
        index = amr_leaf_index_holes[amr_leaf_index_nrholes]
        amr_leaf_index_nrholes = amr_leaf_index_nrholes - 1
    end

    if amr_index_to_leaf[index] isa AMRBranch
        throw(ArgumentError("Leaf index already in use"))
    end

    amr_index_to_leaf[index] = a
end

function amr_find_and_link_all_neighbors_branch(a::Union{AMRBranch,Nothing})
    if a.child isa AMRBranch
        throw(ArgumentError("Cannot find neighbors for a child branch"))
    end

    for idir = 1:3
        if amr_xyzdim[idir] > 0
            neigh = amr_find_neighbor_branch()
            a.neighbor[1, idir] = neigh
            if neigh isa AMRBranch
                amr_link_neighbors_back()
            end
            neigh = amr_find_neighbor_branch()
            a.neighbor[2, idir] = neigh
            if neigh isa AMRBranch
                amr_link_neighbors_back()
            end
        end
    end
end

function amr_find_neighbor_branch()
end

function amr_link_neighbors_back()
end



# Define a monotonically increasing xi, yi, zi
xi = collect(range(0, stop=1, length=6))
yi = collect(range(0, stop=1, length=6))
zi = collect(range(0, stop=1, length=6))

amr_initialize(
    true,
    true,
    true,
    5, 5, 5,
    xi, yi, zi,
    5,
    5,
    5,
    true,
    true,
    true,
    true,
    5,
)