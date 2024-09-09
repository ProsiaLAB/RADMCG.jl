mutable struct AMRBranch
    level::Int
    id::Int
    branch_index::Int
    leaf_index::Int
    leaf::Bool
    neighbor::Array{Int,2}
    child::Array{Int,3}
    parent::AMRBranch
    parent_slot::Array{Int,1}
    ixyz::Array{Int,1}
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
            amr_branch_index_holes = zeros(
                Int,
                nr_branches_max,
            )
            amr_leaf_index_holes = zeros(
                Int,
                nr_branches_max,
            )
            amr_branchindex_nrholes = 0
            amr_leafindex_nrholes = 0
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

    println(amr_finegrid_xi, amr_finegrid_xc)


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