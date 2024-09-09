mutable struct AMRBranch
    level::Int
    id::Int
    branch_index::Int
    leaf_index::Int
    leaf::Bool
    neighbor::Array{Int,2}
    child::Array{Int,3}
    parent::Ptr{AMRBranch}
    parent_slot::Array{Int,1}
    ixyz::Array{Int,1}
end


mutable struct AMRBranchLink
    link::Ptr{AMRBranch}
end


function amr_initialize(
    include_x::Bool,
    include_y::Bool,
    include_z::Bool,
    nx::Int, ny::Int, nz::Int,
    xi::Float64, yi::Float64, zi::Float64,
    level_max::Int,
    fill_base_grid::Int, nr_leaves_max::Int,
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
            amr_index_to_branch = zeros(
                AMRBranchLink,
                nr_branches_max,
            )
            amr_index_to_leaf = zeros(
                AMRBranchLink,
                nr_branches_max,
            )
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
            println(amr_index_to_branch)
            # for i in 1:nr_branches_max
            #     amr_index_to_branch[i].link = Nothing
            #     amr_index_to_leaf[i].link = Nothing
            # end
        end


    end
end

amr_initialize(
    true,
    true,
    true,
    5, 5, 5,
    20.0, 20.0, 20.0,
    100,
    5,
    100,
    100,
    true,
    true,
    true,
    true,
)