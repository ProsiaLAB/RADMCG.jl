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

amr = AMRBranch(1, 1, 1, 1, true, [0 0; 0 0], [0 0 0; 0 0 0; 0 0 0], C_NULL, [0, 0, 0], [0, 0, 0])