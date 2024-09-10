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

    # Inner constructor with default values
    function AMRBranch(level::Int, id::Int; branch_index::Int=0, leaf_index::Int=0, leaf::Bool=false,
        parent::Union{AMRBranch,Nothing}=nothing, parent_slot::Array{Int,1}=Int[],
        ixyzf::Array{Int,1}=Int[])

        # Creating arrays with dimensions based on your specific use-case
        neighbor = Array{AMRBranch,2}(undef, 0, 0)  # Replace with correct dimensions
        child = Array{AMRBranch,3}(undef, 0, 0, 0)  # Replace with correct dimensions

        # Return an instance
        new(level, id, branch_index, leaf_index, leaf, neighbor, child, parent, parent_slot, ixyzf)
    end
end


# Creating an AMRBranch instance with default values for optional fields
branch = AMRBranch(1, 42)

# Creating an AMRBranch instance with custom values for optional fields
parent_slot = [1, 2, 3]
ixyzf = [0, 1, 0]
branch_with_custom_values = AMRBranch(2, 101; branch_index=5, leaf=true, parent_slot=parent_slot, ixyzf=ixyzf)
