mutable struct TreeNode
    value::Int
    children::Vector{TreeNode}

    function TreeNode(value::Int)
        new(value, TreeNode[])
    end
end

function add_child(parent::TreeNode, child::TreeNode)
    push!(parent.children, child)
end

node = TreeNode(1)

add_child(node, TreeNode(2))
add_child(node, TreeNode(3))

println(node)