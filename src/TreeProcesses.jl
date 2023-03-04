module TreeProcesses

export A!, C!, treevalues!, moran, moran2, birthdeath, weighted_coalescent, coalescent,
        maximally_balanced, maximally_unbalanced, to_mean_dataframe, yule

import AbstractTrees
import AbstractTrees: children, descendleft, Leaves, nodevalue, parent, print_tree, PreOrderDFS, PostOrderDFS
import Base: empty!
using BinaryTrees
using DataFrames
import StatsBase: mean, sample
import WeightedSampling: WeightedSampler, adjust_weight!, sample as ws_sample, weight

include("utilities.jl")

"Distance between a node `i` and its ancestor `j` on a tree."
function dist(i, j)
    d = 0
    while i != j
        i = parent(i)
        isnothing(i) && throw(ErrorException("$i is not a descendent of $j."))
        d += 1
    end
    return d
end

## Phylogenetic Observables

function empty!(P::BinaryTree{Vector{T}}) where {T}
    for v::BinaryTree{Vector{T}} in PreOrderDFS(P)
        v.val::Vector{T} .= zero(T)
    end
    nothing
end

"""
    traverse_left_right!(P::BinaryTree{Vector}, i, agg)

Traverse tree "in-order". Apply function `agg` to a node whenever
it is ascended to from the right. Store result in i-th component of 
that nodes value field.

Return number of nodes in tree.
"""
function traverse_left_right!(P::BinaryTree{Vector{T}}, i, agg) where T
    ## Setup
    # 1. find left-most leaf
    cursor = descendleft(P)
    # 2. Set A=1
    cursor.val[i] = 1
    ## Start loop
    ## count nodes on the way
    k = 1
    while cursor!==P
        # 2. move one up if possible,
        #    if not return P
        p = parent(cursor)
        isnothing(p) && break
        k += 1
        # 3. If cursor is a right child, calculate A for the parent node
        #    and move one level up (continue)
        if isrightchild(cursor) || length(children(p))==1
            p.val[i] = agg(p)
            cursor = p
            continue
        end
        cursor = p
        # 4. traverse to the leftmost leaf of the right subtree if it exists, and start over
        if !isnothing(cursor.right)
            cursor = descendleft(cursor.right)
            cursor.val[i] = 1
        end
    end

    return k
end

A!(P::BinaryTree{Vector{T}}) where T = traverse_left_right!(P, 1, p->sum(x->x.val[1], children(p)) + 1)
function C!(P::BinaryTree{Vector{T}}) where T
    A!(P)
    k = traverse_left_right!(P, 2, p->sum(x->x.val[2], children(p)) + p.val[1])
    return k
end

function treevalues!(P::BinaryTree{Vector{T}}) where {T}
    k = C!(P)
    A = Vector{Int}(undef, k)
    C = Vector{Int}(undef, k)
    for (i,x) in enumerate(PreOrderDFS(P))
        @inbounds A[i], C[i] = nodevalue(x)
    end
    return A,C
end

"""
Simulate the coalescent process for n genes.

Return a directed graph with edges pointing towards the root.
"""
function coalescent(n; default_value=[0,0])
    P = [BinaryTree(copy(default_value)) for _ in 1:n] # start with n vertices; store A,C
    while n > 1
        i, j = sample(1:n, 2, replace=false, ordered=true) #sample two distinct vertices to coalesce
        @inbounds l, r = P[i], P[j]
        v = BinaryTree(copy(default_value)) # newly added parental node
        v.left = l # connect sampled nodes to parental node
        v.right = r
        l.parent = v
        r.parent = v

        # replace nodes i,j in active nodes list with
        # new node and current last node.
        # shorten list of admissible draws by one.
        P[i] = v
        P[j] = P[n]
        n -= 1
    end
    return P[1]
end

## Forwards Moran process ##

"""
  Simulate a forward Moran process.

  Start with a population 1...n.

  Return a population (i_1, i_2,...,i_n), and
  a list of moves (i,j) indicating individual at i 
  replaced the one at j.
"""
function moran_moves(n, T)
    t = 1
    moves = Vector{Tuple{Int,Int}}(undef, T)
    v = collect(1:n)
    while t <= T
        i, j = sample(eachindex(v), 2, replace=false)
        v[j] = v[i]
        moves[t] = (i, j)
        t += 1
    end

    return v, moves
end

"""
  Construct a Yule tree with `n` tips.

  At each time step a bifurcation event happens at any of the
  current tips with equal probability.
"""
function yule(n; default_value=[0,0])
    P = [BinaryTree(copy(default_value))]
    k = 1
    while length(P) - k + 1 < n
        i = rand(k:lastindex(P))
        v = P[i]
        l = left!(v, copy(default_value))
        r = right!(v, copy(default_value))
        append!(P, (l, r))
        P[i] = P[k]

        k += 1
    end
    P[1]
end

"""
  birthdeath(n, T, d, b=1.0; N=0)

  Start with a population of size `n`.
  At each timestep an individual duplicates with 
  probability `b`, and one dies with probability `d`.

  Run for `T` timesteps if `N=0`, or until population has
  reached size `N`; whatever happens first.

  Lineages that die out are pruned automatically.

  Return an ancestral tree, or a vector of trees if the process hasn't coalesced.
"""
function birthdeath(n, T, d, b=1.0; N=0, default_value=[0,0])
    P = [BinaryTree(copy(default_value)) for _ in 1:n]
    t = 1
    N > 0 && sizehint!(P, N)
    k = 1
    while t <= T && (N == 0 || length(P) - k + 1 < N)
        k > length(P) && break
        if rand() < b
            i = rand(k:lastindex(P))
            v = P[i]
            l = left!(v, copy(default_value))
            r = right!(v, copy(default_value))
            push!(P, r)
            P[i] = l
        end
        if rand() < d # die
            j = rand(k:lastindex(P))
            w = P[j]
            P[j] = P[k]
            k += 1

            p = parent(w)
            if !isnothing(p)
                pp = parent(p)
                sib = sibling(w)
                if !isnothing(pp)
                    if pp.left === p
                        pp.left = sib
                    else
                        pp.right = sib
                    end
                    sib.parent = pp
                else
                    sib.parent = nothing
                end
            end
        end
        t += 1
    end

    return union(AbstractTrees.getroot.(P[k:end]))
end

"""
  Directly construct an ancestral tree for a Moran process of
    `n` individuals.

  Lineages that die out are pruned automatically.

  See also: @ref(`birth_death`)
"""
moran(n, T; kwargs...) = birthdeath(n, T, 1.0; kwargs...)

"""
  Return a fully imbalanced binary tree of given height.
"""
function maximally_unbalanced(height; default_value=[0,0])
    P = BinaryTree(copy(default_value))
    h = 0
    attach_to = P
    while h < height
        right!(attach_to, copy(default_value))
        attach_to = left!(attach_to, copy(default_value))
        h += 1
    end

    return P
end

"""
  Return a fully balanced binary tree of given height.
"""
function maximally_balanced(height; default_value=[0,0])
    P = [BinaryTree(copy(default_value))]
    root = P[1]
    h = 1
    n = 1
    while h < height
        attach_to = popfirst!(P)
        push!(P, left!(attach_to, copy(default_value)))
        push!(P, right!(attach_to, copy(default_value)))

        n += 2
        h = log2(n + 1)
    end
    return root
end

"""
    weighted_coalescent(n, w=randn(n).^2; default_value=[0, 0], fuse=max)

Simulate a coalescent process for `n` genes, in which coalescing nodes are not selected uniformly, but according to the weight vector `w`.

Starting nodes have recombination rates `w`, which upon coalescence are combined via `fuse`.
Node values are set to `default_value`.

Return a binary tree.
"""
function weighted_coalescent(n, w=randn(n).^2; default_value=[0, 0], fuse=max)
    P = [BinaryTree(copy(default_value)) for _ in 1:n]
    ws = WeightedSampler(w)
    d = ws.d
    while n  > 1
        i, j = ws_sample(ws, 2; ordered=true)
        l, r = P[i], P[j]
        v = BinaryTree(copy(default_value))
        v.left = l
        v.right = r
        l.parent = v
        r.parent = v

        P[i] = v
        adjust_weight!(ws, i, fuse(ws.heap[d+i], ws.heap[d+j]))
        adjust_weight!(ws, j, 0.0)
        n -= 1
    end

    return P[1]
end

end # MODULE
