module TreeProcesses

export  ACD!, moran, moran2, birthdeath, weighted_coalescent, coalescent,
        maximally_balanced, maximally_unbalanced, nichemodel, to_mean_dataframe, yule

import AbstractTrees
using AbstractTrees: children, descendleft, Leaves, nodevalue, parent, print_tree, PreOrderDFS, PostOrderDFS
import Base: empty!
using BinaryTrees
using DataFrames
using DataStructures: PriorityQueue, dequeue!, enqueue!
import Statistics: mean
import StatsBase: sample
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
    traverse_left_right!(P::BinaryTree, agg!; init!)

Traverse tree "in-order". Apply function `agg!` to a node whenever
it is ascended to from the right. Apply function `init!` to leafs. Default is to
set each element of the leaf nodes' value field to `1`. 

Return number of nodes in tree.
"""
function traverse_left_right!(P::BinaryTree, agg!; init! = P->P.val.=1)
    ## Setup
    # 1. find left-most leaf
    cursor = descendleft(P)
    # 2. Set A=1
    init!(cursor)
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
            agg!(p)
            cursor = p
            continue
        end
        cursor = p
        # 4. traverse to the leftmost leaf of the right subtree if it exists, and start over
        if !isnothing(cursor.right)
            cursor = descendleft(cursor.right)
            init!(cursor)
        end
    end

    return k
end

## ---- Specific observables ---- ##

function ACD(p)
    A = 1
    C = 1
    D = 0
    if !isnothing(p.left)
        A += p.left.val[1]
        C += p.left.val[2]
        D += p.left.val[3]
    end
    if !isnothing(p.right)
        A += p.right.val[1]
        C += p.right.val[2]
        D += p.right.val[3]
    end
    # D = sum(x->x.val[2], children(p)) + A
    A, C+A, D+A-1
end

function ACD!(P::BinaryTree, slots=1:3)
    A = Int[]
    C = Int[]
    D = Int[]
    function agg(P)
        P.val[slots] .= ACD(P)
        push!(A, P.val[1])
        push!(C, P.val[2])
        push!(D, P.val[3])
    end
    function init_leaf!(P)
        P.val[slots] .= (1, 1, 0)
        push!(A, 1)
        push!(C, 1)
        push!(D, 0)
    end

    k = traverse_left_right!(P, agg; init! = init_leaf!)
    return k, A, C, D
end


## ---- Processes --- ##

"""
Simulate the coalescent process for n genes.

Return a directed graph with edges pointing towards the root.
"""
function coalescent(n; default_value=()->[0,0])
    P = [BinaryTree(default_value()) for _ in 1:n] # start with n vertices; store A,C
    while n > 1
        i, j = sample(1:n, 2, replace=false, ordered=true) #sample two distinct vertices to coalesce
        @inbounds l, r = P[i], P[j]
        v = BinaryTree(default_value()) # newly added parental node
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
function yule(n; default_value=()->[0,0])
    P = [BinaryTree(default_value())]
    k = 1
    while length(P) - k + 1 < n
        i = rand(k:lastindex(P))
        v = P[i]
        l = left!(v, default_value())
        r = right!(v, default_value())
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
function birthdeath(n, T, d, b=1.0; N=0, default_value=()->[0,0])
    P = [BinaryTree(default_value()) for _ in 1:n]
    t = 1
    N > 0 && sizehint!(P, N)
    k = 1
    while t <= T && (N == 0 || length(P) - k + 1 < N)
        k > length(P) && break
        if rand() < b
            i = rand(k:lastindex(P))
            v = P[i]
            l = left!(v, default_value())
            r = right!(v, default_value())
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
function maximally_unbalanced(height; default_value=()->[0,0])
    P = BinaryTree(default_value())
    h = 0
    attach_to = P
    while h < height
        right!(attach_to, default_value())
        attach_to = left!(attach_to, default_value())
        h += 1
    end

    return P
end

"""
  Return a fully balanced binary tree of given height.
"""
function maximally_balanced(height; default_value=()->[0,0])
    P = [BinaryTree(default_value())]
    root = P[1]
    h = 1
    n = 1
    while h < height
        attach_to = popfirst!(P)
        push!(P, left!(attach_to, default_value()))
        push!(P, right!(attach_to, default_value()))

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
function weighted_coalescent(n, w=randn(n).^2; default_value::F, fuse=max) where F
    P = [BinaryTree(default_value()) for _ in 1:n]
    ws = WeightedSampler(w)
    d = ws.d
    while n  > 1
        i, j = ws_sample(ws, 2; ordered=true)
        l, r = P[i], P[j]
        v = BinaryTree(default_value())
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

## -- Niche model from Goldenfeld (2020) -- ##

e(r, R0) = r/(r+R0) # extinction probability
r(n, ϵ=0.0) = ifelse(n>0, n, ϵ) # speciation rate

"""
    nichemodel(n, σ, R0=10.0; n0=1.0)

Minimal model of a speciation process coupled to fluctuating ecological niches.

Return root node and total number of nodes in the final tree. The latter can be used
to check if the tree generation process completed, or terminated prematurely.

Reference: https://doi.org/10.1073/pnas.1915088117
"""
function nichemodel(n, σ, R0=10.0; n0=1.0, ϵ=0.0, default_value=()->Float64[0.0, 0.0])
    _default_value = () -> vcat(default_value(), [n0, 0.0])
    clock = 0.0
    P = [BinaryTree(_default_value())] # root node
    nnodes = 1
    spec_queue = PriorityQueue(Base.Order.Forward, P[1] => 0.0)
    while nnodes < n && !isempty(spec_queue)
        # Sample a node, speciate, and calculate niche sizes, spec. rates, ext. prob. 
        # of potential child nodes.
        # Advance clock to spec. event.
        p = dequeue!(spec_queue)
        clock = p.val[end]
        np = p.val[end-1]
        nl, nr = np .+ np*σ*randn(2)
        rl, rr = r(nl, ϵ), r(nr, ϵ)
        tl = clock - log(rand())/rl
        tr = clock - log(rand())/rr
        el, er = e(rl, R0), e(rr, R0)

        # for each child roll dice against ext. prob. 
        # if it passes, attach to parent node
        n_children = 0
        if rand() > el
            n_children += 1
            nnodes += 1
            left = child!(p, _default_value(), :left)
            left.val[end-1] = nl
            left.val[end] = tl
            !isinf(tl) && enqueue!(spec_queue, left => tl)
        end
        if nnodes<n && rand() > er
            n_children += 1
            location = n_children == 2 ? :right : :left
            nnodes += 1
            right = child!(p, _default_value(), location)
            right.val[end-1] = nr
            right.val[end] = tr
            !isinf(tr) && enqueue!(spec_queue, right => tr)
        end

        # pruning
        pa = parent(p)

        if n_children == 0
            isnothing(pa) && continue
            sib = sibling(p)
            papa = parent(pa)
            if isnothing(papa)
                # pa was the root node
                # replace with sibling
                pa.left = nothing
                pa.right = nothing
                sib.parent = nothing
                P[1] = sib
                continue
            end
            if isleftchild(pa)
                papa.left = sib
            else
                papa.right = sib
            end
            sib.parent = papa
            pa.left = nothing
            pa.right = nothing
            pa.parent = nothing
            # nnodes -= 2
        end

        if n_children == 1
            chld = p.left
            if isnothing(pa)
                P[1] = chld
                chld.parent = nothing
            end
            chld.parent = pa
            if isleftchild(p)
                pa.left = chld
            elseif isrightchild(p)
                pa.right = chld
                chld.parent = pa
            end
            p.left = p.right = p.parent = nothing
            # nnodes -= 1
        end
    end
    # returning the number of nodes serves as a check whether the tree was constructed fully,
    # or the process halted prematurely.
    return P[1], nnodes
end

end # MODULE
