# TreeProcesses.jl

Some (common) forward and backward processes to generate binary trees as seen
for example in population genetics.

Based on [`BinaryTrees.jl`](https://github.com/skleinbo/BinaryTrees.jl), which in turn implements the interface from [`AbstractTrees.jl`](https://github.com/JuliaCollections/AbstractTrees.jl/)

## Installation

The package is unregistered, therefore install it directly from GitHub.

```julia
julia> ]add https://github.com/skleinbo/BinaryTrees.jl https://github.com/skleinbo/TreeProcesses.jl
```

## Currently implemented

__Forward:__

* `birth_death(n, T, d, b=1.0; N=0)`: Starting from `n` nodes, a randomly selected one splits with probability `b`. Independently, a node dies with probability `d`. Run for `T` timesteps (a birth plus death event are one timestep), or until `N` nodes are present unless `N==0` (default). Lineages that die out are automatically pruned, i.e. after sufficiently (`O(n^2)`) time steps, a most recent common ancestor will be found. Returns all surviving root nodes.
* `maximally_balanced(n)`
* `maximally_unbalanced(n)`
* `moran(n,T)`: `birth_death` with b/d probabilities both equal to `1` to maintain constant population size.
* `yule(n)`: Starting from a root node, split leaves uniformly `n` times.

__Backward:__

* `coalescent(n)`: Neutral coalescent process. Starting from `n` nodes, pick two at random and merge until only one is left.
* `fluctuating_coalescent(n, w; fuse=max)`: Like `coalescent`, but each node is assigned the respective weight from the vector `w` with which it partakes in a coalescent event.
Ancestral nodes' weights are computed from their children's weights by applying `fuse`. The standard coalescent is recovered by choosing all weights equal and setting `fuse=first`.

## Observables

* `C!(t)`: Annotate each node of a tree with 
  * $A$: Number of nodes in the subtree, including itself.
  * $C$: Cumulative number of nodes in the subtree, i.e. $\sum A(i)$ for $i$ in the subtree, including the node itself.  
  
  __Note__: The `val` field of the nodes must be a `Vector` that can store `Int`s, which is the default.


## Example

```julia
julia> using TreeProcesses

julia> T = fluctuating_coalescent(2^4, rand(2^4))
BinaryTree{Vector{Int64}}([0, 0]) 11038476695863080329 with 2 children and no parent.

# calculate and return A & C
julia> treevalues!(T)
([31, 29, 27, 23, 21, 7, 3, 1, 1, 3  …  3, 1, 1, 1, 1, 3, 1, 1, 1, 1], [203, 171, 141, 109, 85, 17, 5, 1, 1, 5  …  5, 1, 1, 1, 1, 5, 1, 1, 1, 1])

# the tree has been annotated with the observable values
julia> T
BinaryTree{Vector{Int64}}([31, 203]) 11038476695863080329 with 2 children and no parent.
```
