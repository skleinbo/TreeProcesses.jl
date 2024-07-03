# TreeProcesses.jl

Some (common) forward and backward processes to generate binary trees as seen
for example in population genetics.

Based on [`BinaryTrees.jl`](https://github.com/skleinbo/BinaryTrees.jl), which in turn implements the interface from [`AbstractTrees.jl`](https://github.com/JuliaCollections/AbstractTrees.jl/)

## Installation

The package is unregistered. Install it directly from GitHub.

```
julia> ]add https://github.com/skleinbo/BinaryTrees.jl https://github.com/skleinbo/TreeProcesses.jl
```

## Currently implemented

__Forward:__

* `birth_death(n, T, d, b=1.0; N=0)`: Starting from `n` nodes, a randomly selected one splits with probability `b`. Independently, a node dies with probability `d`. Run for `T` time steps (a birth plus death event are one time step), or until `N` nodes are present unless `N==0` (default). Lineages that die out are automatically pruned, i.e. after sufficiently many (`O(n^2)`) time steps, a most recent common ancestor will be found. Returns all surviving root nodes.
* `maximally_balanced(n)`
* `maximally_unbalanced(n)`
* `moran(n,T)`: `birth_death` with b/d probabilities both equal to `1` to maintain constant population size.
* `yule(n)`: Starting from a root node, split leafs uniformly `n` times.

__Backward:__

* `coalescent(n)`: Neutral coalescent process. Starting from `n` nodes, pick two at random and merge until only one is left.
* `preferential_coalescent(n, w; fuse=max)`: Like `coalescent`, but each node is assigned the respective weight from the vector `w` with which it coalesces.
Ancestral nodes' weights are computed from their children's weights by applying `fuse`. The standard coalescent is recovered by choosing all weights equal and setting `fuse=first`.

## Node types

Every node is of type `BinaryTree{T}` where its internal state is stored in a variable of type `T`. The state is accessible through the `val` field. The appropriate `T` depends on the process and the observables one wishes to store. Required fields are detailed in the process's docstring. Sensible defaults are provided.

All included processes take a keyword argument `nodevalue` which is a callable that provides the default value when a node is first created.

## Observables

* `ACD!(t)`: Annotate each node of a tree with
  * $A$: Number of nodes in the subtree, including itself.
  * $C$: Cumulative number of nodes in the subtree, i.e. $\sum A(i)$ for $i$ in the subtree, including the node itself.  
  * $D$: Cumulative distance to all nodes in the subtree.
  
  Returns number of nodes, and vectors of `A`,`C` and `D`. The latter are in [post-order](https://en.wikipedia.org/wiki/Tree_traversal).

  __Note__: The node storage `T` (see above) must have a field `observables<:AbstractVector` with at least three elements. The values `A,C,D` are stored in the first three elements.

## Example

```
julia> using TreeProcesses

julia> T = preferential_coalescent(2^4, rand(2^4)) |> first
((w=0.9772233107834184, t=15, observables=[0, 0, 0])) 983bc8d17e328bad with 2 children and no parent.

# calculate and return number of nodes, A & C
julia> ACD!(T)
(31, [1, 1, 3, 1, 5, 1, 7, 1, 1, 3  …  9, 23, 1, 1, 3, 27, 1, 29, 1, 31], [1, 1, 5, 1, 11, 1, 19, 1, 1, 5  …  25, 91, 1, 1, 5, 123, 1, 153, 1, 185], [0, 0, 2, 0, 6, 0, 12, 0, 0, 2  …  16, 68, 0, 0, 2, 96, 0, 124, 0, 154])

# the tree has been annotated with the observables
julia> T
((w=0.9772233107834184, t=15, observables=[31, 185, 154])) 983bc8d17e328bad with 2 children and no parent

julia> using AbstractTrees

julia> print_tree(T, maxdepth=16)
(w=0.977223, t=15, observables=[31, 185, 154])
├─ (w=0.977223, t=14, observables=[29, 153, 124])
│  ├─ (w=0.977223, t=13, observables=[27, 123, 96])
│  │  ├─ (w=0.977223, t=11, observables=[23, 91, 68])
│  │  │  ├─ (w=0.977223, t=9, observables=[13, 43, 30])
│  │  │  │  ├─ (w=0.977223, t=7, observables=[7, 19, 12])
│  │  │  │  │  ├─ (w=0.914926, t=5, observables=[5, 11, 6])
│  │  │  │  │  │  ├─ (w=0.914926, t=3, observables=[3, 5, 2])
│  │  │  │  │  │  │  ├─ (w=0.876, t=0, observables=[1, 1, 0])
│  │  │  │  │  │  │  └─ (w=0.914926, t=0, observables=[1, 1, 0])
│  │  │  │  │  │  └─ (w=0.397427, t=0, observables=[1, 1, 0])
│  │  │  │  │  └─ (w=0.977223, t=0, observables=[1, 1, 0])
│  │  │  │  └─ (w=0.333037, t=4, observables=[5, 11, 6])
│  │  │  │     ├─ (w=0.333037, t=1, observables=[3, 5, 2])
│  │  │  │     │  ├─ (w=0.333037, t=0, observables=[1, 1, 0])
│  │  │  │     │  └─ (w=0.215572, t=0, observables=[1, 1, 0])
│  │  │  │     └─ (w=0.1671, t=0, observables=[1, 1, 0])
│  │  │  └─ (w=0.846859, t=10, observables=[9, 25, 16])
│  │  │     ├─ (w=0.495586, t=6, observables=[3, 5, 2])
│  │  │     │  ├─ (w=0.495586, t=0, observables=[1, 1, 0])
│  │  │     │  └─ (w=0.19761, t=0, observables=[1, 1, 0])
│  │  │     └─ (w=0.846859, t=8, observables=[5, 11, 6])
│  │  │        ├─ (w=0.846859, t=2, observables=[3, 5, 2])
│  │  │        │  ├─ (w=0.122918, t=0, observables=[1, 1, 0])
│  │  │        │  └─ (w=0.846859, t=0, observables=[1, 1, 0])
│  │  │        └─ (w=0.321122, t=0, observables=[1, 1, 0])
│  │  └─ (w=0.549756, t=12, observables=[3, 5, 2])
│  │     ├─ (w=0.216285, t=0, observables=[1, 1, 0])
│  │     └─ (w=0.549756, t=0, observables=[1, 1, 0])
│  └─ (w=0.282131, t=0, observables=[1, 1, 0])
└─ (w=0.219596, t=0, observables=[1, 1, 0])
```
