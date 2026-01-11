import CategoricalArrays: cut
using DataStructures: Queue, dequeue!, enqueue!
import Statistics
using Intervals

"""
  Combine vectors with values for A,C and D into a dataframe,
  averaging C/A, D/A for every unique A.
"""
function to_mean_dataframe(A, C, D)
    df = DataFrame(; A, C, D)
    df.covera = df.C ./ df.A
    df.dovera = df.D ./ df.A
    combine(groupby(df, :A),
      :covera => Statistics.mean => :covera,
      :dovera => Statistics.mean => :dovera,
      nrow
    )
end
function to_mean_dataframe(;gb, obs_f...)
    vals_obs_f = collect(values(obs_f))
    keys_obs_f = keys(obs_f)
    obs = first.(vals_obs_f)
    fs = last.(vals_obs_f)
    @show keys_obs_f, fs
    df = DataFrame(; zip(keys_obs_f, obs)...)
    # df.covera = df.C ./ df.A
    # df.dovera = df.D ./ df.A
    combine(groupby(df, gb),
      Not(gb) .=> fs[2:end] .=> Not(gb),
        # Statistics.mean,
      nrow
    )
end
function to_binned_dataframe(;gb, bins, obs_f...)
    vals_obs_f = collect(values(obs_f))
    keys_obs_f = keys(obs_f)
    obs = first.(vals_obs_f)
    fs = last.(vals_obs_f)
    # @show keys_obs_f, fs
    df = DataFrame(; zip(keys_obs_f, obs)...)
    # subset!(df, gb => (x->x.<=last(bins)))

    binned_col = Symbol("$(gb)_binned")
    int_col = Symbol("$(gb)_int")

    # fmt(l,r,i; leftclosed, rightclosed) = Interval{Closed, Open}(l,r)
    df[!, binned_col] .= cut(df[!, gb], bins, extend=true)
    # df[!, int_col] = parse.(Interval{Float64}, string.(df[:, binned_col]))
    df = combine(groupby(df, binned_col),
      Not(binned_col) .=> fs .=> Not(binned_col),
      binned_col => (x->parse.(Interval{Float64}, string.(first(x)))) => int_col,
      # int_col => first => int_col,
      nrow
    )
    select!(df, Not(binned_col))
    df
end

function time_slice(P::BinaryTree{T}, t) where T
    nodes = BinaryTree{T}[]
    Q = Queue{BinaryTree{T}}()
    enqueue!(Q, P)
    while !isempty(Q)
        P = dequeue!(Q)
        if P.val.t > t
            isnothing(P.left) || enqueue!(Q, P.left)
            isnothing(P.right) || enqueue!(Q, P.right)
        else
            push!(nodes, P)
        end
    end
    return nodes
end

function get_observables(P)
    O = Matrix{eltype(P[1].val.observables)}(undef, length(P), size(P[1].val.observables, 1))
    for i in eachindex(P)
        O[i, :] .= P[i].val.observables
    end
    return O
end