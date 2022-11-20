module WeightedSamplers

export WeightedSampler, adjust_weight!, shorten!

import Random
import StatsBase: sample, sample!

mutable struct WeightedSampler{V,T<:Number}
    v::V
    last::Int
    sum::T
    WeightedSampler(v) = new{typeof(v), eltype(v)}(v./sum(v), length(v), 1.0)
end

sample(ws::WeightedSampler) = sample(Random.GLOBAL_RNG, ws)
function sample(rng, ws::WeightedSampler)
    u = rand(rng)*ws.sum
    i = 1
    s = ws.v[i]
    while s<u
        i += 1
        s += ws.v[i]
    end
    w = ws.v[i]
    adjust_weight!(ws, i, 0)
    
    v = rand(rng)*ws.sum
    j = 1
    s = ws.v[j]
    while s<v
        j += 1
        s += ws.v[j]
    end
    adjust_weight!(ws, i, w)

    return i<j ? (i,j) : (j,i)
end

function adjust_weight!(ws::WeightedSampler, j, w)
    ws.sum = ws.sum - ws.v[j] + w
    ws.v[j] = w
end

function shorten!(ws::WeightedSampler, k=1)
    ws.sum -= sum(ws.v[ws.last-k+1:ws.last])
    ws.last -= k
end

end
