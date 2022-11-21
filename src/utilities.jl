"""
  Combine vectors with values for A and C into a dataframe, averaging C/A for every A.
"""
function to_mean_dataframe(A, C)
    df = DataFrame(; A, C)
    combine(groupby(df, :A), :A => mean => :a, [:A, :C] => ((a, c) -> mean(c ./ a)) => :covera, nrow)
end
