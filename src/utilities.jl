import Statistics

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
