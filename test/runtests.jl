using Test
using TreeProcesses

@testset begin
    d = 8
    tree = weighted_coalescent(2^d)
    A, C = treevalues!(tree)
    @test A[1] == 2^(d+1)-1
end
