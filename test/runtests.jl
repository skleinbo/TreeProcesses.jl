using BinaryTrees
using Test
using TreeProcesses

@testset begin
    d = 8
    tree = preferential_coalescent(2^d)
    k, A, C, D = ACD!(tree[1])
    @test A[end] == 2^(d+1)-1

    @test coalescent(d) isa BinaryTree
    @test yule(d) isa BinaryTree
    @test nichemodel(d, 2.0)[1] isa BinaryTree
end
