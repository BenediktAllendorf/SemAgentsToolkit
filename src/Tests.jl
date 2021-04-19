using SafeTestsets, Test

@safetestset "Basic SEM" begin
    using SemAgentsToolkit, DataFrames, LightGraphs, MetaGraphs
    coef = DataFrame(A = [.0, 0], B = [.5, 0])
    sem = create_sem(coef)

    @test length(edges(sem.g)) == 1
    @test first(edges(sem.g)) == Edge(1, 2)
    @test get_prop(sem.g, Edge(1,2), :weight) == 0.5

    @test sem.g[:A, :name] == 1
    @test sem.g[:B, :name] == 2
    @test sem.g[1, :name] == :A
    @test sem.g[2, :name] == :B
end

@safetestset "Basic InterconnectedSem" begin
    using SemAgentsToolkit, DataFrames, LightGraphs, MetaGraphs

    coef = DataFrame(A = [.0, 0], B = [.5, 0])
    sem = create_sem(coef)
    iSem = InterconnectedSEM(sem, :B, :A, .9)
    
    @test iSem.base == sem
    @test length(edges(iSem.interconnectivity.g)) == 4
    @test size(iSem.interconnectivity.coef) == (4, 4)
    @test get_prop(iSem.interconnectivity.g, Edge(3, 1), :weight) == .45
    @test get_prop(iSem.interconnectivity.g, Edge(3, 2), :weight) == .225
    @test get_prop(iSem.interconnectivity.g, Edge(4, 1), :weight) == .9
    @test get_prop(iSem.interconnectivity.g, Edge(4, 2), :weight) == .45
end

@safetestset "Basic HumanBeing" begin
    using SemAgentsToolkit, DataFrames, LightGraphs, MetaGraphs

    coef = DataFrame(A = [.0, 0], B = [.5, 0])
    sem = create_sem(coef)
    iSem = InterconnectedSEM(sem, :B, :A, .9)

    a = SemAgentsToolkit.HumanBeing(iSem)

    @test a.constructvalues == [.0, 0]
end

@safetestset "Basic HumanBeing with Neighbor" begin
    using SemAgentsToolkit, DataFrames, LightGraphs, MetaGraphs

    coef = DataFrame(A = [.0, 0], B = [.5, 0])
    sem = create_sem(coef)
    iSem = InterconnectedSEM(sem, :B, :A, .9)

    a = SemAgentsToolkit.HumanBeing(iSem)

    @test SemAgentsToolkit.change_construct_value_neighbor_based(a, SemAgentsToolkit.neighborize(:A, Symbol), .123) == [.45 * .123, .225 * .123]
    @test SemAgentsToolkit.change_construct_value_neighbor_based(a, SemAgentsToolkit.neighborize(:B, Symbol), .999) == [.9 * .999, .45 * .999]
end