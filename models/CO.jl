using SemAgentsToolkit, DataFrames

function covid_sem()
    constructs = ["TBS", "RB", "TB", "AM", "PMO", "SN", "BI"]

    coef = DataFrame(zero(Array{Float64}(undef, length(constructs), length(constructs))),
    constructs)

    coef[index_of(constructs, :TBS), :RB] = .552
    coef[index_of(constructs, :TBS), :TB] = .202
    coef[index_of(constructs, :RB), :BI] = .564
    coef[index_of(constructs, :RB), :AM] = -.090
    coef[index_of(constructs, :TB), :AM] = -.472
    coef[index_of(constructs, :TB), :BI] = .049
    coef[index_of(constructs, :PMO), :BI] = .202
    coef[index_of(constructs, :AM), :BI] = -.073
    coef[index_of(constructs, :SN), :BI] = -.037

    r² = DataFrame(zero(Array{Float64}(undef, 2, length(constructs))), constructs)
    r²[1, :BI] = .516
    r²[1, :AM] = .254
    r²[1, :RB] = .300
    r²[1, :TB] = .043

    return create_sem(coef, r²)
end