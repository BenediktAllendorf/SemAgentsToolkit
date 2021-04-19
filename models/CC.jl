using SemAgentsToolkit, DataFrames

function climate_sem()
    constructs = ["K", "TBS", "RB", "TB", "AM", "PMO", "SN", "BI"]

    coef = DataFrame(zero(Array{Float64}(undef, length(constructs), length(constructs))),
    constructs)

    coef[index_of(constructs, :K), :RB] = .262
    coef[index_of(constructs, :K), :TB] = .203
    coef[index_of(constructs, :TBS), :RB] = .351
    coef[index_of(constructs, :TBS), :TB] = .298
    coef[index_of(constructs, :TB), :AM] = -.370
    coef[index_of(constructs, :TB), :BI] = .108
    coef[index_of(constructs, :RB), :AM] = -.145
    coef[index_of(constructs, :RB), :BI] = .322
    coef[index_of(constructs, :AM), :BI] = -.021
    coef[index_of(constructs, :SN), :BI] = .056
    coef[index_of(constructs, :PMO), :BI] = .280

    r² = DataFrame(zero(Array{Float64}(undef, 2, length(constructs))), constructs)
    r²[1, :BI] = .387
    r²[1, :AM] = .195
    r²[1, :RB] = .123
    r²[1, :TB] = .262

    mySem = create_sem(coef, r²)
end