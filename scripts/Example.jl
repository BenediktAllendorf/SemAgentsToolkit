using DrWatson

@quickactivate()

using DataFrames, SemAgentsToolkit

constructs = ["A", "B", "C"]

coef = DataFrame(zeros(Float64, length(constructs), length(constructs)),
                 constructs)

r² = DataFrame(zeros(Float64, 1, length(constructs)), constructs)

coef[index_of(constructs, :A), :C] = .123
coef[index_of(constructs, :B), :C] = .456

r²[1, :C] = .789

mySem = create_sem(coef, r²)

allparams = Dict(
    :f => [minimum, maximum],
    :jitter => [.0, .1, .5],
    :from => [:C],
    :to => [:A],
    :weight => [.1, .5],
    :steps => 100,
    :runs => 3,
    :zscorception => [false, true],
    :dims => (25, 25)
)

dicts = dict_list(allparams)

for dict in dicts
    c = Configuration(;dict...)

    myISem = InterconnectedSEM(mySem, c)

    for data in simulate(myISem, c, [:A, :C])
        safesave(datadir("ttt", savename(forfilename(c), "csv")), data)
    end
end