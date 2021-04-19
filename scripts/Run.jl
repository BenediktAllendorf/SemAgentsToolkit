using DrWatson

@quickactivate()

using SemAgentsToolkit, Statistics, ProgressBars

include(joinpath(@__DIR__, "../models/CO.jl"))
include(joinpath(@__DIR__, "../models/CC.jl"))

if length(ARGS) == 1 && ARGS[1] == "CC"
    model = Dict("CC" => climate_sem())
elseif length(ARGS) == 1 && ARGS[1] == "CO" 
    model = Dict("CO" => covid_sem())
else
    print("Model not found!")
    exit()
end

for (model_name, mySem) in model
    allparams = Dict(
        :f => [minimum, mean, median, maximum],
        :jitter => [.0, .1, .5],
        :from => :BI,
        :to => [:SN, :PMO, :TB, :RB],
        :weight => [-5, -.1, -.01, .01, .1, .5],
        :steps => 100,
        :runs => 3,
        :zscorception => [false, true],
        :dims => (25, 25)
    )

    dicts = dict_list(allparams)

    for dict in ProgressBar(dicts)
        c = Configuration(;dict...)

        myISem = InterconnectedSEM(mySem, c)
     
        for data in simulate(myISem, c, [:SN, :PMO, :TB, :RB, :BI])
            safesave(datadir(model_name, savename(forfilename(c), "csv")), data)
        end
    end
end

