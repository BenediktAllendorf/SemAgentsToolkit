export cacheplots, plot_zscores, plot_heatmap, plot_distributions, plot_column, generateNotebook, load, datapath, loadplot

function loadplot(name, step)
    if backend() isa Plots.HDF5Backend
        throw(DomainError("Select a different backend than HDF5"))
    end

    path = joinpath(plotspath(name, false))
    if !isdir(path)
        throw(DomainError("Nothing cached there"))
    end
    
    step = mod1(step, length(readdir(path)))
        
    Plots.hdf5plot_read(joinpath(path, "$step.hdf5"))
end

function generateBetterSelecter(params)

    sliders_str = []
    filename_str = []

    pad = maximum(length.(keys(params)))
    for name in keys(params)
        values = (x->collect(x))(params[name])
        push!(sliders_str, name *
            raw""": $(@bind""" * " $(varnamesPrefix())"
            * name * raw""" Select(""" * "$(sort((x->collect(x))(params["$(name)"])))" * "))\\\n")

        push!(filename_str, "$(name)=" * "\$($(varnamesPrefix())$name)")
    end

    raw"""md\"\"\" """ * join(sliders_str) * "\"\"\"", "selected_run = \"$(join(filename_str, '_'))\""
end

function get_params_from_filenames(path)
    params = OrderedDict{String, Set{String}}()

    for filename in sort(readdir(path))
        filename = replace(filename, r".csv$"=>"")
        attributes = split.(split(filename, '_'), '=')
        filter!(e -> length(e) == 2, attributes)
        for (x,y) in attributes
            if !(x in keys(params))
                params[x] = Set()
            end

            push!(params[x], [y][1])
        end
    end
    
    params
end

function load(selected_run)
    path = joinpath(datapath(), selected_run)
    groups = []
    if isdir(path)
        for f in filter(x->endswith(x, ".csv"), readdir(path))
            X =  CSV.File(joinpath(path, f)) |> DataFrame
            
            push!(groups, groupby(X, :step, sort = true))
        end
    else
        spl = splitpath(path)
        parent = joinpath(spl[1:(length(spl) - 1)]...)

        for f in filter(x->startswith(last(splitpath(x)), last(spl)), readdir(parent))
            X =  CSV.File(joinpath(parent, f)) |> DataFrame
            
            push!(groups, groupby(X, :step, sort = true))
        end
    end
    
    
    steps = length(keys(groups[1]))
    columns = filter(x -> x ∉["id", "step"], names(groups[1][1]))
    runs = length(groups)
    ids = length(unique(groups[1][1][!, :id]))
    
    
    m = AxisArray(Array{Float64}(undef, steps, length(columns), runs, ids);
        step = 1:steps, column = columns, run = 1:runs, id = 1:ids)
    
    for (run, g) in enumerate(groups)
        for k in keys(g)
            for (ci, c) in enumerate(columns)
                m[k[1]+1, ci, run, 1:end] = g[k][!, c]
            end
        end
    end
    
    m
end

function plotwithribbons(p, data, column)
    means = mean(data, dims=1)[1:end]
    mins = minimum(data, dims=1)[1:end]
    maxs = maximum(data, dims=1)[1:end]

    Plots.plot!(p, means, label = column, ribbon = (means .- mins, maxs .- means))
    
    p
end

let _datapath
    global datapath() = _datapath

    global function datapath(datapath)
        _datapath = datapath
    end

    global function plotspath(name, create = true)
        path = joinpath(datapath(), name)
        if create
            mkpath(path)
        end
        path
    end
end

function cacheplots(M, pfs::AbstractDict)
    backend_before = backend()
    hdf5()
    for step in AxisArrays.axisvalues(M)[1]
        for (plotname, plotfunction) in pfs
            p = plotfunction(M, step)
                
            fullpath = plotspath(plotname)
            Plots.hdf5plot_write(p, joinpath(fullpath, "$step.hdf5"))
        end
    end
    
    Plots.backend(backend_before)
end

cacheplots(M, pfs::AbstractArray) = cacheplots(M, Dict([String(Symbol(f)) for f in pfs] .=> pfs))

function plot_zscores(M::AxisArray, step::Int, columns::AbstractArray; sorted = true)
    p = Plots.plot(title = "#$step", ylabel = "Z-Score", xlabel = "Agent",
            legend = :topleft)
    
    for column in columns
        
        scores = Array(M[step = step, column = column])
        
        if sorted
            sort!(scores, dims=2)
        end
        
        plotwithribbons(p, scores, column)
    end
    
    p
end

plot_zscores(M::AxisArray, step::Int; sorted = true) = plot_zscores(M, step, AxisArrays.axisvalues(M)[2]; sorted = sorted)
plot_zscores(M::AxisArray, step::Int, column::String; sorted = true) = plot_zscores(M, step, [column]; sorted = sorted)

function plot_heatmap(M::AxisArray, step::Int, column, clims = (-Inf, Inf); f::Function = mean, run = 1)
    heatmap(title = "Agents on Grid at step #$(step)",
            reshape(M[column = column, step = step, run = run],
            (floor(Int, (√(size(M)[4]))), floor(Int, (√(size(M)[4]))))), clims = clims)
end

function plot_distributions(M::AxisArray, step::Int, columns::AbstractArray; f::Function = mean)
    p = Plots.plot(title = "#$step", ylabel = "Probability",
        xlabel = "Z-score", legend = :topleft)
    
    for column in columns
        d = Distributions.fit(Normal, f(sort(M[step = step, column = column], dims=2), dims=1))
        Plots.plot!(p, d, label = column * " (μ = $(round(d.μ, digits=2)), σ = $(round(d.σ, digits=2)))")
    end
    
    p
end

plot_distributions(M::AxisArray, step::Int) = plot_distributions(M, step, AxisArrays.axisvalues(M)[2])
plot_distributions(M::AxisArray, step::Int, column::String) = plot_distributions(M, step, [column])

function plot_column(M::AxisArray, column::String)
    p = Plots.plot(title = "$column over Time", xlabel = "Step", ylabel = "Z-Score", legend = :topleft)
    Plots.plot!(minimum(Array(M[column = column]), dims=(2, 3))[1:end], label = "Min")
    Plots.plot!(mean(Array(M[column = column]), dims=(2, 3))[1:end], label = "Mean")
    Plots.plot!(median(Array(M[column = column]), dims =(2,3))[1:end], label = "Median")
    Plots.plot!(maximum(Array(M[column = column]), dims=(2, 3))[1:end], label = "Max")
end

function varnamesPrefix()
    "SemAgentsToolkit_"
end

function foldecCell(code::String)
    c = Pluto.Cell(code)
    c.code_folded = true

    c
end

function generateNotebook(filename::String, datadir::String, betterSelecter::Bool = true; autoactivate::Bool = true)

    cells = [Pluto.Cell("using Markdown, InteractiveUtils, PlutoUI, SemAgentsToolkit, AxisArrays")]
    
    if autoactivate
        push!(cells, Pluto.Cell("""begin
        using Pkg
        Pkg.activate(raw\"\"\"$(Pkg.project().path)\"\"\")
        end"""))
    end

    push!(cells,

        foldecCell("""datapath(raw\"\"\"$datadir\"\"\")"""))
    

    if betterSelecter
        a, b = generateBetterSelecter(get_params_from_filenames(datadir))
        push!(cells, foldecCell(a))
        push!(cells, foldecCell(b))
    else
        push!(cells, foldecCell(raw"""@bind selected_run Select(sort(readdir(datapath())))"""))
    end

    push!(cells, Pluto.Cell("data = load(selected_run)"))

    push!(cells, Pluto.Cell("""
    # usage: plot_column(data, "YOUR_COLUMN_NAME")
    plot_column(data, AxisArrays.axes(data)[2][1])""" ) )

    push!(cells,
        Pluto.Cell(raw"""# Chose one of those two:
            # @bind n Clock(1)
            md\"\"\"Step: $(@bind n Slider(1:size(data)[1], show_value=true))\"\"\"""")
        )

    push!(cells, Pluto.Cell("plot_zscores(data, n)"))

    push!(cells, Pluto.Cell("plot_distributions(data, n)"))

    push!(cells, Pluto.Cell("""
            # usage: plot_heatmap(data, n, "YOUR_COLUMN_NAME")
            plot_heatmap(data, n, AxisArrays.axes(data)[2][1])"""))

    nb =  Pluto.Notebook(cells)

    nb.path = filename
    Pluto.save_notebook(nb)
end