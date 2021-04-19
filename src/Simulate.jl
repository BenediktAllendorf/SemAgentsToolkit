export index_of, create_sem, SemAgent, InterconnectedSEM, setup_adata, init, simulate, plot, forfilename, Configuration

struct SEM
	g::MetaDiGraph
	coef::DataFrame
	r²::DataFrame
end

@with_kw struct Configuration
    f::Function
    from::Symbol
    to::Symbol
    weight::Float64
    steps::Int
	dims::Tuple{Int64, Int64}
	populate::Float64 = 1.
    jitter::Float64 = 0.
    zscorception::Bool = false
	runs::Int = 1
end

forfilename(x::Configuration) = Dict(fn=>String(Symbol(getfield(x, fn))) for fn ∈ fieldnames(Configuration) if fn != :runs)

struct InterconnectedSEM
	interconnectivity::SEM
	base::SEM
	network_edge::Edge
	network_edge_weight::Float64

	function InterconnectedSEM(base::SEM, from, to, weight::Float64)
		columns = vcat(names(base.coef), [neighborize(n) for n in names(base.coef)])
	
		df = DataFrame(fill(.0, ntuple(x->length(columns), 2)), columns)
		
		for (neighbor_construct, v) in total_effect(base, from)
			for (self_construct, v2) in total_influence(base, to)
				df[index_of(columns, neighborize(neighbor_construct)), self_construct] = v * weight * v2
			end
			df[index_of(columns, neighborize(neighbor_construct)), to] = v * weight
		end
	
		df[index_of(columns, neighborize(from)), to] = weight
	
		for (self_construct, v2) in total_influence(base, to)
			df[index_of(columns, neighborize(from)), self_construct] = weight * v2
		end

		new(
			create_sem(df),
			base,
			Edge(index_of(columns, from), index_of(columns, to)),
			weight)
	end

	function InterconnectedSEM(base::SEM, conf::Configuration)
		InterconnectedSEM(base, conf.from, conf.to, conf.weight)
	end

end

mutable struct HumanBeing
	iSem::InterconnectedSEM
	constructvalues::Array{Float64}
end
HumanBeing(iSem::InterconnectedSEM) = HumanBeing(iSem, default_constructs(iSem))

mutable struct SemAgent <: AbstractAgent
	id::Int
    pos::NTuple{2, Int}
	h::HumanBeing
    lastneightborsvalues::Array{Float64}
end

function setup_adata(relevant_constructs::AbstractArray{Symbol}, iSem::InterconnectedSEM)
    adata = []
    for c in relevant_constructs
        idx = index_of(names(iSem.base.coef), c)

        name = Symbol("fn_", c)
        
        eval(quote function $(name)(x) x.h.constructvalues[$(QuoteNode(idx))] end end)
        
        push!(adata, eval(name))
    end

    adata
end

function model_step!(model)
    newvalues = fill(.0, (length(allagents(model)), size(first(allagents(model)).h.iSem.base.coef)[1]))

    for agent in allagents(model)
        if isempty(nearby_agents(agent, model))
            continue
        end

        neighborvalues = model.f(reduce(hcat, [a.h.constructvalues for a in nearby_agents(agent, model)]), dims=2)

        diff = neighborvalues - agent.lastneightborsvalues

        # agent.lastneightborsvalues = neighborvalues

        changes = rand(Normal(0, model.jitter), cnt_constructs(agent.h.iSem), 1)

        for (i, cn) in enumerate(Symbol.(construct_names(agent.h.iSem)))
            changes .+= change_construct_value_neighbor_based(agent.h, neighborize(cn, Symbol), diff[i])
         end

        newvalues[agent.id, :] = agent.h.constructvalues + changes
    end

    if model.zscorception
        for dim in 1:size(newvalues)[2]
            newvalues[:, dim] .= zscore(newvalues[:, dim])
        end
        replace!(newvalues, NaN => 0)
    end

    for agent in allagents(model)
        agent.h.constructvalues .= newvalues[agent.id, :]
    end
end

function simulate(myISem::InterconnectedSEM, conf::Configuration, relevant_constructs::AbstractArray{Symbol})
	simulate(myISem, conf.dims, conf.steps, relevant_constructs, conf.f, conf.runs, conf.jitter, conf.zscorception, conf.populate)
end

function simulate(myISem::InterconnectedSEM, dims::Tuple{Int, Int}, steps::Int, relevant_constructs::AbstractArray{Symbol},
	f::Function,
	runs = 1,
	jitter = .0,
	zscorception = false,
	populate = 1)

    adata = setup_adata(relevant_constructs, myISem)

	agent_df = missing

	chnl = Channel(0, spawn = true) do c
		for run in 1:runs
			model = init(dims, myISem, populate; f = f, jitter = jitter, zscorception = zscorception)
			agent_df, _ = Agents.run!(model, dummystep, SemAgentsToolkit.model_step!, steps; adata = adata)

			for (i, c) in enumerate(adata)
				rename!(agent_df, [string(c) => relevant_constructs[i]])
			end

			put!(c, agent_df)
		end
	end

	# return
	if runs == 1
		take!(chnl)
	else
		chnl
	end
end

function init(grid, myISem, populate; f, jitter, zscorception)

	properties = Dict{Symbol,Any}()
	@pack! properties =  f, jitter, zscorception

    space = GridSpace(grid)
    model = ABM(SemAgent, space; properties = properties)

    for n in 1:populate * grid[1] * grid[2]
        h = generate!(HumanBeing(myISem, fill(.0, size(myISem.base.coef)[1])))
        agent = SemAgent(n, (1, 1), h, default_constructs(myISem))
        add_agent_single!(agent, model)
    end

    model
end




Base.deepcopy(s::SEM) = SEM(deepcopy(s.g), deepcopy(s.coef), deepcopy(s.r²))

function create_sem(coef, r² = DataFrame())
	# build graph with nodes and edges
	adj_m = convert(Matrix, coef)
	tmp_g = SimpleWeightedDiGraph(adj_m)
	
	# create metagraph to store weights
	g = MetaDiGraph(size(adj_m)[1])

	# set names for each node
	for n in vertices(g)
		set_prop!(g, n, :name, Symbol(names(coef)[n]))
	end

	# use names as index
	set_indexing_prop!(g, :name)
	
	# add edges with weights to (meta)graph
	for ind in findall(!iszero, tmp_g.weights)
		e = Edge(ind[2], ind[1])
		add_edge!(g, e)
		set_prop!(g, e, :weight, tmp_g.weights[ind[1], ind[2]])
	end
	
	# set missing r² values to one
	for leaf in [x for x in names(coef) if x ∉ names(r²)]
		insertcols!(r², Symbol(leaf) => 0)
	end
	
	SEM(g, coef, r²)
end

"""
Create SEM-structure from Seminr-Exports like this:

write.csv(pls[["rSquared"]], "r2.csv")
write.csv(pls[["path_coef"]], "coef.csv")

"""
function read_from_seminr(coef_filename, r²_filename)
	# read coef values
	coef = CSV.File(coef_filename) |> DataFrame
	select!(coef, Not(:Column1))

	# read r² values
	r² = CSV.File(r²_filename) |> DataFrame
	delete!(r², 1)
	
	coef, r²
end

function plot(sem::SEM; kwargs...)
	edgelabel = Dict{Any, Any}(
		(src(e), dst(e)) => round(get_prop(sem.g, e, :weight), digits = 3)
		for e in edges(sem.g))

	graphplot(sem.g,
	edgelabel=edgelabel,
	names = names(sem.coef)[1:end],
	arrow = true,
	nodeshape = :circle; kwargs...)
end

function init_coef(sem::SEM, all)	
	if !all
		return OrderedDict{Symbol, Float64}()
	end
	
	target_coef = OrderedDict((Symbol(c), .0) for c in names(sem.coef))
end

function total_effect(sem::SEM, target, all = false)
	g = copy(sem.g)
	
	# dict of nodes affecting the target node
	target_coef = init_coef(sem, all)

	# recursively walk through the graph
	function walk(node, coe)
		
		# init with coefficient of one if origin is new
		if !haskey(target_coef, node)
			target_coef[node] = 0
		end
		
		# multiply additional coefficient
		target_coef[node] += coe

		# get all incoming neighbors
		ns = copy(inneighbors(g, g[node, :name]))
		
		# walk into each neighbor
		for neighborid in ns
			edge = Edge(neighborid, g[node, :name])
			tmp = coe * get_prop(g, edge, :weight)
			rem_edge!(g, edge)
			walk(g[neighborid, :name], tmp)
		end
	end
	
	# start walking
	walk(target, 1)
	
	# delete (temporarily) self-effect and return dict
	delete!(target_coef, target)
end

function total_influence(sem::SEM, origin, all = false)
	# dict of nodes affected by origin node
	origin_coef = init_coef(sem, all)

	# for each outgoing neighbor
	ns = outneighbors(sem.g, sem.g[origin, :name])
	for n in ns
		# get the coefficient to that neighbor
		origin_coef[sem.g[n, :name]] =
			get_prop(sem.g, Edge(sem.g[origin, :name], n), :weight)

		# calculate all neighbors of the neighbor (recursively)
		sub_coefs = total_influence(sem::SEM, sem.g[n, :name])

		# for each neighbors' coefficient
		for (node, sub_coef) in sub_coefs
			if !haskey(origin_coef, node)
				origin_coef[node] = 0
			end
			
			# add the combined effect from the origin node to that neighbor's neighbor
			origin_coef[node] += origin_coef[sem.g[n, :name]] * sub_coef
		end
	end
	
	origin_coef
end


index_of(arr::AbstractArray, construct) = findfirst(x -> x == String(construct), arr)

index_of(sem::SEM, construct) = index_of(names(sem.coef), construct)

neighborize(s) = "[$(s)]"
neighborize(s, T) = T(neighborize(s))

construct_names(sem::SEM) = names(sem.coef)
construct_names(iSem::InterconnectedSEM) = construct_names(iSem.base)

function change_construct_value(h::HumanBeing, construct, change)
	coef = [y for (x, y) in
			merge(total_influence(h.iSem.base, construct, true),
				Dict(construct => 1))]
	
	change * coef
end

function change_construct_value_neighbor_based(h::HumanBeing, construct, change)
	coef = [y for (x, y) in total_influence(h.iSem.interconnectivity, construct, true)]

	(change * coef)[1:size(h.iSem.base.coef)[1]]
end

function plot(asem::InterconnectedSEM; kwargs...)
	edgelabel = Dict{Any, Any}(
		(src(e), dst(e)) => round(get_prop(asem.base.g, e, :weight), digits = 3)
		for e in edges(asem.base.g))
	
	g = deepcopy(asem.base.g)
	# check if edge already exists
	add_edge!(g, asem.network_edge)
	edgelabel[(src(asem.network_edge), dst(asem.network_edge))] =
											"[ $(asem.network_edge_weight) ]"
	
	graphplot(g,
	edgelabel=edgelabel,
	names = names(asem.base.coef)[1:end],
	arrows = true,
	nodeshape = :circle; kwargs...)
end

function generate!(human::HumanBeing, construct)
	idx = index_of(human.iSem.base, construct)

	if .0 != human.constructvalues[idx]
		return human
	end

	for (y, v)  in enumerate(human.iSem.base.coef[!, construct])
		if v == 0 # value does not affect selected construct
			continue
		end

		generate!(human::HumanBeing, names(human.iSem.base.coef)[y])
	end

	# calculate z score based on influences of other construct
	z = human.constructvalues ⋅ human.iSem.base.coef[!, construct]

	# add unexplained variance
	human.constructvalues[idx] = z + rand(Normal(0, 1 -  √(human.iSem.base.r²[!, construct][1])))
end

function generate!(human::HumanBeing)
	for construct in names(human.iSem.base.coef)
		generate!(human, construct)
	end

	human
end

function default_constructs(sem::SEM)
	zeros(size(sem.coef)[1])
end

default_constructs(iSem::InterconnectedSEM) = default_constructs(iSem.base)

cnt_constructs(sem::SEM) = size(sem.coef)[1]
cnt_constructs(iSem::InterconnectedSEM) = cnt_constructs(iSem.base)