module MarkovBlocks

using DataFrames, LinearAlgebra
using StatsBase: mean, ordinalrank, quantile
using NearestNeighbors
using NaNStatistics

export Datap

# Misc Functions


function indicatormat(X::Matrix{Float64}, k::Vector)
    dt = k[2]-k[1]
    km = [k k.+dt]
    nk, nX = size(km, 1), size(X,1)
    A = falses(nX, nk) 
    for i in 1:nk
         A[:,i]  = km[i,1] .≤ view(X,:,1) .< km[i,2] 
    end
    A *= 1.0
    return A    
end


function bindata(Y::Matrix{Float64}, k::Vector)
    A = indicatormat(Y, k)
    n = sum(A, dims=1)
    return (permutedims(A)*view(Y,:,2))./vec(n) 
end


function cor_pairwise(X::AbstractMatrix, Y::AbstractMatrix)
    X[ismissing.(X)] .= 0.0
    Y[ismissing.(Y)] .= 0.0
    numer = X'Y
    denom = sqrt.(sum(abs2, X, dims=1)'sum(abs2, Y, dims=1))
    return numer./denom
end

# Markov Block Bootstrap

fn1(i::Bool, inds::Vector{Int}, rinds::Vector{Int}, k::Int) = i ? inds.+1 : rand(rinds[Not(k)], 6)


function preboot(d::Matrix; opts::NamedTuple=(;))
    y = d[:,1]
    ys = div.(y, opts.blocklength)*opts.blocklength
    dy = y-ys
    uys = unique(ys)
    blocks = [ [dy d[:,2]][y.==ys,:] for y in uys]
# block nearest neighbors
    nbl = length(blocks)
#    @show extrema(y)
    μblock = [mean(b[:,2]) for b in blocks[Not(end)]]
    rank = ordinalrank(μblock)
    kdtree = KDTree(permutedims(1.0*rank))
    blockinds = inrange(kdtree, permutedims(rank), 3)
    blocki = 1:nbl-1
    rinds = collect(1:nbl)
    hasneighbor = diff(uys).==opts.blocklength
    dbli = Dict(k=>fn1(i,v,rinds,k) for (k,i,v) in zip(blocki, hasneighbor, blockinds))
    push!(dbli, nbl=>[rand(rinds[Not(nbl)])])
    return (uys=uys, blocks=blocks, n=length(blocks), dbli=dbli)
end


function boot(x::NamedTuple)
    r = Vector{Int}(undef, x.n)
    r[1] = rand(1:x.n)
    for j in 2:x.n
        r[j] = rand(x.dbli[r[j-1]])
    end
    return reduce(vcat, [b.+[y 0] for (b, y) in zip(x.blocks[r], x.uys)])
end


function bints(X::SubArray, kv::Vector)
    dt = 0.5*abs(kv[2]-kv[1])
    y = bindata(vcat(X...), kv)
    goodi = isfinite.(y)
    return (t=kv[goodi].+dt, O=y[goodi].-nanmean(y))
end

# Unstack and eigenvalue decomposition
function unstak(d::SubDataFrame, opts::NamedTuple; verbose=false)
    d2 = unstack(d, :tsname, :O)[:,2:end]
    j =  [count(ismissing.(x)) for x in eachcol(d2) ] .< (opts.percentmissing*size(d2,1))
    j[1] = false
    p = count(j)
    X = Matrix(d2[:,j])
    C = cor_pairwise(X,X)
    verbose && @show tr(C)
    irange = max(p-opts.neigenvals+1,1):p
    ev = eigvals(Symmetric(C), irange)
    return (Rank=reverse(axes(irange,1)), ev=100*ev./tr(C))
end

# Altogether now


function bootdf(Ds::DataFrame, opts::NamedTuple; nboot::Int=10)
    Dboot = reduce(vcat, [bstat(Ds, opts) for j in 1:nboot])
    Dbq = combine(groupby(Dboot, [:isoint1_vG, :Rank]), :ev=>quantilef=>AsTable, :binrange=>first=>:binrange)
    return Dbq
end


function bstat(Ds::DataFrame, opts::NamedTuple)
    D2 = select(Ds, :tsname, :isoint1_vG, :pb=>ByRow(boot)=>:dO)
    D3 = eigvalf(D2, opts)
    return D3
end


function eigvalf(Ds::DataFrame, opts::NamedTuple; verbose=false)
    Dy = combine(groupby(Ds, [:tsname, :isoint1_vG]), :dO=>(x->bints(x, opts.binv))=>AsTable)
    Dev = combine(x->unstak(x, opts, verbose=verbose), groupby(Dy, :isoint1_vG))
    Dev[!,:binrange] .= "$(opts.binrange) CE"
    return Dev
end


function quantilef(d::SubArray)
   q =  quantile(d, [0.1, 0.99])
   return (qlo=q[1], qup=q[2])
end


include("datap.jl")


function initdf(opts, a; datayears=[0.0, 2019.0])
    D1, dO = Datap.readf(a, datayears)
    pboot = preboot.(dO, opts=opts)

    # Make a reduced dataframe, add a multi-label variable: tsname
    Dsite = D1[:, [:dataSetName, :TSid, :lon, :lat, :isoint1_vG, :atype]]
    Dsite[!,[1,2]] = string.(Dsite[:,[1,2]])
    icet = map(in(["IC", "IW", "GI"]), D1.atype)
    Dsite.isoint1_vG[icet] .= "P_isotope"
    Dsite.tsname = .*(D1.dataSetName, "_", D1.vtype, "_", D1.measMat)
    @show size(Dsite)
    
    ireg = Datap.regioni(opts, Dsite)
    Ds = copy(Dsite)
    Ds.dO = dO
    Ds.pb = pboot
    Ds = Ds[ireg,:]
    return Ds
end







end
