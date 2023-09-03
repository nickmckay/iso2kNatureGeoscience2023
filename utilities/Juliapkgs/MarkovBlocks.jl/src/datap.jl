

module Datap

using DataFrames, NaNStatistics

function span(x::Matrix, k=[0.0, 2019])
    i = k[1].<x.<k[2]
    isempty(x[i]) && return NaN
    mn, mx = nanextrema(x[i])
    return mx-mn
end


function noNaN(d::Matrix, k=[0.0, 2019])
    goodi = isfinite.(d[:,2]) .& (k[1].<d[:,1].<k[2])
    return sortslices(d[goodi,:], dims=1)
end


function readf(a::Dict, datayears= [0.0, 2019]; var="sTS",
    fields=[["dataSetName", "geo_longitude", "geo_latitude"]; "paleoData_".*["variableName","description","iso2kPrimaryTimeseries", "inferredMaterial", "measurementMaterial", "TSid"]; "isotopeInterpretation1_".*["variable","variableGroup","seasonality", "inferredMaterial", "direction"]; "inferredMaterial"])

    println(unique(a[var]["archiveType"]))
    year = vec(a[var]["year"])
    pDv = vec(a[var]["paleoData_values"])
    tspan = span.(year, [datayears])

    
    d1 = filter(f->in(f[1], fields), a[var])
    d2 = Dict(k=>vec(v)   for (k, v) in d1)
    D1 = DataFrame(d2)

    rename!(D1, :geo_longitude=>:lon, :geo_latitude=>:lat,:isotopeInterpretation1_variable=>:isoint1_variable,
        :isotopeInterpretation1_variableGroup=>:isoint1_vG, :isotopeInterpretation1_direction=>:isoint1_dir,
        :paleoData_inferredMaterial=>:infMat,
        :isotopeInterpretation1_inferredMaterial=>:isoint1_infMat,
        :isotopeInterpretation1_seasonality=>:isoint1_seas, 
        :paleoData_variableName=>:vtype, :paleoData_description=>:dscrptn, :paleoData_iso2kPrimaryTimeseries=>:i2kPT,
        :paleoData_measurementMaterial=>:measMat, :paleoData_TSid=>:TSid)

    [x[isempty.(x)].="" for x in eachcol(D1)]
    D1[!,[:lat,:lon]] .*= 1.0
    D1.i2kPT[ismissing.(D1.i2kPT)] .= ""
    D1.atype = [x[1:2] for x in D1.dataSetName]
    D1.hem = [(l>0 ? "NH" : "SH")  for l in D1.lat]


    tspanmin = (first(diff(datayears))>1000.0) ? 750.0 : 75.0
    i1 = (length.(year).>0.0) .& map(in(["d2H","d18O"]), D1.vtype) .&  (tspanmin .< tspan)
    i1 = i1 .& (D1.dataSetName.≠"TR13SITA") .& (D1.dataSetName.≠"IC08VI77")
    i1 = i1 .& map(in(["EffectiveMoisture","P_isotope","Temperature"]), D1.isoint1_vG)

    D1 = D1[i1,:]

    sgndict = Dict("positive"=>1, "negative"=>-1, "postive"=>1, ""=>1, "decrease"=>-1, "increase"=>1, "positve"=>1)
    sgn = [sgndict[x] for x in D1.isoint1_dir]

    dO = hcat.(year[i1], sgn.*pDv[i1])
    dO = noNaN.(dO, [datayears])
    i2 = .!isempty.(dO)
    D1, dO = D1[i2,:], dO[i2]

    return D1, dO
end


function optionf(;region::String="Global", binwidth::Float64=30.0, blocklength::Float64=10.0, 
        percentmissing::Float64=0.15, neigenvals::Int=15, years=[0.0,1950.0])
    binv = collect(years[1]:binwidth:years[2])
    binrange = string(Int(binv[1]),"-",floor(Int,binv[end]+binwidth))
    return (; region, binwidth, binv, binrange, blocklength, percentmissing, neigenvals)
end


function regioni(opts::NamedTuple, Dsite::DataFrame=DataFrame())
    regi = if opts.region=="Global"
        trues(size(Dsite,1))
    elseif opts.region=="Pacific"
        (-50 .< Dsite.lat .< 53) .& ((58 .< Dsite.lon) .| (Dsite.lon .< -79))
    end
    return regi
end






end



