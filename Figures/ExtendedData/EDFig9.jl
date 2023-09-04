
# Extended Data Fig. 9


import Cairo
using DataFrames, Gadfly
using MAT
import Downloads as Dl

# For MarkovBlocks, see the front README for this repo
using MarkovBlocks: Datap
import MarkovBlocks as Mb


f1 = Dl.download("https://lipdverse.org/iso2k/1_0_0/iso2k1_0_0.mat", tempname()*".mat")
a = matread(f1)

D1, dO = Datap.readf(a)
@show size(D1)

# For the three rows in ED Fig. 9:

nboot = 1000

o1 = Datap.optionf(years=[850.0, 1810.0])
Ds1 = Mb.initdf(o1, a)
Ds1b = copy(Ds1)
Dev1 = Mb.eigvalf(Ds1, o1, verbose=true)
@time Dbq1 = Mb.bootdf(Ds1b, o1, nboot=nboot)

o2 = Datap.optionf()
Ds2 = Mb.initdf(o2, a)
Ds2b = copy(Ds2)
Dev2 = Mb.eigvalf(Ds2, o2, verbose=true)
@time Dbq2 = Mb.bootdf(Ds2b, o2, nboot=nboot)

o3 = Datap.optionf(binwidth=10.0, blocklength=10.0, years=[1850.0, 1990.0], percentmissing=0.1)
Ds3 = Mb.initdf(o3, a, datayears=[1800.0,2019.0])
Ds3b = copy(Ds3)
Dev3 = Mb.eigvalf(Ds3, o3, verbose=true)
@time Dbq3 = Mb.bootdf(Ds3b, o3, nboot=nboot)

Dev = vcat(Dev1, Dev2, Dev3)
Dbq = vcat(Dbq1, Dbq2, Dbq3)
Devq = innerjoin(Dev, Dbq, on=[:binrange,:isoint1_vG, :Rank])
Devq.sig = Devq.ev .> Devq.qup


p = plot(Devq, x=:Rank, xgroup=:isoint1_vG, ygroup=:binrange, 
        Geom.subplot_grid(
            layer(ymin=:qlo, ymax=:qup, color=[colorant"deepskyblue"], Geom.errorbar),
            layer(y=:ev, color=[colorant"yellow3"], Geom.point, alpha=:sig),
            free_y_axis=true),
        Scale.xgroup(levels=["Temperature", "EffectiveMoisture", "P_isotope"]),
        Scale.ygroup(levels=["850-1840 CE", "0-1980 CE", "1850-2000 CE"]),
        Guide.xlabel("Eigenvalue Rank"), Guide.ylabel("% Variance"),
        Theme(errorbar_cap_length=0mm, discrete_highlight_color=identity,
        line_width=3pt, alphas=[0, 0.5])
)

draw(PDF("EDfig9.pdf", 6.6inch, 6.6inch), p)
