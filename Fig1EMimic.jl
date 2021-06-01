using StatsPlots
#include("FunctionsModel30.jl")

function Fig1E(solc::compartments,sols::compartments)
        #calculate the motile state percentages
        Ds = NetDirected(sols)
        Fs = NetFree(sols)
        Cs = NetCaged(sols)
        Dc = NetDirected(solc)
        Fc = NetFree(solc)
        Cc = NetCaged(solc)

        data = [Dc, Fc, Cc, Ds, Fs, Cs]
        category = repeat(["Control","Stimulation"],inner = 3)
        std = repeat([0.5],inner = 6) #we don't actually have this...but could.
        name = ["Directed", "Free", "Caged"]

        groupedbar(name,data,yerr = std, group = category, ylabel = "
        Percentage of Secretory Vesicles", xlabel = "Behaviour of Vesicles")

end
