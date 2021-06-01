
using Plots

function RawVis(sol::compartments,class)

    #plot initial bar (the R1 population bar)
    #################################################################
    point1 = [0.05,0.95]
    point2 = [sol.R1,sol.R1]

    #Set up axes and title (dependent on class of solution)
    if class == 0
        x_labels = string.(["Central","Subcortical","Peripheral"])
        plot(point1,point2,xlims = (0,3.6), ylims = (0,10.0),
        linewidth = 5,xticks = (0.5:1:2.5, x_labels),color = :blue,
        label = "Right Directed", title = "Control Cell")
    else
        x_labels = string.(["Central","Subcortical","Peripheral"])
        plot(point1,point2,xlims = (0,3.6),ylims = (0,6.0),
        linewidth = 5,xticks = (0.5:1:2.5, x_labels),color = :blue,
        label = "Right Directed", title = "Stimulated Cell")
    end

    #add the Subcortical and Peripheral right directed
    #Subcortical
    point1 = [1.05,1.95]
    point2 = [sol.R2,sol.R2]
    plot!(point1,point2,color = :blue,linewidth = 5,label = false)
    #Peripheral
    point1 = [2.05,2.95]
    point2 = [sol.R3,sol.R3]
    plot!(point1,point2,color = :blue,linewidth = 5,label = false)
    #################################################################

    #Free populations
    #################################################################
    point1 = [0.05,0.95]
    point2 = [sol.F1,sol.F1]
    plot!(point1,point2,color = :green,linewidth = 5,label = "Free")
    point1 = [1.05,1.95]
    point2 = [sol.F2,sol.F2]
    plot!(point1,point2,color = :green,linewidth = 5,label = false)
    point1 = [2.05,2.95]
    point2 = [sol.F3,sol.F3]
    plot!(point1,point2,color = :green,linewidth = 5,label = false)
    #################################################################

    #Caged Populations
    #################################################################
    point1 = [0.05,0.95]
    point2 = [sol.C1,sol.C1]
    plot!(point1,point2,color = :red,linewidth = 5,label = "Caged")
    point1 = [1.05,1.95]
    point2 = [sol.C2,sol.C2]
    plot!(point1,point2,color = :red,linewidth = 5,label = false)
    point1 = [2.05,2.95]
    point2 = [sol.C3,sol.C3]
    plot!(point1,point2,color = :red,linewidth = 5,label = false)
    #################################################################

    #Left directed populations
    #################################################################
    point1 = [0.05,0.95]
    point2 = [sol.L1,sol.L1]
    plot!(point1,point2,color = :yellow,linewidth = 3,label = "Left Directed")
    point1 = [1.05,1.95]
    point2 = [sol.L2,sol.L2]
    plot!(point1,point2,color = :yellow,linewidth = 3,label = false)
    point1 = [2.05,2.95]
    point2 = [sol.L3,sol.L3]
    plot!(point1,point2,color = :yellow,linewidth = 3,label = false,
    xlabel = "Cellular Region",ylabel = " % of Steady State Population", show=true)

end
