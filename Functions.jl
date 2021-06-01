#Useful functions to be used in our investigations of the system

#calculate the effective on-rates (use for Hill model)
function HillRatesFree(sol::compartments)

    #for the free-to-right in Central
    rateFR1 = ρ1*κfd1/(1+(sol.R1-R1_stim)/R1hat)
    #for the free-to-left in Central
    rateFL1 = (1-ρ1)*κfd1/(1+(sol.L1-L1_stim)/L1hat)

    #for the free-to-right in Subcortical
    rateFR2 = ρ2*κfd2/(1+(sol.R2-R2_stim)/R2hat)
    #for the free-to-left in Subcortical
    rateFL2 = (1-ρ2)*κfd2/(1+(sol.L2-L2_stim)/L2hat)

    #For the free-to-right in Peripheral
    rateFR3 = ρ3*κfd3/(1+(sol.R3-R3_stim)/R3hat)
    #For the free-to-left in Peripheral
    rateFL3 = (1-ρ3)*κfd3/(1+(sol.L3-L3_stim)/L3hat)

    return [rateFR1,rateFL1,rateFR2,rateFL2,rateFR3,rateFL3]
end

function HillRatesCaged(sol::compartments)
    #for the free-to-right in Central
    rateCR1 = ρ1*κcd1/(1+(sol.R1-R1_stim)/R1hat)
    #for the free-to-left in Central
    rateCL1 = (1-ρ1)*κcd1/(1+(sol.L1-L1_stim)/L1hat)
    #for the free-to-right in Subcortical
    rateCR2 = ρ2*κcd2/(1+(sol.R2-R2_stim)/R2hat)
    #for the free-to-left in Subcortical
    rateCL2 = (1-ρ2)*κcd2/(1+(sol.L2-L2_stim)/L2hat)

    return [rateCR1,rateCL1,rateCR2,rateCL2]
end

#just calculate effective on-rate (Free to Directed (both R and L)) for Central
function HillRateFree1(sol::compartments)
    rates = HillRatesFree(sol)
    return rates[1]+rates[2]
end

function HillRateFree2(sol::compartments)
    rates = HillRatesFree(sol)
    return rates[3]+rates[4]
end

function HillRateFree3(sol::compartments)
    rates = HillRatesFree(sol)
    return rates[5]+rates[6]
end

#Calculate the effective on-rates (caged to directed)

function HillRateCaged1(sol::compartments)
    rates= HillRatesCaged(sol)
    return rates[1]+rates[2]
end

function HillRateCaged2(sol::compartments)
    rates= HillRatesCaged(sol)
    return rates[3]+rates[4]
end

#calculate the average free-to-directed transition rate for the system (Hill)
function AvgF2D_Hill(sol::compartments)
    return (sol.F1*HillRateFree1(sol) + sol.F2*HillRateFree2(sol) +
  sol.F3*HillRateFree3(sol))/(sol.F1+sol.F2+sol.F3)
end

#calculate the average free-to-caged transition rate for the system (Hill)
function AvgC2D_Hill_Stim(sol::compartments)
    return (sol.C1*HillRateCaged1(sol) + sol.C2*HillRateCaged2(sol) +
  sol.C3*κcd3stim)/(sol.C1+sol.C2+sol.C3)
end

function AvgC2D_Hill_Cont(sol::compartments)
    return (sol.C1*HillRateCaged1(sol) + sol.C2*HillRateCaged2(sol) +
  sol.C3*κcd3cont)/(sol.C1+sol.C2+sol.C3)
end

#get the total number of directed vesicles
function NetDirected(sol::compartments)
    return sol.R1+sol.R2+sol.R3+sol.L1+sol.L2+sol.L3
end

#calculate total number of free vesicles
function NetFree(sol::compartments)
    return sol.F1+sol.F2+sol.F3
end

#calculate total number of caged vesicles
function NetCaged(sol::compartments)
    return sol.C1+sol.C2+sol.C3
end

#calculate Φ
function Φ(sol::compartments)
    Φ = (NetDirected(sol))/(NetFree(sol)+NetCaged(sol))
    return Φ
end

#calculate Ψ
function Ψ(sol::compartments)
    Ψ = (sol.R1+sol.R2+sol.R3)/(sol.L1+sol.L2+sol.L3)
end
