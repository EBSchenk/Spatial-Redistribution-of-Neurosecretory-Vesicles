#Include all relevant libraries
include("Data-Converter.jl")
include("Visualise-Raw-Data.jl")
include("Functions.jl")
include("Hill-System-of-Equations.jl")
include("Parameters.jl")
include("Fig1EMimic.jl")
include("Parameter-Fitting.jl")


###########################################################
#Solve the Hill30 system

#use fitted values
myfittedparams=[0.8078292108704153, 1.2490750191006295, 0.8851945456804439, 0.29582710227751585, 0.4135296092103562, 0.01546590136875035]
myfittedparams2=[0.7202877255738009, 0.602927567061379, 0.19542701667102558, 1.0512433725788306, 0.8928071585870612, 0.049028584546544046]

#p_cont, p_stim = fittedparams(myfittedparams2)

#For Stimulation
u0=[1.0; zeros(length(species(Hill30))-1) ]
prob = SteadyStateProblem(Hill30, u0, p_stim)
sols = solve(prob, maxiters = 1e5,DynamicSS(Tsit5()))
#@time sols = solve(prob)
solS = reassign(sols.u)
#For Control
prob = SteadyStateProblem(Hill30, u0, p_cont)
solc = solve(prob, maxiters = 1e5,DynamicSS(Tsit5()))
solC = reassign(solc.u)
###########################################################

#Assess model solutions

#Average Transition Rates Data#
###################################################
#Calculate average transition rates (For F->D)#

#For Stimulation (Should be about 0.3):
avgf2ds=AvgF2D_Hill(solS)
#For Control (Should be about 0.1):
avgf2dc=AvgF2D_Hill(solC)
println("rate F->D (stim., 0.3)=", avgf2ds, "; rate F->D (control, 0.1)=", avgf2dc)

#Calculate average transition rates (For C->D)#

#For Stimulation (Should be about 0.15):
avgc2ds=AvgC2D_Hill_Stim(solS)
#For Control (Should be about 0.1):
avgc2dc=AvgC2D_Hill_Cont(solC)
println("rate C->D (stim., 0.15)=", avgc2ds, "; rate C->D (control, 0.1)=", avgc2dc)
###################################################


#Calculate Φ#
###################################################

#For Simulation (Should be about 1.0):
Φs=Φ(solS)
#For Control (Should be about 0.60-0.66):
Φc=Φ(solC)
println("Total transport vs non-transport: Φs (1.0)=", Φs, " Φc(0.6)=", Φc)
#Calculate Ψ#

#For Stimulation (Should be about 2.0):
Ψs=Ψ(solS)
#For Control (Should be about 1.0):
Ψc=Ψ(solC)
println("Transport right vs left: Ψs(2.0)=", Ψs, " Ψc(1.0)=", Ψc)
###################################################


#Calculate the effective on-rates
###################################################

#For F->D

#For stimulation (Should be the same as original model):
rates_S = HillRatesFree(solS)
RateFR1_S = rates_S[1]  #Should be 0.3*ρ1 = 0.3*0.75 = 0.225
RateFL1_S = rates_S[2]  #Should be 0.3*(1-ρ1) = 0.3*0.25 = 0.075
RateFR2_S = rates_S[3]  #Should be 0.4*ρ2 = 0.4*0.6 = 0.24
RateFL2_S = rates_S[4]  #Should be 0.4*(1-ρ2) = 0.4*0.4 = 0.16
RateFR3_S = rates_S[5]  #Should be 0.2*ρ3 = 0.2*0.4 = 0.08
RateFL3_S = rates_S[6]  #Should be 0.2*(1-ρ3) = 0.2*0.6 = 0.12

#For control (Should be lower than original model):
rates_C = HillRatesFree(solC)
RateFR1_C = rates_C[1]
RateFL1_C = rates_C[2]
RateFR2_C = rates_C[3]
RateFL2_C = rates_C[4]
RateFR3_C = rates_C[5]
RateFL3_C = rates_C[6]

##########

#For C->D

#For stimulation (Should be the same as original model):
rates_S2 = HillRatesCaged(solS)
RateCR1_S = rates_S2[1]  #Should be 0.1*ρ1 = 0.1*0.75 = 0.075
RateCL1_S = rates_S2[2]  #Should be 0.1*(1-ρ1) = 0.1*0.25 = 0.025
RateCR2_S = rates_S2[3]  #Should be 0.2*ρ2 = 0.2*0.6 = 0.12
RateCL2_S = rates_S2[4]  #Should be 0.2*(1-ρ2) = 0.2*0.4 = 0.08

#For control (Should be lower than original model):
rates_C2 = HillRatesCaged(solC)
RateCR1_S = rates_C2[1]
RateCL1_S = rates_C2[2]
RateCR2_S = rates_C2[3]
RateCL2_S = rates_C2[4]

#Visualise the control and stimulation distributions
###################################################
RawVis(solS,1)
RawVis(solC,0)

#Fig1E Mimic
#Fig1E(solC,solS)
