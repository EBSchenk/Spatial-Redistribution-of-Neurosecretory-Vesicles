using Optim 
include("ParametersModel30.jl")

negative(x)=-min(0.0, x)

function fittedparams(params::Vector)
    myL1hat, myL2hat, myL3hat, myR1hat, myR2hat, myR3hat = abs.(params)

    p_stim = [κrc1,κlc1,κrf1,κlf1,κfd1,R1_stim,L1_stim,myR1hat,myL1hat,κcd1,κcf1,κfc1,
    κrc2,κlc2,κrf2,κlf2,κfd2,R2_stim,L2_stim,myR2hat,myL2hat,κcd2,κcf2,κfc2,
    κrc3,κlc3,κrf3,κlf3,κfd3,R3_stim,L3_stim,myR3hat,myL3hat,κcd3stim,κcf3,κfc3,
    τ,δ,ζ,
    ρ1,ρ2,ρ3]

    p_cont = [κrc1,κlc1,κrf1,κlf1,κfd1,R1_stim,L1_stim,myR1hat,myL1hat,κcd1,κcf1,κfc1,
     κrc2,κlc2,κrf2,κlf2,κfd2,R2_stim,L2_stim,myR2hat,myL2hat,κcd2,κcf2,κfc2,
     κrc3,κlc3,κrf3,κlf3,κfd3,R3_stim,L3_stim,myR3hat,myL3hat,κcd3cont,κcf3,κfc3,
     τ,δ,0.0,
     ρ1,ρ2,ρ3]
    return p_cont, p_stim
end

function errorfunctional(params::Vector)
    p_cont, p_stim = fittedparams(params)

    probstim = SteadyStateProblem(Hill30, u0, p_stim)
    sols = solve(probstim)
    solS = reassign(sols.u)
    
    probcont = SteadyStateProblem(Hill30, u0, p_cont)
    solc = solve(probcont)
    solC = reassign(solc.u)

    
    #For Stimulation (Should be about 0.3):
    avgf2ds=AvgF2D_Hill(solS)
    #For Control (Should be about 0.1):
    avgf2dc=AvgF2D_Hill(solC)
    #println("rate F->D (stim., 0.3)=", avgf2ds, "; rate F->D (control, 0.1)=", avgf2dc)
    ret=(avgf2ds-0.3)^2+20.0*(avgf2dc-0.1)^2

    #Calculate average transition rates (For C->D)#

    #For Stimulation (Should be about 0.15):
    avgc2ds=AvgC2D_Hill_Stim(solS)
    #For Control (Should be about 0.1):
    avgc2dc=AvgC2D_Hill_Cont(solC)
    #println("rate C->D (stim., 0.15)=", avgc2ds, "; rate C->D (control, 0.1)=", avgc2dc)
    ret+=(avgc2ds-0.15)^2+(avgc2dc-0.1)^2
    ###################################################


    #Calculate Φ#
    ###################################################

    #For Simulation (Should be about 1.0):
    Φs=Φ(solS)
    #For Control (Should be about 0.60-0.66):
    Φc=Φ(solC)
    #println("Total transport vs non-transport: Φs (1.0)=", Φs, " Φc(0.6)=", Φc)
    ret+=(Φs-1.0)^2+(Φc-0.6)^2

    #Calculate Ψ#
    #For Stimulation (Should be about 2.0):
    Ψs=Ψ(solS)
    #For Control (Should be about 1.0):
    Ψc=Ψ(solC)
    #println("Transport right vs left: Ψs(2.0)=", Ψs, " Ψc(1.0)=", Ψc)
    ret+=(Ψs-2.0)^2+(Ψc-1.0)^2
    ###################################################
    return ret
end

function errorandregularistion(params)
    return errorfunctional(params)-0.001*sum(abs.(params))+1000*sum(negative.(params))
end

function fitparams()
    result = optimize(errorfunctional, [0.5,0.5,0.5,0.5,0.5,0.5], show_trace=true, time_limit=100)
    return result
end
