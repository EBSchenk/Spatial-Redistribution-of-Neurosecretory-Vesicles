#Set up default parameters

###BASIC PARAMETERS (Base package)
##################################################
const τ,δ,ζ = 0.15,0.04,0.15

#rhos
const ρ1,ρ2,ρ3 = 0.75,0.60,0.40
#Central
const κrc1, κlc1 = 0.1,0.1
const κrf1, κlf1 = 0.04, 0.1
const κfd1 = 0.3
const κcd1 = 0.1
const κcf1, κfc1 = 0.1,0.04

#Subcortical
const κrc2, κlc2 = 0.1,0.1
const κrf2, κlf2 = 0.04, 0.04
const κfd2 = 0.4
const κcd2 = 0.2
const κcf2, κfc2 = 0.03,0.08

#Peripheral
const κrc3, κlc3 = 0.3,0.3
const κrf3, κlf3 = 0.1, 0.04
const κfd3 = 0.2
const κcd3stim = 0.1
const κcd3cont = 0.1

const κcf3, κfc3 = 0.03,0.2
##################################################

###Hill Model Parameters (Extension Package)
##################################################

#Stimulation population numbers (for fibres)
const R1_stim, L1_stim = 0.65, 0.42
const R2_stim, L2_stim = 0.58, 0.22
const R3_stim, L3_stim = 0.4, 0.1

#Hat parameters
const R1hat, L1hat = 0.5, 0.5
const R2hat, L2hat = 0.5, 0.5
const R3hat, L3hat = 0.5, 0.5
##################################################

###Set up the parameter collections for solvers
##################################################

###Hill30 Parameters
p_stim = [κrc1,κlc1,κrf1,κlf1,κfd1,R1_stim,L1_stim,R1hat,L1hat,κcd1,κcf1,κfc1,
     κrc2,κlc2,κrf2,κlf2,κfd2,R2_stim,L2_stim,R2hat,L2hat,κcd2,κcf2,κfc2,
     κrc3,κlc3,κrf3,κlf3,κfd3,R3_stim,L3_stim,R3hat,L3hat,κcd3stim,κcf3,κfc3,
     τ,δ,ζ,
     ρ1,ρ2,ρ3]

 p_cont = [κrc1,κlc1,κrf1,κlf1,κfd1,R1_stim,L1_stim,R1hat,L1hat,κcd1,κcf1,κfc1,
      κrc2,κlc2,κrf2,κlf2,κfd2,R2_stim,L2_stim,R2hat,L2hat,κcd2,κcf2,κfc2,
      κrc3,κlc3,κrf3,κlf3,κfd3,R3_stim,L3_stim,R3hat,L3hat,κcd3cont,κcf3,κfc3,
      τ,δ,0.0,
      ρ1,ρ2,ρ3]
