#Reaction network for:
#Hill Carrying Capacities implemented for F->D transitions (all compartments)
#Hill Carrying Capacities implemented for C->D transitions (central and subcortical)

#We don't alter the peripheral C->D rate for control as we suspect the decreased
#activity of the actin network (i.e. it is no longer actively recruiting vesicles)
#increases the propensity of vesicles to detach in control. So we would expect
#an increase in the C->D transition as a result of this effect. We shall assume
#however that this effect is perfectly (perhaps unrealistic??) countered by the
#crowding effects that would otherwise drive down the C->D rate in control.


#relevant packages
using Catalyst
using DifferentialEquations

cut(x)=max(x, 0.0)

#set up network
Hill30 = @reaction_network begin
  #actin on-rates
  κrc1, R1-->C1
  κlc1, L1-->C1
  #tubule off-rates
  κrf1, R1-->F1
  κlf1, L1-->F1
  #tubule on-rates
  #κfd1*ρ1, F1-->R1
  #κfd1*(1-ρ1), F1-->L1
  ρ1*κfd1/(1+cut(R1-R1_stim)/R1hat), F1-->R1
  (1-ρ1)*κfd1/(1+cut(L1-L1_stim)/L1hat), F1-->L1
  #actin off-rates
  ρ1*κcd1/(1+cut(R1-R1_stim)/R1hat), C1-->R1
  (1-ρ1)*κcd1/(1+cut(L1-L1_stim)/L1hat), C1-->L1
  #actin-tubule exchange rates
  κcf1, C1 --> F1
  κfc1, F1 --> C1

  #actin on-rates
  κrc2, R2-->C2
  κlc2, L2-->C2
  #tubule off-rates
  κrf2, R2-->F2
  κlf2, L2-->F2
  #tubule on-rates
  ρ2*κfd2/(1+cut(R2-R2_stim)/R2hat), F2-->R2
  (1-ρ2)*κfd2/(1+cut(L2-L2_stim)/L2hat), F2-->L2
  #actin off-rates
  ρ2*κcd2/(1+cut(R2-R2_stim)/R2hat), C2-->R2
  (1-ρ2)*κcd2/(1+cut(L2-L2_stim)/L2hat), C2-->L2
  #actin-tubule exchange rates
  κcf2, C2 --> F2
  κfc2, F2 --> C2

  #actin on-rates
  κrc3, R3-->C3
  κlc3, L3-->C3
  #tubule off-rates
  κrf3, R3-->F3
  κlf3, L3-->F3
  #tubule on-rates
  #κfd3*ρ3, F3-->R3
  #κfd3*(1-ρ3), F3-->L3
  ρ3*κfd3/(1+cut(R3-R3_stim)/R3hat), F3-->R3
  (1-ρ3)*κfd3/(1+cut(L3-L3_stim)/L3hat), F3-->L3
  #actin off-rates
  ρ3*κcd3, C3-->R3
  (1-ρ3)*κcd3, C3-->L3
  #actin-tubule exchange rates
  κcf3, C3 --> F3
  κfc3, F3 --> C3

  #transport between compartments via fibres
  τ, R1-->R2
  τ, R2-->R3
  τ, L3-->L2
  τ, L2-->L1

  #transport between compartments via diffusion
  δ, F1 --> F2
  δ, F2 --> F1
  δ, F2 --> F3
  δ, F3 --> F2

  #inflow
  100.0*(1.0-C1), 0-->C1
  #100.0*(5.0-C1), 0-->C1
  #outflow
  ζ, C3 --> 0
end κrc1 κlc1 κrf1 κlf1 κfd1 R1_stim L1_stim R1hat L1hat κcd1 κcf1 κfc1 κrc2 κlc2 κrf2 κlf2 κfd2 R2_stim L2_stim R2hat L2hat κcd2 κcf2 κfc2 κrc3 κlc3 κrf3 κlf3 κfd3 R3_stim L3_stim R3hat L3hat κcd3 κcf3 κfc3 τ δ ζ ρ1 ρ2 ρ3;
