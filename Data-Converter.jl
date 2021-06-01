#define a struct mocking the output of the reaction network
#i.e. compares with species(rn)

struct compartments
  R1::Float64
  C1::Float64
  L1::Float64
  F1::Float64
  R2::Float64
  C2::Float64
  L2::Float64
  F2::Float64
  R3::Float64
  C3::Float64
  L3::Float64
  F3::Float64
end

#convert steady state output into compartments struct
function reassign(mysol::Vector{Float64})
  return compartments(mysol...)
end
