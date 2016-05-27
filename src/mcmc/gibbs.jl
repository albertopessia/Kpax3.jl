# This file is part of Kpax3. License is MIT.

function gibbs!(data::Matrix{UInt8},
                priorR::PriorRowPartition,
                priorC::PriorColPartition,
                settings::KSettings,
                support::MCMCSupport,
                state::AminoAcidState)
  shuffle!(support.u)

  for a in support.u
  end

  nothing
end
