# PSP_reaction_julia

A selection of drivers and plotting tools for the TurbulentMixingParticleTrackingReactions module

Dependencies:   Julia 1.6 is required.
                TurbulentMixingParticleTrackingReactions can be found at https://github.com/Asoppitt/TurbulentMixingParticleTrackingReactions.jl . Install using \`import Pkg; Pkg.add(url="https://github.com/Asoppitt/TurbulentMixingParticleTrackingReactions.jl")\`
Additionally for the drivers: Random
Additionally for the plotting code: Plots

for use, ensure  PSP_Particletracking_module.jl is in the same folder as the drivers

run any of the drivers in julia, it will save the data as a binary file
(PSP_Particletracking_reaction_driver.jl is currently the easiest to use, as is plug and play, the rest require path and file size specifications)

use Plotting_for_PSP.jl to produce the required plots
