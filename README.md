# interPhaseChangeBubbleFoam

Hi folks,

interPhaseChangeBubbleFoam is a multi-phase solver that can consider bubbles in both the Lagrangian and Eulerian frameworks.

How to use it:

Please source OpenFOAM version 2406, and then run the script ./Allrun. This will compile the solver and execute all tutorial cases. In each case study directory, there is a ParaView state file that can be used to import the pipelines used to visualize the results in this work.

Explanation:

The interPhaseChangeBubbleFoam solver is designed for simulating multiphase flows involving bubbles. It incorporates both the Lagrangian and Eulerian frameworks, which allow for detailed tracking of individual bubbles (Lagrangian) while also capturing the overall fluid flow behavior (Eulerian). This flexibility makes it suitable for a wide range of applications in bubble dynamics and multiphase simulations. The provided tutorial cases and ParaView state files facilitate the visualization and analysis of simulation results, offering a comprehensive workflow for studying bubble behavior in various conditions. 

Please cite the work as: Cite as: Lavari, M.: interPhaseChangeBubbleFoam: A hybrid Eulerian-Lagrangian solver to capture nuclei effects on hydrodynamic cavitation. In Proceedings of CFD with OpenSource Software, 2024, Edited by Nilsson. H., https://dx.doi.org/10.17196/OS_CFD#YEAR_2024.

Thanks,
ML
