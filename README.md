# singlecell
This collection of codes is the implementation of the theoretical model for creating a dynamic cytoskeletal network within a prescribed cell geometry

# Getting Started
1. Download all singlecell Matlab codes
2. Run parameter_gen_units_linked_v2.m to generate the parameter file
3. Run the approrpriate single_cell_units_linked_v3_{CASE}.m to implement the model for different cases. These include when you're interested in studying the impact of nucleus placement (NucPlacement) and when cell aspect ratio is changed (ARvaried)
When prompted, you will identify the cell geometry, initial condition, whether the nucleus should be classified as an obstruction that requires an energetic cost to interact with, and (if needed) where the nucleus will be located
4. Results will be saved with the file marker "_results_store.mat".
5. To analyze results, download the Matlab codes in the PostProcessing folder or the associated Figure folder.

Consult the UserGuide for more information.
