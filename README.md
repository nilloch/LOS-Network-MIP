# LOS-Network-MIP
Julia Implementation of a Mixed-Integer Program with constraints that ensure agents respect obstacles, maintain a connected network,  and produce viable paths with the objective of minimizing the  number of UAS agents in the network. The formulation is implemented in the function ```findStep()``` defined in ```main.jl```. The file ```workspace.jl``` defines a module for use in ```main.jl``` for generating 2-dimensional workspaces with polygonal obstacles.

Created by Collin Hudson as a final project for ASEN 6519 "Optimization: Applications and Algorithms" Fall 2024
