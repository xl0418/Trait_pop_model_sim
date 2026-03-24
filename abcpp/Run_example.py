# Import necessary libraries and functions
from dvtraitsim_shared import DVTreeData, DVParamLiang
from Dvtraitsim_TVP_local import DVSimTVP
import numpy as np

# --- Configuration ---
# Set the directory path where your data is located.
# For this example, we assume 'treedata' is in the same directory as the script.
dir_path = ''
files = dir_path + 'treedata/'

# --- Load Tree Data ---
# DVTreeData loads and processes the phylogenetic tree data from the specified path.
# The 'scalar' argument scales the simulation time. Make it small to reduce running time.
try:
    td = DVTreeData(path=files, scalar=100)
except FileNotFoundError:
    print(f"Error: The directory '{files}' or its contents were not found.")
    print("Please ensure the 'treedata' directory exists and contains the required CSV files.")
    exit()


# --- Set Simulation Parameters ---
# DVParamLiang packages the parameters for the simulation model.
# These are example parameters for a toy simulation.
# gamma: strength of stabilizing selection
# a: strength of competition
# K: carrying capacity
# h: heritability (square root of h^2)
# nu: plasticity
# r: intrinsic rate of increase
# theta: trait optimum
# V00, V01: initial trait variances for the two starting species
# Vmax: maximum trait variance
# inittrait: initial trait value
# initpop: initial population size
# initpop_sigma: standard deviation of initial population size
# break_on_mu: if True, stops simulation if mean population size is < 1
param = DVParamLiang(
    gamma=0.001,
    a=0.001,
    K=100000,
    h=0.5,
    nu=0.001,
    r=1.0,
    theta=0.0,
    V00=1.0,
    V01=1.0,
    Vmax=10.0,
    inittrait=0.0,
    initpop=10000,
    initpop_sigma=10,
    break_on_mu=True
)

# --- Run Simulation ---
# DVSimTVP runs the trait-population co-evolution model.
print("Running DVSimTVP simulation...")
# The simulation can be computationally intensive.
# We can add a progress indicator or just a message that it's running.
result = DVSimTVP(td, param)
print("Simulation finished.")

# --- Display Results ---
# The result is a dictionary containing the state of the simulation at its end.
print("\n--- Simulation Results ---")
for key, value in result.items():
    print(f"{key}:")
    # To avoid printing very large arrays, we can print their shape or a summary
    if isinstance(value, np.ndarray):
        print(f"  Shape: {value.shape}")
        # Optionally, print a small snippet of the array
        print(f"  Snippet: {value[:5] if value.ndim == 1 else value[:2, :5]}")
    else:
        print(f"  {value}")
print("--------------------------")
