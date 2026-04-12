import numpy as np
import matplotlib.pyplot as plt
import os
from lcr_engine import LCRSimulator

# ==========================================
# 1. Load the Legacy Fortran Data
# ==========================================
# UPDATE THIS PATH to point to your actual Fortran output file
fortran_file_path = "v1_legacy_fortran/data/data_complete_ud.dat" 

try:
    fortran_data = np.loadtxt(fortran_file_path)
    t_fortran = fortran_data[:, 1]       # Column 1 is Time
    current_fortran = fortran_data[:, 6] # Column 2 is Current whereas Column 6 is Current in the Underdamped case complete response
except FileNotFoundError:
    print(f"Error: Could not find the file at {fortran_file_path}")
    print("Please update the path in the script to point to your legacy .dat file.")
    exit()

# ==========================================
# 2. Run the New Python Engine
# ==========================================
sim = LCRSimulator()

# We only test the Underdamped case (R=10) to match the legacy file
R_under = np.array([10.0]) 

# Match the exact time span and step size from your Fortran code
t_span = (0, 5.0) 
h = 0.01          
# =====================================================================
# NOTE ON INITIAL CONDITIONS (Transient Matching)
# In the legacy Fortran V1 code, the "Complete Response" (Column 6) 
# was calculated via the superposition of the Natural and Forced responses.
# Because both independent states had an initial velocity of 0.1 A/s, 
# the true initial boundary condition for the Complete Response is the sum:
# di/dt(0)_complete = di/dt(0)_natural + di/dt(0)_forced = 0.1 + 0.1 = 0.2
# =====================================================================
initial_state = np.array([[0.0], [0.2]])

print("Running Python RK4 Engine for comparison...")
t_python, data_python = sim.solve_rk4(t_span, h, initial_state, R_under)
current_python = data_python[:, 0, 0]

# ==========================================
# 3. The Cross-Validation Plot
# ==========================================
plt.figure(figsize=(10, 6))

# Plot Python as a thick solid blue line
plt.plot(t_python, current_python, label="Python Engine (V2)", color='blue', linewidth=2)

# Plot Fortran as red dots. 
# We use 'markevery=5' so it doesn't get too crowded.
plt.plot(t_fortran, current_fortran, 'ro', label="Fortran Engine (V1)", markersize=4, markevery=5)

plt.title("Cross-Validation: Fortran V1 vs. Python V2 (Underdamped)")
plt.xlabel("Time (s)")
plt.ylabel("Current (A)")
plt.legend()
plt.grid(True)

plt.savefig("v2_modern_python/plots/cross_validation_plot.png", dpi=300, bbox_inches='tight')
print("Comparison plot successfully saved as 'cross_validation_plot.png'")