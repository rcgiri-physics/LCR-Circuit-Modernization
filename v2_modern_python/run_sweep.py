# run_sweep.py
import numpy as np
import matplotlib.pyplot as plt
import os
from lcr_engine import LCRSimulator  # Importing custom engine!

# 1. Setup the Experiment
sim = LCRSimulator()
R_values = np.array([10.0, 40.0, 60.0]) # Under, Critical, Overdamped

t_span = (0, 10.0)
h = 0.01

# Initial conditions: i(0)=0, i'(0)=0.1 across all 3 R values
initial_state = np.array([
    [0.0, 0.0, 0.0],  
    [0.1, 0.1, 0.1]   
])

# 2. Run the Engine
print("Running Vectorized RK4 Sweep...")
t_val, data = sim.solve_rk4(t_span, h, initial_state, R_values)
currents = data[:, 0, :]  

# 3. Plot the Results
import matplotlib.pyplot as plt

plt.figure(figsize=(10, 6))
plt.plot(t_val, currents[:, 0], label=f"Underdamped (R={R_values[0]} Ω)")
plt.plot(t_val, currents[:, 1], label=f"Critical (R={R_values[1]} Ω)", linestyle="--")
plt.plot(t_val, currents[:, 2], label=f"Overdamped (R={R_values[2]} Ω)", linestyle=":")

plt.title("LCR Current vs Time: Vectorized Parameter Sweep")
plt.xlabel("Time (s)")
plt.ylabel("Current (A)")
plt.legend()
plt.grid(True)
plt.savefig("v2_modern_python/plots/vectorized_sweep_plot.png", dpi=300, bbox_inches='tight')
print("Plot successfully saved as 'vectorized_sweep_plot.png' in your folder.")