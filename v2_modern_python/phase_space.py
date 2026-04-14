import numpy as np
import matplotlib.pyplot as plt
from lcr_engine import LCRSimulator

# 1. Setup the Experiment
sim = LCRSimulator()
R_values = np.array([10.0, 40.0, 60.0]) # Under, Critical, Overdamped

# We run this for a bit longer (15 seconds) to clearly see the stable limit cycle
t_span = (0, 15.0)
h = 0.01

initial_state = np.array([
    [0.0, 0.0, 0.0],  # i(0)
    [0.2, 0.2, 0.2]   # di/dt(0)
])

# 2. Run the Engine
print("Generating Phase Space Data...")
t_val, data = sim.solve_rk4(t_span, h, initial_state, R_values)

# Extract Current (i) and its derivative (di/dt)
currents = data[:, 0, :]   # X-axis
di_dts = data[:, 1, :]     # Y-axis

# 3. Plot the Phase Space
plt.figure(figsize=(10, 8))

# We plot the Underdamped case as the primary focus
plt.plot(currents[:, 0], di_dts[:, 0], color='blue', linewidth=1.5, label=f"Underdamped (R={R_values[0]} Ω)")
plt.plot(currents[:, 1], di_dts[:, 1], color='green', linewidth=1.0, linestyle="--", alpha=0.7, label=f"Critical (R={R_values[1]} Ω)")
plt.plot(currents[:, 2], di_dts[:, 2], color='red', linewidth=1.0, linestyle=":", alpha=0.7, label=f"Overdamped (R={R_values[2]} Ω)")

# Mark the starting point
plt.plot(currents[0, 0], di_dts[0, 0], 'ko', markersize=8, label="Initial State (t=0)")

plt.title("Phase Space Analysis: Current vs. Rate of Change")
plt.xlabel("Current, i (A)")
plt.ylabel("di/dt (A/s)")
plt.axhline(0, color='black', linewidth=0.5)
plt.axvline(0, color='black', linewidth=0.5)
plt.legend()
plt.grid(True)

plt.savefig("v2_modern_python/plots/phase_space_plot.png", dpi=300, bbox_inches='tight')
print("Phase space mapping complete. Saved as 'phase_space_plot.png'")