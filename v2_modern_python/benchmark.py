import numpy as np
import time
from lcr_engine import LCRSimulator

print("Starting Python Vectorized RK4 Benchmark (1,000,000 steps)...")

# 1. Setup
sim = LCRSimulator()
R_values = np.array([10.0]) # Just one R value to keep the race 1-to-1

t_span = (0, 10.0)
n_steps = 1000000
h = (t_span[1] - t_span[0]) / n_steps

initial_state = np.array([[0.0], [0.2]])

# 2. Race
start_time = time.perf_counter()

# Python Vectorized RK4
t_val, data = sim.solve_rk4(t_span, h, initial_state, R_values)

end_time = time.perf_counter()

# 3. Results
print("Python Benchmark Complete!")
print(f"Execution Time: {end_time - start_time:.4f} seconds")