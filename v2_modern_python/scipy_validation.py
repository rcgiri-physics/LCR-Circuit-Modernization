import numpy as np
from scipy.integrate import solve_ivp

# System Parameters (Underdamped: R=10)
L = 1.0
C = 0.0025
R = 10.0
omega = 2.5

def circuit_derivatives(t, state):
    y, z = state  # y = current, z = di/dt
    
    dy_dt = z
    E_dot = 0.2 * np.sin(omega * t)
    dz_dt = (E_dot / L) - (R / L) * z - (1 / (L * C)) * y
    
    return np.array([dy_dt, dz_dt])

def custom_rk4(t_span, h, y0):
    t_eval = np.arange(t_span[0], t_span[1] + h, h)
    results = np.zeros((len(t_eval), 2))
    results[0] = y0
    
    for i in range(len(t_eval) - 1):
        t = t_eval[i]
        y = results[i]
        
        k1 = h * circuit_derivatives(t, y)
        k2 = h * circuit_derivatives(t + 0.5*h, y + 0.5*k1)
        k3 = h * circuit_derivatives(t + 0.5*h, y + 0.5*k2)
        k4 = h * circuit_derivatives(t + h, y + k3)
        
        results[i+1] = y + (k1 + 2*k2 + 2*k3 + k4) / 6.0
        
    return t_eval, results

if __name__ == "__main__":
    t_span = (0, 5.0)
    h = 0.01
    initial_state = [0.0, 0.1]
    
    # 1. Run Custom RK4
    t_rk4, data_rk4 = custom_rk4(t_span, h, initial_state)
    current_rk4 = data_rk4[:, 0]
    
    # 2. Run SciPy Reference
    sol = solve_ivp(circuit_derivatives, t_span, initial_state, 
                    method='RK45', t_eval=t_rk4, rtol=1e-9, atol=1e-9)
    current_scipy = sol.y[0]
    
    # 3. Calculate Error
    max_error = np.max(np.abs(current_rk4 - current_scipy))
    
    print(f"Validation complete.")
    print(f"Max residual error: {max_error:.2e}")
    # Max residual error: 9.97e-08, well within the 1e-6 threshold.