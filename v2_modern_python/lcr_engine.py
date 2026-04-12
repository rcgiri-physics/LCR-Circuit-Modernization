# lcr_engine.py
import numpy as np

class LCRSimulator:
    def __init__(self, L=1.0, C=0.0025, omega=2.5):
        self.L = L
        self.C = C
        self.omega = omega

    def vectorized_derivatives(self, t, state, R_array):
        y, z = state 
        dy_dt = z
        E_dot = 0.2 * np.sin(self.omega * t)
        
        # NumPy automatically calculates this for all R values
        dz_dt = (E_dot / self.L) - (R_array / self.L) * z - (1 / (self.L * self.C)) * y
        
        return np.array([dy_dt, dz_dt])

    def solve_rk4(self, t_span, h, y0_matrix, R_array):
        t_eval = np.arange(t_span[0], t_span[1] + h, h)
        results = np.zeros((len(t_eval), 2, len(R_array)))
        results[0] = y0_matrix
        
        for i in range(len(t_eval) - 1):
            t = t_eval[i]
            state = results[i]
            
            k1 = h * self.vectorized_derivatives(t, state, R_array)
            k2 = h * self.vectorized_derivatives(t + 0.5*h, state + 0.5*k1, R_array)
            k3 = h * self.vectorized_derivatives(t + 0.5*h, state + 0.5*k2, R_array)
            k4 = h * self.vectorized_derivatives(t + h, state + k3, R_array)
            
            results[i+1] = state + (k1 + 2*k2 + 2*k3 + k4) / 6.0
            
        return t_eval, results