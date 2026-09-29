"""Navion matrices printed in Nelson (1998), 2nd ed., pp. 158 and 199.

Entries are transcribed as printed (rounded). Responses are recomputed,
not digitized textbook curves or reconstructed from published eigenvalues.
"""
import numpy as np
from navion_linear_models import integrate_rk4, DT, T_END

NELSON_PRINTED_MATRIX = {
    "longitudinal": np.array([
        [-0.045, 0.036, 0.0, -32.2],
        [-0.369, -2.02, 176.0, 0.0],
        [0.0019, -0.0396, -2.948, 0.0],
        [0.0, 0.0, 1.0, 0.0],
    ]),
    "lateral_directional": np.array([
        [-0.254, 0.0, -1.0, 0.182],
        [-16.02, -8.40, 2.19, 0.0],
        [4.488, -0.350, -0.760, 0.0],
        [0.0, 1.0, 0.0, 0.0],
    ]),
}


def simulate_nelson_response(model, x0, dt=DT, t_end=T_END):
    """Free response with the same state order and units as the Python model."""
    return integrate_rk4(NELSON_PRINTED_MATRIX[model], x0, dt, t_end)
