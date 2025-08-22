import numpy as np
import matplotlib.pyplot as plt
from typing import Tuple, Optional, Callable

# --------- Core (scalar) translations ---------

def _sqrt_nonneg(x: float) -> float:
    """Mimic R's sqrt behavior: return NaN for negatives instead of raising."""
    return float(np.sqrt(x)) if x >= 0 else float("nan")


def lequil(p1: float, p2: float, L: float, m: float, M: float, t: float, R: float) -> Tuple[float, float]:
    """
    Equilibrium Rates (R: lequil). Returns (l1, l2).
    Scalar implementation that mirrors the original branching exactly.
    """
    K1 = (M*(1 - t) - L*(R - p1))**2 + 4*(M*(1 + t) + L*(R - p2))*L*(p1 - p2)
    K2 = (M*(1 - t) + L*(R - p2))**2 - 4*(M*(1 + t) + L*(R - p1))*L*(p1 - p2)

    if (p1 >= R - (R - p2)*t*M/(M + (R - p2)*L)) and (p1 <= R):
        l1 = 0.0
        l2 = (R - p2)*L*m/(M + (R - p2)*L)

    elif (
        (p1 <= R - (R - p2)*t*M/(M - (R - p2)*L) and (p2 >= R - R*M/(t*M + R*L)) and (p2 <= R - (1 - t)*M/L))
        or ((p1 <= p2) and (p2 >= R - (1 - t)*M/L))
    ):
        l1 = (R - p1)*L*m/(t*M + (R - p1)*L)
        l2 = 0.0

    elif (p1 > p2) and (p1 < R - (R - p2)*t*M/(M + (R - p2)*L)):
        sqrtK1 = _sqrt_nonneg(K1)
        l2 = m*(L*(R - p1) - M*(1 - t) + sqrtK1) / (2*(M*(1 + t) + L*(R - p2)))
        l1 = ((R - p1)*m*L - t*M*l2) / (t*M + L*(R - p1))

    elif (
        ((p1 <= p2) and (p2 <= R - R*M/(t*M + R*L)))
        or ((p1 >= R - (R - p2)*t*M/(M - (R - p2)*L)) and (p1 <= p2) and (p2 >= R - R*M/(t*M + R*L)) and (p2 <= R - (1 - t)*M/L))
    ):
        sqrtK2 = _sqrt_nonneg(K2)
        l1 = m*(L*(R - p2) + M*(1 - t) + sqrtK2) / (2*(M*(1 + t) + L*(R - p1)))
        l2 = ((R - p2)*m*L - M*l1) / (M + L*(R - p2))

    elif (p1 == p2) and (p2 < R - (1 - t)*M/L):
        l1 = m*(L*(R - p1) + M*(1 - t)) / (L*(R - p1) + M*(1 + t))
        l2 = m*(L*(R - p2) - M*(1 - t)) / (L*(R - p2) + M*(1 + t))

    else:
        l1 = -1.0
        l2 = -1.0

    return float(l1), float(l2)


def pequil(l1: float, l2: float, L: float, m: float, M: float, t: float, R: float) -> Tuple[float, float]:
    """
    Equilibrium Prices (R: pequil). Returns (p1, p2).
    """
    if l1 == 0:
        p1 = R - t*M*l2/(L*m)
        p2 = R - M*l2/((m - l2)*L)

    elif l2 == 0:
        p1 = R - t*M*l1/((m - l1)*L)
        p2 = R - M*l1/(L*m)

    elif l1 < m - t*(m - l2):
        p1 = R - M*t*(l1 + l2)/((m - l1)*L)
        p2 = p1 - M*l2*(1/(m - l2) - t/(m - l1))/L

    elif l1 > m - t*(m - l2):
        p2 = R - M*(l1 + l2)/((m - l2)*L)
        p1 = p2 + M*l1*(1/(m - l2) - t/(m - l1))/L

    else:
        # Preserved as in the R code.
        p1 = M*(l1 + l1)/((m - l2)*L)
        p2 = p1

    return float(p1), float(p2)


def p20(p1: float, l2_target: float, L: float, m: float, M: float, t: float, R: float, eps: float = 1e-6) -> Tuple[float, Optional[float]]:
    """
    Given p1, find p2 such that l2_e(p1, p2) ≈ l2_target using bisection.
    Returns (p2_star, l1_at_solution). If infeasible at p2=0, returns (-1, l1_at_p2_0).
    """
    l10_at0, l20_at0 = lequil(p1, 0.0, L, m, M, t, R)
    if l20_at0 < l2_target:
        return -1.0, l10_at0

    x0, x1 = 0.0, float(R)
    x = 0.5*(x0 + x1)
    l1x, l2x = lequil(p1, x, L, m, M, t, R)

    # Bisection loop
    while abs(l2x - l2_target) > eps:
        if l2x > l2_target:
            x0 = x
        else:
            x1 = x
        x = 0.5*(x0 + x1)
        l1x, l2x = lequil(p1, x, L, m, M, t, R)

    return float(x), float(l1x)


# --------- NumPy-friendly wrappers ---------

# Vectorized versions (accept numpy arrays; broadcast over p1, p2 etc.)
_lequil_vec: Callable = np.vectorize(lambda p1, p2, L, m, M, t, R: lequil(float(p1), float(p2), float(L), float(m), float(M), float(t), float(R)),
                                     otypes=[float, float])

_pequil_vec: Callable = np.vectorize(lambda l1, l2, L, m, M, t, R: pequil(float(l1), float(l2), float(L), float(m), float(M), float(t), float(R)),
                                     otypes=[float, float])

_p20_vec: Callable = np.vectorize(lambda p1, l2, L, m, M, t, R, eps=1e-6: p20(float(p1), float(l2), float(L), float(m), float(M), float(t), float(R), float(eps)),
                                  otypes=[float, float])


def lequil_np(p1, p2, L, m, M, t, R) -> Tuple[np.ndarray, np.ndarray]:
    """NumPy-friendly: returns arrays l1, l2 with broadcasting."""
    l1, l2 = _lequil_vec(p1, p2, L, m, M, t, R)
    return np.asarray(l1), np.asarray(l2)


def pequil_np(l1, l2, L, m, M, t, R) -> Tuple[np.ndarray, np.ndarray]:
    """NumPy-friendly: returns arrays p1, p2 with broadcasting."""
    p1, p2 = _pequil_vec(l1, l2, L, m, M, t, R)
    return np.asarray(p1), np.asarray(p2)


def p20_np(p1, l2_target, L, m, M, t, R, eps=1e-6) -> Tuple[np.ndarray, np.ndarray]:
    """NumPy-friendly bisection wrapper (vectorized over inputs)."""
    p2_star, l1_at_sol = _p20_vec(p1, l2_target, L, m, M, t, R, eps)
    return np.asarray(p2_star), np.asarray(l1_at_sol)


# --------- Example Matplotlib usage ---------

if __name__ == "__main__":
    # Example parameters (adjust to your model)
    L, m, M, t, R = 1.0, 0.8, 2.0, 0.3, 5.0
    p1_fixed = 2.0

    # Sweep p2 to see l2(p1_fixed, p2)
    p2_grid = np.linspace(0, R, 400)
    _, l2_vals = lequil_np(p1_fixed, p2_grid, L, m, M, t, R)

    plt.figure()
    plt.plot(p2_grid, l2_vals, label=r"$\ell_2(p_1{=}%.2f, p_2)$" % p1_fixed)
    plt.xlabel(r"$p_2$")
    plt.ylabel(r"$\ell_2$")
    plt.title("Equilibrium l2 vs p2")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()

    # Example: find p2* such that l2 = 0.25 for this p1
    p2_star, l1_at_star = p20(p1_fixed, l2_target=0.25, L=L, m=m, M=M, t=t, R=R, eps=1e-6)
    print(f"p2* ≈ {p2_star:.6f}, with l1 at solution ≈ {l1_at_star:.6f}")
