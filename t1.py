import numpy as np 
import matplotlib.pyplot as plt
import equil as eq 
#
# Example parameters
L, m, M, t, R = 1.0, 0.8, 2.0, 0.3, 5.0
p1_fixed = 2.0

    # Sweep p2 to see l2(p1_fixed, p2)
p2_grid = np.linspace(0, R, 400)
_, l2_vals = eq.lequil_np(p1_fixed, p2_grid, L, m, M, t, R)

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
p2_star, l1_at_star = eq.p20(p1_fixed, l2_target=0.25, L=L, m=m, M=M, t=t, R=R, eps=1e-6)
print(f"p2* ≈ {p2_star:.6f}, with l1 at solution ≈ {l1_at_star:.6f}")

