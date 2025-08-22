import numpy as np
import matplotlib.pyplot as plt
#
# Unobservable M/M/1 Edelson plots
#
def mm1_unobs(L, m, R, C, p):
    # L = arrival rate
    # m = service rate
    # R = service reward
    # C = delay cost
    # p = entrance fee
    #
    # Equilibrium analysis
    if p >= R:
        l_eq = 0.0
    else:
        l_eq = min( max(0, m - C/(R - p)), L)
    # Socially optimal analysis
    l_opt = min( max(0, m - np.sqrt(C*m/R)), L)
    U_opt = l_opt * (R - C/(m - l_opt))
    #
    U_cust = l_eq * (R - p - C/(m - l_eq) )
    Profit = l_eq * p
    S= U_cust + Profit
    return l_eq, l_opt, U_cust, Profit, S, U_opt  
#
if __name__ == "__main__":
    m, R, C = 1, 5, 1    
#     # Example parameters 
  # range of L to sweep
    L_vals = np.linspace(0.0, 2.0, 400)

    # compute equilibrium arrival rates for p=1 and p=2, and socially optimal l_opt
    l_eq_p0 = [mm1_unobs(L, m, R, C, 0)[0] for L in L_vals]
    l_eq_p2 = [mm1_unobs(L, m, R, C, 2)[0] for L in L_vals]
    l_eq_p3 = [mm1_unobs(L, m, R, C, 3.5)[0] for L in L_vals]
    l_opt_vals = [mm1_unobs(L, m, R, C, 0)[1] for L in L_vals]  # p unused for l_opt

    # plotting 
    # vs. L
    plt.figure(figsize=(8, 5))
    plt.plot(L_vals, l_eq_p0, label='$\lambda^{eq}$ (p=0)', linewidth=2)
    plt.plot(L_vals, l_eq_p2, label='$\lambda^{eq}$ (p=2)', linewidth=2)
    plt.plot(L_vals, l_eq_p3, label='$\lambda^{eq}$ (p=3.5)', linewidth=2)
    plt.plot(L_vals, l_opt_vals, label='$\lambda^*$ (social optimum)', linestyle='--', linewidth=2)
    plt.xlabel(r'$\Lambda$ (maximum arrival rate)')
    plt.ylabel('Arrival rate $\lambda$')
    plt.title('Equilibrium and socially optimal arrival rates vs $\Lambda$')
    plt.legend()
    #plt.grid(True)
    plt.ylim(0, max(L_vals)*1.05)
    plt.tight_layout()
    plt.savefig('mm1_unobs_arrival_rates.png', dpi=300)
    plt.show() 
    # vs. p for fixed L
    p_vals = np.linspace(0.0, R, 400)
    L_fixed = 0.7
    l_eq_p = [mm1_unobs(L_fixed, m, R, C, p)[0] for p in p_vals]
    l_opt_p = [mm1_unobs(L_fixed, m, R, C, p)[1] for p in p_vals]
    plt.figure(figsize=(8, 5))
    plt.plot(p_vals, l_eq_p, label='$\lambda^{eq}$', linewidth=2)
    plt.plot(p_vals, l_opt_p, label='$\lambda^*$ (social optimum)', linestyle='--', linewidth=2)
    plt.xlabel('Entrance fee $p$')
    plt.ylabel('Arrival rate $\lambda$')
    plt.title(f'Equilibrium and socially optimal arrival rates vs entrance fee $p$ ($\Lambda$={L_fixed})')
    plt.legend()
    #plt.grid(True)
    plt.ylim(0, L_fixed*1.05)
    plt.tight_layout()
    plt.savefig('mm1_unobs_arrival_rates_vs_p.png', dpi=300)
    plt.show()    
    # Utility and profit vs p for fixed L  
    U_cust_p = [mm1_unobs(L_fixed, m, R, C, p)[2] for p in p_vals]
    Profit_p = [mm1_unobs(L_fixed, m, R, C, p)[3] for p in p_vals]
    S_p = [mm1_unobs(L_fixed, m, R, C, p)[4] for p in p_vals]
    U_opt_p = [mm1_unobs(L_fixed, m, R, C, p)[5] for p in p_vals]
    plt.figure(figsize=(8, 5))
    plt.plot(p_vals, U_cust_p, label='Customer utility', linewidth=2)
    plt.plot(p_vals, Profit_p, label='Provider profit', linewidth=2)
    plt.plot(p_vals, S_p, label='Social welfare', linestyle='--', linewidth=2)
    plt.plot(p_vals, U_opt_p, label='Optimal social welfare', linestyle=':', linewidth=2)
    plt.xlabel('Entrance fee $p$')
    plt.ylabel('Value')
    plt.title(f'Utility, profit, and social welfare vs entrance fee $p$ ($\Lambda$={L_fixed})')
    plt.legend()
    #plt.grid(True)
    plt.ylim(0, max(U_opt_p)*1.05)
    plt.xlim(0, R) 
    plt.axvline(x=R, color='red', linestyle=':', label='$p=R$')
    plt.tight_layout()
    plt.savefig('mm1_unobs_util_profit_vs_p.png', dpi=300)
    plt.show()
   


