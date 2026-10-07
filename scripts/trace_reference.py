"""Cas de référence du bloc A : trace la trajectoire et les vitesses (figure dans outputs/)."""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from aerolargage.trajectoire import simuler, cxs

# Cas de référence (hypothèses de démonstration, voir journal de décisions)
P_REF = dict(m=5000.0, t_ouv=2.0, t_rempl=3.0, cxs_charge=3.0, cxs_max=1250.0, n=2.0,
             w_ref=0.0, psi_from=0.0, z_sol=0.0)
Z0, V_AVION = 400.0, 65.0

sol = simuler(P_REF, Z0, V_AVION)
t, (x, y, z, vx, vy, vz) = sol.t, sol.y
te, Xe = sol.t_events[0][0], sol.y_events[0][0]
print(f"t_impact = {te:.1f} s | x = {Xe[0]:.1f} m | y = {Xe[1]:.2e} m | v_impact = {Xe[3:]}")

fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(8, 10), constrained_layout=True)

ax1.plot(x, z)
ax1.set(xlabel="x, axe de vol [m]", ylabel="z [m]", title="Trajectoire (plan vertical)")
ax1.set_aspect("equal")
ax1.grid(True)

ax2.plot(t, vx, label="vx (horizontale)")
ax2.plot(t, vz, label="vz (verticale)")
for tt in (P_REF["t_ouv"], P_REF["t_ouv"] + P_REF["t_rempl"]):
    ax2.axvline(tt, color="grey", ls=":")
ax2.set(xlabel="t [s]", ylabel="vitesse sol [m/s]", title="Vitesses")
ax2.legend()
ax2.grid(True)

ax3.plot(t, [cxs(ti, P_REF["t_ouv"], P_REF["t_rempl"], P_REF["cxs_charge"], P_REF["cxs_max"], P_REF["n"]) for ti in t])
ax3.set(xlabel="t [s]", ylabel="Cx·S [m²]", title="Surface de traînée", xlim=(0, 15))
ax3.grid(True)

out = Path("outputs")
out.mkdir(exist_ok=True)
fig.savefig(out / "A_reference.png", dpi=120)
print(f"Figure : {out / 'A_reference.png'}")
