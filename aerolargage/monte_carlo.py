"""Monte Carlo de dispersion : empreinte d'impact, contour de confinement, zone de sauvegarde, sensibilités.

Usage (depuis la racine du projet) :
    python -m aerolargage.monte_carlo --config config/reference.toml
    python -m aerolargage.monte_carlo --config config/reference.toml --n 25000 --graine 2 --workers 8

Sorties dans outputs/ : C_tirages.csv (entrées + impacts), C_empreinte.png, C_sensibilite.png.
"""

import argparse
import os
import time
import tomllib
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from scipy.spatial import ConvexHull
from scipy.stats import binomtest, spearmanr

from aerolargage.trajectoire import simuler

MODES = ("nominal", "partiel", "non-ouverture")

# Entrées continues dispersées (clés du dict de tirages), dans l'ordre d'affichage
VARIABLES = ("x0", "y0", "z0", "V_avion", "w_ref", "psi_from_deg", "cxs_max", "t_ouv", "t_rempl")


# ---------------------------------------------------------------- tirages


def charger_config(chemin):
    with open(chemin, "rb") as fic:
        return tomllib.load(fic)


def tirer(cfg, n, rng):
    """Tire n jeux d'entrées. Retourne un dict de tableaux de taille n.

    Lois : normales (écarts-types de [dispersion]), bornées par écrêtage (np.clip) là où
    une valeur négative n'a pas de sens. Mode de fonctionnement tiré selon [pannes] :
    0 = nominal, 1 = gonflement partiel, 2 = non-ouverture.
    """
    nom, disp, pan = cfg["nominal"], cfg["dispersion"], cfg["pannes"]
    t = {
        "x0": rng.normal(0.0, disp["x0"], n),
        "y0": rng.normal(0.0, disp["y0"], n),
        "z0": rng.normal(nom["z0"], disp["z0"], n),
        "V_avion": rng.normal(nom["V_avion"], disp["V_avion"], n),
        "w_ref": np.clip(rng.normal(nom["w_ref"], disp["w_ref"], n), 0.0, None),
        "psi_from_deg": rng.normal(nom["psi_from_deg"], disp["psi_from_deg"], n),
        "cxs_max": nom["cxs_max"] * (1.0 + rng.normal(0.0, disp["cxs_max_rel"], n)),
        "t_ouv": np.clip(rng.normal(nom["t_ouv"], disp["t_ouv"], n), disp["t_ouv_min"], None),
        "t_rempl": np.clip(rng.normal(nom["t_rempl"], disp["t_rempl"], n), disp["t_rempl_min"], None),
    }
    u = rng.random(n)
    mode = np.zeros(n, dtype=int)
    mode[u < pan["p_non_ouverture"] + pan["p_gonflement_partiel"]] = 1
    mode[u < pan["p_non_ouverture"]] = 2
    t["mode"] = mode
    facteur = rng.uniform(pan["facteur_partiel_min"], pan["facteur_partiel_max"], n)
    t["facteur_partiel"] = np.where(mode == 1, facteur, 1.0)
    return t


def parametres(cfg, t, i):
    """Construit (p, z0, V_avion) pour le tirage i, au format attendu par trajectoire.simuler."""
    nom = cfg["nominal"]
    mode = t["mode"][i]
    p = dict(
        m=nom["m"], cxs_charge=nom["cxs_charge"], n=nom["n"], z_sol=nom["z_sol"],
        t_ouv=np.inf if mode == 2 else t["t_ouv"][i],          # non-ouverture : la voilure ne s'ouvre jamais
        t_rempl=t["t_rempl"][i],
        cxs_max=max(t["cxs_max"][i] * t["facteur_partiel"][i], nom["cxs_charge"]),
        w_ref=t["w_ref"][i],
        psi_from=np.radians(t["psi_from_deg"][i]),
    )
    return p, t["z0"][i], t["V_avion"][i]


# ---------------------------------------------------------------- exécution


def _trajectoire(args):
    """Une trajectoire -> (t_impact, x, y, z, vx, vy, vz) à l'impact. Fonction de module pour multiprocessing."""
    p, z0, V = args
    sol = simuler(p, z0, V)
    if sol.status != 1:   # 1 = arrêt sur événement (impact) ; sinon t_max atteint ou échec du solveur
        raise RuntimeError(f"Pas d'impact détecté (status={sol.status}) pour p={p}, z0={z0}, V={V}")
    return (sol.t_events[0][0], *sol.y_events[0][0])


def executer(cfg, t, workers=1):
    """Lance toutes les trajectoires. Retourne un dict de tableaux d'impacts.

    Le point de largage dispersé (x0, y0) est ajouté APRÈS intégration : le modèle est invariant
    par translation horizontale (vent et atmosphère ne dépendent que de z), donc décaler le point
    de largage décale l'impact d'autant.
    """
    n = len(t["mode"])
    taches = [parametres(cfg, t, i) for i in range(n)]
    if workers > 1:
        with ProcessPoolExecutor(max_workers=workers) as ex:
            res = list(ex.map(_trajectoire, taches, chunksize=max(1, n // (8 * workers))))
    else:
        res = [_trajectoire(a) for a in taches]
    res = np.array(res)
    return dict(t_imp=res[:, 0], x=res[:, 1] + t["x0"], y=res[:, 2] + t["y0"],
                vx=res[:, 4], vy=res[:, 5], vz=res[:, 6])


# ---------------------------------------------------------------- analyse


def mahalanobis(xy, centre, cov):
    e = xy - centre
    return np.sqrt(np.einsum("ij,jk,ik->i", e, np.linalg.inv(cov), e))


def contour_confinement(xy, probabilite, masque_forme):
    """Ellipse empirique contenant la fraction `probabilite` des impacts.

    Forme (centre, orientation, rapport d'axes) : médiane et covariance des impacts `masque_forme`
    (mode nominal). Taille k : quantile empirique de la distance de Mahalanobis sur TOUS les impacts,
    pannes comprises -> le niveau de probabilité ne repose pas sur une hypothèse gaussienne.
    Ellipse = {p : (p - centre)^T cov^-1 (p - centre) <= k^2}.
    """
    centre = np.median(xy[masque_forme], axis=0)
    cov = np.cov(xy[masque_forme].T)
    k = np.quantile(mahalanobis(xy, centre, cov), probabilite)
    return centre, cov, k


def points_ellipse(centre, cov, k, npts=200):
    vals, vecs = np.linalg.eigh(cov)
    a = np.linspace(0.0, 2.0 * np.pi, npts, endpoint=False)
    cercle = np.stack([np.cos(a), np.sin(a)])
    return (centre[:, None] + vecs @ (k * np.sqrt(vals)[:, None] * cercle)).T


def zone_sauvegarde(contour, marge, npts=72):
    """Contour (+) disque de rayon `marge` (somme de Minkowski). Retourne (sommets anti-horaires, aire)."""
    a = np.linspace(0.0, 2.0 * np.pi, npts, endpoint=False)
    disque = marge * np.stack([np.cos(a), np.sin(a)], axis=1)
    nuage = (contour[:, None, :] + disque[None, :, :]).reshape(-1, 2)
    hull = ConvexHull(nuage)
    return nuage[hull.vertices], hull.volume      # en 2D, .volume = aire ; sommets anti-horaires


def dans_polygone_convexe(xy, sommets):
    """Booléen par point : à l'intérieur d'un polygone convexe aux sommets anti-horaires."""
    a, b = sommets, np.roll(sommets, -1, axis=0)
    cross = ((b[:, 0] - a[:, 0])[None, :] * (xy[:, 1:2] - a[None, :, 1])
             - (b[:, 1] - a[:, 1])[None, :] * (xy[:, 0:1] - a[None, :, 0]))
    return np.all(cross >= 0.0, axis=1)


def sensibilites(t, imp, masque):
    """Corrélation de rang de Spearman entrée -> impact (x axe de vol, y travers), sur le mode nominal.

    Classement par max(|rho_x|, |rho_y|). Limite : Spearman ne capte que les effets monotones.
    """
    res = [(v, spearmanr(t[v][masque], imp["x"][masque]).statistic,
            spearmanr(t[v][masque], imp["y"][masque]).statistic) for v in VARIABLES]
    return sorted(res, key=lambda r: -max(abs(r[1]), abs(r[2])))


def intervalle_cp(k, n, niveau=0.95):
    """Intervalle de Clopper-Pearson (exact) sur une proportion k/n."""
    ci = binomtest(k, n).proportion_ci(confidence_level=niveau, method="exact")
    return ci.low, ci.high


# ---------------------------------------------------------------- sorties


def sauver_csv(chemin, t, imp):
    cols = list(VARIABLES) + ["mode", "facteur_partiel"] + ["t_imp", "x", "y", "vx", "vy", "vz"]
    data = np.column_stack([t[c] for c in cols[:len(VARIABLES) + 2]] + [imp[c] for c in cols[len(VARIABLES) + 2:]])
    np.savetxt(chemin, data, delimiter=",", header=",".join(cols), comments="", fmt="%.6g")


def tracer(dossier, cfg, t, imp, contour, zone, sens):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # Palette de référence (skill dataviz) : 3 premières couleurs catégorielles, encre, surface
    serie = ("#2a78d6", "#eb6834", "#1baf7a")
    encre, encre2, surface = "#0b0b0b", "#52514e", "#fcfcfb"
    plt.rcParams.update({"figure.facecolor": surface, "axes.facecolor": surface, "axes.edgecolor": encre2,
                         "axes.labelcolor": encre, "xtick.color": encre2, "ytick.color": encre2,
                         "text.color": encre, "grid.color": "#e3e2dd", "grid.linewidth": 0.6})
    mc = cfg["monte_carlo"]

    # --- Empreinte
    fig, ax = plt.subplots(figsize=(9, 8), constrained_layout=True)
    zp = np.vstack([zone, zone[:1]])
    ax.fill(zp[:, 0], zp[:, 1], color="#e9e8e3", zorder=0, label=f"zone de sauvegarde (contour + {mc['marge']:.0f} m)")
    ax.plot(zp[:, 0], zp[:, 1], color=encre2, lw=1.2, ls="--", zorder=1)
    marqueurs = ("o", "s", "^")
    for k, nom in enumerate(MODES):
        sel = t["mode"] == k
        if sel.any():
            ax.scatter(imp["x"][sel], imp["y"][sel], s=6 if k == 0 else 22, marker=marqueurs[k], color=serie[k],
                       alpha=0.35 if k == 0 else 0.9, linewidths=0, zorder=2 + k, label=f"{nom} ({sel.sum()})")
    cp = np.vstack([contour, contour[:1]])
    ax.plot(cp[:, 0], cp[:, 1], color=encre, lw=2, zorder=6,
            label=f"contour de confinement {100 * mc['probabilite_confinement']:g} %")
    ax.plot(0, 0, marker="X", ms=11, color=encre, ls="none", zorder=7, label="point de largage nominal")
    # Rose des directions dans le coin haut gauche (coordonnées relatives aux axes)
    psi = np.radians(cfg["nominal"]["psi_from_deg"])
    for (dx, dy), txt, y0 in (((1.0, 0.0), "axe de vol", 0.95),
                              ((-np.cos(psi), -np.sin(psi)), f"vent nominal {cfg['nominal']['w_ref']:g} m/s (souffle vers)", 0.89)):
        o = np.array([0.09, y0])
        ax.annotate("", xy=o + 0.05 * np.array([dx, dy]), xytext=o - 0.05 * np.array([dx, dy]), xycoords="axes fraction",
                    arrowprops=dict(arrowstyle="-|>", color=encre, lw=1.5))
        ax.text(0.16, y0, txt, transform=ax.transAxes, va="center", color=encre2, fontsize=9)
    ax.set(xlabel="x, axe de vol [m]", ylabel="y, travers [m]", aspect="equal",
           title=f"Empreinte d'impact - {len(t['mode'])} tirages (graine {mc['graine']})")
    ax.grid(True)
    ax.legend(loc="best", fontsize=9, framealpha=0.95)
    fig.savefig(dossier / "C_empreinte.png", dpi=130)
    plt.close(fig)

    # --- Sensibilités + vitesse verticale d'impact
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    noms = [s[0] for s in sens][::-1]
    yy = np.arange(len(noms))
    a1.barh(yy + 0.2, [s[1] for s in sens][::-1], height=0.38, color=serie[0], label="impact x (axe de vol)")
    a1.barh(yy - 0.2, [s[2] for s in sens][::-1], height=0.38, color=serie[1], label="impact y (travers)")
    a1.set_yticks(yy, noms)
    a1.axvline(0, color=encre2, lw=0.8)
    a1.set(xlim=(-1, 1), xlabel="corrélation de rang de Spearman [-]", title="Sensibilités (mode nominal)")
    a1.grid(True, axis="x")
    a1.legend(fontsize=9, loc="lower right")

    sel = t["mode"] == 0
    a2.hist(-imp["vz"][sel], bins=40, color=serie[0])
    a2.set(xlabel="vitesse verticale d'impact [m/s]", ylabel="nombre de tirages",
           title="Vitesse verticale d'impact (mode nominal)")
    a2.grid(True, axis="y")
    fig.savefig(dossier / "C_sensibilite.png", dpi=130)
    plt.close(fig)


# ---------------------------------------------------------------- programme principal


def main(argv=None):
    ap = argparse.ArgumentParser(description="Monte Carlo de dispersion d'un largage de charge lourde.")
    ap.add_argument("--config", default="config/reference.toml", help="fichier TOML de configuration")
    ap.add_argument("--n", type=int, help="nombre de tirages (remplace [monte_carlo].n_tirages)")
    ap.add_argument("--graine", type=int, help="graine aléatoire (remplace [monte_carlo].graine)")
    ap.add_argument("--workers", type=int, default=os.cpu_count(), help="nombre de processus (1 = séquentiel)")
    ap.add_argument("--sortie", default="outputs", help="dossier de sortie")
    args = ap.parse_args(argv)

    cfg = charger_config(args.config)
    mc = cfg["monte_carlo"]
    if args.n is not None:
        mc["n_tirages"] = args.n
    if args.graine is not None:
        mc["graine"] = args.graine
    n, P = mc["n_tirages"], mc["probabilite_confinement"]

    rng = np.random.default_rng(mc["graine"])
    t = tirer(cfg, n, rng)
    t0 = time.perf_counter()
    imp = executer(cfg, t, workers=args.workers)
    duree = time.perf_counter() - t0

    xy = np.column_stack([imp["x"], imp["y"]])
    nominal = t["mode"] == 0
    centre, cov, k = contour_confinement(xy, P, nominal)
    contour = points_ellipse(centre, cov, k)
    zone, aire = zone_sauvegarde(contour, mc["marge"])
    dedans = dans_polygone_convexe(xy, zone)
    n_dehors = int((~dedans).sum())
    lo, hi = intervalle_cp(n_dehors, n)

    # Vérification de la couverture sur données indépendantes : contour calé sur la 1re moitié,
    # fraction contenue mesurée sur la 2e moitié (doit être proche de P).
    h = n // 2
    c1, cov1, k1 = contour_confinement(xy[:h], P, nominal[:h])
    couverture_test = float(np.mean(mahalanobis(xy[h:], c1, cov1) <= k1))

    sens = sensibilites(t, imp, nominal)
    vh = np.hypot(imp["vx"], imp["vy"])

    dossier = Path(args.sortie)
    dossier.mkdir(exist_ok=True)
    sauver_csv(dossier / "C_tirages.csv", t, imp)
    tracer(dossier, cfg, t, imp, contour, zone, sens)

    demi_axes = k * np.sqrt(np.linalg.eigvalsh(cov))
    print(f"=== Monte Carlo : {n} tirages, graine {mc['graine']}, {args.workers} processus, {duree:.1f} s ===")
    print(f"Modes : " + ", ".join(f"{m} {np.sum(t['mode'] == i)}" for i, m in enumerate(MODES)))
    print(f"Centre de l'empreinte nominale : x = {centre[0]:.1f} m, y = {centre[1]:.1f} m (largage en 0,0)")
    print(f"Contour {100 * P:g} % : ellipse de demi-axes {demi_axes[1]:.0f} m x {demi_axes[0]:.0f} m")
    print(f"  couverture vérifiée sur la 2e moitié des tirages : {100 * couverture_test:.2f} % (visé {100 * P:g} %)")
    print(f"Zone de sauvegarde (contour + {mc['marge']:g} m) : aire {aire / 1e6:.3f} km²")
    print(f"  impacts hors zone : {n_dehors}/{n} = {n_dehors / n:.2e}, IC 95 % Clopper-Pearson [{lo:.1e}, {hi:.1e}]")
    for i, m in enumerate(MODES[1:], start=1):
        sel = t["mode"] == i
        if sel.any():
            print(f"  dont mode {m} : {np.sum(sel & ~dedans)}/{sel.sum()} hors zone")
    print("Vitesses d'impact (pires cas, pour le bloc D) :")
    for i, m in enumerate(MODES):
        sel = t["mode"] == i
        if sel.any():
            print(f"  {m:14s} |vz| max {-imp['vz'][sel].min():6.1f} m/s, médiane {-np.median(imp['vz'][sel]):5.1f} m/s"
                  f" | v horizontale max {vh[sel].max():5.1f} m/s")
    print("Sensibilités (Spearman, mode nominal), de la plus forte à la plus faible :")
    for v, rx, ry in sens:
        print(f"  {v:13s} rho_x = {rx:+.2f}   rho_y = {ry:+.2f}")
    print(f"Sorties : {dossier / 'C_tirages.csv'}, {dossier / 'C_empreinte.png'}, {dossier / 'C_sensibilite.png'}")


if __name__ == "__main__":
    main()
