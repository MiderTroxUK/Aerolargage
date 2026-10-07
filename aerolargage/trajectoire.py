"""Trajectoire 3-DDL d'une charge sous parachute (masse ponctuelle).

Conventions (unités SI partout) :
    - repère terrestre supposé galiléen, x dans l'axe de vol, y à gauche, z vers le haut ;
    - z = altitude au-dessus du niveau de la mer [m] (requis par l'atmosphère ISA) ;
    - temps t = 0 à la sortie de rampe.
"""

import numpy as np

# Atmosphère standard ISA, troposphère (ISO 2533 / US Standard Atmosphere 1976)
T0 = 288.15      # température au niveau de la mer [K]
P0 = 101325.0    # pression au niveau de la mer [Pa]
L_T = 0.0065     # gradient thermique [K/m]
R_S = 287.05     # constante spécifique de l'air sec [J/(kg.K)]
G0 = 9.80665     # pesanteur standard [m/s²]
H_TROPO = 11000.0 # altitude limite de la troposphère [m]
Z_MIN = -2000.0   # altitude minimale [m]


def rho(z):
    """Masse volumique de l'air ISA [kg/m³] à l'altitude z [m] au-dessus du niveau de la mer.

    Valide dans la troposphère uniquement : {Z_MIN} <= z <= {H_TROPO} m.
    """

    if z <Z_MIN or z > H_TROPO:
        raise ValueError(f"Altitude {z} m hors de la troposphère ({Z_MIN}..{H_TROPO} m)")

    T = T0 - L_T * z
    P = P0 * (T / T0) ** (G0 / (R_S * L_T))
    return P / (R_S * T)


def vent(z, w_ref, psi_from, z_sol=0.0, z_ref=10.0, alpha=1.0 / 7.0):
    """Vecteur vent [wx, wy, 0] [m/s] à l'altitude z [m], profil en loi puissance.

    Profil : |w| = w_ref * ((z - z_sol) / z_ref) ** alpha, nul au sol.

    Paramètres
        z        : altitude au-dessus du niveau de la mer [m]
        w_ref    : vitesse du vent mesurée à z_ref au-dessus du sol [m/s]
        psi_from : direction D'OÙ VIENT le vent (convention météo) [rad],
                   angle mesuré depuis +x, positif vers +y.
                   psi_from = 0 -> vent venant de l'avant (vent de face), il souffle vers -x.
        z_sol    : altitude du sol au-dessus du niveau de la mer [m]
        z_ref    : hauteur de mesure au-dessus du sol [m] (10 m = standard OMM)
        alpha    : exposant du profil [-] (1/7 ~ terrain dégagé, atmosphère neutre ; hypothèse)
    """

    if z < z_sol:
        return np.array([0.0, 0.0, 0.0])

    w = w_ref * ((z - z_sol) / z_ref) ** alpha
    return np.array([-w * np.cos(psi_from), -w * np.sin(psi_from), 0.0])


def cxs(t, t_ouv, t_rempl, cxs_charge, cxs_max, n=2.0):
    """Surface de traînée totale Cx·S(t) [m²].

    Phases :
        t < t_ouv                       -> déploiement : seule la plateforme traîne
        t_ouv <= t < t_ouv + t_rempl    -> gonflement : croissance en puissance n
        t >= t_ouv + t_rempl            -> voilure pleinement gonflée

    Paramètres
        t          : temps depuis la sortie de rampe [s]
        t_ouv      : délai d'ouverture (début du gonflement) [s]
        t_rempl    : temps de remplissage [s], > 0
        cxs_charge : Cx·S de la plateforme seule [m²]
        cxs_max    : Cx·S voilure gonflée (plateforme incluse) [m²]
        n          : exposant de la loi de gonflement [-] (hypothèse, voir journal)
    """

    if t < t_ouv:
        return cxs_charge
    elif t < t_ouv + t_rempl:
        return cxs_charge + (cxs_max - cxs_charge) * ((t - t_ouv) / t_rempl) ** n
    else:
        return cxs_max


def f(t, X, p):
    """Second membre dX/dt = f(t, X) du modèle 3-DDL (masse ponctuelle).

    X = [x, y, z, vx, vy, vz] : position [m] et vitesse sol [m/s].
    p = dict des paramètres :
        m [kg], t_ouv [s], t_rempl [s], cxs_charge [m²], cxs_max [m²], n [-],
        w_ref [m/s], psi_from [rad], z_sol [m]

    Retourne dX/dt = [vx, vy, vz, ax, ay, az] (np.ndarray de taille 6).
    """
    
    x, y, z, vx, vy, vz = X
    m = p["m"]
    t_ouv = p["t_ouv"]
    t_rempl = p["t_rempl"]
    cxs_charge = p["cxs_charge"]
    cxs_max = p["cxs_max"]
    n = p["n"]
    w_ref = p["w_ref"]
    psi_from = p["psi_from"]
    z_sol = p.get("z_sol", 0.0)

    # Masse volumique de l'air
    rho_air = rho(z)

    # Vecteur vent
    w_vec = vent(z, w_ref, psi_from, z_sol=z_sol)

    # Vitesse relative à l'air
    v_rel = np.array([vx, vy, vz]) - w_vec
    v_rel_norm = np.linalg.norm(v_rel)

    # Surface de traînée totale Cx·S(t)
    cxs_t = cxs(t, t_ouv, t_rempl, cxs_charge, cxs_max, n=n)

    # Accélération aérodynamique (traînée)
    a_drag = -0.5 * rho_air * cxs_t / m * v_rel_norm * v_rel

    # Accélération gravitationnelle
    a_gravity = np.array([0.0, 0.0, -G0])

    # Accélération totale
    a_total = a_drag + a_gravity

    return np.array([vx, vy, vz, a_total[0], a_total[1], a_total[2]])
