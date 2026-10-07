# Carnet de projet — Démonstrateur zone de sauvegarde & choc à l'impact

Documents courants : journal de décisions, glossaire FR/EN, banque de questions d'entretien.
Tout est non classifié et basé sur des sources publiques ou des hypothèses explicitement marquées.

---

## 1. Journal de décisions

Format : **hypothèse** — pourquoi — limite — avec plus de temps / des données d'essai.

### Bloc A — Trajectoire

| # | Hypothèse | Pourquoi | Limite | Avec plus de temps / données |
|---|---|---|---|---|
| A1 | Début à la sortie de rampe ; extraction et basculement de la plateforme non modélisés | Hors périmètre MVP | Ignore la dynamique d'extraction (parachute extracteur, tip-off) | Modèle d'extraction ou conditions initiales issues d'essais |
| A2 | Voilure équivalente unique (Cx·S total) au lieu d'une grappe | Simplicité | Perte de rendement de grappe, ouvertures non simultanées | Facteur d'efficacité de grappe (Knacke), modèle par voilure |
| A3 | Masse ponctuelle 3-DDL (translations seules) | Suffisant pour le point d'impact | Pas d'oscillation pendulaire, pas d'attitude à l'impact | 6-DDL ou multi-corps si l'attitude à l'impact compte |
| A4 | Négligés : masse du parachute, masse ajoutée, poussée d'Archimède (≈ 0,25 % du poids), Coriolis ; repère terrestre galiléen | Effets faibles sur une descente de ~1 min | Masse ajoutée non négligeable pendant le gonflement (pic de choc) | Modèle de masse ajoutée (Knacke) |
| A5 | Loi de gonflement Cx·S(t) = Cx·S_charge + (Cx·S_max − Cx·S_charge)·τⁿ, n = 2 | Croissance lente puis rapide, continue aux transitions | Le pic du choc d'ouverture dépend fortement de n — source à vérifier | Identifier n et t_rempl sur essais (bloc F) |
| A6 | Cx·S plateforme seule = 3 m² | Ordre de grandeur | Hypothèse non sourcée | Mesure / données constructeur |
| A7 | Atmosphère ISA troposphère, ρ(z) analytique ; domaine étendu à −2 000 m | ISO 2533 ; l'intégrateur évalue f sous le sol avant de localiser l'impact | Atmosphère standard ≠ atmosphère du jour | Radiosondage du jour (T, p, humidité) |
| A8 | Vent en loi puissance α = 1/7, direction constante, nul sous le sol | Modèle courant de couche limite | Pas de rotation du vent avec l'altitude (Ekman), pas de rafales | Profil de vent mesuré (ballon, radiosondage) |
| A9 | Convention vent « d'où il vient » (météo), angle depuis +x vers +y | Convention opérationnelle | Source classique d'erreur de signe | Test unitaire dédié |
| A10 | Intégration RK45 (`solve_ivp`), max_step = 1 s, arrêt sur événement z = z_sol | Pas adaptatif + localisation de l'impact | `max_step` sert à ne pas sauter le gonflement, pas à la précision | Convergence vérifiée à la main : 6 cm sur l'impact, 0,4 % sur v_z (à mettre en test, bloc B) |

### Bloc C — Monte Carlo

| # | Hypothèse | Pourquoi | Limite | Avec plus de temps / données |
|---|---|---|---|---|
| C1 | Valeurs de dispersion (point de largage σ 40/15 m, vent σ 2 m/s et 20°, Cx·S ±5 %, t_ouv σ 0,3 s…) | Démonstration | Aucune donnée réelle | Statistiques d'essais et de prévision météo |
| C2 | Lois normales écrêtées (`np.clip`) | Simplicité | Accumulation de probabilité à la borne (≠ vraie troncature) | Lois tronquées ou ajustées sur données |
| C3 | Pannes : non-ouverture 1 %, gonflement partiel 3 % (Cx·S × U[0,4 ; 0,8]) | Démonstration | Probabilités inventées | Données de fiabilité (non publiques) |
| C4 | Point de largage dispersé ajouté après intégration | Modèle invariant par translation horizontale | Faux si vent/atmosphère dépendaient de (x, y) | — |
| C5 | Contour = ellipse : forme du nuage nominal, taille = quantile empirique sur tous les tirages | Niveau de probabilité sans hypothèse gaussienne | Forme elliptique imposée (zone bimodale mal décrite) | Région de plus haute densité (KDE) |
| C6 | Zone de sauvegarde = contour 99 % ⊕ marge 100 m | Envelope + marge | Marge arbitraire ; le niveau 99 % est piloté par les pannes (1 %) | Justifier marge et niveau par le risque accepté |
| C7 | Sensibilités par corrélation de rang de Spearman sur le mode nominal | Simple, robuste | Ne capte que les effets monotones | Indices de Sobol |

---

## 2. Glossaire FR / EN

| Français | English |
|---|---|
| point de largage calculé | CARP (Computed Air Release Point) |
| zone de sauvegarde | safety footprint |
| charge lourde / largage de charge lourde | heavy load / heavy drop |
| surface de traînée (Cx·S) | drag area (C_D·S) |
| surface nominale | nominal area |
| voilure / grappe | canopy / cluster |
| vitesse de descente (stabilisée) | rate of descent |
| vitesse limite | terminal velocity |
| équilibre dynamique, régime stabilisé | steady state |
| pression dynamique | dynamic pressure |
| atmosphère type internationale | ISA (International Standard Atmosphere) |
| équilibre hydrostatique | hydrostatic equilibrium |
| couche limite atmosphérique | atmospheric boundary layer |
| loi puissance | power law |
| vent traversier / cisaillement de vent | crosswind / wind shear |
| portée | forward throw |
| temps de remplissage | filling time |
| choc à l'ouverture | opening shock |
| masse ajoutée | added mass |
| loi binomiale | binomial distribution |
| intervalle de Clopper-Pearson (exact) | Clopper-Pearson interval |
| échantillonnage préférentiel | importance sampling |
| contour de confinement | containment contour |
| corrélation de rang | rank correlation (Spearman) |
| raideur tangente | tangent stiffness |
| différences centrées | central difference |
| masse concentrée | lumped mass |
| condition de Courant | CFL condition |
| pas de temps critique | critical time step |
| mise à l'échelle de masse | mass scaling |
| test aux limites | boundary testing |
| sensibilité par compensation | cancellation |
| recette | acceptance testing |
| MCO (maintien en condition opérationnelle) | sustainment / maintenance |

---

## 3. Banque de questions d'entretien

Format : question — meilleure réponse après correction.

1. **Quelle surface de traînée pour 5 t descendant à 8 m/s au niveau de la mer ?**
   Équilibre dynamique (a = 0) : ½ρv²·Cx·S = mg ⇒ Cx·S = 2mg/(ρv²) = 98 100 / 78,4 ≈ 1 250 m². Avec Cx ≈ 0,8 : S₀ ≈ 1 560 m², soit 2 à 3 voilures de ~30 m de diamètre → grappe.

2. **3 impacts hors zone sur 1 000 tirages : quelle confiance dans p̂ = 0,003 ?**
   Loi binomiale, σ = √(p(1−p)/N) ≈ 0,0017 → ±58 % à 1σ, ≈ ±100 % à 95 %. Règle : erreur relative ≈ 1/√k (k = nombre d'événements). Peu d'événements → intervalle exact de Clopper-Pearson.

3. **Combien de tirages pour estimer p = 10⁻³ à ±20 % (1σ) ?**
   1/√k = 0,2 ⇒ k = 25 événements ⇒ N = k/p = 25 000.

4. **Implicite ou explicite pour un impact ?**
   Implicite : peu de pas, chacun coûteux (factorisation de la raideur tangente, itérations de Newton), inconditionnellement stable — adapté au quasi-statique. Explicite : différences centrées, masse concentrée, a = M⁻¹(F_ext − F_int) sans inversion — pas très bon marché mais conditionnellement stable (Δt ≤ L_min/c). Impact court et très non linéaire (contact, écrasement, grandes déformations) → explicite.

5. **Pas de temps critique d'un élément acier de 10 mm ?**
   c = √(E/ρ) ≈ 5 170 m/s ⇒ Δt ≈ 1,9 µs ⇒ ~52 000 pas pour 100 ms (~57 000 avec le facteur de sécurité 0,9 de LS-DYNA). Attention : les éléments écrasés du pad réduisent Δt → mass scaling, à contrôler.

6. **Dérive au vent depuis 300 m, descente 8 m/s, vent 10 m/s ?**
   37,5 s × 10 m/s = 375 m. Négligé : portée (vitesse avion pendant déploiement/gonflement), cisaillement du vent. Archimède ≈ 0,25 % du poids : négligeable (à chiffrer, pas « ça dépend »).

7. **Pourquoi écrire la traînée −½ρ·Cx·S·‖v_r‖·v_r et pas −½ρ·Cx·S·v_z² ?**
   (1) Le signe : le carré efface le signe, ‖v_r‖·v_r garde la direction → la traînée s'oppose toujours au mouvement relatif. (2) Le couplage : l'intensité dépend de la vitesse totale, pas de chaque composante séparément.

8. **Pourquoi la vitesse relative à l'air dans la traînée ?**
   La traînée naît du mouvement par rapport à l'air : v_r = v − w(z). C'est par là — et seulement par là — que le vent agit sur la charge.

9. **Pourquoi un modèle 3D et pas 2D ?**
   La zone de sauvegarde est une surface au sol (x, y) ; la direction du vent est dispersée → vent traversier → trajectoire hors du plan. La 2D sert de cas de vérification (sans vent : y = 0 exactement).

10. **Accélération verticale à 1 000 m, voilure gonflée, à l'arrêt puis à 8 m/s ?**
    À l'arrêt : traînée nulle ⇒ a_z = −g (vitesse nulle ≠ accélération nulle). À 8 m/s : traînée 8,88 m/s² ⇒ a_z ≈ −0,93 m/s² ; l'air est moins dense en altitude, la vitesse limite y est plus élevée (≈ 8,4 m/s) : v_lim(z) = √(2mg/(ρ(z)·Cx·S)) varie pendant la descente.

11. **Pourquoi 46,1 s et pas 400/8 = 50 s ?**
    Pendant le déploiement et le début du gonflement (n = 2), quasi chute libre jusqu'à ~21 m/s : 66 m perdus en 5 s (~3 s gagnées) ; puis vitesse limite un peu > 8 m/s en altitude (~1 s). Conséquence : la voilure s'ouvre à ~62 m/s horizontal et ~21 m/s vertical → c'est là que se joue le choc à l'ouverture.

12. **Pourquoi une fonction appelée par l'intégrateur ne doit-elle pas lever d'erreur juste sous le sol ?**
    `solve_ivp` évalue f à des étages intermédiaires qui peuvent passer sous le sol (démo : z = −1,08 m) avant de localiser l'événement. On borne (vent nul, ρ ISA prolongée) ; on teste la borne et juste au-delà.

13. **Comment choisir le niveau de probabilité d'une zone de sauvegarde ?** *(à approfondir demain)*
    Constat du Monte Carlo : avec 1 % de non-ouvertures, un contour à 99 % se place à la frontière de ce mode de panne (18/44 non-ouvertures hors zone). Le niveau se fixe d'abord par le risque accepté face aux modes de défaillance, puis on dimensionne N pour estimer la probabilité de sortie avec un intervalle de confiance.
