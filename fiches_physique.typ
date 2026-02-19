// --- CONFIGURATION DU DOCUMENT ---
#set page(
  paper: "a4", 
  margin: (x: 1cm, y: 1cm), 
  numbering: "1",
)
#set text(font: "Linux Libertine", lang: "fr", size: 9pt)
#set par(justify: true)

// --- STYLE DES TITRES POUR L'INDEX ---
// On définit comment les titres de chapitres doivent apparaître visuellement
#show heading.where(level: 1): it => {
  set align(center)
  v(10pt)
  text(1.2em, weight: "bold", fill: blue)[#it.body]
  line(length: 100%, stroke: 1pt + blue)
  v(5pt)
}

// --- MACROS ---
// La macro chapitre appelle maintenant une "heading" pour être indexée
#let chapitre(nom) = heading(level: 1, outlined: true, nom)

#let loi(nom, corps) = block(stroke: 1pt + red, inset: 8pt, radius: 2pt, width: 100%, [*#nom* : \ #corps])
#let def(mot, corps) = block(fill: rgb("#e6f4ea"), inset: 8pt, radius: 4pt, width: 100%, [*#mot* \ #corps])

// --- PAGE DE TITRE ET INDEX ---
#align(center)[
  #text(2.5em, weight: "bold")[FICHES DE RÉVISION - PHYSIQUE MP2I] \
  #text(1.5em, style: "italic")[Lycée La Martinière Monplaisir]
]

#v(20pt)

// On personnalise l'apparence de l'index
#show outline.entry.where(level: 1): it => {
  v(6pt)
  strong(it)
}

#outline(title: "Sommaire des révisions", indent: auto)

#pagebreak() // On commence les fiches sur une nouvelle page

// --- DÉBUT DU CONTENU EN 2 COLONNES ---
#show: rest => columns(2, rest)

// ==========================================
// MÉCANIQUE DU POINT
// ==========================================
#align(center)[#text(1.5em, weight: "bold", fill: luma(100))[I. MÉCANIQUE DU POINT]]

#chapitre("M1. Cinématique")

*Vecteurs de base :*
Position $arrow(O P)$, Vitesse $arrow(v) = (d arrow(O P)) / (d t)$, Accélération $arrow(a) = (d arrow(v)) / (d t)$

*1. Coordonnées Cartésiennes $(x, y, z)$ :*
$ arrow(v) = dot(x) arrow(u_x) + dot(y) arrow(u_y) + dot(z) arrow(u_z) $
$ arrow(a) = (d^2 x)/(d t^2) arrow(u_x) + (d^2 y)/(d t^2) arrow(u_y) + (d^2 z)/(d t^2) arrow(u_z) $

*2. Coordonnées Cylindriques $(r, theta, z)$ :*
Position : $arrow(O P) = r arrow(u_r) + z arrow(u_z)$
$ arrow(v) = dot(r) arrow(u_r) + r dot(theta) arrow(u_theta) + dot(z) arrow(u_z) $
$ arrow(a) = ( (d^2 r)/(d t^2) - r dot(theta)^2 ) arrow(u_r) + ( 2 dot(r) dot(theta) + r (d^2 theta)/(d t^2) ) arrow(u_theta) + (d^2 z)/(d t^2) arrow(u_z) $

*3. Coordonnées Sphériques $(rho, theta, phi)$ :*
$ arrow(v) = dot(rho) arrow(u_r) + rho dot(theta) arrow(u_theta) + rho sin(theta) dot(phi) arrow(u_phi) $

*Repère de Frenet (Mouvements circulaires) :*
$ arrow(a) = v^2/R arrow(N) + (d v)/(d t) arrow(T) $

#chapitre("M2. Dynamique")

#loi("Principe Fondamental de la Dynamique (PFD)", $ m arrow(a) = sum arrow(F)_i $)
*Quantité de mouvement :* $arrow(p) = m arrow(v)$

*Forces usuelles :*
- *Gravitation :* $arrow(g)(A) = - (G m) / r^2 arrow(u_r)$
- *Force électrique :* $arrow(F) = q arrow(E)$ avec $arrow(E) = - (d V)/(d x) arrow(u_x)$
- *Rappel d'un ressort :* $arrow(F) = -k(l - l_0) arrow(u_x)$
- *Frottement fluide visqueux (Stokes) :* $arrow(f) = - 6 pi eta R arrow(v)$
- *Poussée d'Archimède :* $arrow(cal(A)) = - rho V arrow(g)$

#def("Mouvements Spécifiques", [
  *1. Balistique :*
  Trajectoire : $z = - g / (2 v_0^2 cos^2(alpha)) x^2 + tan(alpha) x + h$
  Portée max (à 45°) : $D = (v_0^2 sin(2 alpha)) / g$
  
  *2. Chute verticale freinée :*
  Vitesse limite : $v_oo = (2 R^2 g (rho_f - rho_b)) / (9 eta)$
  Temps caractéristique : $tau = (2 rho_b R^2) / (9 eta)$

  *3. Oscillateurs (Ressort / Pendule) :*
  Équation type : $(d^2 X)/(d t^2) + omega_0^2 X = 0$
  Ressort : $omega_0 = sqrt(k/m)$ | Pendule : $omega_0 = sqrt(g/l)$
  Période propre : $T_0 = (2 pi) / omega_0$

  *4. Particule dans un champ magnétique :*
  Mouvement circulaire uniforme de rayon : $R = (m v_0) / (q B)$
])

#chapitre("M3. Énergétique")

*Travail et Puissance :*
Travail élémentaire : $delta W = arrow(F) dot arrow(d l) = arrow(F) dot arrow(v) d t$
Puissance : $cal(P) = (delta W) / (d t) = arrow(F) dot arrow(v)$

#def("Forces conservatives & Énergie Potentielle", [
  Une force est conservative si son travail ne dépend pas du chemin suivi. Elle dérive d'une énergie potentielle :
  $ arrow(F) = - (d E_p) / (d x) arrow(u_x) $
  - *Pesanteur :* $E_p = m g z + E_0$
  - *Élastique :* $E_p = 1/2 k (Delta l)^2$
  - *Électrique :* $E_p = q V + E_0$
])

#loi("Théorèmes Énergétiques", [
  *Théorème de l'Énergie Cinétique (TEC) :*
  $ Delta E_c = 1/2 m v_f^2 - 1/2 m v_i^2 = sum W(arrow(F)_i) $
  $ (d E_c) / (d t) = sum cal(P)_i $

  *Théorème de l'Énergie Mécanique (TEM) :*
  Avec $E_m = E_c + E_p$.
  $ Delta E_m = sum W(arrow(F)_"non conservatives") $
  Si seules des forces conservatives travaillent, $E_m$ est constante.
])

*Équilibre et Stabilité :*
Position d'équilibre $x_"eq"$ atteinte si $(d E_p)/(d x) = 0$.
- *Stable (creux) :* $(d^2 E_p)/(d x^2) > 0$
- *Instable (bosse) :* $(d^2 E_p)/(d x^2) < 0$

Approximation harmonique autour d'un équilibre stable :
$ E_p approx E_0 + k/2 (x - x_"eq")^2 $ avec $k = (d^2 E_p)/(d x^2)_(x_"eq")$

// ==========================================
// ÉLECTRICITÉ
// ==========================================
#align(center)[#text(1.5em, weight: "bold", fill: luma(100))[II. ÉLECTRICITÉ]]

#chapitre("EL1. Lois de Base & Théorèmes")

*ARQS (Approximation des Régimes Quasi-Stationnaires) :*
Valable si la taille du circuit $d$ est très inférieure à la longueur d'onde $lambda = c/f$.

#loi("Lois de Kirchhoff", [
  - *Loi des Nœuds :* La somme des intensités entrant dans un nœud est égale à la somme des intensités qui en sortent ($sum i_k = 0$ en algébrique).
  - *Loi des Mailles :* La somme algébrique des tensions le long d'une maille fermée est nulle ($sum u_k = 0$).
])

#def("Composants & Relations Courant-Tension", [
  - *Résistance* $R$ : $u = R i$ (Loi d'Ohm) ; Puissance dissipée : $cal(P) = R i^2$
  - *Condensateur* $C$ : $i = C (d u_C)/(d t)$ ; Continuité de $u_C$
  - *Bobine* $L$ : $u_L = L (d i)/(d t)$ ; Continuité de $i$
])

*Diviseurs et Théorème de Millman :*
- *Diviseur de tension (série) :* $u_k = E (R_k) / (sum R_i)$
- *Diviseur de courant (parallèle) :* $i_k = I_0 (G_k) / (sum G_i)$ avec la conductance $G = 1/R$
- *Théorème de Millman :* Pour un nœud $M$ relié à des branches $k$ (potentiels $V_k$, résistances $R_k$) :
  $ V_M = (sum (V_k) / (R_k) + sum eta_j) / (sum 1 / (R_k)) $
  (avec $eta_j$ les courants des sources idéales de courant).

#chapitre("EL2. Systèmes du 1er Ordre (RC, RL)")

Pour un circuit RC ou RL soumis à un échelon de tension $E$, l'équation différentielle s'écrit toujours sous la forme :
#loi("Équation canonique du 1er ordre", $ tau (d y)/(d t) + y = y_oo $)
Avec :
- *Circuit RC :* $y = u_C$ et temps caractéristique $tau = R C$
- *Circuit RL :* $y = i$ et temps caractéristique $tau = L / R$

*Solution générale (charge / établissement) :*
$ y(t) = y_oo + (y(0) - y_oo) e^(-t / tau) $

#def("Bilans Énergétiques", [
  - *Énergie stockée par un condensateur :* $E_C = 1/2 C u_C^2$
  - *Énergie stockée par une bobine :* $E_L = 1/2 L i^2$
  - En chargeant un condensateur avec un échelon $E$, la moitié de l'énergie fournie par le générateur est inévitablement dissipée par effet Joule dans la résistance (Rendement = 50%).
])

#chapitre("EL3. Systèmes du 2nd Ordre (RLC)")

#loi("Équation canonique du 2nd ordre", $ (d^2 x)/(d t^2) + omega_0/Q (d x)/(d t) + omega_0^2 x = f(t) $)
Pour un *RLC série* ($x = u_C$) : 
- Pulsation propre : $omega_0 = 1 / sqrt(L C)$
- Facteur de qualité : $Q = 1/R sqrt(L/C)$

*Régimes libres (sans second membre) :*
Dépendent du discriminant $Delta$ de l'équation caractéristique $r^2 + (omega_0/Q) r + omega_0^2 = 0$.

#def("Les 3 Régimes de l'Oscillateur Amorti", [
  *1. Régime Pseudo-Périodique ($Q > 1/2$) :* Faible amortissement. Oscillations sinusoïdales amorties par une exponentielle.
  - Pseudo-pulsation : $omega = omega_0 sqrt(1 - 1/(4Q^2))$
  - Pseudo-période : $T = (2 pi) / omega$
  - Décrément logarithmique : $delta = ln( (x(t)) / (x(t+T)) ) = (pi) / (Q sqrt(1 - 1/(4Q^2)))$
  
  *2. Régime Critique ($Q = 1/2$) :* Retour le plus rapide possible à l'équilibre sans oscillation. Solution de la forme $x(t) = (A t + B) e^(-omega_0 t)$.
  
  *3. Régime Apériodique ($Q < 1/2$) :* Fort amortissement. Retour lent à l'équilibre sans oscillation. Deux racines réelles négatives.
])

// ==========================================
// OPTIQUE
// ==========================================
#align(center)[#text(1.5em, weight: "bold", fill: luma(100))[III. OPTIQUE GÉOMÉTRIQUE]]

#chapitre("OG1. Bases & Réfraction")

*Nature de la lumière :*
Dans le vide, la lumière se propage à $c approx 3.00 times 10^8 "m/s"$.
Dans un milieu matériel transparent, homogène et isotrope (THI) :
#def("Indice de réfraction", $ n = c / v >= 1 $)

#loi("Lois de Snell-Descartes", [
  *1. Réflexion :* Le rayon réfléchi est dans le plan d'incidence et $i_r = -i_1$.
  *2. Réfraction :* $n_1 sin(i_1) = n_2 sin(i_2)$
])

*Réflexion totale :*
Possible uniquement si $n_1 > n_2$ (passage vers un milieu *moins* réfringent).
Angle limite d'incidence : $i_"lim" = arcsin(n_2 / n_1)$. Si $i_1 > i_"lim"$, toute la lumière est réfléchie.

#chapitre("OG2. Systèmes Optiques & Lentilles")

*Stigmatisme et Aplanétisme :*
- *Stigmatisme rigoureux :* Tout rayon issu d'un point objet $A$ émerge en passant par un unique point image $A'$. (Rare en pratique, sauf pour le miroir plan).
- *Conditions de Gauss :* Rayons paraxiaux (proches de l'axe optique et peu inclinés par rapport à lui). Assurent un stigmatisme et un aplanétisme approchés.

#def("Lentilles Minces Sphériques", [
  - *Centre optique* $O$ : Un rayon passant par $O$ n'est pas dévié.
  - *Foyers* $F$ (objet) et $F'$ (image) : Symétriques par rapport à $O$.
  - *Distance focale image* : $f' = overline(O F') = - overline(O F)$.
  - *Vergence* : $V = 1/f'$ (en dioptries $delta$, ou $m^(-1)$).
])

#loi("Relations de conjugaison et Grandissement", [
  *Origine au centre (Relation de Descartes) :*
  $ 1/overline(O A') - 1/overline(O A) = 1/f' $
  $ gamma = overline(A' B') / overline(A B) = overline(O A') / overline(O A) $

  *Origine aux foyers (Relation de Newton) :*
  $ overline(F A) dot overline(F' A') = - f'^2 $
  $ gamma = - f' / overline(F A) = - overline(F' A') / f' $
])

#chapitre("OG3. Instruments d'Optique")

*L'Œil :*
Modélisé par une lentille convergente (cristallin) de vergence variable, et un écran (rétine) à distance fixe.
- *Punctum Remotum (PR) :* Point vu net sans accommoder (à l'infini pour un œil normal).
- *Punctum Proximum (PP) :* Point vu net avec accommodation maximale (environ $25 "cm"$).

#def("Lunette Astronomique (Système Afocal)", [
  Composée de deux lentilles convergentes : l'Objectif ($L_1$, grande focale $f'_1$) et l'Oculaire ($L_2$, petite focale $f'_2$).
  - *Afocal :* Un objet à l'infini donne une image à l'infini (l'œil de l'observateur n'accommode pas, ce qui évite la fatigue).
  - *Condition d'afocalité :* Le foyer image de l'objectif $F'_1$ et le foyer objet de l'oculaire $F_2$ sont confondus.
  - *Encombrement (longueur du tube) :* $L = overline(O_1 O_2) = f'_1 + f'_2$.
])

#loi("Grossissement de la Lunette", [
  Le grossissement $G$ est le rapport de l'angle sous lequel on voit l'image à travers l'instrument ($theta'$) sur l'angle sous lequel on voit l'objet à l'œil nu ($theta$).
  $ G = theta' / theta = - f'_1 / f'_2 $
])

// ==========================================
// ONDES
// ==========================================
#align(center)[#text(1.5em, weight: "bold", fill: luma(100))[IV. ONDES ET PROPAGATION]]

#chapitre("OP1. Ondes & Propagation")

*1. Définitions et Célérité* 
- *Onde :* Perturbation locale d'une grandeur physique se propageant de proche en proche.
- *Célérité :* $c = (M_1 M_2) / tau$ où $tau$ est le retard entre deux points.
- *Forme mathématique (1D) :* 
  $ psi(x, t) = f(x - c t) + g(x + c t) quad "ou" quad psi(x, t) = F(t - x/c) + G(t + x/c) $

*2. Ondes Progressives Planes Harmoniques (OPPH)* 
#loi("Expression de l'OPPH", $ psi(x, t) = A cos(omega t - k x + phi) $)
- *Périodicités :* $omega = (2 pi) / T$ (temporelle) et $k = (2 pi) / lambda$ (spatiale).
- *Vitesse de phase :* $v_phi = omega / k = lambda / T$.
- *Déphasage :* $Delta phi = omega tau$ (temporel) et $Delta phi = -k Delta x$ (spatial, sens $x$ croissant).
- *États vibratoires :* En phase si $Delta x = p lambda$ ; en opposition si $Delta x = (2p+1) lambda / 2$.

*3. Effet Doppler (non relativiste)* 
#loi("Fréquence apparente", $ f = f_0 (1 - (arrow(v) dot arrow(u)_x) / c) $)
- *Rapprochement :* $f = f_0 (1 + v/c) > f_0$ (plus aigu/bleu).
- *Éloignement :* $f = f_0 (1 - v/c) < f_0$ (plus grave/rouge).

*4. Ondes Stationnaires* 
#loi("Forme découplée", $ psi(x, t) = A cos(omega t) cos(k x + phi) $)
- *Nœuds :* Amplitude nulle à tout instant.
- *Ventres :* Amplitude maximale.
- *Modes (corde de longueur L) :* $lambda_n = (2 L) / n$ et $f_n = n f_1$.

#chapitre("OP2. Diffraction")

*1. Caractérisation* 
- Se produit si la taille de l'ouverture $a approx lambda$.
- *Étalement angulaire :* $theta approx lambda / a$.
- *Conditions de Fraunhofer :* Observation en champ lointain ($O E >> a^2 / lambda$).

*2. Résolution des instruments* 
- *Tâche d'Airy :* Rayon angulaire du premier zéro $theta approx 1.22 lambda / D_0$.
- *Critère de Rayleigh :* Deux points sont résolus si $Delta theta >= 1.22 lambda / D_0$.

#chapitre("OP3. Interférences")

*1. Interférences à deux ondes* 
- *Différence de phase :* $Delta phi = (2 pi delta) / lambda$ avec $delta$ la différence de marche.
#loi("Relation de Fresnel", $ I = I_1 + I_2 + 2 sqrt(I_1 I_2) cos((2 pi delta) / lambda) $)
- *Contraste :* $C = (I_max - I_min) / (I_max + I_min) = (2 sqrt(I_1 I_2)) / (I_1 + I_2)$.

*2. Interférences lumineuses* 
- *Chemin optique :* $delta(A B) = c tau = integral_A^B n(M) dif l$.
- *Principe de Fermat :* La lumière suit le chemin optique extremum (minimal).

*3. Dispositif des Trous d'Young* 
- *Différence de marche :* $delta = (n a x) / D$.
- *Intensité sur l'écran :* $I = 2 I_0 (1 + cos((2 pi n a x) / (lambda D)))$.
- *Interfrange :* $i = (lambda D) / (n a)$.


// ==========================================
// THERMODYNAMIQUE
// ==========================================
#align(center)[#text(1.5em, weight: "bold", fill: luma(100))[V. THERMODYNAMIQUE]]

#chapitre("TH1. Introduction & Gaz Parfait")

*1. Vocabulaire & Variables d'état*
- *Systèmes :* Ouvert (matière + énergie), Fermé (énergie seule), Isolé (rien).
- *Variables :* *Extensives* (proportionnelles à $n$ : $V, m, U$) | *Intensives* (indépendantes de $n$ : $P, T$).
- *Transformations :* Isobare ($P=$ cste), Isochore ($V=$ cste), Isotherme ($T=$ cste), Adiabatique ($Q=0$), Quasistatique (suite d'états d'équilibre interne).

*2. Grandeurs Physiques*
- *Température :* $T(K) = theta(degree C) + 273,15$.
- *Pression :* Force par unité de surface $P = F/A$ (1 bar = $10^5$ Pa).
- *Masse volumique :* $mu = m/V$.

*3. Modèle du Gaz Parfait (GP)*
#loi("Équation d'état du GP", $ P V = n R T $)
- *Maxwell (Modèle cinétique) :* $P = 1/3 m mu_N v^2$ où $v$ est la vitesse quadratique moyenne.
- *Lien microscopique :* $1/2 m v^2 = d_"lib"/2 k_B T$ et $R = k_B N_A$.

*4. Énergie Interne & Capacité Thermique*
- *1ère loi de Joule :* $U$ du GP ne dépend que de $T$. Pour un GP monoatomique : $U = 3/2 n R T$.
#def("Capacité à volume constant", [
  $ C_V = (partial U) / (partial T)|_V $ 
  - GP monoatomique : $C_V = 3/2 n R$
  - GP diatomique : $C_V = 5/2 n R$ 
])
- *Loi de Dulong-Petit :* Pour les solides à haute température, $C_"mol" approx 3 R$.

#chapitre("TH2. Premier Principe")

*1. Énoncé du Premier Principe*
#loi("Conservation de l'énergie", $ Delta E_c + Delta U = W + Q $)
- *Forme infinitésimale :* $d U = delta W + delta Q$.
- *Travail des forces de pression :* $delta W = -P_"ext" dif V$.

*2. Modes de Transfert Thermique*
#def("Conduction (Loi de Fourier)", [
  $ arrow(j)_(t h) = -lambda arrow("grad")(T) $
  $ Phi_(t h) = -lambda A (d T) / (d x) $ 
  *Résistance thermique :* $R_(t h) = Delta T / Phi_(t h) = e / (lambda A)$ 
])

#def("Convection (Loi de Newton)", [
  $ |Phi_(t h)| = h A |T_p - T_f| $ 
  *Évolution de la température* (sphère dans thermostat) : 
  $ T(t) = T_"ext" + (T_0 - T_"ext") exp(-t/tau) $ avec $tau = (mu c a) / (3 h)$ 
])

#def("Rayonnement", [
  - *Stefan-Boltzmann :* $Phi_(t h) = sigma T^4 A$ 
  - *Loi de Wien :* $lambda_"max" T approx 2,898 times 10^(-3) m dot K$ 
])

*3. Enthalpie & Lois de Laplace*
- *Enthalpie :* $H = U + P V$. Pour une transf. monobare : $Delta H = Q$.
- *Relation de Mayer :* $C_P - C_V = n R$ (pour un GP).
- *Coefficient de Laplace :* $gamma = C_P / C_V$.

#loi("Lois de Laplace (GP, adiabatique, quasistatique)", [
  - $P V^gamma =$ cste
  - $T V^(gamma-1) =$ cste
  - $T^gamma P^(1-gamma) =$ cste 
])


*4. Bilans pour le Gaz Parfait ($n$ moles)*
#table(
  columns: (1fr, 1.2fr, 1.2fr),
  inset: 7pt,
  fill: (x, y) => if y == 0 { luma(230) },
  align: center + horizon,
  [*Transformation*], [*Travail $W$ reçu*], [*Transfert $Q$ reçu*],
  
  [Isochore \ ($V = "cste"$)], 
  [$W = 0$], 
  [$Q = Delta U = C_V Delta T$],
  
  [Isobare \ ($P = "cste"$)], 
  [$W = -n R Delta T$], 
  [$Q = Delta H = C_P Delta T$],
  
  [Isotherme \ ($T = "cste"$)], 
  [$W = n R T ln(V_i / V_f)$], 
  [$Q = -W$],
  
  [Adiabatique \ ($Q = 0$)], 
  [$W = Delta U = C_V Delta T$], 
  [$Q = 0$]
)

// Insère l'image centrée dans ta colonne
#align(center)[
  #image("transfert_thermique.jpeg", width: 90%)
  #text(size: 8pt, style: "italic")[Schéma des modes de transfert : conduction, convection, rayonnement]
]

// ==========================================
// MÉCANIQUE SEMESTRE 2
// ==========================================
#align(center)[#text(1.5em, weight: "bold", fill: luma(100))[VI. MÉCANIQUE DU SOLIDE & GRAVITATION]]

#chapitre("M4. Mouvement de Rotation")

*1. Moment d'une force*
caractérise l'effet d'une force sur la rotation.
#loi("Moment en O", $ arrow(M)_O = arrow(O P) times arrow(F) $)
- *Changement de point :* $arrow(M)_B = arrow(M)_A + arrow(B A) times arrow(F)$
- *Bras de levier (b) :* $M_Delta = b dot F$ (distance minimale à l'axe).

*2. Moment cinétique*
#loi("Moment cinétique", $ arrow(L)_O = arrow(O P) times arrow(p) = arrow(O P) times (m arrow(v)) $)
- *Rotation axe fixe :* $L_Delta = I_Delta omega$ avec $I_Delta$ le moment d'inertie.
- *Moment d'inertie (point) :* $I_Delta = m r^2$.

#loi("Théorème du Moment Cinétique (TMC)", [
  Dans un référentiel galiléen, pour un point O fixe :
  $ (d arrow(L)_O) / (d t) = sum arrow(M)_(O, i) $
  *Conservation :* Si le système est isolé ou soumis à une force centrale, $arrow(L)_O$ est constant.
])

#chapitre("M5. Mouvement d'un Solide")

*1. Moments d'inertie classiques ($I_Delta$)*
#table(
  columns: (1fr, 1fr), fill: (x, y) => if y == 0 { luma(230) },
  [*Corps (masse m)*], [*Moment d'inertie*],
  [Tige (longueur L, centre)], [$1/12 m L^2$],
  [Tige (extrémité)], [$1/3 m L^2$],
  [Cylindre plein / Disque], [$1/2 m R^2$],
  [Sphère pleine], [$2/5 m R^2$]
)

*2. Énergie et Puissance*
- *Énergie cinétique de rotation :* $E_c = 1/2 I_Delta omega^2$.
- *Puissance d'un couple :* $cal(P) = arrow(Gamma) dot arrow(omega)$.

#loi("Analogie Translation / Rotation", [
  #table(
    columns: (1fr, 1fr), stroke: none,
    [*Translation*], [*Rotation (axe fixe)*],
    [Masse $m$], [Moment inertie $I_Delta$],
    [Vitesse $v$], [Vitesse angulaire $omega$],
    [Force $F$], [Moment $M_Delta$],
    [Quantité mouv. $p$], [Moment cinétique $L_Delta$],
    [$P = F v$], [$P = M_Delta omega$]
  )
])

#chapitre("M6. Forces Centrales & Kepler")

*1. Propriétés Fondamentales*
Une force est centrale si $arrow(F) = f(r) arrow(u)_r$.
- *Conservation :* $arrow(L)_O$ est constant $=>$ Mouvement plan.
#loi("Loi des Aires", [
  Le rayon vecteur balaye des surfaces égales en des temps égaux.
  Vitesse aréolaire : $ (d A) / (d t) = C / 2 = (r^2 dot(theta)) / 2 = L_O / (2 m) $
])

*2. Potentiel Effectif & Trajectoires*
$ E_(p, "eff") = (m C^2) / (2 r^2) + E_p (r) $
- *État lié ($E_m < 0$) :* Trajectoire elliptique.
- *État de diffusion ($E_m >= 0$) :* Parabole ou Hyperbole.

*3. Lois de Kepler*
#loi("3ème Loi de Kepler", $ T^2 / a^3 = (4 pi^2) / (G M_"att") $)

*4. Vitesses Cosmiques (Attracteur $M, R$)*
- *1ère vitesse cosmique :* $v_1 = sqrt(g R)$ (orbite circulaire rasante).
- *Vitesse de libération :* $v_2 = sqrt(2 g R)$ (échapper à l'attraction).
