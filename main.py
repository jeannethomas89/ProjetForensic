from __future__ import annotations
from concurrent.futures import ProcessPoolExecutor, as_completed
import matplotlib
matplotlib.use('TkAgg')
import random
import numpy as np
import matplotlib.pyplot as plt
from math import prod

from matplotlib.widgets import Button, Slider

from Allele import Allele
from Individu import Individu
from ListValueMarker import get_marker_model, marker_names


N = 200 # Correspond à la taille de la population
# On a besoin de N grand, pour que les résultats mathématiques soient respectés
GENERATIONS = 10
N_RUNS = 100 #nombre de fois que l'on répète la simulation pour un paramètre de coalescence donné
THETA_MIN, THETA_MAX = 0.0000001, 1.0 # valeur maximale et minimale que peut prendre theta
N_THETA_PTS = 60 #nombre de subdivision totale de l'intervalle [theta_min, theta_max]


# on stocke la liste des valeurs des allèles avec leur proba pour optimiser l'algorithme
MARKER_CACHE = {m: get_marker_model(m) for m in marker_names()}


def sample_allele_value(marker: str) -> int:
    """
    Réalise un tirage multinomial dont les poids sont les probabilités des différents allèles dans le MARKER_CACHE[marker]
    :param marker: str : marqueur sur lequel on réalise le tirage
    :return: int : un allèle parmi le MARKER_CACHE[marker]
    """
    model = MARKER_CACHE[marker]
    return int(random.choices(model.alleles, weights=model.probs, k=1)[0])


def make_generation0(n: int, markers: list[str]) -> list[Individu]:
    """
    Simule la génération 0
    :param n: int : le nombre d'individus
    :param markers: list[str] : liste des marqueurs considérés
    :return: list[Individu] : génération 0
    """
    generation = []
    for _ in range(n):
        genotypes = {}
        for m in markers:
            a1 = Allele(sample_allele_value(m), m, parent=None)
            a2 = Allele(sample_allele_value(m), m, parent=None)
            genotypes[m] = (a1, a2)
        generation.append(Individu(genotypes))
    return generation


def maybe_mutate(marker: str, allele: Allele, p_mut: float) -> Allele:
    """
    Réalise un test de mutation et si le test est positif, on renvoie l'allèle modifié
    :param marker: str : nom du marqueur sur lequel est l'allèle
    :param allele: Allele : allèle à laquel, on veut faire subir le test
    :param p_mut: float : coefficient de mutation
    :return: Allele : allèle ayant muté ou non
    """
    if random.random() < p_mut:
        return Allele(sample_allele_value(marker), marker, parent=allele, mutated=True)
    return allele


def make_child(p1: Individu, p2: Individu, markers: list[str], p_mut: float) -> Individu:
    """
    Créé un individu enfant associé à ses deux parents en considérent la mutation
    :param p1: Individu : premier parent considéré
    :param p2: Individu : deuxième parent considéré
    :param markers: list[str] : liste des marqueurs considérés
    :param p_mut: float : paramètre de mutation
    :return: Individu : individu enfant
    """
    genotype = {}
    for m in markers:
        g1, g2 = p1.genotypes[m], p2.genotypes[m]

        child_a1_val = random.choice(g1)
        child_a2_val = random.choice(g2)

        a1 = Allele(child_a1_val.value, m, parent=child_a1_val)
        a2 = Allele(child_a2_val.value, m, parent=child_a2_val)

        if p_mut > 0:
            a1 = maybe_mutate(m, a1, p_mut)
            a2 = maybe_mutate(m, a2, p_mut)
        genotype[m] = (a1, a2)
    return Individu(genotype)


def simulate(n: int, nbr_generations: int,
             markers: list[str], p_mut: float) -> list[list[Individu]]:
    """
    Simule des générations de taille n
    :param n: int : taille de la génération
    :param nbr_generations: int : nombre de générations à simuler
    :param markers: list[str] : liste de marqueurs considérés
    :param p_mut: float : coefficient de mutation
    :return: list[list[Individu]] : liste de toutes les générations
    """
    history: list[list[Individu]] = []
    current = make_generation0(n, markers)
    history.append(current)
    for _ in range(nbr_generations - 1):
        next_gen: list[Individu] = []
        for _ in range(n):
            p1 = random.choice(current)
            p2 = random.choice(current)
            if len(current) > 1:
                while p2 is p1:
                    p2 = random.choice(current)
            next_gen.append(make_child(p1, p2, markers, p_mut))
        current = next_gen
        history.append(current)
    return history

def get_freq_map(pop: list[Individu], markers: list[str]) -> dict:
    """
    Stock toutes les fréquences alléliques d'une génération afin d'y accéder plus tard en 0(1) opérations
    :param pop: list[Individu] : population considérée
    :param markers: list[str] : liste des marqueurs
    :return: dict : dictionnaire des fréquences alléliques
    """
    freq_map = {m: {} for m in markers}
    n = len(pop)
    for ind in pop:
        for m in markers:
            a1, a2 = ind.genotypes[m]
            v1, v2 = a1.value, a2.value
            freq_map[m][v1] = freq_map[m].get(v1, 0) + 1
            freq_map[m][v2] = freq_map[m].get(v2, 0) + 1
    for m in markers:
        for val in freq_map[m]:
            freq_map[m][val] /= 2*n
    return freq_map


def num_calcul(marker: str, guilty: Individu, trace: Individu) -> float:
    """
    Calcul le numérateur du LR, on suppose que DoParameter et DiParameter sont constants par marqueur et donc
    pour optimiser, on les stocke une fois
    :param marker: str : le nom du marqueur considéré
    :param guilty: Individu : la personne coupable
    :param trace: Individu : trace sur la scène de crime
    :return: float : le numérateur du LR
    """

    res = 1.0
    g_geno = guilty.genotypes[marker]
    t_geno = trace.genotypes[marker]

    for i in range(2):
        a_g = g_geno[i]
        a_t = t_geno[i]
        x, y = a_g.DoParameter, a_g.DiParameter

        if a_g.value != a_t.value:
            prob = x
            if a_t.value != 0:
                prob *= y / len(MARKER_CACHE[marker].alleles)
            else:
                prob *= (1 - y)
        else:
            prob = (1 - x) * (1 - y)
        res *= prob
    return res


def den_calcul(marker: str, trace: Individu, freq_map_m: dict) -> float:
    """
    Calcul le dénominateur du LR
    :param marker: str : le nom du marqueur
    :param trace: Individu : la trace de l'individu
    :param freq_map_m: dict : fréquence allélique de la population dont fait partie l'individu qui a produit la trace
    :return: float : le dénominateur du LR
    """
    a1_val = trace.genotypes[marker][0].value
    a2_val = trace.genotypes[marker][1].value

    #On suppose que DoParameter et DiParameter sont constants par marqueur et donc pour optimiser on les stocks une fois
    temp_allele = trace.genotypes[marker][0]
    x, y = temp_allele.DoParameter, temp_allele.DiParameter

    f1 = freq_map_m.get(a1_val, 0)
    f2 = freq_map_m.get(a2_val, 0)

    if a1_val == 0 and a2_val == 0:
        return x ** 2 * (1 - y) ** 2
    if a1_val == 0:
        return 2 * x * (1 - y) * ((1 - x) * f2 + x * y * f2)
    if a2_val == 0:
        return 2 * x * (1 - y) * ((1 - x) * f1 + x * y * f1)
    if a1_val != a2_val:
        return 2 * ((1 - x) * f1 + x * y * f1) * ((1 - x) * f2 + x * y * f2)

    return ((1 - x) * f1 + x * y * f1) ** 2


def LR_global(markers: list[str], trace: Individu, suspect: Individu, freq_map: dict) -> float:
    """
    Calcul le produit du LR
    :param markers: list[str] : liste des marqueurs considérés
    :param trace: Individu : la trace de l'individu
    :param suspect: Individu : suspect considéré
    :param freq_map: dict : fréquence allélique de la population dont fait partie l'individu qui a produit la trace
    :return: float : le produit du LR
    """
    lrs = []
    for m in markers:
        d = den_calcul(m, trace, freq_map[m])
        if d > 0:
            lrs.append(num_calcul(m, suspect, trace) / d)
    return prod(lrs) if lrs else 1.0


def _compute_one_theta(args) -> float:
    """
    Calcul la moyenne du log(LR) pour une valeur donnée de Theta en simulant N_RUNS fois
    :param args: list : k : le nom du thread, theta : paramètre de coascendance, markers : liste des marqueurs,
     n : la taille d'une population, gens_count : nombre de générations à simuler, n_runs : nombre de fois qu'on calcule
     le LR
    :return: list : k : nom du thread, float : LR du coupable, float : LR de l'innocent
    """
    k, theta, markers, n, gens_count, n_runs = args
    log_g_runs, log_i_runs = [], []
    mu = ((1/theta) -1)/(2 * N)
    if mu > 1 : mu = 1

    for _ in range(n_runs):
        hist = simulate(n, gens_count, markers, mu)
        pop = hist[-1]
        freq_map = get_freq_map(pop, markers)

        guilty = random.choice(pop)
        trace = guilty.clone()
        for m in markers:
            for a in trace.genotypes[m]:
                if a.DropOut(): a.DropIn(sample_allele_value(m))

        innocent = random.choice([ind for ind in pop if ind is not guilty])

        lr_g = LR_global(markers, trace, guilty, freq_map)
        lr_i = LR_global(markers, trace, innocent, freq_map)

        if lr_g > 0: log_g_runs.append(np.log10(lr_g))
        if lr_i > 0: log_i_runs.append(np.log10(lr_i))

    return k, np.mean(log_g_runs) if log_g_runs else 0.0, np.mean(log_i_runs) if log_i_runs else 0.0


def compute_curves(markers: list[str]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Calcul de manière parallélisé les différentes valeurs du LR en fonction de theta
    :param markers: list[str] : liste des marqueurs
    :return: tuple[np.ndarray, np.ndarray, np.ndarray] : la valeur de theta, les valeurs du log(LR) pour le coupable
    puis pour l'innocent
    """
    theta_vals     = np.linspace(THETA_MIN, THETA_MAX, N_THETA_PTS)
    log10_guilty   = np.zeros(N_THETA_PTS)
    log10_innocent = np.zeros(N_THETA_PTS)
    done           = [0]

    args_list = [
        (k, float(theta), markers, N, GENERATIONS, N_RUNS)
        for k, theta in enumerate(theta_vals)
    ]

    print(f"Calcul parallèle : {N_THETA_PTS} θ × {N_RUNS} runs …")
    with ProcessPoolExecutor() as executor:
        #stockage des références des threads
        futures = {executor.submit(_compute_one_theta, a): a[0] for a in args_list}
        for future in as_completed(futures):
            k, val_g, val_i = future.result()
            log10_guilty[k]   = val_g
            log10_innocent[k] = val_i
            #affichage graphique dans la console
            done[0] += 1
            pct = int(done[0] / N_THETA_PTS * 40)
            bar = "█" * pct + "░" * (40 - pct)
            print(f"  [{bar}]  {done[0]}/{N_THETA_PTS} terminés"
                  f"  (dernier log₁₀(LR_c)={val_g:.2f})", end="\r")

    print("\nTerminé.")
    return theta_vals, log10_guilty, log10_innocent

def build_ui(theta_vals, log10_guilty, log10_innocent):
    """
    Construit l'interface graphique
    :param theta_vals: np.ndarray : valeurs de theta
    :param log10_guilty: np.ndarray : valeurs du log(LR) du coupable
    :param log10_innocent: np.ndarray : valeurs du log(LR) de l'innocent
    """


    BG = "#fdf6f7"
    CARD = "#f0ece9"
    TEXT = "#000000"
    GRID = "#dcd7d4"

    GUILTY_C = "#b55454"
    INNOC_C = "#5d7293"
    SAGE_C = "#7da88d"
    SLIDER_C = "#8da0b6"

    plt.rcParams.update({
        "figure.facecolor": BG,
        "axes.facecolor": CARD,
        "axes.edgecolor": "#bcbcbc",
        "axes.labelcolor": TEXT,
        "axes.titlecolor": TEXT,

        "text.color": TEXT,
        "xtick.color": TEXT,
        "ytick.color": TEXT,
        "font.family": "monospace",

        "legend.facecolor": CARD,
        "legend.edgecolor": GRID,
        "legend.labelcolor": TEXT,

        "grid.color": GRID,
        "grid.linestyle": "-",
        "grid.alpha": 0.6,

        "axes.spines.top": False,
        "axes.spines.right": False,
    })

    fig = plt.figure(figsize=(11, 7))
    fig.patch.set_facecolor(BG)
    ax = fig.add_axes([0.11, 0.22, 0.84, 0.68])

    state = {"view": "guilty"}

    views = {
        "guilty":   (log10_guilty,   GUILTY_C, "LR coupable"),
        "innocent": (log10_innocent, INNOC_C,  "LR innocent"),
    }

    def ylim_for(data):
        """
        Enlève les 2% des valeurs les plus extrêmes de chaque côté
        :param data: np.ndarray : valeur du log(LR)
        :return: renvoie les deux valeurs les plus extrêmes une fois les deux pourcents exclus, plus une petite marge
        pour l'effet visuel
        """
        finite = data[np.isfinite(data)]
        if len(finite) == 0:
            return -1, 1
        lo = np.percentile(finite, 2)
        hi = np.percentile(finite, 98)
        margin = max((hi - lo) * 0.12, 0.5)
        return lo - margin, hi + margin

    def draw_view(name):
        """
        Trace le graphique et affiche les boutons
        :param name: str : nom du graphique
        :return: None
        """
        data, color, label = views[name]
        ax.cla()
        ax.set_facecolor(CARD)
        ax.grid(True)

        ax.plot(theta_vals, data, color=color, lw=2.0, label=label)

        g = slider.val if hasattr(slider, "val") else theta_vals[0]
        ax.axvline(g, color="#f0d060", lw=1.5, linestyle="-.", alpha=0.85)

        # valeur log10(LR) au curseur
        idx = min(np.searchsorted(theta_vals, g), len(theta_vals) - 1)
        val = data[idx]
        ylo, yhi = ylim_for(data)
        ax.annotate(f"log₁₀(LR) = {val:.2f}",
                    xy=(g, val), xytext=(12, 0), textcoords="offset points",
                    color="#f0d060", fontsize=9,
                    bbox=dict(boxstyle="round,pad=0.3", fc=CARD, ec="#f0d060", alpha=0.9))

        ax.set_xlabel("Coefficient de coascendance  θ", labelpad=8, fontsize=11)
        ax.set_ylabel("log₁₀(LR)", labelpad=8, fontsize=11)
        ax.set_title(f"log₁₀(LR) en fonction de θ  —  {label}",
                     fontsize=13, pad=12, color=TEXT)
        ax.legend(facecolor=CARD, edgecolor=GRID, labelcolor=TEXT, fontsize=10)
        ax.set_xlim(theta_vals[0], theta_vals[-1])
        ax.set_ylim(ylo, yhi)
        fig.canvas.draw_idle()

    ax.set_facecolor(CARD)
    
    ax_slider = fig.add_axes([0.15, 0.08, 0.60, 0.03], facecolor=SLIDER_C)
    slider = Slider(ax=ax_slider, label="θ",
                    valmin=THETA_MIN, valmax=THETA_MAX,
                    valinit=THETA_MIN,
                    color="#f0d060", initcolor="none")
    slider.label.set_color(TEXT)
    slider.valtext.set_color(TEXT)
    slider.on_changed(lambda val: draw_view(state["view"]))
    
    btn_guilty_ax = fig.add_axes([0.15, 0.01, 0.28, 0.045])
    btn_innoc_ax  = fig.add_axes([0.47, 0.01, 0.28, 0.045])

    btn_guilty = Button(btn_guilty_ax, "● Coupable",
                        color=GUILTY_C + "99", hovercolor=GUILTY_C + "cc")
    btn_innoc  = Button(btn_innoc_ax,  "● Innocent",
                        color=INNOC_C  + "33", hovercolor=INNOC_C  + "66")
    btn_guilty.label.set_color(TEXT)
    btn_innoc.label.set_color(TEXT)

    def select_guilty(event):
        state["view"] = "guilty"
        btn_guilty.ax.set_facecolor(GUILTY_C + "99")
        btn_innoc.ax.set_facecolor(INNOC_C   + "33")
        draw_view("guilty")

    def select_innocent(event):
        state["view"] = "innocent"
        btn_guilty.ax.set_facecolor(GUILTY_C + "33")
        btn_innoc.ax.set_facecolor(INNOC_C   + "99")
        draw_view("innocent")

    btn_guilty.on_clicked(select_guilty)
    btn_innoc.on_clicked(select_innocent)

    draw_view("guilty")
    plt.show()


if __name__ == "__main__":

    from multiprocessing import freeze_support
    freeze_support()

    markers = marker_names()
    theta_vals, log10_g, log10_i = compute_curves(markers)
    build_ui(theta_vals, log10_g, log10_i)
