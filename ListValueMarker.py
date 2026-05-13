from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List


@dataclass(frozen=True)
class MarkerModel:
    """
    Initialise la structure de données sur un marqueur stockant les allèles avec leurs probabilités
    """

    alleles: List[int]
    probs: List[float]


MARKER_MODELS: Dict[str, MarkerModel] = {
    "CSF1PO": MarkerModel(
        alleles=[6, 7, 8, 9, 10, 11, 12, 13, 14, 15],
        probs=  [0.002, 0.005, 0.004, 0.019, 0.276, 0.321, 0.305, 0.05, 0.017, 0.002],
    ),
    "D2S1338": MarkerModel(
        alleles=[15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27],
        probs=  [0.001, 0.043, 0.174, 0.076, 0.14, 0.134, 0.051, 0.063, 0.111, 0.096, 0.086, 0.022, 0.003],
    ),
    "D3S1358": MarkerModel(
        alleles=[11, 12, 13, 14, 15, 16, 17, 18, 19, 20],
        probs=  [0.001, 0.001, 0.005, 0.102, 0.3, 0.259, 0.201, 0.12, 0.011, 0.001],
    ),
    "D5S818": MarkerModel(
        alleles=[7, 8, 9, 10, 11, 12, 13, 14, 15],
        probs=  [0.002, 0.004, 0.032, 0.063, 0.349, 0.367, 0.171, 0.012, 0.001],
    ),
    "D7S820": MarkerModel(
        alleles=[7, 8, 9, 10, 11, 12, 13, 14],
        probs=  [0.02, 0.178, 0.159, 0.256, 0.21, 0.136, 0.034, 0.006],
    ),
    "D8S1179": MarkerModel(
        alleles=[8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18],
        probs=  [0.01, 0.006, 0.086, 0.068, 0.137, 0.288, 0.228, 0.13, 0.04, 0.005, 0.001],
    ),
    "D10S1248": MarkerModel(
        alleles=[8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19],
        probs=  [0.001, 0.001, 0.003, 0.012, 0.065, 0.269, 0.303, 0.201, 0.116, 0.028, 0.002, 0.001],
    ),
    "D13S317": MarkerModel(
        alleles=[8, 9, 10, 11, 12, 13, 14, 15],
        probs=  [0.121, 0.074, 0.065, 0.302, 0.274, 0.118, 0.043, 0.002],
    ),
    "D16S539": MarkerModel(
        alleles=[5, 8, 9, 10, 11, 12, 13, 14, 15],
        probs=  [0.001, 0.018, 0.154, 0.098, 0.291, 0.266, 0.147, 0.024, 0.001],
    ),
    "D22S1045": MarkerModel(
        alleles=[8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19],
        probs=  [0.001, 0.001, 0.012, 0.129, 0.02, 0.005, 0.048, 0.334, 0.31, 0.129, 0.009, 0.002],
    ),
    "TPOX": MarkerModel(
        alleles=[6, 7, 8, 9, 10, 11, 12],
        probs=  [0.002, 0.003, 0.558, 0.097, 0.066, 0.255, 0.02],
    ),
    "VWA": MarkerModel(
        alleles=[11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21],
        probs=  [0.001, 0.001, 0.003, 0.101, 0.122, 0.229, 0.262, 0.184, 0.084, 0.011, 0.002],
    ),
}

def marker_names() -> List[str]:
    """
    Retourne la liste des noms des marqueurs
    :return: List[str] : liste des noms des marqueurs
    """
    return list(MARKER_MODELS.keys())


def get_marker_model(marker: str) -> MarkerModel:
    """
    Retourne sur un marker donné les allèles avec leur probabilité
    :param marker: str : le nom du marqueur
    :return: MarkerModel : les allèles avec leur probabilité sur le marqueur
    """
    return MARKER_MODELS[marker]
