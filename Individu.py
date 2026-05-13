from __future__ import annotations

from typing import Dict, Tuple

from Allele import Allele


class Individu:
    """
    Représente un profil génétique composé de génotypes sur chaque marqueur que l'on considère
    """
    def __init__(self, genotypes : Dict[str, Tuple[Allele, Allele]]):
        """
        initialise la valeur du génotype
        :param genotypes: Dict[str, Tuple[Allele, Allele]] : dictionnaire des génotypes prenant en clé le nom du
        marqueur et en valeur un tuple d'allèle qui est le génotype.
        """
        self.__genotypes = genotypes

    @property
    def genotypes(self) -> Dict[str, Tuple[Allele, Allele]]:
        """
        Retourne le dictionnaire des génotypes
        :return: Dict[str, Tuple[Allele, Allele]] : dictionnaire des génotypes
        """
        return self.__genotypes

    def clone(self) -> Individu:
        """
        Copie l'instance actuelle en la déliant de la référence
        :return: Individu : copie de l'instance actuelle
        """
        temp = []
        for gen in self.genotypes.values():
            temp.append((gen[0].clone(), gen[1].clone()))
        zip_obj = zip(self.genotypes.keys(), temp)
        return Individu(dict(zip_obj))
