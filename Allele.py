from __future__ import annotations

import random
import ListValueMarker

class Allele:
    """
        Représente un allèle
    """
    def __init__(self, value, marker, parent, DoParameter=0.002, DiParameter=0.002, mutated: bool = False):
        """
        Initialise la valeur de l'allèle, le marqueur où se situe l'allèle, le parent de l'allèle, le paramètre de
        Drop-out qui est initialement fixé à 0.002, le paramètre de Drop-in qui est fixé initialement à 0.002, et le
        s'il y a eu ou non mutation (non par défaut).

        :param value: int
        :param marker: str
        :param parent: Allele
        :param DoParameter: float
        :param DiParameter: float
        :param mutated: bool
        """
        self.__value = value
        self.__marker = marker
        self.__listAlleleMarker = ListValueMarker.get_marker_model(marker).alleles
        self.__parent = parent
        self.__DoParameter = DoParameter
        self.__DiParameter = DiParameter
        self.__mutated = mutated

    @property
    def value(self):
        """
        Retourne la valeur de l'allèle
        :return: valeur de l'allèle
        """
        return self.__value

    @value.setter
    def value(self, newValue):
        """
        Modifie la valeur de l'allèle
        :param newValue: nouvelle valeur
        :return: None
        """
        self.__value = newValue

    @property
    def marker(self) -> str:
        """
        Retourne le nom du marqueur
        :return: str : nom du marqueur
        """
        return self.__marker

    @property
    def parent(self):
        """
        Retourne le parent
        :return: Allele : parent
        """
        return self.__parent

    @property
    def mutated(self) -> bool:
        """
        Retourne s'il y a eu mutation
        :return: bool : vrai s'il y a eu mutation, faux sinon
        """
        return self.__mutated

    @property
    def DoParameter(self):
        """
        Retourne le parametre de Drop-out
        :return: float : parametre de Drop-out
        """
        return self.__DoParameter

    @DoParameter.setter
    def DoParameter(self, newDoParameter):
        """
        Modifie le paramètre de Drop-out
        :param newDoParameter: float : nouveau paramètre de Drop-out
        :return: float : parametre de Drop-out
        """
        self.__DoParameter = newDoParameter

    @property
    def DiParameter(self):
        """
        Retourne le parametre de Drop-in
        :return: float : parametre de Drop-in
        """
        return self.__DiParameter

    @DiParameter.setter
    def DiParameter(self, newDiParameter):
        """
        Modifie le parametre de Drop-in
        :param newDiParameter: float : nouveau parametre de Drop-in
        :return: float : parametre de Drop-in
        """
        self.__DiParameter = newDiParameter

    def DropOut(self):
        """
        Retourne s'il y a eu un événement de Drop-out
        :return: bool : vrai s'il y a eu Drop-out, faux sinon
        """
        U = random.random()
        if self.DoParameter > U:
            self.value = 0
            return True
        return False

    def DropIn(self, value : int):
        """
        Retourne s'il y a eu un événement de Drop-in
        :param value: int : valeur de remplacement à considérer s'il y a eu Drop-in
        :return: int : retourne la nouvelle valeur de l'allèle
        """
        U = random.random()
        if self.DiParameter > U:
            self.value = value
        return self.value

    def clone(self) -> Allele:
        """
        Fait une copie de l'instance actuelle en la déliant de la référence
        :return: Allele : copie de l'instance actuelle
        """
        return Allele(self.value, self.marker, self.parent, DoParameter=self.DoParameter, DiParameter=self.DiParameter, mutated=self.mutated)

    def __eq__(self, other) -> bool:
        """
        Vérifie l'égalité entre deux instances d'allèle
        :param other: Allele : allèle à comparer avec l'instance actuelle
        :return: bool : Vrai si elles possèdent la même valeur, faux sinon
        """
        if isinstance(other, self.__class__):
            return self.value == other.value
        else :
            return False

    def __nq__(self, other) -> bool:
        """
        Vérifie l'inégalité entre deux instances d'allèle
        :param other: Allele : allèle à comparer avec l'instance actuelle
        :return: bool : Faux si elles possèdent la même valeur, vrai sinon
        """
        return not self.__eq__(other)