"""
The module contains abstract interfaces of the classes.
The interfaces are aimed to be schemas for further classes implementations.
Following these schemas will ensure compatibility of the code with the entire project.
"""
from abc import ABCMeta, abstractmethod


class AbstractStructureDesigner(object):
    """
    The class builds the atomic structure represented by a list of atoms
    and their neighbouring.
    """

    __metaclass__ = ABCMeta

    @abstractmethod
    def get_neighbours(self, query):
        pass

    @property
    @abstractmethod
    def atom_list(self):
        pass

    @property
    @abstractmethod
    def num_of_nodes(self):
        pass

    @property
    @abstractmethod
    def num_of_species(self):
        pass


class AbstractBasis(object):
    """
    The class contains information about sets of quantum numbers and
    dimensionality of the Hilbert space.
    It is also equipped with the member functions translating quantum numbers
    into a raw index and vise versa.
    """

    __metaclass__ = ABCMeta

    @abstractmethod
    def qn2ind(self, qn):
        """
        The member function trasform a dictionary of quantum numbers into a matrix index

        :param qn:
        :type qn:

        :return ind:   index
        :rtype:        int
        """
        pass

    @abstractmethod
    def ind2qn(self, ind):
        """
        The member function trasform a dictionary of quantum numbers into a matrix index

        :param ind:    index
        :type ind:     int

        :return qn:
        :rtype:
        """
        pass

    @property
    @abstractmethod
    def orbitals_dict(self):
        pass
