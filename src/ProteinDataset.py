from abc import ABCMeta, abstractmethod
from typing import Tuple, List
from src.utils.Types import ProteinDataset
from src.utils.Protein import Protein
import os
import CDPL.Chem as Chem

PROTEIN_BENCHMARK = 'data/benchmark5.5/structures/'


class BaseProteinDataset(metaclass=ABCMeta):

    def __init__(self):
        pass

    def __call__(self, *args, **kwargs) -> ProteinDataset:
        """
        Load dataset.
        """
        return self.load(*args, **kwargs)

    def __iter__(self):
        raise NotImplementedError

    @abstractmethod
    def load(self, *args, **kwargs) -> ProteinDataset:
        """
        Abstract method to implement load a dataset. This function is called when calling self.
        :param args:
        :param kwargs:
        :return: Proteins
        """
        raise NotImplementedError


class BenchmarkDataset(BaseProteinDataset):

    def __init__(self, path: str = None, bound: bool = True, randomize: bool = True):
        self.path = PROTEIN_BENCHMARK if path is None else path
        self.bindingType = 'b' if bound else 'u'
        self.randomize = randomize
        super(BenchmarkDataset, self).__init__()

    def __iter__(self):
        pdbCodes = self.getPDBCodes()
        for pdbCode in pdbCodes:
            p1 = self.loadProtein('{}{}_l_{}.pdb'.format(self.path, pdbCode, self.bindingType))
            p2 = self.loadProtein('{}{}_r_{}.pdb'.format(self.path, pdbCode, self.bindingType))

            if self.randomize:
                randomized1 = p1.makeRandomRotation().makeRandomTranslation()
                randomized2 = p2.makeRandomRotation().makeRandomTranslation()

            else:
                randomized1 = Protein(p1)
                randomized2 = Protein(p2)

            yield (randomized1, randomized2), (p1, p2)

    def load(self, path: str = None, bound: bool = True, **kwargs) -> ProteinDataset:
        if path is None:
            path = self.path

        bindingType = 'b' if bound else 'u'
        proteinPairs = []
        pdbCodes = self.getPDBCodes(**kwargs)
        for pdbCode in pdbCodes:
            p1 = self.loadProtein('{}{}_l_{}.pdb'.format(path, pdbCode, bindingType), **kwargs)
            p2 = self.loadProtein('{}{}_r_{}.pdb'.format(path, pdbCode, bindingType), **kwargs)

            randomized1 = p1.makeRandomRotation().makeRandomTranslation()
            randomized2 = p2.makeRandomRotation().makeRandomTranslation()

            proteinPairs.append(((randomized1, randomized2), (p1, p2)))

        return proteinPairs

    def loadProtein(self, path: str, **kwargs) -> Protein:
        protein = Protein()
        protein.fromFile(path).prepare()
        return protein

    def getPDBCodes(self, path: str = None, **kwargs) -> List[str]:
        if path is None:
            path = self.path

        pdbCodes = set()
        for fileName in os.listdir(path):
            pdbCode = fileName.split('_')[0]
            if pdbCode in pdbCodes:
                continue
            else:
                pdbCodes.add(pdbCode)

        return list(pdbCodes)
