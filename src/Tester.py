"""
Testing interface to easily test modules like loss/scoring functions, alignment, orientation, ...

Loads the dataset and  evaluates it based on found interaction sites.

Note: maybe add evaluation based on interaction maps --> predict whether proteins are interacting with each, but
disregard bindinig site.
"""

# just some default functions
from ProteinDataset import BenchmarkDataset
from BindingSite import ShapeMatching
from Scoring import ShapeScoring

from typing import Callable, Iterable

from src.utils.Protein import Protein


class Tester:

    def __init__(self,
                 dataset: Callable[[], Iterable] = None,
                 bindingSiteFinder: Callable[[Protein, Protein], Iterable] = None,
                 scorer: Callable[[Protein, Protein], float] = None,
                 ):
        self.dataset = dataset if dataset is not None else BenchmarkDataset
        self.bindingSiteFinder = bindingSiteFinder if bindingSiteFinder is not None else ShapeMatching
        self.scorer = scorer if scorer is not None else ShapeScoring()

    def run(self, randomize=True, maxSolutions: int = 10, *args, **kwargs):
        proteins = self.dataset(randomize=randomize, *args, **kwargs)
        solutions = []
        for inputProteins, outputProteins in proteins:
            currentSolutions = []
            p1, p2, = inputProteins

            baselineScore = self.getBaseline(*outputProteins)
            currentSolutions.append((*outputProteins, baselineScore))

            solutionNr = 0
            bindingSolutions = self.bindingSiteFinder(p1, p2, *args, **kwargs)
            for alignedP1, alignedP2 in bindingSolutions:
                score = self.scorer(alignedP1, alignedP2, *args, **kwargs)
                currentSolutions.append((alignedP1, alignedP2, score))
                solutionNr += 1
                if solutionNr == maxSolutions:
                    break

            solutions.append(currentSolutions)
        return solutions

    def getBaseline(self, trueProtein1: Protein, trueProtein2: Protein):
        raise NotImplementedError
