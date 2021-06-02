from abc import ABCMeta, abstractmethod
from src.utils.Protein import Protein


class BaseScoringFunction(metaclass=ABCMeta):

    def __init__(self):
        pass

    def __call__(self, p1: Protein, p2: Protein, *args, **kwargs) -> float:
        return self.score(p1, p2, *args, **kwargs)

    @abstractmethod
    def score(self, p1: Protein, p2: Protein, *args, **kwargs) -> float:
        """

        :param p2:
        :param p1:
        :param args:
        :param kwargs:
        :return:
        """
        raise NotImplementedError


class ShapeScoring(BaseScoringFunction):

    def __init__(self):
        super(ShapeScoring, self).__init__()

    def score(self, p1: Protein, p2: Protein, *args, **kwargs) -> float:
        raise NotImplementedError


class PharmacophoreScoring(BaseScoringFunction):  # TODO: maybe add clash functions etc in the future

    def __init__(self):
        super(PharmacophoreScoring, self).__init__()

    def score(self, p1: Protein, p2: Protein, *args, **kwargs) -> float:
        raise NotImplementedError

# TODO: force field scoring
