import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from Error import InvalidModelError

class Model:
    """
    """
    def __init__(self, model: list, x0: list, start=0, stop=10, step=10000):
        self.model = model
        self.x_0 = x0
        self.time = np.linspace(start, stop, num=step)
        self.simulation = odeint(self.model, self.x_0, self.time)
        self.__Validate()

    def __Validate(self):
        if len(self.model) != len(self.x_0):
            raise InvalidModelError(f"InvalidModelError: Model length: {len(self.model)} != x0 length: {len(self.x_0)}")

    def PhasePlane(self):
        # TODO
        return

    def Nullcline(self):
        # TODO
        return

    def Summary(self):
        # TODO
        return
        