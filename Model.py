import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from typing import Callable

class Model:
    """
    """
    def __init__(self, matrix: np.ndarray, rates: np.ndarray):
        self.matrix = matrix
        self.rates = rates
        self.expressions = np.dot(self.matrix, self.rates)
    
    def Normalize(self, x, y):
        return x/(np.sqrt(x**2 + y**2)), y/(np.sqrt(x**2 + y**2))
    
    def PhasePortrait(self, xFunc: Callable, yFunc: Callable,
                      xArgs: tuple, yArgs: tuple,
                      X: list, Y: list,
                      scale: int=10, nullcline: bool=True):
        
        DX, DY = self.Normalize(xFunc(*xArgs), yFunc(*yArgs))
        XScale, YScale = np.array([row[::scale] for row in X[::scale]]), np.array([row[::scale] for row in Y[::scale]])
        DXScale, DYScale = np.array([row[::scale] for row in DX[::scale]]), np.array([row[::scale] for row in DY[::scale]])
        
        plt.figure()
        plt.xlabel("X Concentration")
        plt.ylabel("Y Concentration")
        plt.title("X-Y Phase Plane with Nullcline")
        plt.quiver(XScale, YScale, DXScale, DYScale, color='grey')
        if nullcline:
            plt.contour(X, Y, DX, levels=[0], linewidths=1.5, colors='C0')
            plt.contour(X, Y, DY, levels=[0], linewidths=1.5, colors='C0')
            
        plt.show()
        return
    
    def Behavior(self, model: Callable, x0: list,
                 start: float, end: float, step: float):
        
        time = np.arange(start, end, step)
        simulation = odeint(model, x0, time)
        
        plt.figure()
        plt.xlabel("Time")
        plt.ylabel("Concentration")
        plt.title("System Behavior")
        for i in range(len(simulation[:,])):
            plt.plot(time, simulation[:, i], linewidth=1.5, label=f"S{i+1}")
        return
    
    def Nullcline(self):
        pass
    
    def Bifurcation(self):
        pass
            