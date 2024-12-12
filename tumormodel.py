from math import exp
from random import random

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

class TumorGrowthModel:
    """ Tumor Growth Model"""
    def __init__(self, t_data, V_data):
        self.t_data = t_data
        self.V_data = V_data

    def __str__(self):
        """ String representation"""
        return (f"TumorGrowthModel:\n"
                f"numbers of time points  : {len(self.t_data)},\n"
                f"numbers of observations : {len(self.V_data)}")

    def __repr__(self):
        """ Technical representation"""
        return f"{self.__class__.__name__}: (t_data={self.t_data}, V_data={self.V_data})"

    # Generate a fake data
    @staticmethod
    def generate_fake_data():
        """ Generate fake data"""
        t_data = [random() for _ in range(20)]
        V_data = [exp(data_t) - 1.0 for data_t in t_data]
        return np.array(t_data), np.array(V_data)

    @staticmethod
    def Runge_method(model_growth, t, V0, dt, *params):
        """
        Compute ODE-Solver Runge-Kutta method
        :param model_growth: The growth model
        :param t: Time (in days)
        :param V0: Begin condition of the tumor volume (mm3)
        :param dt: difference in time
        :param params: Parameters
        :return: An array of tumorvolume (mm3)
        """
        V = [V0]
        for i in range(1, len(t)):
            t_current = t[i - 1]
            V_current = V[-1]
            y1 = dt * model_growth(t_current, V_current, *params)
            y2 = dt * model_growth(t_current + dt / 2, V_current + y1 / 2, *params)
            y3 = dt * model_growth(t_current + dt / 2, V_current + y2 / 2, *params)
            y4 = dt * model_growth(t_current + dt, V_current + y3, *params)
            V_new = V_current + (y1 + 2 * y2 + 2 * y3 + y4) / 6
            V.append(V_new)
        return np.array(V)

    @staticmethod
    def fit_model_brute_force(model_instance, t_data, V_data, p0, num_iterations=10000, step_size=0.01):
        """
        Fit model using Brute Force method
        :param model_instance: Model wrapper to be fitted
        :param t_data: Time data
        :param V_data: Tumor volume data
        :param p0: Initial guess for the parameters
        :param num_iterations: Number of iterations
        :param step_size: Number of step size
        :return: The best fitted parameters for the model
        """
        def model_wrapper(t, *params):

            return model_instance.wrapper(t, *params, V0=V_data[0])  # Include V_data as an argument

        best_params = np.array(p0)
        best_cost = np.sum((model_wrapper(t_data, *best_params) - V_data) ** 2)

        for _ in range(num_iterations):
            new_params = best_params + np.random.uniform(-step_size, step_size, len(p0))
            cost = np.sum((model_wrapper(t_data, *new_params) - V_data) ** 2)

            if cost < best_cost:
                best_params = new_params
                best_cost = cost

        return best_params

    @staticmethod
    def calculate_residuals(V_data, V_sim):
        """
        Calculate residuals
        :param V_data: observed tumor volume (mm3)
        :param V_sim: Simulated tumorvolume (mm3)
        :return:The sum of residuals
        """
        residuals = V_data - V_sim
        rss = np.sum(residuals ** 2)
        return rss

    @staticmethod
    def calculate_aic(n, rss, k):
        """
        Calculate AIC
        :param n: length of V_data
        :param rss: sum of residuals
        :param k: length of the parameters
        :return: Calculation AIC
        """
        return n * np.log(rss / n) + 2 * k

    @staticmethod
    def calculate_bic(n, rss, k):
        """
        Calculate BIC
        :param n: length of V_data
        :param rss: sum of residuals
        :param k: length of the parameters
        :return: Calculation BIC
        """
        return n * np.log(rss / n) + k * np.log(n)


class LogisticModel(TumorGrowthModel):
    """ Class of Logistic Model"""
    def __str__(self):
        """ String representation """
        return (f"dV/dt = c * V * (V_max - V);,\n"
                f"V0=V0")

    def __repr__(self):
        """ Technical representation"""
        return f"LogisticModel(c=?, V_max=?, V0=?)"

    @staticmethod
    def growth(t, V, c, V_max):
        """
        Logistic growth model equation
        :param t: Time (days)
        :param V: Tumor volume (mm3)
        :param c: Growth rate
        :param V_max: Maximum tumor volume (mm3)
        :return: The equation
        """
        return c * V * (V_max - V)

    @staticmethod
    def wrapper(t, *params, V0):
        """
        Wrapper function for model fitting
        :param t: time (days)
        :param V0: Begin condition of the tumor volume (mm3)
        :param c: Growth rate
        :param V_max: Maximum tumor volume (mm3)
        :return: An array of tumorvolume (mm3)
        """
        dt = t[1] - t[0]  # Assuming uniform time steps
        V = [V0]
        for i in range(1, len(t)):
            V_new = V[-1] + dt * LogisticModel.growth(t[i], V[-1], *params)
            V.append(V_new)
        return np.array(V)

    @staticmethod
    def runge(t, V0, c, V_max, dt):
        """
        Runge-Kutta method
        :param t: time (days)
        :param V0: Begin condition of the tumor volume (mm3)
        :param c: Growth rate
        :param V_max: Maximum tumor volume (mm3)
        :param dt: Difference in time
        :return: Runge-Kutta solution of Logistic Model
        """
        return TumorGrowthModel.Runge_method(LogisticModel.growth, t, V0, dt, c, V_max)

# Mendelsohn Model
class MendelsohnModel(TumorGrowthModel):
    """ Class of Mendelsohn equations"""
    def __str__(self):
        """ String representation """
        return (f"dV/dt = c * V^d;\n,"
                f"V0=V0")
    def __repr__(self):
        """ Technical representation"""
        return f"MendelsohnModel(c=?, d=?, V0=?)"

    @staticmethod
    def growth(t, V, c, d):
        """
        Mendelsohn growth model equation
        :param t: Time (days)
        :param V: Tumor volume (mm3)
        :param c: Growth rate
        :param d: Growth speed/curve
        :return: The equation"""
        return c * V * np.exp(d)

    @staticmethod
    def wrapper(t, *params, V0):
        """
        Wrapper function for model fitting
        :param t: time (days)
        :param V0: Begin condition of the tumor volume (mm3)
        :param c: Growth rate
        :param d: Growth speed/curve
        :return: An array of tumorvolume (mm3)
        """
        dt = t[1] - t[0]
        V = [V0]
        for i in range(1, len(t)):
            V_new = V[-1] + dt * MendelsohnModel.growth(t[i], V[-1], *params)
            V.append(V_new)
        return np.array(V)
    @staticmethod
    def runge(t, V0, c, d, dt):
        """
        Runge-Kutta method
        :param t: time (days)
        :param V0: Begin condition of the tumor volume (mm3)
        :param c: Growth rate
        :param d: Growth speed/curve
        :param dt: Difference in time
        :return: Runge-Kutta solution of Mendelsohn Model"""
        return TumorGrowthModel.Runge_method(MendelsohnModel.growth, t, V0, c, d, dt)

# Gompertz model
class GompertzModel(TumorGrowthModel):
    """ Class of Gompertz Model"""
    def __str__(self):
        """ String representation """
        return (f"dV/dt = c * V * ln(V_max/V);\n,"
                f"V0=V0")

    def __repr__(self):
        """ Technical representation"""
        return f"GompertzModel(c=?, V_max=?, V0=?)"

    @staticmethod
    def growth(t, V, c, V_max):
        """
        Gompertz growth model equation
        :param t: Time (days)
        :param V: Tumor volume (mm3)
        :param c: Growth rate
        :param V_max: Maximum tumor volume (mm3)
        :return: The equation"""
        return c * V * np.log(V_max / V)

    @staticmethod
    def wrapper(t, *params, V0):
        """
        Wrapper function for model fitting
        :param t: time (days)
        :param V0: Begin condition of the tumor volume (mm3)
        :param c: Growth rate
        :param V_max: Maximum tumor volume (mm3)
        :return: An array of tumorvolume (mm3)
        """
        dt = t[1] - t[0]
        V = [V0]
        for i in range(1, len(t)):
            V_new = V[-1] + dt * GompertzModel.growth(t[i], V[-1], *params)
            V.append(V_new)
        return np.array(V)

    @staticmethod
    def runge(t, V0, c, V_max, dt):
        """
        Runge-Kutta method
        :param t: time (days)
        :param V0: Begin condition of the tumor volume (mm3)
        :param c: Growth rate
        :param V_max: Maximum tumor volume (mm3)
        :param dt: Difference in time
        :return: Runge-Kutta solution of Gompertz Model"""
        return TumorGrowthModel.Runge_method(GompertzModel.growth, t, V0, c, V_max, dt)

class VonBertalanffyModel(TumorGrowthModel):
    """ Class of Von Bertalanffy Model"""
    def __str__(self):
        """ String representation """
        return (f"dV/dt = c * V^(2/3) - d * V;\n,"
                f"V0=V0")
    def __repr__(self):
        """ Technical representation"""
        return f"VonBertalanffyModel(c=?, d=?, V0=?)"
    @staticmethod
    def growth(t, V, c, d):
        """
        Von Bertalanffy growth model equation
        :param t: Time (days)
        :param V: Tumor volume (mm3)
        :param c: Growth rate
        :param d: Growth speed/curve
        :return: The equation
        """
        return c * V ** (2 / 3) - d * V

    @staticmethod
    def wrapper(t, *params, V0):
        """
        Wrapper function for model fitting
        :param t: time (days)
        :param V0: Begin condition of the tumor volume (mm3)
        :param c: Growth rate
        :param d: Growth speed/curve
        :return: An array of tumorvolume (mm3)
        """
        dt = t[1] - t[0]
        V = [V0]
        for i in range(1, len(t)):
            V_new = V[-1] + dt * VonBertalanffyModel.growth(t[i], V[-1], *params)
            V.append(V_new)
        return np.array(V)
    @staticmethod
    def runge(t, V0, c, d, dt):
        """
        Runge-Kutta method
        :param t: time (days)
        :param V0: Begin condition of the tumor volume (mm3)
        :param c: Growth rate
        :param d: Growth speed/curve
        :param dt: Difference in time
        :return: Runge-Kutta solution of Von Bertalanfyy Model"""
        return TumorGrowthModel.Runge_method(VonBertalanffyModel.growth, t, V0, c, d, dt)

class MontrollModel(TumorGrowthModel):
    """ Class of Montroll Model"""

    def __str__(self):
        """ String representation """
        return (f"dV/dt = c * V * (V_max^d - V^d);\n,"
                f"V0=V0")

    def __repr__(self):
        """ Technical representation"""
        return f"MontrollModel(c=?, V_max=?, d=?, V0=?)"

    @staticmethod
    def growth(t, V, c, V_max, d):
        """"
        Montroll growth model equation
        :param t: Time (days)
        :param V: Tumor volume (mm3)
        :param c: Growth rate
        :param V_max: Maximum tumor volume (mm3)
        :return: The equation
        """
        return c * V * (V_max ** d - V ** d)

    @staticmethod
    def wrapper( t, *params, V0):
        """
        Wrapper function for model fitting
        :param t: time (days)
        :param V0: Begin condition of the tumor volume (mm3)
        :param c: Growth rate
        :param V_max: Maximum tumor volume (mm3)
        :param d: Growth speed/curve
        :return: An array of tumorvolume (mm3)
        """
        dt = t[1] - t[0]
        V = [V0]
        for i in range(1, len(t)):
            V_new = V[-1] + dt * MontrollModel.growth(t[i], V[-1], *params)
            V.append(V_new)
        return np.array(V)
    @staticmethod
    def montroll_Runge(t, V0, c, V_max, d, dt):
        """
        Compute ODE-Solver Runge-Kutta method
        :param t: time (days)
        :param V0: Begin condition of the tumor volume (mm3)
        :param c: Growth rate
        :param V_max: Maximum tumor volume (mm3)
        :param d: Growth speed/curve
        :param dt: Difference in time
        :return: Runge-Kutta solution of Montroll Model
        """
        V = [V0]  # Beginwaarde van het volume
        for i in range(1, len(t)):
            t_current = t[i - 1]
            V_current = V[-1]

            # Bereken de groeisnelheid voor Montroll's model
            y1 = dt * c * V_current * (V_max**d - V_current**d)
            y2 = dt * c * (V_current + y1 / 2) * (V_max**d - (V_current + y1 / 2)**d)
            y3 = dt * c * (V_current + y2 / 2) * (V_max**d - (V_current + y2 / 2)**d)
            y4 = dt * c * (V_current + y3) * (V_max**d - (V_current + y3)**d)

            # Bereken de nieuwe waarde van V
            V_new = V_current + (y1 + 2 * y2 + 2 * y3 + y4) / 6
            V.append(V_new)

        return np.array(V)
    @staticmethod
    def runge(t, V0, c, V_max, d, dt):
        """
        Runge-Kutta method
        :param t: time (days)
        :param V0: Begin condition of the tumor volume (mm3)
        :param c: Growth rate
        :param V_max: Maximum tumor volume (mm3)
        :param d: Growth speed/curve
        :param dt: Difference in time
        :return: Runge-Kutta solution of Montroll Model
        """
        return MontrollModel.montroll_Runge(t, V0, c, V_max, d, dt)


class AlleeEffectModel(TumorGrowthModel):
    """ Class of Allee Effect Model"""

    def __str__(self):
        """ String representation """
        return (f"dV/dt = c * (V - V_min) * (V_max - V);"
                f"\nV0=V0")

    def __repr__(self):
        """ Technical representation"""
        return f"AlleeEffectModel(c=?, V_max=?, V_min=?, V0=?)"

    @staticmethod
    def growth(t, V, c, V_min, V_max):
        """
        Allee Effect growth model equation
        :param t: Time (days)
        :param V: Tumor volume (mm3)
        :param c: Growth rate
        :param V_min: Minimum tumor volume (mm3)
        :param V_max: Maximum tumor volume (mm3)
        :return: The equation
        """
        # Allee effect: growth rate depends on V and the boundaries V_min and V_max
        if V <= V_min or V >= V_max:
            return 0
        return c * (V - V_min) * (V_max - V)

    @staticmethod
    def wrapper(t, *params, V0):
        """
        Wrapper function for model fitting
        :param t: time (days)
        :param V0: Begin condition of the tumor volume (mm3)
        :param c: Growth rate
        :param V_min: Minimum tumor volume (mm3)
        :param V_max: Maximum tumor volume (mm3)
        :return: An array of tumorvolume (mm3)
        """
        dt = t[1] - t[0]
        V = [V0]
        for i in range(1, len(t)):
            V_new = V[-1] + dt * AlleeEffectModel.growth(t[i], V[-1], *params)
            V.append(V_new)
        return np.array(V)

    @staticmethod
    def allee_Runge(t, V0, c, V_min, V_max, dt):
        """
        Compute ODE-Solver Runge-Kutta method
        :param t: time (days)
        :param V0: Begin condition of the tumor volume (mm3)
        :param c: Growth rate
        :param V_min: Minimum tumor volume (mm3)
        :param V_max: Maximum tumor volume (mm3)
        :param dt: Difference in time
        :return: Runge-Kutta solution of Allee Effect Model
        """
        V = [V0]  # Beginwaarde van het volume
        for i in range(1, len(t)):
            t_current = t[i - 1]
            V_current = V[-1]

            # Bereken de groeisnelheid voor het Allee-effect model (direct in de Runge methode)
            y1 = dt * c * (V_current - V_min) * (V_max - V_current)
            y2 = dt * c * (V_current + y1 / 2 - V_min) * (V_max - (V_current + y1 / 2))
            y3 = dt * c * (V_current + y2 / 2 - V_min) * (V_max - (V_current + y2 / 2))
            y4 = dt * c * (V_current + y3 - V_min) * (V_max - (V_current + y3))

            # Bereken de nieuwe waarde van V
            V_new = V_current + (y1 + 2 * y2 + 2 * y3 + y4) / 6
            V.append(V_new)

        return np.array(V)

    @staticmethod
    def runge(t, V0, c, V_min, V_max, dt):
        """
        Runge-Kutta method
        :param t: time (days)
        :param V0: Begin condition of the tumor volume (mm3)
        :param c: Growth rate
        :param V_min: Minimum tumor volume (mm3)
        :param V_max: Maximum tumor volume (mm3)
        :param dt: Difference in time
        :return: Runge-Kutta solution of Allee Effect Model
        """
        return AlleeEffectModel.allee_Runge(t, V0, c, V_min, V_max, dt)


class ModelEvaluator:
    """ Class of Model Evaluator """
    def __init__(self, t_data, V_data):
        self.t_data = t_data
        self.V_data = V_data
        self.models = []

    def add_model(self, model_class, **kwargs):
        """
        :param model_class: The class of the model to be added.
        :param kwargs: arguments passed to model_class.
        """
        self.models.append(model_class(self.t_data, self.V_data))

    def evaluate(self, t_forward=None):
        """
        Evaluate the model
        :param t_forward: Forward time for which model should be evaluated
        :return: A list of each model and their formula, BIC and AIC.
        """
        if t_forward is None:
            t_forward = np.linspace(0, 120, 100)  # Standard time range for prediction

        # Initialize result list
        results = []
        # Dictionary of the initial parameters
        p0_dict = {
            "LogisticModel": [0.01, 8],
            "MendelsohnModel": [0.01, 0.1],
            "GompertzModel": [0.01, 8],
            "AlleeEffectModel": [0.01, 0, 8],
            "VonBertalanffyModel": [0.1, 0.01],
            "MontrollModel": [0.01, 8, 0.1]
        }

        for model in self.models:
            model_name = model.__class__.__name__
            if model_name not in p0_dict:
                print(f"Initial parameters for {model_name} not defined!")
                continue

            # Get initial parameters for this model
            p0 = p0_dict[model_name]
            print(f"Fitting {model_name} with initial parameters {p0}...")

            # Fit the model for the class instances
            params = TumorGrowthModel.fit_model_brute_force(
                model, self.t_data, self.V_data, p0=p0, num_iterations=10000, step_size=0.01)

            # Simulate the model
            V_sim = model.runge(self.t_data, self.V_data[0], *params, self.t_data[1] - self.t_data[0])

            rss = model.calculate_residuals(self.V_data, V_sim)
            aic = model.calculate_aic(len(self.V_data), rss, len(params))
            bic = model.calculate_bic(len(self.V_data), rss, len(params))

            # Tabel of model, formula, AIC and BIC
            results.append({
                'Model': repr(model),
                'Formula': str(model),
                'AIC': aic,
                'BIC': bic})

            # Visualization
            plt.figure(figsize=(10, 6))
            plt.scatter(self.t_data, self.V_data, color='r', label='Data')
            plt.plot(self.t_data, V_sim, label=f"{model.__class__.__name__}", color='b')
            plt.title("Tumorgroei Modellen vs. Data")
            plt.xlabel("Tijd (dagen)")
            plt.ylabel("Tumorvolume (mm³)")
            plt.legend()
            plt.grid(True)
            plt.show()


        return results



if __name__ == "__main__":
    # Stel zelf tijd en volume data in
    tdata = [
        3.46, 4.58, 5.67, 6.64, 7.63, 8.41, 9.32, 10.27, 11.19,
        12.39, 13.42, 15.19, 16.24, 17.23, 18.18, 19.29, 21.23, 21.99,
        24.33, 25.58, 26.43, 27.44, 28.43, 30.49, 31.34, 32.34, 33.00,
        35.20, 36.34, 37.29, 38.50, 39.67, 41.37, 42.58, 45.39, 46.38,
        48.29, 49.24, 50.19, 51.14, 52.10, 54.00, 56.33, 57.33, 59.38,
    ]
    Vdata = [
        0.0158, 0.0264, 0.0326, 0.0445, 0.0646, 0.0933, 0.1454, 0.2183, 0.2842,
        0.4977, 0.6033, 0.8441, 1.2163, 1.4470, 2.3298, 2.5342, 3.0064, 3.4044,
        3.2046, 4.5241, 4.3459, 5.1374, 5.5376, 4.8946, 5.0660, 6.1494, 6.8548,
        5.9668, 6.6945, 6.6395, 6.8971, 7.2966, 7.2268, 6.8815, 8.0993, 7.2112,
        7.0694, 7.4971, 6.9974, 6.7219, 7.0523, 7.1095, 7.0694, 8.0562, 7.2268,
    ]

    t_forward = np.linspace(0.1, 60, 500)
    model = ModelEvaluator(tdata, Vdata)
    model.add_model(LogisticModel)
    model.add_model(GompertzModel)
    model.add_model(AlleeEffectModel)
    model.add_model(MontrollModel)
    model.add_model(MendelsohnModel)
    model.add_model(VonBertalanffyModel)

    result = model.evaluate(t_forward)
    # model.visualize(t_vooruit)

    df = pd.DataFrame(result).set_index('Model').sort_values("AIC")
    df1 = pd.DataFrame(result).set_index('Model').sort_values("BIC")
    print(df)
    print(df1)
