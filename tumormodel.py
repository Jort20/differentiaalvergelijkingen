import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

class GrowthModel:
    def __init__(self, t_data, V_data):
        self.t_data = t_data
        self.V_data = V_data


    def __str__(self):
        pass

    def __repr__(self):
        pass

    def fit_model_brute_force(self, model_wrapper, t_data, V_data, p0, num_iterations=10000, step_size=0.01):
        def model(t, *params):
            return model_wrapper(t, *params, V_data)  # Include V_data as an argument

        best_params = p0
        best_cost = np.sum((model(t_data, *best_params) - V_data) ** 2)

        for _ in range(num_iterations):
            new_params = best_params + np.random.uniform(-step_size, step_size, len(p0))
            cost = np.sum((model(t_data, *new_params) - V_data) ** 2)

            if cost < best_cost:
                best_params = new_params
                best_cost = cost

        return best_params
    def simulate(self):
        # raise value error of implementatie
        pass

    @staticmethod
    def calculate_aic(n, rss, k):
        return n * np.log(rss / n) + 2 * k

    @staticmethod
    def calculate_bic(n, rss, k):
        return n * np.log(rss / n) + k * np.log(n)

    @staticmethod
    def calculate_residuals(V_data, V_sim):
        return np.sum((V_data - V_sim) ** 2)

# Logistic Model
class Logistic_model(GrowthModel):
    @staticmethod
    def growth(t, V, c, V_max):
        return c * V * (V_max - V)

    @staticmethod
    def logistic_wrapper(self, t, c, V_max, V_data):
        V0 = V_data[0]
        dt = t[1] - t[0]
        V = [V0]
        for i in range(1, len(t)):
            V_new = V[-1] + dt * self.growth(t[i], V[-1], c, V_max)
            V.append(V_new)
        return np.array(V)
class

class ModelEvaluator:
    def __init__(self, t_data, V_data):
        self.t_data = t_data
        self.V_data = V_data
        self.models = []

    def add_model(self, model_class, **kwargs):
        self.models.append(model_class(self.t_data, self.V_data))

    def evaluate(self, t_predict):
        results = []

        for model in self.models:
            params = model.fit(p0=[0.01] * 2, num_iterations=10000)

            V_sim = model.simulate(t_predict, *params)
            rss = model.calculate_residuals(self.V_data, V_sim)
            aic = model.calculate_aic(len(self.V_data), rss, len(params))
            bic = model.calculate_bic(len(self.V_data), rss, len(params))

            results.append({
                'model': model.__class__.__name__,
                'params': params,
                'rss': rss,
                'aic': aic,
                'bic': bic
            })

        return pd.DataFrame(results)

    def visualize(self, t_predict):
        plt.figure(figsize=(10, 6))
        plt.scatter(self.t_data, self.V_data, color="red", label="Data")





