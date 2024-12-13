import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Basis klasse
class BaseModel:
    """
    Basis klasse voor groeimodellen. 
    Geeft methoden voor numerieke integratie, residuals en modelselectie metrics (AIC en BIC).
    """
    
    def __init__(self, t_data, V_data):
        """
        Initialisatie van de basis klasse met data.

        Parameters:
        t_data (np.array): Array met tijddata.
        V_data (np.array): Array met volumedata.
        """

        self.t_data = t_data
        self.V_data = V_data

    
    def runga_method(self, model_growth, t, V0, dt, *params):
        """
        Numerieke integratie van een modelgroei-functie met behulp van de Runge-Kutta 4e orde methode.

        Parameters:
        model_growth (function): Groeifunctie die de groeirate berekent.
        t (np.array): Array met tijddata.
        V0 (float): Startwaarde van het volume.
        dt (float): Tijdstapgrootte.
        *params: Overige parameters van het modelgroei.

        Returns:
        np.array: Array met gesimuleerd volume over de tijd.
        """
        V = [V0]  # Startwaarde
        for i in range(1, len(t)):
            t_current = t[i - 1]
            V_current = V[-1]
            # Berekeningen volgens Runge-Kutta 4e orde
            y1 = dt * model_growth(t_current, V_current, *params)
            y2 = dt * model_growth(t_current + dt / 2, V_current + y1 / 2, *params)
            y3 = dt * model_growth(t_current + dt / 2, V_current + y2 / 2, *params)
            y4 = dt * model_growth(t_current + dt, V_current + y3, *params)
            V_new = V_current + (y1 + 2 * y2 + 2 * y3 + y4) / 6
            V.append(V_new)
        return np.array(V)


    @staticmethod
    def calculate_residuals(V_data, V_sim):
        """
        Berekent de residuals tussen de gemeten data en gesimuleerde data.

        Parameters:
        V_data (np.array): Array met gemeten volumedata.
        V_sim (np.array): Array met gesimuleerde volumedata.

        Returns:
        float: Residual Sum of Squares (RSS).
        """
        residuals = V_data - V_sim
        rss = np.sum(residuals**2)  # Residual Sum of Squares
        return rss


    @staticmethod
    def calculate_aic(n, rss, k):
        """
        Berekent de Akaike Information Criterion (AIC) voor modelselectie.

        Parameters:
        n (int): Aantal datapunten.
        rss (float): Residual Sum of Squares.
        k (int): Aantal parameters in het model.

        Returns:
        float: AIC waarde.
        """
        return n * np.log(rss / n) + 2 * k

    @staticmethod
    def calculate_bic(n, rss, k):
        """
        Berekent de Bayesian Information Criterion (BIC) voor modelselectie.

        Parameters:
        n (int): Aantal datapunten.
        rss (float): Residual Sum of Squares.
        k (int): Aantal parameters in het model.

        Returns:
        float: BIC waarde.
        """
        return n * np.log(rss / n) + k * np.log(n)

# Specifieke groeimodellen
class LogisticModel(BaseModel):
    """
    Logistische groeimodel.
    """
    @staticmethod
    def growth(t, V, c, V_max):
        """
        Logistische groeifunctie.

        Parameters:
        t (float): Tijd.
        V (float): Volume op tijd t.
        c (float): Groeiparameter.
        V_max (float): Maximum volume.

        Returns:
        float: Groei snelheid.
        """
        return c * V * (V_max - V)

    # Simulatie met Runge-Kutta
    def simulate(self, t, V0, c, V_max, dt):
        """
        Simulatie van het logistische groeimodel met Runge-Kutta methode.

        Parameters:
        t (np.array): Array met tijddata voor de simulatie.
        V0 (float): Startwaarde van het volume.
        c (float): Groeiparameter.
        V_max (float): Maximum volume.
        dt (float): Tijdstapgrootte.

        Returns:
        np.array: Array met gesimuleerd volume over de tijd.
        """
        return self.runga_method(self.growth, t, V0, dt, c, V_max)

class GompertzModel(BaseModel):
    """
    Gompertz groeimodel.
    """

    @staticmethod
    def growth(t, V, c, V_max):
        """
        Gompertz groeifunctie.

        Parameters:
        t (float): Tijd.
        V (float): Volume op tijd t.
        c (float): Groeiparameter.
        V_max (float): Maximum volume.

        Returns:
        float: Groei snelheid.
        """
        return c * V * np.log(V_max / V)

    # Simulatie met Runge-Kutta
    def simulate(self, t, V0, c, V_max, dt):
        """
        Simulatie van het Gompertz groeimodel met Runge-Kutta methode.

        Parameters:
        t (np.array): Array met tijddata voor de simulatie.
        V0 (float): Startwaarde van het volume.
        c (float): Groeiparameter.
        V_max (float): Maximum volume.
        dt (float): Tijdstapgrootte.

        Returns:
        np.array: Array met gesimuleerd volume over de tijd.
        """
        return self.runga_method(self.growth, t, V0, dt, c, V_max)

class VonBertalanffyModel(BaseModel):
    """
    Von Bertalanffy groeimodel.
    """
    @staticmethod
    def growth(t, V, c, d):
        """
        Von Bertalanffy groeifunctie.

        Parameters:
        t (float): Tijd.
        V (float): Volume op tijd t.
        c (float): Groeiparameter.
        d (float): Groeiparameter.

        Returns:
        float: Groei snelheid.
        """
        return c * V**(2/3) - d * V

    # Simulatie met Runge-Kutta
    def simulate(self, t, V0, c, d, dt):
        """
        Simulatie van het Von Bertalanffy groeimodel met Runge-Kutta methode.

        Parameters:
        t (np.array): Array met tijddata voor de simulatie.
        V0 (float): Startwaarde van het volume.
        c (float): Groeiparameter.
        d (float): Groeiparameter.
        dt (float): Tijdstapgrootte.

        Returns:
        np.array: Array met gesimuleerd volume over de tijd.
        """
        return self.runga_method(self.growth, t, V0, dt, c, d)

class MendelsohnModel(BaseModel):
    """
    Mendelsohn groeimodel.
    """
    @staticmethod
    def growth(t, V, c, D):
        """
        Mendelsohn groeifunctie.

        Parameters:
        t (float): Tijd.
        V (float): Volume op tijd t.
        c (float): Groeiparameter.
        D (float): Groeiparameter.

        Returns:
        float: Groei snelheid.
        """
        return c * V**D

    # Simulatie met Runge-Kutta
    def simulate(self, t, V0, c, D, dt):
        """
        Simulatie van het Mendelsohn groeimodel met Runge-Kutta methode.

        Parameters:
        t (np.array): Array met tijddata voor de simulatie.
        V0 (float): Startwaarde van het volume.
        c (float): Groeiparameter.
        D (float): Groeiparameter.
        dt (float): Tijdstapgrootte.

        Returns:
        np.array: Array met gesimuleerd volume over de tijd.
        """
        return self.runga_method(self.growth, t, V0, dt, c, D)

class MontrollModel(BaseModel):
    """
    Montroll groeimodel.
    """
    @staticmethod
    def growth(t, V, c, V_max, d):
        """
        Montroll groeifunctie.

        Parameters:
        t (float): Tijd.
        V (float): Volume op tijd t.
        c (float): Groeiparameter.
        V_max (float): Maximum volume.
        d (float): Groeiparameter.

        Returns:
        float: Groei snelheid.
        """
        return c * V * (V_max**d - V**d)


    def simulate(self, t, V0, c, V_max, d, dt):
        """
        Simulatie van het Montroll groeimodel met Runge-Kutta methode.

        Parameters:
        t (np.array): Array met tijddata voor de simulatie.
        V0 (float): Startwaarde van het volume.
        c (float): Groeiparameter.
        V_max (float): Maximum volume.
        d (float): Groeiparameter.
        dt (float): Tijdstapgrootte.

        Returns:
        np.array: Array met gesimuleerd volume over de tijd.
        """
        return self.runga_method(self.growth, t, V0, dt, c, V_max, d)

class AlleeModel(BaseModel):
    """
    Allee groeifunctie met minimum en maximum drempelwaarden.
    """
    @staticmethod
    def growth(t, V, c, V_min, V_max):
        """
        Allee groeifunctie.

        Parameters:
        t (float): Tijd.
        V (float): Volume op tijd t.
        c (float): Groeiparameter.
        V_min (float): Minimum drempelwaarde van volume.
        V_max (float): Maximum drempelwaarde van volume.

        Returns:
        float: Groei snelheid.
        """
        if V <= V_min or V >= V_max:
            return 0
        return c * (V - V_min) * (V_max - V)

    # Simulatie met Runge-Kutta
    def simulate(self, t, V0, c, V_min, V_max, dt):
        """
        Simulatie van het Allee groeimodel met Runge-Kutta methode.

        Parameters:
        t (np.array): Array met tijddata voor de simulatie.
        V0 (float): Startwaarde van het volume.
        c (float): Groeiparameter.
        V_min (float): Minimum drempelwaarde van volume.
        V_max (float): Maximum drempelwaarde van volume.
        dt (float): Tijdstapgrootte.

        Returns:
        np.array: Array met gesimuleerd volume over de tijd.
        """
        return self.runga_method(self.growth, t, V0, dt, c, V_min, V_max)

# Data generatie
class DataHandler:
    @staticmethod
    def generate_fake_data():
        """
        Generatie van gesimuleerde tijd en volume data.

        Returns:
        np.array: Gesimuleerde tijddata.
        np.array: Gesimuleerde volumedata.
        """
        t_data = [
            3.46,  4.58,  5.67,  6.64,  7.63,  8.41,  9.32, 10.27, 11.19,
            12.39, 13.42, 15.19, 16.24, 17.23, 18.18, 19.29, 21.23, 21.99,
            24.33, 25.58, 26.43, 27.44, 28.43, 30.49, 31.34, 32.34, 33.00,
            35.20, 36.34, 37.29, 38.50, 39.67, 41.37, 42.58, 45.39, 46.38,
            48.29, 49.24, 50.19, 51.14, 52.10, 54.00, 56.33, 57.33, 59.38,
        ]
        V_data = [
            0.0158, 0.0264, 0.0326, 0.0445, 0.0646, 0.0933, 0.1454, 0.2183, 0.2842,
            0.4977, 0.6033, 0.8441, 1.2163, 1.4470, 2.3298, 2.5342, 3.0064, 3.4044,
            3.2046, 4.5241, 4.3459, 5.1374, 5.5376, 4.8946, 5.0660, 6.1494, 6.8548,
            5.9668, 6.6945, 6.6395, 6.8971, 7.2966, 7.2268, 6.8815, 8.0993, 7.2112,
            7.0694, 7.4971, 6.9974, 6.7219, 7.0523, 7.1095, 7.0694, 8.0562, 7.2268, 
        ]
        return np.array(t_data), np.array(V_data)

# Model evaluatie
class Evaluator:
    def __init__(self, t_data, V_data):
        self.t_data = t_data
        self.V_data = V_data

    # Fitting en evalueren van een model
    def fit_and_evaluate(self, model_class, p0, t_forecast, num_iterations=10000, step_size=0.01):
        def model(t, *params):
            """
            Fit een model aan de data en evalueer met AIC en BIC.

            Parameters:
            model_class (class): Groeimodelklasse.
            t_data (np.array): Array met tijddata.
            V_data (np.array): Array met volumedata.
            param_guesses (tuple): Initiale schattingen van de modelparameters.
            dt (float): Tijdstapgrootte voor simulatie.

            Returns:
            tuple: Simulatie resultaten, RSS, AIC en BIC.
            """
            instance = model_class(self.t_data, self.V_data)
            return instance.simulate(t, self.V_data[0], *params, t_forecast[1] - t_forecast[0])

        best_params = p0
        best_cost = np.sum((model(self.t_data, *best_params) - self.V_data) ** 2)

        # Optimalisatie van parameters door iteratief te zoeken naar lagere afwijkingen
        for _ in range(num_iterations):
            new_params = best_params + np.random.uniform(-step_size, step_size, len(p0))# Kleine willekeurige stap in parameterruimte
            cost = np.sum((model(self.t_data, *new_params) - self.V_data) ** 2) # Bereken de afwijking voor de nieuwe parameters
            if cost < best_cost:
                best_params = new_params
                best_cost = cost
        # Simuleer met de geoptimaliseerde parameters
        instance = model_class(self.t_data, self.V_data)
        V_sim = instance.simulate(t_forecast, self.V_data[0], *best_params, t_forecast[1] - t_forecast[0])
        return best_params, V_sim


    def compare_models(self, models, t_forecast):
        """
        Vergelijk verschillende groeimodellen op basis van AIC en BIC.

        Parameters:
        models (list): Lijst van modelklassen.
        t_data (np.array): Array met tijddata.
        V_data (np.array): Array met volumedata.
        param_guesses_list (list of tuples): Initiale schattingen van modelparameters voor elk model.
        dt (float): Tijdstapgrootte voor simulatie.

        Returns:
        DataFrame: Vergelijking van modellen met respect tot AIC en BIC.
        """
        results = {}
        for model_name, (model_class, p0) in models.items():
            params, V_sim = self.fit_and_evaluate(model_class, p0, t_forecast)
            rss = BaseModel.calculate_residuals(self.V_data, V_sim[:len(self.t_data)])
            aic = BaseModel.calculate_aic(len(self.V_data), rss, len(params))
            bic = BaseModel.calculate_bic(len(self.V_data), rss, len(params))
            results[model_name] = {
                'params': params,
                'rss': rss,
                'aic': aic,
                'bic': bic,
                'V_sim': V_sim
            }
        return results

# Voorbeeldgebruik
# Genereer gesimuleerde data
t_data, V_data = DataHandler.generate_fake_data()
evaluator = Evaluator(t_data, V_data)
t_forecast = np.linspace(0, 120, 100)  # Tijdreeks voor voorspelling

# Definieer modellen met startparameters
models = {
    'Logistic': (LogisticModel, [0.01, 7]),
    'Gompertz': (GompertzModel, [0.11, 7.7]),
    'Von Bertalanffy': (VonBertalanffyModel, [0.5, 0.2]),
    'Mendelsohn': (MendelsohnModel, [0.01, 0.1]),
    'Montroll': (MontrollModel, [0.01, 8, 0.1]),
    'Allee': (AlleeModel, [0.05, 0, 7.5])
}

# Vergelijk modellen
results = evaluator.compare_models(models, t_forecast)

# Visualisatie van resultaten
plt.figure(figsize=(10, 6))
plt.scatter(t_data, V_data, color="red", label="Data")

for model_name, result in results.items():
    plt.plot(t_forecast, result['V_sim'], label=f"{model_name} (AIC: {result['aic']:.2f})")

plt.title("Tumorgroei Modellen vs. Data")
plt.xlabel("Tijd (dagen)")
plt.ylabel("Tumorvolume (mm³)")
plt.legend()
plt.grid(True)
plt.show()

# Toon resultaten als tabel
df = pd.DataFrame.from_dict(results, orient='index')
df = df[['aic', 'bic']].sort_values(by='aic')
print(df)
