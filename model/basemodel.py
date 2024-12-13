#imports
from math import exp
from random import random
import numpy as np

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
    def __str__(self):
        """ String representatie """
        return f"dV/dt = c * V * (V_max - V)"

    def __repr__(self):
        """ Technische representatie """
        return f"LogisticModel(BaseModel)"

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
    def __str__(self):
        """ String representatie """
        return f"dV/dt = c * V * ln(V_max/V)"

    def __repr__(self):
        """ Technische representatie"""
        return f"GompertzModel(BaseModel)"
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
    def __str__(self):
        """ String representatie """
        return f"dV/dt = c * V^(2/3) - d * V"
    def __repr__(self):
        """ Technische representatie"""
        return f"VonBertalanffyModel(BaseModel)"
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
    def __str__(self):
        """ String representatie """
        return f"dV/dt = c * V^d"

    def __repr__(self):
        """ Technische representatie"""
        return f"MendelsohnModel(BaseModel)"

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
    def __str__(self):
        """ String representation """
        return f"dV/dt = c * V * (V_max^d - V^d)"

    def __repr__(self):
        """ Technical representation"""
        return f"MontrollModel(BaseModel)"

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
    def __str__(self):
        """ String representatie """
        return f"dV/dt = c * (V - V_min) * (V_max - V)"

    def __repr__(self):
        """ Technische representatie"""
        return f"AlleeEffectModel(BaseModel)"

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
    """
    Class DataHandler
    Genereert neppe data set, als het nodig is.
    """
    def generate_data(self, length=45):
        """
        Genereert gesimuleerde tijd en volume nepdataset.

        Parameters: Give your own length of the range of data
        Returns:
        np.array: Gesimuleerde tijddata.
        np.array: Gesimuleerde volumedata.
        """
        t_data = [random() for _ in range(length)]
        V_data = [exp(data_t) - 1.0 for data_t in t_data]
        return np.array(t_data), np.array(V_data)

# Model evaluatie
class Evaluator:
    """
    Klasse Evaluator
    Geeft methoden om het model te fitten en te evalueren.
    """
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

        best_params = p0 # initiele parameters
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
            model_instance = model_class(self.t_data, self.V_data)
            results[model_name] = {
                'Formule': str(model_instance),
                'Model': repr(model_instance),
                'params': params,
                'rss': rss,
                'AIC': aic,
                'BIC': bic,
                'V_sim': V_sim,

            }
        return results





