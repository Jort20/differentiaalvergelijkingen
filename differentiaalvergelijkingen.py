import csv

import numpy as np
import matplotlib.pyplot as plt

class tumor_growth:
    """Simulate a tumor growth model"""
    def __init__(self, tijd, V0, t_data, v_data):
        self.t_forward = tijd
        self.V0 = V0
        self.t_data = t_data
        self.v_data = v_data

    def __str__(self):
        return (f"Tumor growth model: "
                f"time : {self.time},\n"
                f"V0 : {self.V0}")
    def __repr__(self):
        pass


    # Generate fake data for comparison
    def generate_fake_data(self):
        """
        Generate a fake data set
        :return t_data; time  and V_data; observations
        """
        t_data = [0, 13, 20, 32, 42, 55, 65, 75, 85, 88, 95, 98, 107, 115, 120]
        V_data = [250, 255, 550, 575, 576, 800, 1050, 1250, 1750, 2000, 2550, 2750, 3000, 3500, 4000]
        return t_data, V_data

    # Tumor growth models
    def gompertz_growth(self, t, V, c, V_max):
        return c * V * np.log(V_max / V)

    def logistic_growth(self, t, V, c, V_max):
        return c * V * (1 - V / V_max)

    def mendelsohn_growth(self, t, V, c, d):
        return c * V **d

    def montroll_growth(self, t, V, c, V_max, d):
        return c * V * ( 1- V**d/V_max**d)

    def allee_effect_growth(self, t, V, c, V_max, V_min):
        return c * (V-V_min) * (1-V/V_max)

    def von_bertalanffy_growth(self, t, V, c, d):
        return c * V**(2/3) - d * V


    def Hean_method(self, model_growth, t, V0, c, V_max, dt):
        V = [V0]
        for i in range(1, len(t)):
            t_current = t[i - 1]
            V_current = V[-1]

            y1 = model_growth(t_current, V_current, c, V_max)

            V_euler = V_current + dt * y1
            y2 = model_growth(t_current + dt, V_euler, c, V_max)

            V_new = V_current + dt * 0.5 * (y1 + y2)
            V.append(V_new)

        return np.array(V)

    def Runga_method(self, model_growth, t, V0, c, V_max, dt):

        V = [V0]
        for i in range(1, len(t)):
            t_current = t[i - 1]
            V_current = V[-1]

            y1 = dt * model_growth(t_current, V_current, c, V_max)
            y2 = dt * model_growth(t_current + dt / 2, V_current + y1 / 2, c, V_max)
            y3 = dt * model_growth(t_current + dt / 2, V_current + y2 / 2, c, V_max)
            y4 = dt * model_growth(t_current + dt, V_current + y3, c, V_max)

            V_new = V_current + (y1 + 2 * y2 + 2 * y3 + y4) / 6
            V.append(V_new)

        return np.array(V)

    def Runga_method2(self, model_growth, t, V0, c, V_max, d, dt):
        V = [V0]
        for i in range(1, len(t)):
            t_current = t[i - 1]
            V_current = V[-1]
            y1 = dt * model_growth(t_current, V_current, c, V_max, d)
            y2 = dt * model_growth(t_current + dt / 2, V_current + y1 / 2 , c, V_max, d)
            y3 = dt * model_growth(t_current + dt / 2, V_current + y2 / 2 , c, V_max, d)
            y4 = dt * model_growth(t_current + dt, V_current + y3, c, V_max, d)

            V_new = V_current + (y1 + 2 * y2 + 2 * y3 + y4) / 6
            V.append(V_new)
        return np.array(V)

    def mendelsohn_Rangu(self, t, V0, c, d, dt):
        return self.Runga_method(self.mendelsohn_growth, t, V0, c, d, dt)
    def montroll_Rangu(self, t, V0, c, V_max, d, dt):
        return self.Runga_method2(self.montroll_growth, t, V0, c, V_max, d, dt)
    def allee_effect_Rangu(self, t, V0, c, V_max, V_min, dt):
        return self.Runga_method2(self.allee_effect_growth, t, V0, c, V_max, V_min, dt)

    def montroll_wrapper(self, t, c, V_max, d):
        dt = t[1] - t[0]
        V = [self.V0]
        for i in range(1, len(t)):
            if self.V0 > 0:
                V_new = V[-1] + dt * self.montroll_growth([i], V[-1], c, V_max, d)
                V.append(V_new)
            if self.V0 < 0:
                v_new = 0
                V.append(v_new)
        return np.array(V)

    def allee_effect_wrapper(self, t, c, V_max, V_min):
        dt = t[1] - t[0]
        V = [self.V0]
        for i in range(1, len(t)):
            if self.V0 > 0:
                V_new = V[-1] + dt * self.allee_effect_growth([i], V[-1], c, V_max, V_min)
                V.append(V_new)
            if self.V0 < 0:
                v_new = 0
                V.append(v_new)
        return np.array(V)

    def mendelsohn_wrapper(self, t, c, d):
        """Wrapper om mendelsohn-growth werkend te maken met curve_fit."""
        dt = t[1] - t[0]
        V = [self.V0]
        for i in range(1, len(t)):
            if self.V0 > 0 :
                V_new = V[-1] + dt * self.mendelsohn_growth(t[i], V[-1], c, d)
                V.append(V_new)
            if self.V0 < 0:
                v_new = 0
                V.append(v_new)
        return np.array(V)

    @staticmethod
    def fit_model_brute_force(model_wrapper, t_data, V_data, p0, num_iterations=1000, step_size=0.01):
        def model(t, *params):
            return model_wrapper(t, *params)

        best_params = p0
        best_cost = np.sum((model(t_data, *best_params) - V_data) ** 2)

        for _ in range(num_iterations):
            new_params = best_params + np.random.uniform(-step_size, step_size, len(p0))
            cost = np.sum((model(t_data, *new_params) - V_data) ** 2)

            if cost < best_cost:
                best_params = new_params
                best_cost = cost

        return best_params

    def evaluate_models(self):
        if self.t_forward is None:
            self.t_forward = np.linspace(0,120,100)
        # Mendelsohn
        initial_params_mendelsohn = [0.1, 1]
        params_mendelsohn = self.fit_model_brute_force(self.mendelsohn_wrapper, self.t_data, self.v_data, p0=initial_params_mendelsohn)
        dt = self.t_forward[1] - self.t_forward[0]
        V_sim_mendelsohn = self.mendelsohn_Rangu(self.t_forward, self.V0, *params_mendelsohn, dt)
        # Montroll
        initial_params_montroll = [0.1, 4000, 1]
        params_montroll = self.fit_model_brute_force(self.montroll_wrapper, self.t_data, self.v_data, p0=initial_params_montroll)
        V_sim_montroll = self.montroll_Rangu(self.t_forward, self.V0, *params_montroll, dt=dt)
        # Allee
        initial_params_allee = [0.1, 4000, 0]
        params_allee = self.fit_model_brute_force(self.allee_effect_wrapper, self.t_data, self.v_data, p0=initial_params_allee)
        V_sim_allee = self.allee_effect_Rangu(self.t_forward, self.V0, *params_allee, dt=dt)
        plt.figure(figsize=(10,5))
        plt.scatter(self.t_data, self.v_data, color="red", label="Nepdata (observaties)")
        plt.plot(t_forward, V_sim_mendelsohn, label=f"Mendelsohn Model, params={params_mendelsohn}",
                 color="purple")
        plt.plot(t_forward, V_sim_montroll, label=f"Montroll Model, params={params_montroll}",
                 color="blue")
        plt.plot(t_forward, V_sim_allee, label=f"Allee Effect Model, params={params_allee}", color="green")
        plt.title("Tumorgroei Modellen vs. Nepdata")
        plt.xlabel("Tijd (dagen)")
        plt.ylabel("Tumorvolume (mm³)")
        plt.legend()
        plt.grid(True)
        plt.show()

    def read_data(self):
        """
        Read text file data
        :return: sorted_t_data, sorted_v_data
        """
        with open("tumor_growth.txt") as file:
            text = csv.reader(file, delimiter='\t')
            header = next(text)
            t_data = []
            V_data = []
            for line in text:
                t_data.append(float(line[1]))
                V_data.append(float(line[2]))
            t_data = [round(num) for num in t_data]
            V_data = [round(num) for num in V_data]
            sorted_t_data, sorted_V_data = zip(*sorted(zip(t_data, V_data)))

            sorted_t_data = list(sorted_t_data)
            sorted_V_data = list(sorted_V_data)
        return sorted_t_data, sorted_V_data


if __name__ == "__main__":
    """ Main """
    V0 = 250
    t_data = [0, 13, 20, 32, 42, 55, 65, 75, 85, 88, 95, 98, 107, 115, 120]
    V_data = [250, 255, 550, 575, 576, 800, 1050, 1250, 1750, 2000, 2550, 2750, 3000, 3500, 4000]
    params = {"c": 0.1, "V_max": 4000}
    params_montroll = {"c": 0.1, "V_max": 4000, "d": 1}
    params_allee = {"c": 0.1, "V_max": 4000, "V_min": 0}
    initial_params_mendelsohn = {"c": 0.1, "d": 1}

    t_forward = np.linspace(0, 120, 100)
    tg = tumor_growth(t_forward, V0, t_data, V_data)
    tg.evaluate_models()

    #V_sim_gompertz = tg.ODE_solver(tg.gompertz_growth, params)

    #V_sim_logistic = tg.ODE_solver(tg.logistic_growth, params)
    #V_sim_mendelsohn = tg.ODE_solver(tg.mendelsohn_growth, params_mendelsohn)
    #V_sim_montroll = tg.ODE_solver(tg.montroll_growth, params_montroll)
    #V_sim_allee = tg.ODE_solver(tg.allee_effect_growth, params_allee)
    #V_sim_von_bertalanffy = tg.ODE_solver(tg.von_bertalanffy_growth, params_mendelsohn)

    #plt.figure(figsize=(10, 6))

    #plt.plot(tijd, V_sim_gompertz, label="Gompertz Model", color="blue")
    #plt.plot(tijd, V_sim_logistic, label="Logistisch Model", color="green")
    #plt.plot(tijd, V_sim_montroll, label="Morello Model", color="yellow")
    #plt.plot(tijd, V_sim_allee, label="Allee Effect Model", color="red")
    #plt.plot(tijd, V_sim_mendelsohn, label="Mendelsohn Model", color="orange")
    #plt.plot(tijd, V_sim_von_bertalanffy, label="Von Bertalanffy Model", color="purple")


    #plt.scatter(t_data, V_data, color="red", label="Nepdata (observaties)")

    #plt.title("Tumorgroei Modellen vs. Nepdata")
    #plt.xlabel("Tijd (dagen)")
    #plt.ylabel("Tumorvolume (mm³)")
    #plt.legend()
    #plt.grid(True)
    #plt.show()