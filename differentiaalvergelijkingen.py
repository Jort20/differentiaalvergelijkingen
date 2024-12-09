import csv

import numpy as np
import matplotlib.pyplot as plt

class tumor_growth:
    """Simulate a tumor growth model"""
    def __init__(self, tijd, V0):
        self.time = tijd
        self.V0 = V0

    def __str__(self):
        pass
    def __repr__(self):
        pass

    # Models
    def gompertz_growth(self, t, V, c, V_max):
        return c * V * np.log(V_max / V)

    def logistic_growth(self, t, V, c, V_max):
        return c * V * (1 - V / V_max)

    def mendelsohn_growth(self, t, V, c, d):
        return c * V **d

    def montroll_growth(self, t, V, c, V_max, d):
        return c * V * ( 1- V**d/V_max**d)

    def allee_effect_growth(self, t, V, c, V_max, V_min):
        return c * (1-V_min/V) * (1-V/V_max)

    def von_bertalanffy_growth(self, t, V, c, d):
        return c * V**(2/3) - d * V

    # Solver
    def ODE_solver(self, model, params, dt=0.1):
        V = [self.V0]
        for i in range(1, len(self.time)):
            V_new = V[-1] + dt * model(self.time[i], V[-1], **params)
            V.append(V_new)
        return V


    # Generate fake data for comparison
    def generate_fake_data(self):
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

            # Convert back to lists if needed
            sorted_t_data = list(sorted_t_data)
            sorted_V_data = list(sorted_V_data)

            print("Sorted t_data:", sorted_t_data)
            print("Sorted V_data:", sorted_V_data)
    


        #t_data = [0, 13, 20, 32, 42, 55, 65, 75, 85, 88, 95, 98, 107, 115, 120]
        #V_data = [250, 255, 550, 575, 576, 800, 1050, 1250, 1750, 2000, 2550, 2750, 3000, 3500, 4000]
        return sorted_t_data, sorted_V_data


if __name__ == "__main__":


    V0 = 36
    params = {"c": 0.1, "V_max": 2353}
    params_montroll = {"c": 0.1, "V_max": 2353, "d": 1}
    params_allee = {"c": 0.1, "V_max": 2353, "V_min": 0}
    params_mendelsohn = {"c": 0.1, "d": 1}

    tijd = np.linspace(0, 38, 100)
    tg = tumor_growth(tijd, V0)
    t_data, V_data = tg.generate_fake_data()
    V_sim_gompertz = tg.ODE_solver(tg.gompertz_growth, params)

    V_sim_logistic = tg.ODE_solver(tg.logistic_growth, params)
    V_sim_mendelsohn = tg.ODE_solver(tg.mendelsohn_growth, params_mendelsohn)
    V_sim_montroll = tg.ODE_solver(tg.montroll_growth, params_montroll)
    V_sim_allee = tg.ODE_solver(tg.allee_effect_growth, params_allee)
    V_sim_von_bertalanffy = tg.ODE_solver(tg.von_bertalanffy_growth, params_mendelsohn)
    plt.figure(figsize=(10, 6))

    plt.plot(tijd, V_sim_gompertz, label="Gompertz Model", color="blue")
    plt.plot(tijd, V_sim_logistic, label="Logistisch Model", color="green")
    plt.plot(tijd, V_sim_montroll, label="Morello Model", color="yellow")
    plt.plot(tijd, V_sim_allee, label="Allee Effect Model", color="red")
    plt.plot(tijd, V_sim_mendelsohn, label="Mendelsohn Model", color="orange")
    plt.plot(tijd, V_sim_von_bertalanffy, label="Von Bertalanffy Model", color="purple")


    plt.scatter(t_data, V_data, color="red", label="Nepdata (observaties)")

    plt.title("Tumorgroei Modellen vs. Nepdata")
    plt.xlabel("Tijd (dagen)")
    plt.ylabel("Tumorvolume (mm³)")
    plt.legend()
    plt.grid(True)
    plt.show()