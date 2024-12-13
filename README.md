# Tumorgroei Modellen met Python

Dit project implementeert en vergelijkt verschillende wiskundige modellen voor tumorgroei, inclusief het Gompertz-, Logistic-, Mendelsohn-, Montroll-, Von Bertalanffy- en Allee-effect-model. De implementatie omvat methoden om deze modellen te simuleren, te fitten aan data en de resultaten te visualiseren. Daarnaast worden numerieke methoden zoals de Heun- en Runge-Kutta-methode gebruikt om de modellen op te lossen.

---

## Inhoudsopgave

1. [Achtergrond](#achtergrond)
2. [Modellen en Parameters](#modellen-en-parameters)
3. [Numerieke Methoden](#numerieke-methoden)
4. [Functies en Functionaliteit](#functies-en-functionaliteit)
5. [Installatie](#installatie)
6. [Gebruik](#gebruik)
7. [Resultaten](#resultaten)
8. [Problemen](#problemen)
9. [Authors and Acknowledgments](#Author_and_Acknowledgments)
10. [License](#license)
11. [Changelog](#changelog)
12. [Referenties](#referenties)

---

## Achtergrond

De groei van tumoren is een complex proces dat vaak wordt gemodelleerd met behulp van differentiaalvergelijkingen. Deze modellen simuleren hoe een tumor in volume groeit over de tijd en helpen onderzoekers om groeiprofielen te begrijpen en behandelingsstrategieën te ontwikkelen. 

De zes modellen in dit project hebben elk unieke eigenschappen die ze geschikt maken voor verschillende toepassingen:
- **Gompertz-model**: Wordt vaak gebruikt in tumorbiologie vanwege zijn flexibiliteit en realistische groei-limieten.
- **Logistisch model**: Eenvoudig en effectief, vaak gebruikt voor populatiedynamica of groei met beperkende factoren.
- **Von Bertalanffy-model**: Beschrijft zowel groei als afbraak, word gebruikt in biologie en fysiologie.
- **Mendelsohn-model**: Beschrijft een ongelimiteerd groei van tumor, vaak gebruikt voor populatie groei of groei van tumors zonder limieten, biologisch onrealistisch. 
- **Montroll-model**: Het model geeft een beter fit en predictie van het buigpunt, vaak gebruikt voor populatie groei en tumor groei.
- **Allee-effect-model**: Het model neemt het negative en positive groei mee, vaak gebruikt voor ecologisch en biologische doeleinden.
---

## Modellen en Parameters

Hieronder worden de modellen en hun parameters in detail beschreven:

### Gompertz-model

De differentiaalvergelijking voor het Gompertz-model is:
dV/dt = c * V * ln(V_max / V)

**Parameters**:
- `V`: Het tumorvolume op een bepaald tijdstip `t` (mm³).
- `c`: De groeisnelheidsparameter (dag⁻¹). Een hogere waarde geeft een snellere initiële groei aan.
- `V_max`: Het maximale tumorvolume (mm³), dat de asymptotische limiet voor de groei vertegenwoordigt.
- `t`: Tijd (dagen).

### Logistisch model

De logistische groeivergelijking is:
dV/dt = c * V * (1 - V / V_max)

**Parameters**:
- `V`: Het tumorvolume (mm³).
- `c`: De groeisnelheid (dag⁻¹). Deze parameter bepaalt hoe snel de tumor zijn maximale capaciteit nadert.
- `V_max`: Het maximumvolume dat de tumor kan bereiken vanwege fysieke of biologische beperkingen.
- `t`: Tijd (dagen).


### Von Bertalanffy-model

Het Von Bertalanffy-model beschrijft de balans tussen groei en afbraak:
dV/dt = c * V^(2/3) - d * V

**Parameters**:
- `V`: Het tumorvolume (mm³).
- `c`: De groeifactor (dag⁻¹). Een hogere waarde wijst op een snellere initiële groei.
- `d`: De afbraaksnelheid (dag⁻¹). Dit beschrijft hoe snel de tumor krimpt door interne processen.
- `t`: Tijd (dagen).

### Mendelsohn model

Het Mendelsohn groei model:
dV/dt = c * V^d

- `V`: Het tumorvolume (mm³).
- `c`: De groeifactor (dag⁻¹). Het bepaalt de waarde van het initiële groei.
- `d`: De groeisnelheid van het tumorvolume. Dit beschrijft het snelheid hoe snel de tumor groeit.
- `t`: Tijd (dagen).

### Montroll model

Het Montroll groei model:
dV/dt = c * V *(Vmax^d - V^d)

- `V`: Het tumorvolume (mm³)
- `c`: De groeifactor (dag⁻¹). Het bepaalt de waarde van het initiële groei.
- `Vmax`: Het maximumvolume dat het tumor kan bereiken vanwege fysieke of biologische beperkingen.
- `d`: De groeisnelheid van het tumorvolume. Dit beschrijft hoe sterk het groei en afbraak is van het tumorvolume. 
- `t`: Tijd (dagen).


### Allee effect model

Het Allee effect groei model:
dV/dt = c * (V-Vmin) * (Vmax - V)

- `V`: Het tumorvolume (mm³)
- `c`: De groeifactor (dag⁻¹). Het bepaalt de waarde van het initiële groei.
- `Vmin`: Het minimumvolume dat het tumorvolume kan zijn. Dit beschrijft het threshold waarde van het groei.
- `Vmax`: Het maximumvolume dat het tumor kan bereiken vanwege fysieke of biologische beperkingen.
- `t`: Tijd (dagen).
 
---

## Numerieke Methoden

Aangezien de modellen gebaseerd zijn op niet-lineaire differentiaalvergelijkingen, gebruiken we numerieke methoden om ze op te lossen:

### Runge-Kutta Methode (Vierde Orde)

De **Runge-Kutta-methode** van de vierde orde is een veelgebruikte techniek voor numerieke integratie van gewone differentiaalvergelijkingen (ODE's). Deze methode wordt vaak toegepast wanneer een analytische oplossing moeilijk te verkrijgen is. De vierde-orde methode biedt een uitstekende balans tussen nauwkeurigheid en rekenefficiëntie.

#### Algemene Formule

Het algemene formule voor de vierde-orde Runge-Kutta-methode is:

$$
V_{\text{new}} = V + \frac{k_1 + 2k_2 + 2k_3 + k_4}{6}
$$

waarbij:
- $(V_{\text{new}}\)$ de geschatte waarde is na één stap van de integratie.
- $(V\)$ de huidige waarde van de functie is.
- $(\Delta t\)$ de tijdstap is.

#### Berekening van de coëfficiënten

De waarden $(k_1\)$, $(k_2\)$, $(k_3\)$, en $(k_4\)$ worden als volgt berekend:

1. **Eerste coëfficiënt (k1):**
   
$$
k_1 = \Delta t \cdot f(t, V)
$$
   Dit is de snelheid van verandering van de functie op het begin van de tijdstap.

2. **Tweede coëfficiënt (k2):**
   
$$
k_2 = \Delta t \cdot f\left(t + \frac{\Delta t}{2}, V + \frac{k_1}{2}\right)
$$

   Dit is de snelheid van verandering van de functie halverwege de tijdstap, gecorrigeerd door de helft van \(k_1\).

3. **Derde coëfficiënt (k3):**
   
$$
k_3 = \Delta t \cdot f\left(t + \frac{\Delta t}{2}, V + \frac{k_2}{2}\right)
$$

   Dit is de snelheid van verandering op hetzelfde halverwege-tijdstip als \(k_2\), maar met een interim-waarde gecorrigeerd door $(k_2\)$.

4. **Vierde coëfficiënt (k4):**
   
$$
k_4 = \Delta t \cdot f(t + \Delta t, V + k_3)
$$

   Dit is de snelheid van verandering aan het einde van de tijdstap, gecorrigeerd door $(k_3\)$.

#### Toepassing

Na het berekenen van de vier coëfficiënten $(k_1\)$, $(k_2\)$, $(k_3\)$, en $(k_4\)$, wordt de nieuwe waarde $(V_{\text{new}}\)$ van de functie berekend met de gewogen som van deze coëfficiënten:

$$
V_{\text{new}} = V + \frac{k_1 + 2k_2 + 2k_3 + k_4}{6}
$$

Deze gewogen gemiddelde aanpak zorgt ervoor dat de methode nauwkeuriger is dan eenvoudiger methoden zoals de **Euler-methode**.



## Voorbeeld Code (Python)

```python
def runge_kutta_4(f, V, t, dt):
    k1 = dt * f(t, V)
    k2 = dt * f(t + dt / 2, V + k1 / 2)
    k3 = dt * f(t + dt / 2, V + k2 / 2)
    k4 = dt * f(t + dt, V + k3)
    
    V_new = V + (k1 + 2*k2 + 2*k3 + k4) / 6
    return V_new
```

---

## Functies en Functionaliteit

### Functies:
- `runga_method`: Voert numerieke integratie uit voor een gegeven groeifunctie met behulp van de Runge-Kutta 4e orde methode.
- `calculate_residuals`: Berekent de Residual Sum of Squares (RSS) tussen gemeten en gesimuleerde volumedata.
- `calculate_aic`: Berekent de Akaike Information Criterion (AIC) voor modelselectie.
- `calculate_bic`: Berekent de Bayesian Information Criterion (BIC) voor modelselectie.
### Specifieke groeimodellen
#### Elk van de volgende klassen representeert een specifiek groeimodel:
- `LogisticModel`: Voor het logistische groeimodel, met simulatie en groeifunctie.
- `GompertzModel`: Voor het Gompertz groeimodel, met simulatie en groeifunctie.
- `VonBertalanffyModel`: Voor het Von Bertalanffy groeimodel, met simulatie en groeifunctie.
- `MendelsohnModel`: Voor het Mendelsohn groeimodel, met simulatie en groeifunctie.
- `MontrollModel`: Voor het Montroll groeimodel, met simulatie en groeifunctie.
- `AlleeModel`: Voor het Allee groeimodel, met simulatie en groeifunctie.

#### Elk specifiek model klasse bevat:
- `growth`: De groeifunctie die de groeisnelheid berekent.
- `simulate`: Voert een simulatie uit van het model over een gegeven tijdsinterval met behulp van de Runge-Kutta methode.

### Genereren nepdataset
#### Datahandeler
- `generate_data`: Genereert een tijds data en volume data met behulp van random getallen.

### Model evaluatie
#### Evaluator:
- `fit_and_evaluate`: Fitte een model aan de gegevens en evalueer de kwaliteit met behulp van AIC en BIC.
- `compare_models`: Vergelijk meerdere groeimodellen op basis van AIC en BIC en retourneer een DataFrame met resultaten.




---

## Installatie

1. Clone de repository:
   ```bash
   git clone https://github.com/jouw-repo/tumorgroei-modellen.git
   ```
2. Installeer de packages
   ```bash
   pip install numpy pandas matplotlib random math
   ```
---

## Gebruik

### Stap 1: Importeer de Module

De module moet geïmporteerd worden vanuit de Python-bestand waarin de class is opgeslagen (test.py in dit geval).

```python
from basemodel import DataHandler, Evaluator, LogisticModel, GompertzModel, VonBertalanffyModel, MendelsohnModel, MontrollModel, AlleeModel
```

### Stap 2: Genereer of definieer je eigen Data

Je kunt de nepdatasets die door de module worden gegenereerd gebruiken, of je kunt je eigen tijd- en volumegegevens instellen.

```python

# Voorbeeld van nepdataset genereren

# Bij het genereren van t_data neemt hij random getalen van een lengte dat wordt meegegeven.
t_data = [random() for _ in range(length)]

# Bij het genereren van V_data neemt hij een exponentiele waarde -1 voor elke tijds punt in t_data.
V_data = [exp(data_t) - 1.0 for data_t in t_data]

# Of je kunt je eigen data instellen
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
```

### Stap 3: Stel de Tijdspanne in voor Simulatie (Optioneel)

Je kunt de tijdspanne voor de simulatie aanpassen via de t_vooruit parameter. Als deze parameter niet wordt opgegeven, wordt de standaardwaarde van np.linspace(0, 120, 100) gebruikt, wat betekent dat de simulatie loopt van 0 tot 120 dagen met 100 punten.

```python

# Stel een eigen tijdspanne in voor betere resultaten
t_vooruit = np.linspace(0, 150, 150)  # Van 0 tot 150 dagen met 150 punten
```
### Stap 4: Voer Model Evaluatie en Visualisatie uit

```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from test import DataHandler, Evaluator, LogisticModel, GompertzModel, VonBertalanffyModel, MendelsohnModel, MontrollModel, AlleeModel

# Voeg hier je eigen data toe
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
# Genereer data
t_data, V_data = DataHandler.generate_data(length) # Geef lengte mee van de data, die je wilt genereren

# Instantieer de evaluator
evaluator = Evaluator(t_data, V_data)

# Definieer de modellen en beginwaarden voor parameters
t_forecast = np.linspace(0, 60, 45)
models = {
    'Logistic': (LogisticModel, [0.01, 7]),
    'Gompertz': (GompertzModel, [0.11, 7.5]),
    'Von Bertalanffy': (VonBertalanffyModel, [0.5, 0.2]),
    'Mendelsohn': (MendelsohnModel, [0.01, 0.1]),
    'Montroll': (MontrollModel, [0.01, 8, 0.1]),
    'Allee': (AlleeModel, [0.05, 0, 7.5])
}

# Modelvergelijking
results = evaluator.compare_models(models, t_forecast)

# Visualiseer de resultaten
plt.figure(figsize=(10, 6))
plt.scatter(t_data, V_data, color="red", label="Data")

for model_name, result in results.items():
    params = result['params']  # Zorg dat 'params' onderdeel is van result
    param_str = ", ".join([f"{p:.4f}" for p in params])  # Formatteer de parameters
    plt.plot(t_forecast, result['V_sim'], label=f"{model_name} (params: [{param_str}])")

plt.title("Tumorgroei Modellen vs. Data")
plt.xlabel("Tijd (dagen)")
plt.ylabel("Tumorvolume (mm³)")
plt.legend()
plt.grid(True)
plt.show()

# Toon de statistische vergelijking
df = pd.DataFrame.from_dict(results, orient='index')
df = df[['AIC', 'BIC', 'params']].sort_values(by='AIC')  # Voeg 'params' toe
print(df)

```
###
Voor een uitgebreidere uitleg zie toelichting.ipynb
## Resultaten

Na het uitvoeren van de bovenstaande code:

Visualisatie: Er worden grafieken getoond van de tumorgroei volgens de zes modellen (Gompertz, Logistic, Von Bertalanffy, Mendelsohn, Montroll en Allee-effect) in vergelijking met de werkelijke gegevens.

Model Evaluatie: De AIC- en BIC-waarden en optimale parameters worden berekend voor elk model en weergegeven in de console. Deze waarden helpen je bij het kiezen van het beste model. Uitgebreidere uitleg over wat wat is te vinden in de Jupyter Notebook toelichting.ipynb en laat zien hoe het module wordt toegepast p een realistische dataset.

## Problemen

Bij problemen in de module kun je een mail sturen (zie de **Authors** sectie).


## Authors and Acknowledgments

Developed by:
- Jort Gommers: [j.r.gommers@st.hanze.nl](mailto:j.r.gommers@st.hanze.nl)
- Akastia Christo: [m.a.christo@st.hanze.nl](mailto:m.a.christo@st.hanze.nl)

## License
No specific licensing applies 

## Changelog

- **v1.0**: Initial release with basic functionality.


### Referenties
1. *Gompertz function*. [Link](https://www.tmlep.com/clinical-learning/2023-01-23-when-did-this-tumour-start-the-need-for-a-gompertzian-understanding-of-tumour-growth-kinetics)
2. *Logistic regression*. [Link](https://www.spiceworks.com/tech/artificial-intelligence/articles/what-is-logistic-regression/)
3. Derek H. Ogle (2006). *Growth (von Bertalanffy) Notes*. [Link](https://derekogle.com/NCNRS349/modules/PREP/NOTES/Growth)
4. *Runge-Kutta Method*. [link](https://www.sciencedirect.com/topics/mathematics/runge-kutta-method)
5. *Heun's method*. [link](https://en.wikipedia.org/wiki/Heun%27s_method)
6. *Mendelsohn*, [link](https://bmccancer.biomedcentral.com/articles/10.1186/s12885-016-2164-x)
7. *Montroll* [link](https://www.mdpi.com/2227-7390/12/8/1195)
8. *Allee effect* [link] (https://www.sciencedirect.com/science/article/abs/pii/S1476945X16300745)
9. *Realistische data* [link] (https://journals.plos.org/ploscompbiol/article?id=10.1371%2Fjournal.pcbi.1007178#sec002)
