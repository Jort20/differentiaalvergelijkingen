# Tumorgroei Modellen met Python

Dit project implementeert en vergelijkt verschillende wiskundige modellen voor tumorgroei, inclusief het Gompertz-, Logistic- en Von Bertalanffy-model. De implementatie omvat methoden om deze modellen te simuleren, te fitten aan data en de resultaten te visualiseren. Daarnaast worden numerieke methoden zoals de Heun- en Runge-Kutta-methode gebruikt om de modellen op te lossen.

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

De drie modellen in dit project hebben elk unieke eigenschappen die ze geschikt maken voor verschillende toepassingen:
- **Gompertz-model**: Wordt vaak gebruikt in tumorbiologie vanwege zijn flexibiliteit en realistische groei-limieten.
- **Logistisch model**: Eenvoudig en effectief, vaak gebruikt voor populatiedynamica of groei met beperkende factoren.
- **Von Bertalanffy-model**: Beschrijft zowel groei als afbraak, word gebruikt in biologie en fysiologie.

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



---

## Numerieke Methoden

Aangezien de modellen gebaseerd zijn op niet-lineaire differentiaalvergelijkingen, gebruiken we numerieke methoden om ze op te lossen:

### Heun-methode
Een verbeterde versie van de Euler-methode, die een correctiefactor toevoegt om nauwkeuriger te zijn:

V_new = V + Δt * (f(t, V) + f(t + Δt, V_euler)) / 2

Hierbij wordt \( f(t, V) \) gegeven door het groeimodel.

### Runge-Kutta-methode (vierde orde)
Een populaire methode voor numerieke integratie:

V_new = V + (k1 + 2k2 + 2k3 + k4) / 6

waarbij:

k1 = Δt * f(t, V) k2 = Δt * f(t + Δt/2, V + k1/2) k3 = Δt * f(t + Δt/2, V + k2/2) k4 = Δt * f(t + Δt, V + k3)

Enzovoorts, met hogere orde correcties.

---

## Functies en Functionaliteit

- **`generate_fake_data`**: Genereert gesimuleerde data voor tumorvolumes en tijd.
- **`fit_model`**: Past een model aan op de data om de optimale parameters te vinden.
- **`logistic_wrapper`, `gompertz_wrapper`, `von_bertalanffy_wrapper`**: Wrapper-functies voor compatibiliteit met `curve_fit`.
- **Visualisatie**: Grafieken van modelvergelijkingen en data met Matplotlib.

---

## Installatie

1. Clone de repository:
   ```bash
   git clone https://github.com/jouw-repo/tumorgroei-modellen.git
   ```
2. Installeer de packages
   ```bash
   pip install numpy scipy matplotlib
   ```
---

## Gebruik

### Stap 1: Importeer de Module

De module moet geïmporteerd worden vanuit de Python-bestand waarin de class is opgeslagen (test.py in dit geval).

```python
from test import TumorGrowthModels
```

### Stap 2: Genereer of definieer je eigen Data

Je kunt de nepdatasets die door de module worden gegenereerd gebruiken, of je kunt je eigen tijd- en volumegegevens instellen.
```python
# Voorbeeld van nepdataset genereren
t_data, V_data = TumorGrowthModels.generate_fake_data()

# Of je kunt je eigen data instellen
t_data = np.array([0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
V_data = np.array([250, 300, 450, 600, 750, 1000, 1300, 1700, 2200, 2700, 3200])
```

### Stap 3: Maak een Model aan

Maak een object van de TumorGrowthModels class, geef de tijd- en volumegegevens door aan de constructor.

```python
model = TumorGrowthModels(t_data, V_data)
```
### Stap 4: Stel de Tijdspanne in voor Simulatie (Optioneel)

Je kunt de tijdspanne voor de simulatie aanpassen via de t_vooruit parameter. Als deze parameter niet wordt opgegeven, wordt de standaardwaarde van np.linspace(0, 120, 100) gebruikt, wat betekent dat de simulatie loopt van 0 tot 120 dagen met 100 punten.

```python

# Stel een eigen tijdspanne in
t_vooruit = np.linspace(0, 150, 150)  # Van 0 tot 150 dagen met 150 punten
```
### Stap 5: Voer Model Evaluatie en Visualisatie uit

Roep de evaluate_models methode aan om de modellen te fitten, de simulaties uit te voeren, de resultaten te visualiseren en de AIC/BIC-waarden te berekenen.

```python
# Voer model evaluatie uit en visualiseer de resultaten
model.evaluate_models(t_vooruit)
```
### Volledige Voorbeeld

from test import TumorGrowthModels
import numpy as np

```python
# Genereer of stel je eigen tijd- en volumegegevens in
t_data = np.array([0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
V_data = np.array([250, 300, 450, 600, 750, 1000, 1300, 1700, 2200, 2700, 3200])

# Maak een model aan met je eigen data
model = TumorGrowthModels(t_data, V_data)

# Stel een eigen tijdspanne in voor de simulatie
t_vooruit = np.linspace(0, 150, 150)  # Van 0 tot 150 dagen met 150 punten

# Voer model evaluatie uit en visualiseer de resultaten
model.evaluate_models(t_vooruit)
```

## Resultaten

Na het uitvoeren van de bovenstaande code:

Visualisatie: Er worden grafieken getoond van de tumorgroei volgens de drie modellen (Gompertz, Logistic, Von Bertalanffy) in vergelijking met de werkelijke gegevens.

Model Evaluatie: De AIC- en BIC-waarden worden berekend voor elk model en weergegeven in de console. Deze waarden helpen je bij het kiezen van het beste model.
   

## Problemen

Bij problemen in de module kun je een mail sturen (zie de **Authors** sectie).


## Authors and Acknowledgments

Developed by:
- Jort Gommers: [j.r.gommers@st.hanze.nl](mailto:j.r.gommers@st.hanze.nl)

## License
No specific licensing applies 

## Changelog

- **v1.0**: Initial release with basic functionality.


### Referenties
1. *Gompertz function*. [Link](https://en.wikipedia.org/wiki/Gompertz_function)
2. *Logistic regression*. [Link](https://en.wikipedia.org/wiki/Logistic_regression)
3. Derek H. Ogle (2006). *Growth (von Bertalanffy) Notes*. [Link](https://derekogle.com/NCNRS349/modules/PREP/NOTES/Growth)
4. *Runge-Kutta Method*. [link](https://www.sciencedirect.com/topics/mathematics/runge-kutta-method)
5. *Heun's method*. [link](https://en.wikipedia.org/wiki/Heun%27s_method)


