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

Het Allee effect model
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

De algemene formule voor de vierde-orde Runge-Kutta-methode is:

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

Visualisatie: Er worden grafieken getoond van de tumorgroei volgens de zes modellen (Gompertz, Logistic, Von Bertalanffy, Mendelsohn, Montroll en Allee-effect) in vergelijking met de werkelijke gegevens.

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
1. *Gompertz function*. [Link](https://www.tmlep.com/clinical-learning/2023-01-23-when-did-this-tumour-start-the-need-for-a-gompertzian-understanding-of-tumour-growth-kinetics)
2. *Logistic regression*. [Link](https://en.wikipedia.org/wiki/Logistic_regression)
3. Derek H. Ogle (2006). *Growth (von Bertalanffy) Notes*. [Link](https://derekogle.com/NCNRS349/modules/PREP/NOTES/Growth)
4. *Runge-Kutta Method*. [link](https://www.sciencedirect.com/topics/mathematics/runge-kutta-method)
5. *Heun's method*. [link](https://en.wikipedia.org/wiki/Heun%27s_method)
6. *Mendelsohn*, [link](https://bmccancer.biomedcentral.com/articles/10.1186/s12885-016-2164-x)
7. *Montroll* [link](https://www.mdpi.com/2227-7390/12/8/1195)
8. *Allee effect* [link} (https://www.sciencedirect.com/science/article/abs/pii/S1476945X16300745)



