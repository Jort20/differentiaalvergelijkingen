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
7. [Resultaten en Visualisatie](#resultaten-en-visualisatie)
8. [Referenties](#referenties)

---

## Achtergrond

De groei van tumoren is een complex proces dat vaak wordt gemodelleerd met behulp van differentiaalvergelijkingen. Deze modellen simuleren hoe een tumor in volume groeit over de tijd en helpen onderzoekers om groeiprofielen te begrijpen en behandelingsstrategieën te ontwikkelen. 

De drie modellen in dit project hebben elk unieke eigenschappen die ze geschikt maken voor verschillende toepassingen:
- **Gompertz-model**: Wordt vaak gebruikt in tumorbiologie vanwege zijn flexibiliteit en realistische groei-limieten.
- **Logistisch model**: Eenvoudig en effectief, vaak gebruikt voor populatiedynamica of groei met beperkende factoren.
- **Von Bertalanffy-model**: Beschrijft zowel groei als afbraak, nuttig voor biologie en fysiologie.

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

