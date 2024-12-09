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
\[
\frac{dV}{dt} = c \cdot V \cdot \ln\left(\frac{V_{\text{max}}}{V}\right)
\]

**Parameters**:
- \( V \): Het tumorvolume op een bepaald tijdstip \( t \) (mm³).
- \( c \): De groeisnelheidsparameter (\( \text{dag}^{-1} \)). Een hogere waarde geeft een snellere initiële groei aan.
- \( V_{\text{max}} \): Het maximale tumorvolume (mm³), dat de asymptotische limiet voor de groei vertegenwoordigt.
- \( t \): Tijd (dagen).

### Logistisch model

De logistische groeivergelijking is:
\[
\frac{dV}{dt} = c \cdot V \cdot \left(1 - \frac{V}{V_{\text{max}}}\right)
\]

**Parameters**:
- \( V \): Het tumorvolume (mm³).
- \( c \): De groeisnelheid (\( \text{dag}^{-1} \)). Deze parameter bepaalt hoe snel de tumor zijn maximale capaciteit nadert.
- \( V_{\text{max}} \): Het maximumvolume dat de tumor kan bereiken vanwege fysieke of biologische beperkingen.
- \( t \): Tijd (dagen).

**Gedrag**:
- Bij \( V \ll V_{\text{max}} \): Groeit de tumor bijna exponentieel.
- Bij \( V \to V_{\text{max}} \): Neemt de groeisnelheid af tot nul.

### Von Bertalanffy-model

Het Von Bertalanffy-model beschrijft de balans tussen groei en afbraak:
\[
\frac{dV}{dt} = c \cdot V^{\frac{2}{3}} - d \cdot V
\]

**Parameters**:
- \( V \): Het tumorvolume (mm³).
- \( c \): De groeifactor (\( \text{dag}^{-1} \)). Een hogere waarde wijst op een snellere initiële groei.
- \( d \): De afbraaksnelheid (\( \text{dag}^{-1} \)). Dit beschrijft hoe snel de tumor krimpt door interne processen.
- \( t \): Tijd (dagen).

**Gedrag**:
- Het model is geschikt voor systemen waar de groei wordt bepaald door het oppervlak (\( V^{2/3} \)) en de afbraak proportioneel is aan het volume.

---

## Numerieke Methoden

Aangezien de modellen gebaseerd zijn op niet-lineaire differentiaalvergelijkingen, gebruiken we numerieke methoden om ze op te lossen:

### Heun-methode
Een verbeterde versie van de Euler-methode, die een correctiefactor toevoegt om nauwkeuriger te zijn:
\[
V_{n+1} = V_n + \Delta t \cdot \frac{f(t_n, V_n) + f(t_{n+1}, V_e)}{2}
\]
Hierbij wordt \( f(t, V) \) gegeven door het groeimodel.

### Runge-Kutta-methode (vierde orde)
Een populaire methode voor numerieke integratie:
\[
V_{n+1} = V_n + \frac{1}{6}(k_1 + 2k_2 + 2k_3 + k_4)
\]
waarbij:
\[
k_1 = f(t_n, V_n), \quad k_2 = f(t_n + \frac{\Delta t}{2}, V_n + \frac{\Delta t}{2} k_1)
\]
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

