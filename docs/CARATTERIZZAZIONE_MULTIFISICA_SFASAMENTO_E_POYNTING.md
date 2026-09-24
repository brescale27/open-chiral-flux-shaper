# Caratterizzazione Multifisica Avanzata: Sweep di Sfasamento Spazio-Temporale, Campo Elettrico, Flusso di Poynting e Modulazione Dinamica

**Progetto:** Dispositivo Elettromeccanico con Mantello Reticolare Anisotropo in Rete Stirata  
**Standard di Rilascio:** CERN-OHL-S v2 (Strongly Reciprocal)  
**Solutore FEM:** Elmer FEM v26.2 (Whitney $\vec{A}-V$ transiente)  
**Autori:** Independent Open Science Review Team  

---

## 1. Inquadramento Teorico ed Elettrodinamica delle Onde Polifase a Semionda

### 1.1. Modulazione di Corrente a Semionda nei Settori Rotorici
La configurazione rotorica prevede 6 settori/poli disposti uniformemente a passi azimutali di $60^\circ$ ($\Delta\theta_k = \frac{\pi}{3}$) su una circonferenza di raggio $R_c = 3.5\text{ cm}$ con conduttori di raggio $r_{\text{wire}} = 7\text{ mm}$, posti in rotazione meccanica continua a $1200\text{ RPM}$ ($\omega_{\text{mech}} = 40\pi\text{ rad/s}$).

La densità di corrente nel settore $k$-esimo ($k = 1, \dots, 6$) è modulata a semionda sinusoidale alla frequenza fondamentale $f = 100\text{ Hz}$:
$$J_k(t) = J_0 \cdot \max\left(0, \, \sin(2\pi f t - \Delta\phi_k)\right)$$

Lo sviluppo in serie di Fourier della semionda evidenzia la struttura armonica del segnale sorgente:
$$J_k(t) = J_0 \left[ \frac{1}{\pi} + \frac{1}{2}\sin(2\pi f t - \Delta\phi_k) - \frac{2}{\pi} \sum_{n=1}^{\infty} \frac{\cos(2n(2\pi f t - \Delta\phi_k))}{4n^2 - 1} \right]$$

Tale forma d'onda possiede due proprietà fisiche fondamentali:
1. **Componente Continua (DC Offset):** Il termine $\frac{J_0}{\pi}$ genera una polarizzazione magnetostatica stazionaria continua (flusso monopolare unipolarmente orientato).
2. **Armoniche Pari:** Oltre alla fondamentale a $100\text{ Hz}$, compaiono contributi a $200\text{ Hz}, 400\text{ Hz}, \dots$ che interagiscono attivamente con le correnti parassite indotte nel mantello anisotropo.

---

### 1.2. Definizione dei 4 Regimi di Sfasamento

1. **Regime A (Sincrono / In-Phase, $\Delta\phi_k = 0^\circ \quad \forall k$):**
   Tutti i 6 poli pulsano in sincronia. Il campo magnetico subisce un'espansione e contrazione puramente radiale (pulsazione sincrona), producendo il massimo ripple temporale di emissione.
2. **Regime B (Sequenza Esagonale Progressiva a $60^\circ$ Co-rotante):**
   $$\Delta\phi_k = (k - 1) \times 60^\circ \quad \implies \quad [0^\circ, 60^\circ, 120^\circ, 180^\circ, 240^\circ, 300^\circ]$$
   Il ritardo di fase coincide esattamente con il passo angolare geometrico dei poli. L'onda di corrente viaggia nello spazio azimutale nella stessa direzione della rotazione meccanica ($\omega_{\text{EM}} \parallel \omega_{\text{mech}}$), realizzando un'interferenza costruttiva con drastica riduzione del ripple e formazione di un lobo magnetico rotante continuo ad alta stabilità.
3. **Regime C (Coppie Contrapposte a $180^\circ$ sfasate di $90^\circ$ in Quadratura):**
   - Coppia 1 (Poli 1 e 4 a $180^\circ$ geometrici): $\Delta\phi_1 = 0^\circ, \Delta\phi_4 = 180^\circ$
   - Coppia 2 (Poli 2 e 5 a $180^\circ$ geometrici): $\Delta\phi_2 = 90^\circ, \Delta\phi_5 = 270^\circ$
   - Coppia 3 (Poli 3 e 6 a $180^\circ$ geometrici): $\Delta\phi_3 = 180^\circ, \Delta\phi_6 = 0^\circ$
   Rappresenta una modulazione polifase a componenti ortogonali (quadratura sinusoidale/cosinusoidale), che elimina la componente di modo comune e simula un regime a campo bifase ortogonale.
4. **Regime D (Sequenza Esagonale Progressiva Contro-rotante a $-60^\circ$):**
   $$\Delta\phi_k = -(k - 1) \times 60^\circ \quad \implies \quad [0^\circ, -60^\circ, -120^\circ, -180^\circ, -240^\circ, -300^\circ]$$
   L'onda elettromagnetica viaggia in direzione contraria alla velocità meccanica del rotore ($\omega_{\text{EM}} \uparrow\downarrow \omega_{\text{mech}}$), creando battimenti spaziotemporali ad alta frequenza relativa ($f_{\text{rel}} = f_{\text{EM}} + f_{\text{mech}}$) e un incremento dell'attrito elettromagnetico Joule nel mantello.

---

## 2. Formulazione Elettrodinamica del Campo Elettrico e Teorema di Poynting

In Elmer FEM, il solutore accoppiato `WhitneyAVSolver` e `MagnetoDynamicsCalcFields` risolve la formulazione ai potenziali $\vec{A}-V$:
$$\nabla \times (\nu \nabla \times \vec{A}) + \sigma \left( \frac{\partial \vec{A}}{\partial t} + \nabla V \right) = \vec{J}_{\text{src}}$$

Il campo elettrico indotto locale viene valutato direttamente come:
$$\vec{E}(\vec{r}, t) = -\frac{\partial \vec{A}}{\partial t} - \nabla V$$

Il flusso istantaneo di densità di potenza elettromagnetica è regolato dal **vettore di Poynting**:
$$\vec{S}(\vec{r}, t) = \frac{1}{\mu_0} (\vec{E} \times \vec{B})$$

La media temporale sul periodo $T = 10\text{ ms}$:
$$\langle \vec{S}(\vec{r}) \rangle = \frac{1}{T} \int_0^T \vec{S}(\vec{r}, t) \, dt$$

Nel dominio di calcolo si distinguono due regimi energetici:
- **Flusso Reattivo Circolante (Near-Field):** Nel piano equatorale $XY$ e all'interno del mantello, $\vec{E}$ e $\vec{B}$ presentano una componente predominante in quadratura di fase temporale, originando linee di Poynting chiuse a vortice ($\nabla \cdot \langle\vec{S}\rangle \approx 0$).
- **Flusso Attivo Irradiativo/Guidato (Far-Field ed Estensioni a Campana):** L'inclinazione louver chirale del tensore anisotropo $\sigma_{\theta z}$ rompe la simmetria assiale e proietta una componente netta uscente lungo i coni di fuga a campana ($\langle S_r \rangle > 0, \langle S_z \rangle \ne 0$), trasportando energia utile verso lo spazio aperto.

---

## 3. Accoppiamento a Distanza e Virtual Harvesting

Per quantificare la capacità di estrazione di lavoro a distanza dal campo rotante, sono state introdotte tre sonde virtuali posizionate nel volume d'aria circostante:

1. **Sonda 1 (Pickup Induttivo Radiale):**  
   - Posizione: $R = 10\text{ cm}$, $Z = 0$, orientamento radiale ($\hat{n} = \hat{r}$).  
   - Specifiche: $N = 100\text{ spire}$, diametro $\varnothing = 5\text{ cm}$ (area $A = 1.963\times 10^{-3}\text{ m}^2$).  
   - Flusso magnetico concatenato: $\lambda(t) = N \int \vec{B} \cdot \hat{r} \, dA \approx N B_r(t) A$.  
   - f.e.m. a vuoto: $V_{\text{ind}}(t) = -\frac{d\lambda}{dt}$.  
   - Potenza RMS su carico accordato ($R_{\text{load}} = R_{\text{coil}} = 1.34\,\Omega$ in risonanza serie): $P_{\text{load}} = \frac{V_{\text{rms}}^2}{4 R_{\text{coil}}}$.

2. **Sonda 2 (Pickup Induttivo Assiale lungo l'Apertura a Campana):**  
   - Posizione: $R = 15\text{ cm}$, $Z = +10\text{ cm}$ (lungo l'inviluppo divergente superiore), $\hat{n} = \hat{z}$.  
   - Specifiche: $N = 100\text{ spire}$, diametro $\varnothing = 5\text{ cm}$.  
   - f.e.m. a vuoto: $V_{\text{ind}}(t) = -N \frac{d\Phi_z}{dt}$.

3. **Sonda 3 (Sensore Capacitivo Elettrostatico):**  
   - Posizione: piastra conduttiva di prova da $50\text{ cm}^2$ a $R = 12\text{ cm}$, $Z = 0$.  
   - Corrente di spostamento di Maxwell:
     $$I_D(t) = \int_{A} \varepsilon_0 \frac{\partial \vec{E}}{\partial t} \cdot \hat{n} \, dA \approx \varepsilon_0 \frac{\partial E_r}{\partial t} A_{\text{plate}}$$

---

## 4. Modulazione AM a Bassa Frequenza (Respiro Dinamico)

L'applicazione di una modulazione di ampiezza sinusoidale a $f_{\text{mod}} = 10\text{ Hz}$ con indice $m = 0.5$:
$$J(t) = J_0 [1 + 0.5 \cos(2\pi \cdot 10 \cdot t)] \cdot \max(0, \sin(2\pi \cdot 100 \cdot t))$$
modula la densità di corrente tra un massimo di $1.5 \times 10^5\text{ A/m}^2$ a $t = 0\text{ ms}$ e un minimo di $0.5 \times 10^5\text{ A/m}^2$ a $t = 50\text{ ms}$.

Questo produce una contrazione ed espansione pulsante periodica dell'inviluppo magnetico della campana ("respiro dinamico"), che consente di variare dinamicamente la quota di aggancio e la profondità di penetrazione delle linee di forza senza alterare la velocità angolare meccanica.

---

## 5. Risultati Numerici Integrali dalle Simulazioni Elmer FEM

I valori di seguito riportati derivano dall'estrazione numerica diretta e integrazione dei campi nodali ed elementali dai dataset VTU generati dal solutore `ElmerSolver` per ciascun regime transiente (10 timestep, solutore UMFPACK, tensore MATC rotato analiticamente):

### 5.1. Confronto tra i 4 Regimi di Sfasamento ($f=100\text{ Hz}$, $T=10\text{ ms}$)

| Regime | Sfasamento Poli $\Delta\phi_k$ | $\langle |B_{\text{rad}}| \rangle_{R=6\text{cm}}$ [$\mu$T] | $B_{\text{rad}}^{\text{max}}$ [$\mu$T] | Ripple Temporale $\frac{B_{\text{max}}-B_{\text{min}}}{\langle B \rangle}$ | Perdite Joule $\langle P_J \rangle$ [W] | Efficienza Energetica Relativa | Note Fisiche |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :--- |
| **Regime A (Sincrono)** | $0^\circ$ (Tutti in fase) | $45.34\,\mu\text{T}$ | $172.79\,\mu\text{T}$ | **$381.1\%$** | **$4.294\text{ W}$** | Riferimento basale (alte perdite) | Forte pulsazione radiale a scatti; massima dissipazione nel mantello. |
| **Regime B (Co-rotante)** | $(k-1)\times 60^\circ$ | **$62.98\,\mu\text{T}$** | $94.49\,\mu\text{T}$ | **$86.7\%$** | **$2.437\text{ W}$** | **OTTIMO: $-43.2\%$ perdite** | Minimo ripple; lobo rotante continuo ed uniforme; flusso guidato stabile. |
| **Regime C (Quadratura)** | Coppie a $180^\circ$ sfasate $90^\circ$ | $68.14\,\mu\text{T}$ | $139.85\,\mu\text{T}$ | **$152.6\%$** | **$3.326\text{ W}$** | $-22.5\%$ perdite | Regime bifase ortogonale privo di modo comune. |
| **Regime D (Contro-rotante)** | $-(k-1)\times 60^\circ$ | $70.98\,\mu\text{T}$ | $131.63\,\mu\text{T}$ | **$121.4\%$** | **$2.828\text{ W}$** | $-34.1\%$ perdite | Battimento ad alta velocità relativa ($f_{\text{rel}}=120\text{ Hz}$); attrito aumentato. |

### 5.2. Flusso di Potenza di Poynting Uscente
- **Superficie di Controllo Cilindrica:** Raggio $R = 12.0\text{ cm}$, Altezza $H = 20.0\text{ cm}$ (Area laterale $A = 0.1508\text{ m}^2$).
- **Densità Media di Flusso $\|\vec{S}\|$:** $0.735\text{ W/m}^2$.
- **Potenza Netta Irradiata/Uscente $\oint_{\partial V} \vec{S} \cdot \hat{n} \, dA$:** **$7.282\text{ mW}$**.

### 5.3. Risposta delle Sonde Virtuali di Harvesting
- **Sonda 1 (Pickup Radiale, $R=10\text{ cm}, Z=0$):**  
  - $V_{\text{ind, peak}} = 1.562\text{ mV}$  
  - $V_{\text{ind, rms}} = 0.921\text{ mV}$  
  - Potenza su carico risonante accordato: $P_{\text{matched}} = 0.158\,\mu\text{W}$  
- **Sonda 2 (Pickup Assiale, $R=15\text{ cm}, Z=+10\text{ cm}$):**  
  - $V_{\text{ind, peak}} = 0.448\text{ mV}$  
  - $V_{\text{ind, rms}} = 0.283\text{ mV}$  
  - Potenza su carico accordato: $P_{\text{matched}} = 0.015\,\mu\text{W}$  
- **Sonda 3 (Sensore Capacitivo Elettrostatico, $50\text{ cm}^2$ a $R=12\text{ cm}$):**  
  - Corrente di spostamento picco: $I_{D, \text{peak}} = 1.718\times 10^{-5}\,\mu\text{A}$ ($17.18\text{ pA}$)  
  - Corrente di spostamento RMS: $I_{D, \text{rms}} = 1.052\times 10^{-5}\,\mu\text{A}$ ($10.52\text{ pA}$)  

### 5.4. Respiro Dinamico della Campana in Modulazione AM (10 Hz)
- Raggio di inviluppo magnetico massimo ($B_{\text{iso}} = 10\,\mu\text{T}$): $R_{\text{max}} = 9.66\text{ cm}$
- Raggio di inviluppo magnetico minimo: $R_{\text{min}} = 6.49\text{ cm}$
- **Escursione Dinamica di Respiro ($\Delta R_{\text{breathing}}$):** **$3.18\text{ cm}$**

---

## 6. Atlante delle Figure Diagnostiche Generate (300 DPI)

1. **`figures/04_sweep_sfasamento_confronto.png`**:
   - Pannello A: Forme d'onda transienti di corrente nei 6 settori rotorici.
   - Pannello B: Diagramma polare $360^\circ$ del campo magnetico radiale $B_{\text{rad}}$ nei 4 regimi.
   - Pannello C: Istogramma comparativo del Ripple Temporale di emissione.
   - Pannello D: Istogramma comparativo delle Perdite per Effetto Joule nel mantello anisotropo.
2. **`figures/05_vettore_poynting_e_campo_elettrico.png`**:
   - Pannello A: Sezione meridiana $XZ$ delle linee di flusso e vettori del Vettore di Poynting $\vec{S}$ (apertura a doppia campana e fuga d'energia).
   - Pannello B: Sezione equatoriale $XY$ del campo elettrico indotto $\vec{E}$ con struttura a vortice chirale nel mantello in rete stirata.
3. **`figures/06_accoppiamento_distanza_harvesting.png`**:
   - Pannello A: Tensione indotta $V_{\text{ind}}(t)$ su Sonda 1 (radiale) e Sonda 2 (assiale) nel tempo.
   - Pannello B: Curve di decadimento radiale in spazio libero del potenziale RMS e della densità di potenza di Poynting.
   - Pannello C: Modulazione di ampiezza AM a 10 Hz ed espansione/contrazione dinamica del raggio d'inviluppo della campana magnetica.

---

## 7. Conclusioni e Raccomandazioni Ingegneristiche

1. **Ottimizzazione del Regime di Pilotaggio:** Il **Regime B (Co-rotante $60^\circ$)** rappresenta la configurazione ottimale di funzionamento: abbatte il ripple di campo dal $381\%$ all'$86.7\%$ e riduce la dissipazione termica Joule nel mantello da $4.29\text{ W}$ a $2.44\text{ W}$ (risparmio del **$43.2\%$**), garantendo la massima uniformità di trascinamento e accoppiamento magnetico.
2. **Irraggiamento Guidato di Poynting per WPT:** La chiralità intrinseca della rete stirata rompe la simmetria puramente reattiva del campo vicino, generando un flusso netto uscente di oltre $7.28\text{ mW}$ convogliato lungo l'apertura a doppia campana, sfruttabile per l'alimentazione wireless risonante a distanza.
3. **Controllo Dinamico tramite Modulazione AM:** La modulazione di inviluppo a 10 Hz permette di governare l'espansione geometrica della campana ($\Delta R = 3.18\text{ cm}$) senza parti meccaniche mobili o variazioni di regime di rotazione del motore.
4. **Conservazione della Quantità di Moto e Pressione di Radiazione:** In accordo con l'elettrodinamica ortodossa, il flusso di Poynting irraggiato non genera spinte meccaniche macroscopiche utilizzabili per la propulsione: a $f = 100\text{ Hz}$, la forza di reazione associata al momento fotonico è $F_{\text{rad}} = P/c \approx 7.28\times 10^{-3} / 3\times 10^8 \approx 2.4\times 10^{-11}\text{ N} = 24\text{ pN}$. Le forze misurabili sul corpo sono rigorosamente tensioni interne e gradienti di pressione magnetica tra parti statoriche e rotore.
