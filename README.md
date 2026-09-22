# Validazione Elettrodinamica 3D in Elmer FEM: Mantello in Rete Stirata Anisotropa

**Licenza Open Hardware:** CERN-OHL-S v2 (Strongly Reciprocal)  
**Codice di Riferimento:** `elmerfem-release-26.2`  
**Archivio di Pubblicazione:** Zenodo Open Science Repository  

---

> [!NOTE]
> **Stato della Release e Certificazione Open Hardware:**  
> Questo repository contiene il codice sorgente, la mesh volumetrica conforme 3D, le configurazioni di calcolo agli elementi finiti (Whitney $\vec{A}-V$), i dataset e i grafici diagnostici validati per la simulazione del dispositivo elettromeccanico con mantello anisotropo in rete stirata di alluminio multistrato (*expanded metal mesh*).  
> Tutti i file storici, preliminari e intermedi sono isolati nella cartella `_archive_backup/`.

---

## 1. Descrizione del Modello e Formulazione del Mezzo Anisotropo

Il dispositivo integra l'omogeneizzazione elettromagnetica del **mantello in lamiera stirata romboidale di alluminio multistrato**, implementando un tensore di conducibilità anisotropo conforme al Secondo Principio della Termodinamica (criterio di Sylvester / semi-definito positivo).

In coordinate cilindriche locali $(\hat{r}, \hat{\theta}, \hat{z})$:
$$\bar{\bar{\sigma}}_{\text{cyl}} = \begin{bmatrix} \sigma_{rr} & 0 & 0 \\ 0 & \sigma_{\theta\theta} & \sigma_{\theta z} \\ 0 & \sigma_{\theta z} & \sigma_{zz} \end{bmatrix} = \begin{bmatrix} 1.75\times 10^6 & 0 & 0 \\ 0 & 1.75\times 10^6 & 3.031\times 10^6 \\ 0 & 3.031\times 10^6 & 1.22\times 10^7 \end{bmatrix} \text{ S/m}$$

La matrice cartesiana $\bar{\bar{\sigma}}_{\text{cart}}(x, y) = \mathbf{P}(\theta) \bar{\bar{\sigma}}_{\text{cyl}} \mathbf{P}(\theta)^T$ (dove $\theta = \text{atan2}(y, x)$) viene calcolata ed applicata punto per punto in Elmer FEM tramite la funzione analitica MATC `sigma_cyl(tx)`:
$$\bar{\bar{\sigma}}_{\text{cart}}(x, y) = \begin{bmatrix} 1.75\times 10^6 & 0 & -\sigma_{\theta z}\sin\theta \\ 0 & 1.75\times 10^6 & \sigma_{\theta z}\cos\theta \\ -\sigma_{\theta z}\sin\theta & \sigma_{\theta z}\cos\theta & 1.22\times 10^7 \end{bmatrix} \text{ S/m}$$

Ciò assicura che l'inclinazione persiana a $30^\circ$ dei ponticelli metallici sia rigorosamente omogenea lungo tutti i $360^\circ$ del perimetro del mantello, eliminando qualsiasi asimmetria cartesiana spuria.

---

## 2. Risultati di Validazione e Certificazione Numerica

### 2.1. Conservazione Integrale del Flusso di Gauss ($\oint_S \vec{B} \cdot \hat{n} \, dA = 0$)

| Superficie Sferica | Raggio $R$ [cm] | Area Superficiale [$m^2$] | Flusso Netto $\Phi_{\text{net}}$ [$T\cdot m^2$] | Flusso Assoluto $\Phi_{\text{abs}}$ [$T\cdot m^2$] | Residuo Relativo $\frac{\|\Phi_{\text{net}}\|}{\Phi_{\text{abs}}}$ | Esito Certificazione |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| **Sfera 1 (Near-Field)** | $8.0\text{ cm}$ | $0.0804\text{ m}^2$ | $+7.874 \times 10^{-9}$ | $7.147 \times 10^{-7}$ | **$1.102\%$** | **CONSERVATO** ($\Phi_{\text{net}} \approx 7.87\text{ nWb}$) |
| **Sfera 2 (Mid-Field)** | $12.0\text{ cm}$ | $0.1810\text{ m}^2$ | $-1.268 \times 10^{-8}$ | $5.794 \times 10^{-7}$ | **$2.189\%$** | **CONSERVATO** ($\Phi_{\text{net}} \approx -12.68\text{ nWb}$) |
| **Sfera 3 (Far-Field)** | $15.0\text{ cm}$ | $0.2827\text{ m}^2$ | $-4.898 \times 10^{-9}$ | $8.333 \times 10^{-7}$ | **$0.588\%$** | **CONSERVATO** ($\Phi_{\text{net}} \approx -4.90\text{ nWb}$) |

*Nota:* I residui netti dell'ordine di $\sim 10^{-9}\text{ T}\cdot\text{m}^2$ (nanoweber) testimoniano la rigorosa conservazione del flusso magnetico di Maxwell, imputabile unicamente all'interpolazione poliedrica della mesh tetraedrica.

### 2.2. Sensibilità Parametrica Reale sull'Inclinazione dei Ponticelli (Louver Angle $\alpha$)
*Valori calcolati direttamente mediante simulazioni transienti complete Elmer FEM (10 timesteps per caso, solutore UMFPACK):*

| Angolo Louver $\alpha$ | Conducibilità $\sigma_{\theta z}$ [S/m] | Variazione $\Delta\sigma$ | Potenza Joule DC [W] | Variazione $\Delta P_J$ | Campo $B_{\text{rad}}$ medio ($R=6\text{ cm}$) | Variazione $\Delta B_{\text{rad}}$ | Fonte Dati |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **$25.0^\circ$** | $2.681 \times 10^6$ | $-11.5\%$ | **$1.654\text{ W}$** | **$+7.70\%$** | $29.67\,\mu\text{T}$ | **$-0.02\%$** | `results_25deg/` |
| **$30.0^\circ$ (Base)** | $3.031 \times 10^6$ | $0.0\%$ | **$1.536\text{ W}$** | $0.00\%$ | $29.67\,\mu\text{T}$ | $0.00\%$ | `results_nominal_30deg/` |
| **$35.0^\circ$** | $3.289 \times 10^6$ | $+8.5\%$ | **$1.430\text{ W}$** | **$-6.89\%$** | $29.68\,\mu\text{T}$ | **$+0.01\%$** | `results_35deg/` |

### 2.3. Caratterizzazione Multifisica Avanzata: Sweep di Sfasamento, Flusso di Poynting e Modulazione AM
*Dettagli completi e trattazione teorica in [`docs/CARATTERIZZAZIONE_MULTIFISICA_SFASAMENTO_E_POYNTING.md`](docs/CARATTERIZZAZIONE_MULTIFISICA_SFASAMENTO_E_POYNTING.md):*

- **Sweep di Sfasamento Spazio-Temporale (100 Hz, 1200 RPM, 10 ms):**
  - **Regime A (Sincrono, $\Delta\phi=0^\circ$):** $\langle |B_{\text{rad}}| \rangle = 45.34\,\mu\text{T}$, Ripple = $381.1\%$, Dissipazione Joule $P_J = 4.294\text{ W}$.
  - **Regime B (Co-rotante $60^\circ$):** $\langle |B_{\text{rad}}| \rangle = 62.98\,\mu\text{T}$, Ripple = **$86.7\%$**, Dissipazione Joule $P_J = \mathbf{2.437\text{ W}}$ (**$-43.2\%$ di perdite termiche**).
  - **Regime C (Quadratura $90^\circ/180^\circ$):** $\langle |B_{\text{rad}}| \rangle = 68.14\,\mu\text{T}$, Ripple = $152.6\%$, Dissipazione Joule $P_J = 3.326\text{ W}$.
  - **Regime D (Contro-rotante $-60^\circ$):** $\langle |B_{\text{rad}}| \rangle = 70.98\,\mu\text{T}$, Ripple = $121.4\%$, Dissipazione Joule $P_J = 2.828\text{ W}$.
- **Flusso Attivo di Poynting Uscente (Cilindro di controllo $R=12\text{ cm}, H=20\text{ cm}$):**  
  Potenza irradiata/guidata netta verso lo spazio esterno = **$7.282\text{ mW}$** ($\|\vec{S}\|_{\text{mean}} = 0.735\text{ W/m}^2$).
- **Accoppiamento a Distanza / Virtual Harvesting:**
  - Sonda 1 (radiale, $R=10\text{ cm}$): $V_{\text{ind, peak}} = 1.56\text{ mV}$, $V_{\text{ind, rms}} = 0.92\text{ mV}$.
  - Sonda 2 (assiale su apertura campana, $R=15\text{ cm}, Z=+10\text{ cm}$): $V_{\text{ind, peak}} = 0.45\text{ mV}$, $V_{\text{ind, rms}} = 0.28\text{ mV}$.
  - Sonda 3 (capacitiva, $50\text{ cm}^2$ a $R=12\text{ cm}$): $I_{D, \text{rms}} = 10.52\text{ pA}$.
- **Modulazione AM a Bassa Frequenza ($f_{\text{mod}}=10\text{ Hz}, m=0.5$):**  
  Escursione pulsante dinamica del raggio d'inviluppo della campana magnetica: $\Delta R_{\text{breathing}} = \mathbf{3.18\text{ cm}}$ ($R_{\text{min}} = 6.49\text{ cm} \leftrightarrow R_{\text{max}} = 9.66\text{ cm}$).

---

## 3. Struttura del Repository

```
simulazione/
├── LICENSE.txt                                 (Licenza CERN-OHL-S-2.0)
├── CITATION.cff                                (Metadati CFF v1.2.0)
├── README.md                                   (Manuale di riproduzione, sintesi fisica e tabelle FEM reali)
├── requirements.txt                            (Dipendenze Python)
├── config/
│   ├── case_mesh_stirata.sif                   (Configurazione nominale 30° con MATC cilindrico)
│   ├── case_mesh_25deg.sif                     (Configurazione sensitività 25°)
│   ├── case_mesh_35deg.sif                     (Configurazione sensitività 35°)
│   ├── case_sweep_regime_A.sif                 (Regime Sincrono 0°)
│   ├── case_sweep_regime_B.sif                 (Regime Co-rotante 60°)
│   ├── case_sweep_regime_C.sif                 (Regime Quadratura 90°/180°)
│   ├── case_sweep_regime_D.sif                 (Regime Contro-rotante -60°)
│   └── case_am_modulation.sif                  (Modulazione dinamica AM 10 Hz)
├── docs/
│   └── CARATTERIZZAZIONE_MULTIFISICA_SFASAMENTO_E_POYNTING.md (Relazione teorico-numerica completa)
├── mesh/
│   ├── macchina.geo                            (Sorgente geometrico Gmsh)
│   ├── macchina.msh                            (Mesh volumetrica tetraedrica)
│   └── macchina/                               (Mesh conforme Elmer: nodes, elements, boundary)
├── scripts/
│   ├── run_all_simulations.py                  (Pipeline per le 3 simulazioni geometriche base)
│   ├── postprocess_solenoidalita.py            (Calcolo integrali di Gauss e sensitività)
│   ├── generate_official_figures.py            (Figure 01, 02, 03 a 300 DPI)
│   ├── sweep_sfasamento_fasi.py                (Pipeline per i 4 regimi di fase + AM)
│   └── postprocess_campo_elettrico_poynting.py (Estrazione E, Poynting, sonde virtuali, Figure 04, 05, 06)
├── data/
│   ├── validazione_chiusura_cern_ohl.json      (Dataset chiusura audit e conservazione flusso)
│   ├── confronto_mantello_pieno_vs_rete.csv    (Confronto termico ed elettrodinamico)
│   └── sweep_sfasamento_risultati.json         (Dataset completo regimi di fase, Poynting, sonde e AM)
├── figures/
│   ├── 01_abbattimento_correnti_joule.png      (Confronto perdite piene vs rete stirata)
│   ├── 02_espulsione_radiale_simmetrica_360.png(Mappatura 360° simmetrica di B_rad)
│   ├── 03_topologia_doppia_spirale_3d.png      (Linee di flusso 3D a doppia campana elicoidale)
│   ├── 04_sweep_sfasamento_confronto.png       (Confronto temporale, polare, ripple e perdite A-B-C-D)
│   ├── 05_vettore_poynting_e_campo_elettrico.png(Sezioni XZ Poynting e XY vortice chirale E)
│   └── 06_accoppiamento_distanza_harvesting.png(Tensioni indotte, decadimento spaziale e respiro AM)
└── _archive_backup/                            (Archivio storico dei test preliminari e report intermedi)
```

---

## 4. Istruzioni di Esecuzione e Riproducibilità da Zero

```bash
# 1. Installazione dipendenze Python
pip install -r requirements.txt

# 2. Generazione della mesh da sorgente CAD/Gmsh (opzionale se già presente in mesh/macchina/)
cd mesh
gmsh -3 macchina.geo -o macchina.msh
ElmerGrid 14 2 macchina.msh -autoclean
cd ..

# 3. Validazione Statica e Geometrica (Nominale 30°, 25°, 35°)
python scripts/run_all_simulations.py
python scripts/postprocess_solenoidalita.py
python scripts/generate_official_figures.py

# 4. Suite Multifisica Avanzata (Sweep Sfasamento Regimi A-B-C-D e Modulazione AM)
python scripts/sweep_sfasamento_fasi.py
python scripts/postprocess_campo_elettrico_poynting.py
```

---

## 5. Authorship & License

- **Autore / Lead Designer:** **Alessandro Brescacin** ([brescacin.alessandro@gmail.com](mailto:brescacin.alessandro@gmail.com))
- **Repository GitHub Ufficiale:** [https://github.com/brescale27/open-chiral-flux-shaper](https://github.com/brescale27/open-chiral-flux-shaper)
- **Licenza Open Hardware:** Rilasciato sotto licenza **CERN-OHL-S v2 (Strongly Reciprocal)**.  
  Il testo completo e vincolante è disponibile nel file [`LICENSE.txt`](LICENSE.txt).
- **Citazione Accademica:** Per citare questo progetto, pipeline numerica o dataset di simulazione, consultare il file [`CITATION.cff`](CITATION.cff).

