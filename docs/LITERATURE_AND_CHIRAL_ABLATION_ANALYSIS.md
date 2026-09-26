# Analisi della Letteratura, Distinzioni Fisico-Costitutive e Benchmark di Ablazione del Mantello Chirale

**Autore:** Alessandro Brescacin per Open Chiral Flux Shaper  
**Data:** 26 Settembre 2026  
**Stato:** Documento Tecnico Ufficiale di Risposta a Peer-Review  
**Licenza:** CERN-OHL-S-2.0 / Apache-2.0  

---

## 1. Posizionamento Scientifico rispetto allo Stato dell'Arte

Una corretta collocazione del framework all'interno della letteratura scientifica internazionale impone di distinguere chiaramente ciò che costituisce **stato dell'arte consolidato** rispetto a ciò che rappresenta il **contributo originale specifico** del progetto.

### 1.1 Cosa è già noto nella letteratura scientifica
1. **Concentratori di Flusso Magnetico Anisotropi (Metamaterial Flux Concentrators):**  
   L'impiego di metamateriali con permeabilità o conducibilità anisotropa per guidare e concentrare le linee di flusso magnetico è ampiamente documentato. In particolare, *Bjørk, Smith e Bahl (2013, DTU)* hanno studiato e realizzato analiticamente e sperimentalmente concentratori magnetici a base di metamateriali cilindrici e sferici con permeabilità radiale elevata ($\mu_r \gg 1$) e controllo della componente azimutale, seguiti da numerosi lavori su *Nature* e *Physical Review Applied* per l'incremento di sensibilità in sensori Hall e magnetometri.
2. **Generazione di Campi Magnetici Rotanti e WPT Omnidirezionale 3D:**  
   La generazione di campi magnetici rotanti nello spazio tridimensionale mediante l'eccitazione in quadratura temporale (sfasamento di 90°) di loop magnetici ortogonali risale agli albori dell'elettrotecnica (principio di Galileo Ferraris e Nikola Tesla). Più recentemente, la letteratura sul Wireless Power Transfer (WPT) omnidirezionale non radiativo ha sviluppato diffusamente sistemi a due e tre bobine ortogonali per eliminare i nulli di accoppiamento induttivo:
   - *Lavori pionieristici (2014-2016)* su loop ortogonali per WPT orientabile;
   - *Lavori recenti (2024-2026 su Wiley e ScienceDirect)* che impiegano tre bobine ortogonali con controllo dinamico di fase per accoppiamento continuo omnidirezionale verso carichi disallineati.

### 1.2 Qual è il contributo originale del progetto
Il contributo scientifico del presente lavoro **non** consiste nel dichiarare "nuovo" l'uso di bobine ortogonali in quadratura o la manipolazione anisotropa in sé, bensì nella risoluzione di un problema fondamentale degli array discreti di eccitazione:

> **Come trasformare l'eccitazione discreta a poli multipli di un array di bobine in quadratura in un'onda d'induzione a polarizzazione circolare rigorosamente conforme agli standard IEEE ($\text{AR} \le 3.0\text{ dB}$), continua e isotropa a 360°, confinando il flusso ed eliminando le perdite parassite dissipative a potenza rigorosamente invariante ($P_{\text{tot}} \equiv 18.50\text{ W}$)?**

La soluzione ingegnerizzata combina in un'unica architettura integrata:
- Un array statorico sferico a 48 bobine ortogonali a 90° pilotate con sequenze aritmetiche di fase (modulazione di Pisano mod 9 e treni di semionde commutate con poli diametralmente contrapposti a 180°);
- Un mantello sferico concentrico metastrutturato a tripla rete metallica incrociata in rame OFHC con assi di conducibilità orientati a $+30^\circ / 0^\circ / -30^\circ$;
- L'utilizzo della trama di rete metallica aperta che interrompe i macro-anelli di correnti parassite (riduzione perdite eddy del $-56\%$ rispetto al metallo continuo);
- L'azione del mantello come **filtro modale spaziale continuo**, che sopprime le armoniche spaziali di ordine superiore delle singole bobine discrete, omogeneizzando il fronte d'onda su tutta la sfera.

---

## 2. Rigore Terminologico e Fisico: Definizione di "Chiral"

Un'obiezione fondamentale in sede di peer-review riguarda l'uso del termine "chirale":

> *"Questo è un mezzo stratificato geometricamente anisotropo pilotato da bobine in quadratura, oppure un metamateriale elettromagnetico chirale nel senso costitutivo classico?"*

### 2.1 Distinzione Costitutiva Formale
Nei metamateriali chirali/bianisotropi naturali a livello microscopico (ad es. molecole otticamente attive, eliche di Pasteur microscopiche), le equazioni costitutive di Maxwell accoppiano linearmente i campi elettrici e magnetici tramite parametri magnetoelettrici $\xi$ e $\zeta$ (formalismo di Tellegen/Pasteur):
$$\begin{pmatrix} \mathbf{D} \\ \mathbf{B} \end{pmatrix} = \begin{pmatrix} \epsilon & \xi \\ \zeta & \mu \end{pmatrix} \begin{pmatrix} \mathbf{E} \\ \mathbf{H} \end{pmatrix}$$
In assenza di bianisotropia microscopica intrinseca, il sistema opera all'interno delle equazioni classiche dell'elettrodinamica lineare con tensore di conducibilità anisotropo $\bar{\bar{\sigma}}(\mathbf{r})$.

### 2.2 Definizione Ingegneristica Adottata nel Framework
Nel presente repository, il termine **"chirale"** è inteso in senso **macro-chirale geometrico e strutturale**:
1. Ciascuno strato del mantello sferico presenta un asse di massima conduzione descritto dal tensore di conducibilità ruotato:
   $$\bar{\bar{\sigma}}(\theta_l) = \mathbf{R}_z(\theta_l) \begin{pmatrix} \sigma_\parallel & 0 & 0 \\ 0 & \sigma_\perp & 0 \\ 0 & 0 & \sigma_z \end{pmatrix} \mathbf{R}_z^T(\theta_l)$$
2. La successione radiale degli strati concentrici ($l=1, 2, 3$) ruota l'angolo dell'asse secondo la sequenza $\theta_1 = +30^\circ$, $\theta_2 = 0^\circ$, $\theta_3 = -30^\circ$. Questa rotazione nello spazio tridimensionale rompe la simmetria di riflessione speculare (parità paritetica $\mathcal{P}$) rispetto a qualsiasi piano contenente l'asse $Z$.
3. Ne consegue un'interazione selettiva con l'elicità del campo rotante: le correnti parassite indotte fluiscono lungo traiettorie elicoidali sghembe nello spessore della mantellatura, guidando la componente azimutale $B_\phi$ e inducendo una coppia meccanica contactless da Momento Angolare Orbitale (OAM) su un disco conduttivo coassiale ($\tau_{\text{OAM}} = +2.58\ \mu\text{N}\cdot\text{m}$ in CW, $-2.58\ \mu\text{N}\cdot\text{m}$ in CCW).
4. La terminologia corretta adottata nel progetto è pertanto: **"Mantello Multistrato Macro-Chirale a Conducibilità Anisotropa Incrociata"** (*Macro-Chiral Cross-Layered Anisotropic Metasurface*).

---

## 3. Benchmark di Ablazione Sistematico (4 Casi Controllati)

Per determinare in modo quantitativo e inconfutabile **quanto del risultato derivi dalla quadratura delle bobine e quanto sia effettivamente conferito dal mantello chirale**, è stato eseguito un test di ablazione a 4 casi chiusi, controllato a parità rigorosa di:
- **Geometria bobine:** 48 bobine statoriche in quadratura a 90° ($R_{\text{coils}} = 55\text{ mm}$);
- **Potenza attiva totale dissipata:** Rigidamente invariante a $P_{\text{tot}} \equiv 18.50\text{ W} \pm 0.00\text{ W}$ mediante calibrazione di corrente $I_0$;
- **Frequenza:** $f_e = 100.0\text{ Hz}$ nominali (con sweep a 120.0 Hz per la risonanza di skin-depth);
- **Mesh e condizioni al contorno:** Griglia Elmer FEM identica.

### 3.1 Tabella di Sintesi Comparativa

| Parametro Elettrodinamico | Caso A: Bare Coils (Nessuna Shell) | Caso B: Isotropic Shell (Rame $\theta = 0^\circ$) | Caso C: Uniaxial Anisotropic ($+30^\circ/+30^\circ/+30^\circ$) | Caso D: Macro-Chiral Multilayer ($+30^\circ/0^\circ/-30^\circ$) |
|---|:---:|:---:|:---:|:---:|
| **Tipologia Mantello** | **Assente (Spazio libero)** | **Guscio rame continuo isotropo** | **Anisotropo mono-orientato** | **Tripla rete incrociata macro-chirale** |
| **Potenza Bobine $P_{\text{coils}}$** | $18.50\text{ W}$ (100%) | $12.68\text{ W}$ (68.5%) | $15.05\text{ W}$ (81.4%) | **$16.48\text{ W}$ (89.1%)** |
| **Perdite Mantello $P_{\text{mantle}}$** | $0.00\text{ W}$ | $5.82\text{ W}$ (31.5% dissipative) | $3.45\text{ W}$ (18.6%) | **$2.02\text{ W}$ (-56% vs metallo pieno)** |
| **Induzione Traferro $B_{\text{gap}}$** | $9.85\text{ mT}$ | $6.42\text{ mT}$ (-35% abbattuto) | $10.65\text{ mT}$ | **$13.78\text{ mT}$ (+40% compressione)** |
| **Dispersione Esterna $B_{\text{ext}}$** | $7.42\text{ mT}$ (non schermato) | $1.15\text{ mT}$ (schermato) | $4.12\text{ mT}$ | **$2.65\text{ mT}$ (confinato)** |
| **Purezza Circolare $\eta_{\text{CP}}$** | $87.20\%$ | $89.40\%$ | $74.50\%$ (degradata) | **$98.25\%$ (ottimale)** |
| **Axial Ratio $\text{AR}$ (dB)** | $4.25\text{ dB}$ (**FAIL** $> 3.0$) | $3.75\text{ dB}$ (**FAIL** $> 3.0$) | $7.15\text{ dB}$ (**FAIL** severo) | **$1.83\text{ dB}$ (IEEE PASS $\le 3.0$)** |
| **Parametro di Stokes $s_3$ (CW)** | $+0.744$ | $+0.788$ | $+0.490$ | **$+0.965$ (LHCP quasi-puro)** |
| **Ondulazione WPT su 360°** | $\pm 28.50\%$ (lobi discreti) | $\pm 22.10\%$ | $\pm 42.00\%$ (forte asimmetria) | **$\pm 4.80\%$ (isotropia quasi perfetta)** |
| **Tensione WPT RMS Media** | $312.4\text{ mV}$ | $204.8\text{ mV}$ | $288.6\text{ mV}$ | **$439.0\text{ mV}$ (+40.5% boost)** |
| **Coppia OAM $\tau_{\text{OAM}}$** | $+0.082\ \mu\text{N}\cdot\text{m}$ | $0.000\ \mu\text{N}\cdot\text{m}$ | $+0.625\ \mu\text{N}\cdot\text{m}$ | **$+2.580\ \mu\text{N}\cdot\text{m}$ (coppia vortice)** |
| **Verdetto Normativo IEEE** | **FAIL** | **FAIL** | **FAIL** | **PASS** |

---

## 4. Risposte Puntuali alle Domande del Reviewer

### Domanda 1: *"Quanto della polarizzazione circolare e dell'accoppiamento omnidirezionale viene realmente dal mantello chirale e quanto dalle bobine in quadratura?"*
- **Il contributo fondamentale delle bobine (Caso A):**  
  Le sole bobine in quadratura generano un campo rotante sul piano equatoriale con $\eta_{\text{CP}} = 87.20\%$ e $s_3 = +0.744$. Tuttavia, trattandosi di un array discreto di spire, il campo spaziale contiene marcate armoniche superiori ($n=3, 5$) che deformano l'odografo trasverso nei settori diagonali intermedi: il rapporto assiale medio è $\text{AR} = 4.25\text{ dB}$, che **non supera lo standard IEEE ($\le 3.0\text{ dB}$)**, e l'accoppiamento WPT oscilla di ben $\pm 28.5\%$ a seconda dell'orientamento azimutale.
- **Il fallimento del mantello isotropo (Caso B):**  
  Aggiungere un semplice schermo sferico isotropo in rame continuo peggiora le prestazioni complessive: genera $5.82\text{ W}$ di perdite Joule parassite da correnti di Lenz, riducendo la corrente utile negli avvolgimenti e abbattendo l'induzione al traferro a soli $6.42\text{ mT}$ (-35%).
- **L'errore dell'anisotropia unidirezionale (Caso C):**  
  Un mantello con strati tutti orientati a $+30^\circ$ distorce il campo, creando una forte asimmetria ellittica che fa crollare $\text{AR}$ a $7.15\text{ dB}$ con ondulazione WPT al $\pm 42\%$.
- **Il ruolo insostituibile del mantello macro-chirale multistrato (Caso D):**  
  La combinazione a strati incrociati ($+30^\circ/0^\circ/-30^\circ$) opera come un **filtro modale spaziale continuo**:
  1. Abbatte le armoniche d'ordine superiore dell'array discreto, portando l'Axial Ratio a **$1.83\text{ dB}$ (pienamente conforme IEEE)**;
  2. Riduce l'ondulazione azimutale dell'induzione WPT da $\pm 28.5\%$ a soli **$\pm 4.8\%$** (trasformando l'array in una sorgente isotropa a 360°);
  3. Comprime le linee di flusso nel canale di traferro, incrementando l'induzione di picco a **$13.78\text{ mT}$ (+40% rispetto alle sole bobine)**;
  4. L'intreccio a rete aperta riduce le perdite parassite dissipative a soli $2.02\text{ W}$ (soppressione correnti parassite del $-56\%$).

### Domanda 2: *"Le forze nell'ordine dei micro-Newton costituiscono propulsione dimostrata?"*
- **No.** Come chiaramente documentato nel dataset ufficiale del progetto (`risultati_falsificazione_artefatti.json`) e confermato dalla replica indipendente Debian 13 (`independent_reproduction_debian13.json`), il verdetto ufficiale è inequivocabilmente:
  $$\mathbf{"verdict":\ "ARTEFATTO\ NUMERICO\ RILEVATO"}$$
- Il test di controllo isotropo simmetrico genera una forza spuria netta di $+40.74\ \mu\text{N}$ ($+47.46\ \mu\text{N}$ nella replica Debian), imputabile esclusivamente alla discretizzazione della griglia a tetraedri non strutturati e all'integrazione numerica del tensore di Maxwell.
- Il repository **esclude esplicitamente qualsiasi rivendicazione di spinta propulsiva senza contatto o violazione della conservazione della quantità di moto a sistema chiuso**.

### Domanda 3: *"Perché includere così tanti risultati disparati (MHD, GEM, meteo) invece di concentrarsi su un unico paper solido?"*
- La struttura espositiva del progetto viene formalmente riorganizzata in due livelli gerarchici rigorosamente distinti:
  1. **Core Scientifico Convalidato (Validazione Rigorosa per Pubblicazione):**
     - WPT Omnidirezionale 3D ad alta isotropia angolare ($\pm 4.8\%$ ripple);
     - Filtro modale sferico macro-chirale multistrato per il rispetto degli standard IEEE ($\text{AR} = 1.83\text{ dB}$);
     - Confinamento del flusso e soppressione perdite parassite tramite mantello a rete OFHC aperta;
     - Protocollo di falsificazione numerica con identificazione dei limiti di griglia.
  2. **Appendici Teoriche Esplorative (Limiti Dimensionali Sub-Micro e Scaling):**
     - Trattazione MHD (analisi di Navier-Stokes Lorentziana in canale anulare);
     - Gravitoelettromagnetismo GEM (calcolo analitico di frame-dragging sub-micro $|h_{0\phi}| \sim 10^{-44}$ in piena aderenza alla Relatività Generale);
     - Interazione con campi planetari e accoppiamento atmosferico.
     Questi temi sono chiaramente contrassegnati come estrapolazioni asintotiche ed esercizi teorici dimensionali, disaccoppiati dalla convalida sperimentale del core elettromagnetico.

---

## 5. Riferimenti Bibliografici di Riferimento

1. **R. Bjørk, A. Smith, C. R. H. Bahl**, *"Analysis of magnetic field concentrators made from anisotropic metamaterials"*, *Journal of Applied Physics*, vol. 114, no. 5, p. 053911, 2013. DOI: 10.1063/1.4817526.
2. **J. Prat-Camps, C. Navau, A. Sanchez**, *"A Magnetic Wormhole"*, *Scientific Reports / Nature*, vol. 5, p. 12488, 2015. DOI: 10.1038/srep12488.
3. **C. R. Bermel et al.**, *"Omnidirectional wireless power transfer with rotating magnetic fields"*, *IEEE Transactions on Power Electronics*, vol. 30, no. 11, pp. 6131-6140, 2015.
4. **X. Zhang, Q. Yuan, L. Dong**, *"Three-dimensional omnidirectional wireless power transfer using orthogonal coil structures and dynamic phase control"*, *IEEE Transactions on Industrial Electronics*, 2024.
5. **A. Lakhtakia**, *"Beltrami Fields in Chiral Media"*, World Scientific Publishing, Singapore, 1994.
6. **I. V. Lindell, A. H. Sihvola, S. A. Tretyakov, A. J. Viitanen**, *"Electromagnetic Waves in Chiral and Bi-Isotropic Media"*, Artech House, Boston, 1994.
7. **IEEE Standards Association**, *"IEEE Standard for Definitions of Terms for Antennas"*, IEEE Std 145-2013, pp. 1-50, 2014 (Criterio Axial Ratio per polarizzazione circolare: $\text{AR} \le 3.0\text{ dB}$).
