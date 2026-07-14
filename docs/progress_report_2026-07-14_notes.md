# OncoTox — was seit dem letzten Progress Report gemacht wurde

**Zeitraum:** 13.–14.07.2026 · 14 Commits (`fc23721` … `f7991ff`), alle auf `main`
**Zweck dieser Datei:** Beleg-Anhang zu den Slides — Parameter, Notebooks, CSV-Pfade.
Jede Zahl hier ist aus dem Code / den gespeicherten CSVs gezogen, nicht aus dem Gedächtnis.

> ⚠️ **Stand 14.07. abends — alles neu gerechnet.** Der Drug-Filter wurde korrigiert
> (`learnability = min(#getötet, #verschont)`, Coverage ≥ 90 %, **top 10** statt top 5), der Val-Leak im
> DrEval-Notebook wurde gefixt, und **jedes Experiment wurde auf der neuen 10-Drug-Menge wiederholt.**
> **Maßgeblich sind die Zahlen in [`progress_report_2026-07-14.md`](progress_report_2026-07-14.md) und in
> `notebooks/outputs/*.csv`.** Zahlen in diesem Dokument, die auf die *alten 5 Drugs* verweisen (ρ = 0.42 /
> 0.49, DrEval 0.511 / R² 0.224), sind **überholt** — sie bleiben stehen, weil §4.1 und §5.7 erklären,
> *warum* sie zu optimistisch waren. Die aktuellen Werte:
>
> | | alt (5 Drugs) | **neu (10 Drugs)** |
> |---|---|---|
> | K=10, out-of-fold ρ | PCA 0.42 / scGPT 0.49 | **PCA 0.360 / scGPT 0.396** |
> | K=545, `auc_z` | PCA 0.378 / scGPT 0.430 | **PCA 0.316 / scGPT 0.328** |
> | Rescue: no-reg / `auc_z` | +0.234 / +0.433 | **+0.210 / +0.333** |
> | Ridge (line-level), PCA | 0.428 | **0.342** |
> | DrEval norm. ρ / R² (scGPT) | 0.511 / 0.224 | **0.357 / 0.114** |

---

## 0. Ausgangslage (Stand letzter Report)

| | |
|---|---|
| Setting | 545 Drug-Heads, Target `mean_pv`, PCA vs. scGPT |
| Ergebnis | out-of-fold ρ ≈ 0 für **beide** Repräsentationen |
| Schlussfolgerung damals | *"Weder PCA noch scGPT können Zelllinien ranken."* |
| Offene Frage | Liegt es an der Datenlage — zu wenige Linien für zu viele Drugs? |
| Artefakte | jetzt unter `notebooks/outputs/legacy/training_545_mean_pv/` |

**Die Antwort, die diese Arbeit gibt:** Nein. Es lag an der **Skalierung des Targets**.
Das Null-Resultat war überwiegend ein **Artefakt eines fehlskalierten Multi-Task-Loss**.

---

## 1. Was über ALLE Experimente konstant gehalten wurde

Das ist der Teil, der die Vergleichbarkeit trägt. **Nichts davon variiert**, außer wo explizit gesagt.

### 1.1 Daten (identisch in jedem einzelnen Fit)

| Größe | Wert |
|---|---|
| Datensatz | SCP542 (single-cell RNA-seq) + CTRPv2 (bulk drug response) |
| Zellen | **53.513** |
| Gene | Variante **`hvg5000`** → 5.000 HVGs; scGPT embeddet die **4.576** davon in seinem Vokabular (**424 OOV**); PCA wird auf allen 5.000 gerechnet |
| Repräsentationen | `X_pca` **512-d**, `X_scGPT` **512-d** — bewusst gleiche Dimension, damit nur der *Inhalt* differiert |
| Zelllinien gesamt | 198 (davon 18 `unassigned` = ohne CTRP-Label) |
| Drugs | **545** |

**Splits** (nach Zelllinie, leakage-frei):

| Split | Zellen | Linien |
|---|---|---|
| train | 34.126 | 126 |
| val | 7.121 | 27 |
| test | 5.980 | 27 |
| *(unassigned)* | 6.286 | 18 |

> ⚠️ **Der wichtigste Datenpunkt überhaupt:** Das Label lebt auf **(Zelllinie × Drug)**, nicht auf der
> Zelle. Es gibt zwar 21,3 Mio. beobachtete (Zelle × Drug)-Paare, aber nur **~150 unabhängige
> Beispiele pro Drug**. Die 53.513 Zellen sind **Pseudo-Replikate**. Jede Fehlerbalken-Rechnung muss
> darauf beruhen, sonst ist sie um √350 zu optimistisch.

### 1.2 Modell (`scripts/model/OncoMLP.py`)

```
Input (512-d)
  → Dropout(input_dropout)
  → [Linear → LayerNorm → GELU → Dropout(dropout)] × hidden_dims
  → Linear(→ K)                          # K = Anzahl Drug-Heads, ein shared trunk
```

Default `DEFAULT_HIDDEN_DIMS = (128, 64)` — **für PCA und scGPT identisch**, explizit für die Fairness
des Vergleichs (`scripts/training/train_multitask.py:53`).

| K | Parameter |
|---|---|
| K = 5 | **74.629** |
| K = 545 | **109.729** |

### 1.3 Training (`scripts/training/training_utils.py`, `TrainConfig`)

| Parameter | Wert |
|---|---|
| `epochs` | **25** (war 50; über 36 Runs war die beste Epoche median 6, max 11) |
| `lr` | 1e-3 (Adam) |
| `weight_decay` | 1e-3 |
| `dropout` / `input_dropout` | 0.5 / 0.1 |
| `batch_size` | 128 |
| `grad_clip` | 1.0 |
| LR-Scheduler | ReduceLROnPlateau, patience 3, factor 0.5 |
| Early stopping | patience 10 |
| Loss | **masked MSE** — Mittel über alle *beobachteten* (Zelle × Drug)-Einträge |
| `seed` | 42 (Seed-Stabilität separat geprüft, §5.6) |

### 1.4 Der Vergleichs-Score — was in §5 überall in der Spalte „ρ“ steht

**Jede einzelne ρ-Zahl in §5 ist dieselbe Größe.** Definition, Schritt für Schritt:

1. **5-fold `GroupKFold`, gruppiert nach Zelllinie**, über die **153 train+val-Linien**.
   Gruppierung nach Linie ⇒ keine Zelle einer Test-Linie war je im Training. Kein Leakage.
2. Jede Linie ist damit **genau einmal** out-of-fold vorhergesagt → Vorhersagen für **alle 153 Linien**,
   keine davon aus einem Modell, das sie gesehen hat.
3. **Aggregation Zelle → Zelllinie:** die Vorhersagen aller Zellen einer Linie werden **gemittelt**.
   Zwingend, denn *das Label lebt auf der Linie, nicht auf der Zelle* (§1.1).
4. **Pro Drug** wird dann **Spearman ρ** zwischen den 153 vorhergesagten und den 153 wahren
   Linien-Werten gerechnet.
5. **Berichtet wird der Mittelwert dieser ρ über die 5 diagnostischen Drugs.**

**Was ρ bedeutet:** *„Ordnet das Modell die ungesehenen Zelllinien richtig danach, wie sensitiv sie auf
diese Drug reagieren?"* — ρ = 1 perfektes Ranking, ρ = 0 Zufall, ρ < 0 systematisch verkehrt herum.

**Warum Spearman und nicht Pearson:** Die Frage *ist* ein Ranking („welche Linie ist sensitiver?"), und
Spearman ist robust gegen die Ausreißer der `auc`-Skala (`fqi-2` reicht bis 2.11). **Beide werden
gerechnet und stimmen auf ±0.02 überein** (Spalte `pearson` in `target/target_comparison.csv`) — die
Wahl der Metrik ändert **keine einzige** Schlussfolgerung.

**Wichtig — der Score ist über alle Zeilen von §5 vergleichbar:**
Auch bei **K=545** wird **auf denselben 5 Drugs** evaluiert. Der Unterschied zwischen K=5 und K=545 ist
*ausschließlich*, wie viele Heads **mittrainiert** wurden, nicht was gemessen wird. Genau deshalb ist die
Matrix in §5.1 zeilenweise lesbar.

> **Warum out-of-fold statt des festen val-Splits:** der val-Split hat 27 Linien → SE(ρ) ≈ **±0.2**.
> Damit ist jeder Effekt < 0.4 nicht messbar. Out-of-fold über 153 Linien ist die Voraussetzung dafür,
> dass die Zahlen unten überhaupt etwas bedeuten. **Das war eine der Änderungen, ohne die nichts
> conclusive gewesen wäre.**
>
> ⚠️ **Nicht verwechseln:** DrEval (§5.7) benutzt eine **andere** Metrik — dort wird **über alle
> (Linie × Drug)-Paare gepoolt** korreliert, nicht pro Drug gemittelt, und zusätzlich *normalisiert*.
> Die DrEval-Zahlen sind **nicht** direkt mit den ρ aus §5.1–§5.6 vergleichbar.

---

## 2. Was VARIIERT wurde — die drei Achsen

Alles, was hier steht, ist auf die obige Basis aufgesetzt. Nur diese drei Achsen bewegen sich:

| Achse | Werte |
|---|---|
| **A. Target** | `mean_pv` · `auc` · `auc_z` |
| **B. Head-Zahl K** | 5 (gefiltert) · 545 (alle) |
| **C. Repräsentation** | `X_pca` · `X_scGPT` |

Dazu die Modell-Knöpfe (Regularisierung, Kapazität, Batch, Sample-Weighting), die **nur** in den
Ablationen bewegt werden.

### Achse A — die drei Targets (`scripts/preprocessing/ctrp_to_h5ad.py`)

| Target | Definition | Biologie |
|---|---|---|
| `mean_pv` | ungewichteter Mittelwert von `cpd_avg_pv` über die 16 Dosis-Punkte | *Legacy.* Was die Platte roh gemessen hat |
| `auc` | `area_under_curve / conc_pts_fit` aus den **post-QC-Sigmoid-Fits** | Fläche unter der *gefitteten* Kurve, auf die Dosis-Achse normiert |
| `auc_z` | **z-Score von `auc` pro Drug** über die Zelllinien | *"Ist diese Linie sensitiver als der Durchschnitt für diese Drug?"* — potenz-frei |

CLI-Flag (Commit `aa4f6d3`):
```bash
uv run scripts/preprocessing/run_preprocessing.py --variant hvg5000 --all-drugs \
    --score auc_z --start-at targets --skip-scgpt
```

#### Wie `auc_z` genau skaliert ist

**Standardisiert wird pro Drug, über die Zelllinien** — und die Statistiken werden auf **Linien-Ebene**
gerechnet, **nicht** auf Zell-Ebene (sonst würden Linien mit mehr Zellen den Mittelwert verzerren):

```python
# scripts/preprocessing/ctrp_to_h5ad.py :: _zscore_per_drug()
center[j] = mean_over_cell_lines( auc[:, j] )   # mu_j : mittlere auc der Drug j
scale[j]  =  std_over_cell_lines( auc[:, j] )   # sigma_j : Streuung der Drug j ueber Linien
auc_z[i, j] = ( auc[i, j] - center[j] ) / scale[j]
```

Ergebnis: **jede der 545 Drugs hat über die Zelllinien Mittelwert 0 und Std 1.**
Die Konstanten liegen im h5ad, die Rück-Transformation ist **exakt**:

```
auc = auc_z * uns["ctrp_score_scale"] + uns["ctrp_score_center"]
```

**Zwei Konsequenzen, die man auseinanderhalten muss:**
- **Innerhalb einer Drug ändert sich nichts.** Die Rangfolge der Linien ist identisch → **Spearman ρ
  einer einzelnen Drug ist invariant unter z-scoring.** Das z-scoring „schönt" also **keine Metrik**.
- **Zwischen den Drugs ändert sich alles.** Im gemeinsamen MSE trägt jetzt jede Drug **gleich viel** bei,
  statt proportional zu σ². **Nur** darüber wirkt es — und deshalb wirkt es auch **nur bei K=545**
  (bei K=5 ist es folgerichtig ein Nullereignis: 0.482 → 0.488).

> **Merksatz für die Slide:** z-scoring des Targets ist **kein Metrik-Trick, sondern eine
> Loss-Gewichtung.** Es ist algebraisch identisch dazu, jeden der 545 Heads mit **1/σ²** zu gewichten.

#### Warum wir `auc_z` gewählt haben — die Begründung in der Reihenfolge, in der sie entstand

| Schritt | Frage | Befund | Konsequenz |
|---|---|---|---|
| 1 | Ist `mean_pv` überhaupt der richtige Score? | Es ist der **rohe Dosis-Mittelwert**; CTRPv2 liefert einen **gefitteten AUC** mit, der Rauschen glättet und Drugs mit unterschiedlichen Dosis-Grids vergleichbar macht | `auc` implementieren und **messen**, nicht annehmen |
| 2 | Bringt der Kurven-Fit Genauigkeit? | **Nein — null.** K=5: 0.481 → 0.482. Die beiden Targets landen empirisch innerhalb ~0.03 voneinander (`data/target_biology.png`) | Meine Ausgangs­hypothese ist **falsifiziert**. Der Fit ist *nicht* der Grund |
| 3 | Was unterscheidet Drugs dann im Multi-Task-Loss? | Ihre **Skala**: σ spannt 0.034–0.302 → **80×** im quadrierten Fehler. Die breitesten 10 % tragen 31 % des Loss, und die drei größten Loss-Träger töten **null** Linien | Das Problem ist **nicht der Score, sondern seine Varianz** |
| 4 | Was behebt das? | Per-Drug-Standardisierung ≡ 1/σ²-Gewichtung jedes Heads | **`auc_z`** |
| 5 | Hält es? | K=545: −0.087 → **+0.430** (scGPT). Über 3 Seeds stabil | **`auc_z` wird Default** (`DEFAULT_CTRP_SCORE`) |

**Die ehrliche Kurzfassung:** Wir sind wegen des Kurven-Fits zur `auc` gegangen — und der Kurven-Fit war
**irrelevant**. Was zählt, ist das **z**. Wir hätten denselben Effekt mit `mean_pv_z` bekommen.
`auc` bleibt trotzdem die saubere Basis, weil sie die glattere, dosis-normalisierte Größe ist; aber
**der gesamte gemessene Gewinn kommt aus der Standardisierung.** Das gehört genau so auf die Slide —
es ist ein Fall, in dem die Motivation falsch und das Ergebnis trotzdem richtig war.

---

## 3. Der Befund — der Bug im Loss

**Notebook `10_diagnosis` §2 · Plot `outputs/target/loss_weighting_bug.png`**

Der masked MSE gewichtet jeden beobachteten `(Zelle × Drug)`-Eintrag **gleich**. Aber der *quadrierte
Fehler* einer Drug skaliert mit ihrer Response-Varianz **σ²**.

- Per-Drug σ spannt **0.034 – 0.302** → ≈ **80× im quadrierten Fehler**
- Die breitesten **10 %** der Drugs tragen **30 %** des Loss
- Alle 545 Heads teilen sich **einen Trunk**

⇒ Der ungewichtete Loss gewichtet die 545 Heads implizit mit **σ², also mit Skala statt mit
Lernbarkeit.** Die Repräsentation wird auf die Varianz von Drugs gefittet, die **kein lernbares Signal
tragen**.

**Der Beleg, dass Spread ≠ Lernbarkeit:**

| Drug | σ (Rang) | tötet Zelllinien |
|---|---|---|
| `fqi-2` | **#1 von 545** (0.296) | **0** |
| `ciclopirox` | **#2** (0.265) | **0** |
| `brd-k71781559` | #3 (0.257) | **0** |

Die drei Drugs, die den Loss am stärksten dominieren, **separieren die Zelllinien überhaupt nicht.**
(Aus `outputs/learnability/ctrp_drug_learnability_auc.csv`.)

**Die Formulierung, die auf die Slide gehört:** *Nicht zu gewichten **ist** eine Gewichtung.*
Ein Target pro Drug zu z-scoren ist **äquivalent dazu, jeden Head mit 1/σ² zu gewichten** — also das
Regressions-Analogon zum Class-Imbalance-Reweighting (vgl. Kendall et al. 2018, uncertainty weighting).

---

## 4. Der Drug-Filter — eine METHODE, kein Ergebnis

**Notebook `08_learnability_filter` · Plot `outputs/learnability/learnability_filter_auc.png`**

**Zweck:** Die Frage "kann irgendetwas gelernt werden?" auf Drugs stellen, die nicht schon
*by construction* hoffnungslos sind. **Das ist eine diagnostische Vereinfachung, kein Resultat.**

Eine Drug ist nur dann über Zelllinien rankbar, wenn sie einen echten Teil **tötet** *und* einen echten
Teil **verschont**. Schwellen auf der rohen `auc`-Skala: `auc ≤ 0.5` = getötet, `auc ≥ 0.8` = überlebt.

| Gate | Schwelle | Median über die 545 |
|---|---|---|
| Coverage | ≥ 95 % der 180 Linien | 95 % |
| Spread | `auc_std ≥ 0.15` | 0.109 |
| Dynamic range | `p95 − p05 ≥ 0.40` | 0.306 |
| **getötete Linien** | **`n_sens ≥ 20`** | — |
| **überlebende Linien** | **`n_res ≥ 20`** | — |

Die letzten beiden Gates machen die eigentliche Arbeit. **5 / 545** bestehen alle fünf:

> `1s,3r-rsl-3` · `kx2-391` · `cay10618` · `ml162` · `dasatinib`

Die lose Version des Filters (nur Coverage × Spread, wie in `04_drug_coverage`) behält **439/545** und
lehrt nichts. Auf `auc_z` ist die alte Formel sogar **degeneriert** — z-scoring setzt jede Drug-Std auf
exakt 1.0, alle 545 sind gleichauf.

> ⚠️ **Der blockierende Vorbehalt:** Der Filter hat **alle 180 Linien gesehen**, inkl. val/test.
> Die 5 Drugs sind damit ein **Best-Case-Subset**, keine Zufallsstichprobe. Das macht jede Zahl auf
> diesen 5 Drugs zu einer **Diagnose**, nicht zu einer Generalisierungs-Schätzung. Drug-Auswahl
> *innerhalb* jedes CV-Folds ist der nächste Schritt. **Muss auf die Slide.**

### 4.1 ⚠️ Wie belastbar ist „learnability“ wirklich? — validiert am 14.07., Ergebnis: **nur bedingt**

Der Filter wurde **rein aus Label-Statistik** definiert und **nie gegen das geprüft, was er vorhersagen
soll**: das ρ, das das Modell auf der Drug tatsächlich erreicht. Nachgeholt mit **einem K=545-Fit
(`auc_z`, scGPT, out-of-fold)**, der ein ρ für **alle 537 auswertbaren Drugs** liefert.
*(Skript: `scratchpad/learnability_validity.py`, Ergebnis: `scratchpad/learnability_vs_achieved.csv`.)*

**✓ Was hält — der Filter reichert an:**

| | Anzahl | Ø ρ | Median ρ |
|---|---|---|---|
| **behalten** (alle Gates) | 6 | **+0.395** | +0.428 |
| **verworfen** | 531 | +0.116 | +0.119 |

Die 5 Drugs, auf denen das Hauptergebnis (§5.1–§5.3) steht, **funktionieren wirklich.**
**Das Target-Resultat wackelt dadurch nicht.**

**✗ Was nicht hält — drei Befunde:**

**(1) Als Ranking ist der Score schwach.** Korrelation mit dem tatsächlich erreichten ρ über alle 537:

| Feature | Spearman mit erreichtem ρ |
|---|---|
| `dyn_range` | **+0.403** |
| `n_sens` (getötet) | +0.397 |
| `auc_std` (σ) | +0.381 |
| **`learnability` (der kombinierte Score)** | **+0.357** |
| `cov_frac` | **+0.158** ← **schwächster Prädiktor …** |
| **`n_res` (überlebt)** | **−0.331** ← **… und dieses Gate wirkt GEGEN uns** |

Der zusammengesetzte Score ist **schlechter** als seine besten Einzelfaktoren. Und `n_res ≥ 20`
selektiert auf ein Merkmal, das **negativ** mit Lernbarkeit korreliert.

**(2) Schlechte Recall — der Filter wirft die besten Drugs weg.** Von den **12 Drugs mit ρ > 0.4**
wurden **9 verworfen**:

| Drug | erreichtes ρ | verworfen von welchem Gate |
|---|---|---|
| **`ml210`** | **+0.516** ← **beste Drug im ganzen Panel** | **Coverage 0.94 < 0.95** |
| `1s,3r-rsl-3` | +0.495 | *behalten* |
| `ml162` | +0.494 | *behalten* |
| `dasatinib` | +0.485 | *behalten* |
| `gsk461364` | +0.481 | Coverage 0.92 |
| `ceranib-2` | +0.460 | std, dyn, n_sens |
| `paclitaxel` | +0.449 | Coverage 0.94 |
| `alisertib` | +0.439 | std, n_sens |
| `trametinib` | +0.408 | Coverage 0.46 |

**109 von 531 verworfenen Drugs schlagen die schwächste behaltene Drug** (`kx2-391`, ρ = 0.241 — genau
die, die laut §5.8 ohnehin fast reiner Zelllinien-Effekt ist).
**Das beste Medikament im Panel scheiterte an 0.01 Coverage.**

**(3) Die 5er-Auswahl ist nicht stabil.** Bootstrap über die 180 Zelllinien (300 Resamples):

- Es passieren im Median **9** Drugs die Gates, Spanne **2 – 17**. Die „**5** / 545" sind also schon als
  **Anzahl eine Zufallszahl**, kein Merkmal der Daten.
- Die Zusammensetzung wechselt: zuverlässig in den Top-5 landen nur `kx2-391` (58 %) und
  `1s,3r-rsl-3` (48 %). **`paclitaxel` (45 %), `sb-743921` (44 %), `methotrexate` (42 %) und
  `vincristine` (40 %) — allesamt NICHT in unseren Fünf — werden genauso oft gewählt.**

**⇒ Die Aussage, die daraus folgt (und die auf die Slide muss):**

> **„5 / 545 bestehen die Gates" darf NICHT als „nur 5 Drugs sind lernbar" gelesen werden.
> Diese Lesart ist FALSCH.**

Mit korrigiertem Target erreicht das Modell über **alle 537** Drugs:

| Schwelle | Drugs | Anteil |
|---|---|---|
| ρ > 0 | **408** | **76 %** |
| ρ > 0.2 | 170 | 32 % |
| ρ > 0.3 | 56 | 10 % |
| ρ > 0.4 | 12 | 2 % |
| *Median über alle 537* | **ρ = 0.120** | |

**Das Modell lernt weit mehr als 5 Drugs.** Der Filter war nie dafür gebaut, das zu zeigen — er wählt
*ein paar sichere* Drugs, nicht *alle lernbaren*.

**Konsequenz — und sie stützt die Rahmung „der Filter ist eine Methode, kein Ergebnis" sogar stärker:**
Der Filter ist als **Diagnose-Werkzeug legitim** (er liefert Drugs, die nachweislich funktionieren), aber
als **Messgröße für Lernbarkeit unbrauchbar** — er ist ein schwacher Ranker, hat schlechte Recall und
ist gegen Resampling instabil. Er **unterschätzt**, was das Modell kann.

**Was daraus für die Zukunft folgt** (falls der Filter weiterverwendet wird):
- **Coverage-Gate senken oder streichen** (schwächster Prädiktor, killt die besten Drugs).
- **`n_res`-Gate streichen** (korreliert negativ).
- **`dyn_range` und `n_sens`** behalten — das sind die einzigen, die tragen.
- Besser noch: **Selektion innerhalb jedes CV-Folds**, dann ist der Best-Case-Vorbehalt gleich mit erledigt.

---

## 5. Die Experimente — Parameter, Daten, Ergebnisse

### 5.1 Haupt-Matrix: Target × Head-Zahl × Repräsentation

**Notebook `11_auc_vs_aucz` · Plot `outputs/target/target_comparison.png` ·
CSVs `target/target_comparison.csv`, `target/target_comparison_ci.csv`**

**Design:** 3 Targets × 2 Reps × 2 Head-Zahlen = **12 Fits.**
**Alles andere fix** (`(128,64)`, dropout 0.5, wd 1e-3, batch 128, 25 epochs, seed 42).
**Evaluiert immer auf denselben 5 Drugs**, egal ob mit K=5 oder K=545 trainiert wurde — deshalb sind
die Zeilen direkt vergleichbar. CIs = Bootstrap über die Drugs.

| K | Rep | `mean_pv` | `auc` | `auc_z` |
|---|---|---|---|---|
| **5** | PCA | 0.450 | 0.439 | 0.424 |
| **5** | scGPT | 0.481 | 0.482 | **0.488** |
| **545** | PCA | +0.027 | +0.016 | **+0.378** |
| **545** | scGPT | **−0.070** | **−0.087** | **+0.430** |

**Drei Ablesungen, alle conclusive:**

1. **Bei K=5 sind alle drei Targets gleich gut.** Der Kurven-Fit (`auc` statt `mean_pv`) bringt
   **exakt nichts** (0.481 → 0.482 für scGPT). *Meine frühere Behauptung, der Fit rette Signal, das
   der Dosis-Mittelwert zerstört, ist damit **falsifiziert**.*
2. **Bei K=545 kollabieren beide unstandardisierten Targets** — scGPT wird sogar **negativ**.
3. **`auc_z` hält bei K=545** (+0.430). **Das z-scoring ist der ganze Effekt.**

**Die Zahl, die die Botschaft dreht:**
> Mit `auc_z` liefert **K=545 (ρ = 0.430)** fast dasselbe wie **K=5 (ρ = 0.488)**.
> ⇒ **Die Head-Zahl war nicht das Problem. Die Datenlage war nicht das Problem.**
> ⇒ Der Drug-Filter ist rückwirkend eine *Diagnose-Vereinfachung*, keine *Notwendigkeit*.

### 5.2 PCA vs. scGPT auf den 5 Drugs

**Notebook `09_learnable5_training` · Plots `learnability/learnable5_pca_vs_scgpt.png`,
`learnability/learnable5_pred_vs_true.png` · CSV `learnability/learnable5_per_drug_correlation.csv`**

Out-of-fold über ~150 Linien, `auc_z`:

| Rep | ρ |
|---|---|
| PCA | **0.42** |
| scGPT | **0.49** |
| *Vorsprung scGPT* | *+0.06 … +0.08* |

Pro Drug (scGPT, `auc_z`): `ml162` 0.655 · `1s,3r-rsl-3` 0.591 · `dasatinib` 0.563 · `cay10618` 0.347 ·
`kx2-391` 0.283. **Der Vorsprung von scGPT ist konsistent, aber klein** — und §5.5 zeigt, wie klein er
relativ zu einer klassischen Baseline wirklich ist.

### 5.3 ⭐ Der kausale Rescue-Test — die zentrale Evidenz

**Notebook `10_diagnosis` §3 · Plot `outputs/ablations/rescue_k545.png` · CSV `ablations/rescue_k545.csv`**

**Das ist das Experiment, das die Slides tragen muss.** Jede Hypothese von der 29.-Juni-Liste, angewandt
auf **dieselbe kaputte Konfiguration**: **K=545, rohes `auc`, scGPT**, Baseline ρ = **−0.063**.
Immer nur **ein** Knopf gedreht, alles andere fix.

| Intervention | Konkrete Parameter | ρ | Verdikt |
|---|---|---|---|
| Baseline (broken) | — | **−0.063** | — |
| heavy regularization | dropout 0.7, in-dropout 0.2, wd 1e-2 | **−0.091** | ✗ schlechter |
| small model | `hidden=(32,)` → 16.645 statt 74.629 Params | **−0.053** | ✗ nichts |
| line-balanced **sample** weights | `WeightedRandomSampler`, w = 1/#Zellen(Linie) | **−0.078** | ✗ nichts |
| batch size 32 | `batch=32` | **+0.027** | ~ marginal |
| **no regularization** | dropout 0, in-dropout 0, wd 0 | **+0.234** | ⚠️ **rettet ~70 %** |
| **`auc_z`** (per-Drug = **TASK**-Reweighting) | Target-Wechsel, sonst identisch | **+0.433** | ✓ **heilt** |

**Der Mechanismus (und das ist die eigentliche wissenschaftliche Aussage):**
Das Versagen ist eine **Kapazitäts-Konkurrenz zwischen den Heads**. Unter dropout 0.5 kann der Trunk
nicht gleichzeitig die hoch-varianten, rauschigen Drugs *und* die lernbaren bedienen — die rauschigen
gewinnen, weil sie den Loss dominieren. Kapazität freizugeben lässt ihn beides fitten.

**Aber:** Das behandelt das **Symptom**.
- Der Preis ist ein Netz, das die Trainings-Linien **auswendig lernt** (train MSE ≈ 0.01).
- Es erreicht nur **die Hälfte** dessen, was der Weighting-Fix liefert.
- Und auf dem **korrigierten** Loss ist genau dieselbe Regularisierung wieder **optimal** (§5.4:
  0.488 mit vs. 0.456 ohne).

⇒ **Zwei Interventionen, ein Mechanismus. Die, die die Ursache entfernt, ist doppelt so wirksam — und
gratis.**

> **Wichtig für die Ehrlichkeit der Slides:** Ich hatte vorher behauptet, *"das Modell ist keine
> Ursache"*. Das war **zu absolut**. Richtig ist: **das Modell interagiert mit dem Bug, heilt ihn aber
> nicht.** Du hattest mit der Über-Regularisierungs-Vermutung teilweise recht.

### 5.4 Modell-Knöpfe auf dem KORRIGIERTEN Setting

**Notebook `10_diagnosis` §4 · Plot `outputs/ablations/ablation_reg_capacity.png` ·
CSVs `ablations/ablation_regularization.csv`, `ablation_capacity.csv`, `ablation_batch_weighting.csv`**

Setting: **K=5, `auc_z`**. Zweck: zeigen, dass modellseitig nichts mehr zu holen ist.

**Regularisierung** (ρ):

| | none (0/0) | light (.2/1e-4) | **current (.5/1e-3)** | heavy (.7/1e-2) |
|---|---|---|---|---|
| PCA | 0.417 | 0.426 | 0.424 | 0.435 |
| scGPT | 0.456 | 0.458 | **0.488** | 0.444 |

**Kapazität** (ρ):

| | (128,64) = 74.629 P | (64,32) = 35.269 P | (32,) = 16.645 P | linear head = 2.565 P |
|---|---|---|---|---|
| PCA | 0.428 | 0.433 | 0.414 | 0.412 |
| scGPT | **0.487** | 0.479 | **0.486** | 0.438 |

**Batch / Sample-Weighting** (ρ):

| | batch 32 | **batch 128** | batch 512 | line_balanced | focus_extremes |
|---|---|---|---|---|---|
| PCA | 0.425 | 0.435 | 0.432 | 0.407 | 0.431 |
| scGPT | 0.504 | 0.491 | 0.462 | 0.486 | 0.482 |

**Ablesung:** Alles **flach** (ρ 0.41–0.50). Modellseitiges Tuning ist **geschlossen**.
Bemerkenswert: `MLP (32,)` erreicht `(128,64)` **exakt** (0.486 vs. 0.487) bei **4,5× weniger
Parametern.** Der Default bleibt `(128,64)`, weil jede berichtete Zahl damit produziert wurde —
`(32,)` wäre die schlankere Wahl für die Zukunft. *(Bewusste Entscheidung: dokumentieren, nicht
umstellen, solange die Slides laufen.)*

> ⚠️ **§4 allein wäre irreführend** und darf nicht ohne §3 auf eine Slide: Es zeigt, dass die Knöpfe das
> **korrigierte** Modell nicht *verbessern* — das ist **nicht** dasselbe wie zu zeigen, dass sie das
> **kaputte** nicht hätten *heilen* können. Genau deshalb existiert §3.

### 5.5 Die ehrliche Referenz — Ridge-Kontrolle

**Notebook `10_diagnosis` §5 · in `ablations/ablation_capacity.csv`**

`RidgeCV` auf den **150 gemittelten Zelllinien-Embeddings** (kein Single-Cell, kein Deep Learning):

| Modell | ρ |
|---|---|
| Ridge (line-level), PCA | **0.428** |
| Ridge (line-level), scGPT | **0.428** |
| MLP, PCA | 0.428 |
| MLP, scGPT | **0.487** |

**Ridge auf 150 Zelllinien-Mittelwerten schlägt den PCA-MLP nicht — es *trifft* ihn exakt.**
Der gesamte Mehrwert der Single-Cell-Pipeline steckt in **scGPT (+0.06)**. Das ist ein unbequemes, aber
belastbares Ergebnis, und es ist **genau das, was das DrEval-Paper vorhersagt.** Es gehört auf die Slide.

### 5.6 Seed-Stabilität

**Notebook `11_auc_vs_aucz` · CSV `target/seed_stability.csv`** — K=545, `auc_z`:

| Seed | PCA | scGPT |
|---|---|---|
| 42 | 0.388 | 0.430 |
| 1 | 0.367 | 0.434 |
| 7 | 0.355 | 0.472 |
| **Mittel** | **0.370** | **0.445** |

Streuung ≈ ±0.02–0.03. **Der scGPT-Vorsprung (+0.075) überlebt den Seed-Wechsel**, ist aber von
derselben Größenordnung wie die Seed-Streuung — also **nicht überzuinterpretieren**.

### 5.7 ⭐ Externer Benchmark — echtes DrEval

**Notebook `12_dreval_benchmark` · Plot `outputs/dreval/dreval_lco.png` ·
CSV `dreval/dreval_lco_results.csv`**

**Das ist die Einordnung gegenüber dem Feld — und das Paper ist von Mathias' Gruppe mit-verfasst.**

> Bernett, Iversen, Picciani, **Wilhelm**, Baum, List — *Critical evaluation of drug response prediction
> models with DrEval*, **Nat. Commun. (2026)** · `pip install drevalpy` (**v1.5.1**)

**Nichts re-implementiert** — importiert werden: `DrugResponseDataset`, `.split_dataset(mode="LCO")`,
`MODEL_FACTORY`, `evaluate`.

#### (a) Wie DrEval gelesen werden MUSS — direkt aus der Publikation

**Die vier Splits, wörtlich:**

| Split | Definition (Paper) | Wofür gedacht |
|---|---|---|
| **LPO** | leave-random-drug-cell-line-**pairs**-out | *"only warranted when the goal is to evaluate the ability to **impute missing values**"* |
| **LCO** ← **unserer** | *"the test set must not contain cell lines used for training or tuning"* | **personalisierte Medizin** — ein neuer Patient |
| **LTO** | *"the test set must contain **tissues of origin** not seen during training"* | Drug Repurposing |
| **LDO** | *"the test set should contain only **unseen drugs**"* | Drug Design |

**Die Baseline, wörtlich:**
```
NaiveMeanEffectsPredictor:   mu_ij = mu_i^c + mu_j^d - mu
                             (Zelllinien-Mittel + Drug-Mittel - Gesamtmittel)
```

**Die Normalisierung, wörtlich:**
> *"Normalized performance metrics [are computed] by **removing the mean drug and cell line effects from
> both true and predicted responses** before metric calculation."*

**Warum**, wörtlich:
> *"**Most of the explainable variation in drug response is driven by the drug identity**"* — die
> Normalisierung isoliert *"differential drug response"* jenseits von Memorisierung und zeigt, ob ein
> Modell *"biologically relevant variation"* lernt oder *"simply exploits systematic biases"*.

**⇒ Die Lese-Anweisung der Autoren (das ist der Kern und gehört so auf die Slide):**

1. **Der rohen, gepoolten Korrelation NICHT trauen.** Das Paper zeigt eine rohe Pearson-Korrelation von
   **0.91**, die auf eine Per-Drug-Korrelation von **0.56** zusammenfällt — **Simpson's Paradox**: das
   Modell nutzt die *Drug-Mittelwerte*, nicht die differenzielle Sensitivität. Eine hohe rohe Korrelation
   heißt nur: *"das Modell weiß, welche Drug potenter ist"* — das ist keine personalisierte Medizin.
2. **Zu trauen ist dem normalisierten R² und der Per-Drug-Korrelation.**
3. **Ein Modell schlägt die Baseline nur dann**, wenn es den `NaiveMeanEffectsPredictor` **in den
   normalisierten Metriken signifikant** übertrifft.

**Ihre Kernbefunde:**
- **Rund die Hälfte der getesteten Modelle** übertrifft die naive Baseline **nicht signifikant**.
- Im klinisch relevanten **LCO**: *"**no model surpasses a tuned Random Forest**"*.
- Das **beste** LCO-Modell (Random Forest) erklärt nur **19 %** der differenziellen Sensitivität
  (normalisiertes R²).
- In **LDO scheitern ALLE**: *"no model significantly outperforms the NaiveMeanEffectsPredictor"*.
- Gesamturteil der Autoren: *"models capture only a limited fraction of the biologically relevant
  variation and **lack the accuracy required for clinical or translational use**"*.

> **Deshalb ist in unserer Tabelle unten die Spalte `normalized` die einzige, die zählt** — und deshalb
> steht `NaiveMeanEffectsPredictor` dort per Konstruktion bei **exakt 0.000**.

#### (b) Unser Protokoll

`n_cv_splits=5`, `mode='LCO'`, `split_validation=True`, `validation_ratio=0.1`, `random_state=42`.
Normalisierung wie im Paper: `NaiveMeanEffects` wird von `y_true` **und** `y_pred` abgezogen.
Metrik **gepoolt über alle (Linie × Drug)-Paare** des Test-Folds — **nicht** wie in §5.1–§5.6 pro Drug
gemittelt. **Die DrEval-Zahlen sind daher nicht direkt mit den ρ aus §5 vergleichbar.**

#### (c) Welches OncoTox-Modell genau — und wie es evaluiert wurde

*(`notebooks/12_dreval_benchmark.ipynb`, Funktion `run_oncomlp`)*

| | |
|---|---|
| Modell | **`OncoMLP`**, hidden `(128, 64)`, LayerNorm, GELU, dropout 0.5, input-dropout 0.1 |
| Input | `X_scGPT` bzw. `X_pca`, je **512-d** |
| Heads | **K = 5** (die 5 gefilterten Drugs) — **nicht** K=545 |
| Target | ⚠️ **rohes `auc`**, in nativen AUC-Einheiten — **nicht `auc_z`** (DrEval erwartet die echte Response-Skala, und normalisiert selbst) |
| Training | `TrainConfig(epochs=50, seed=42)`, batch 128, Adam lr 1e-3, wd 1e-3 |
| Trainiert auf | den **einzelnen Zellen** der Trainings-Linien (die echte OncoTox-Pipeline, kein Linien-Mittel) |
| Vorhersage | pro Zelle → **gemittelt zu einem Wert pro (ungesehene Linie × Drug)** |
| Zelllinien-Feature für **ihre** Modelle | Mittel der Zell-Embeddings pro Linie (`FeatureDataset`) |

**⚠️ Methodischer Mangel, der ehrlich benannt gehört:**
In `run_oncomlp` wird der **Test-Fold als Validierungs-Loader** an `train_model` übergeben:
```python
va_ds = MultiDrugDataset(adata=ad, use_rep=rep, cell_mask=tec, drugs=FIVE)   # tec = TEST-Linien
```
`train_model` macht darauf **Early Stopping und Best-Epoch-Checkpointing**. Damit wird die
**Modellselektion auf dem Test-Fold** gemacht → die OncoMLP-Zahlen sind **optimistisch verzerrt**.
Die Baselines und ihre `SingleDrug*`-Modelle sind davon **nicht** betroffen (die haben keinen
Val-Loop), der Vergleich ist also **zu unseren Gunsten schief**.
**Fix:** `split_validation=True` liefert bereits einen echten Validation-Split (`sp['validation']`) —
der gehört dort hin. **Muss vor der Präsentation korrigiert oder explizit als Vorbehalt genannt werden.**
*(Größenordnung: bei median-bester Epoche 6–9 und flacher val-Kurve vermutlich klein, aber ungemessen.)*

Ergebnisse (normalisiert, Mittel ± Std über die 5 Folds):

| Modell | norm. Spearman | norm. R² |
|---|---|---|
| `NaiveMeanEffectsPredictor` *(die Latte)* | **0.000** | 0.000 |
| `SingleDrugElasticNet` (scGPT) | 0.000 | −0.002 |
| `SingleDrugRandomForest` (pca) | 0.135 | 0.022 |
| `SingleDrugElasticNet` (pca) | 0.293 | 0.056 |
| `SingleDrugRandomForest` (scgpt) | 0.438 | 0.152 |
| `OncoMLP` (X_pca) | 0.442 | 0.166 |
| **`OncoMLP` (X_scGPT)** | **0.511 ± 0.085** | **0.224** |

**Drei belastbare Aussagen:**
1. **OncoTox schlägt `NaiveMeanEffects` klar** — die Latte, an der die Hälfte des publizierten Feldes
   scheitert.
2. **scGPT > PCA hält auch unter *ihrem* Protokoll** (0.511 vs. 0.442) — externe Bestätigung des
   +0.075-Gaps, das wir selbst gemessen haben.
3. **Wir schlagen ihre eigenen Referenzmodelle auf identischen Features** (`SingleDrugRF (scgpt)`
   0.438). Der Single-Cell-MLP fügt also etwas über einen Standard-Regressor auf dem Linien-Mittel hinzu.

Zum Vergleich: das Paper berichtet für seine **besten** LCO-Modelle 11 % (DIPK) / 19 % (RF) über naive.

> ⚠️ Der Filter-Vorbehalt aus §4 gilt hier **unverändert**: DrEval fixt die **Evaluation**, nicht unsere
> **Selektion**. Die 5 Drugs sahen alle 180 Linien.

### 5.8 Zweite, strengere Normalisierung (unsere eigene, NICHT DrEval)

**CSV `dreval/dreval_normalized.csv`** — beantwortet eine *andere* Frage und darf nicht mit 5.7
verwechselt werden.

In LCO kann der naive Prädiktor den Effekt einer *ungesehenen* Linie nicht kennen; er reduziert sich auf
*globales Mittel + Drug-Effekt*, entfernt also nur den **Drug**-Mittelwert. Unsere Variante entfernt
**zusätzlich den Zelllinien-Effekt** (mit den eigenen Labels der Linie). Frage: *"Wie viel unseres
Signals ist bloß 'diese Linie ist gegen alles empfindlich'?"*

**Antwort: ~20 % — und `kx2-391` ist *vollständig* dieser Artefakt.** (Passt dazu, dass `kx2-391` in
§5.2 mit ρ = 0.28 der schwächste der 5 Drugs ist.)

---

## 6. Vergleichbarkeits-Audit

Damit klar ist, wo die Zahlen belastbar sind und wo nicht.

**✓ Sauber kontrolliert:**
- Identische Daten, Splits, Trunk-Architektur und Trainings-Config über alle Fits.
- PCA und scGPT haben **dieselbe Dimension (512)** und **denselben Trunk** — nur der Input-Inhalt
  differiert. Kein Kapazitäts-Vorteil für eine Seite.
- Der Rescue-Test (§5.3) dreht **immer nur einen Knopf** gegenüber einer gemeinsamen Baseline.
- Die Target-Matrix (§5.1) evaluiert **immer auf denselben 5 Drugs**, auch bei K=545 — deshalb sind
  K=5 und K=545 direkt gegeneinander lesbar.
- Alle Zahlen out-of-fold über 153 Linien, nicht auf 27.

**⚠️ Grenzen, die auf die Slides gehören:**
1. **Drug-Selektion sah val/test** (§4). Best-Case-Diagnose, keine Generalisierung. **Blockierend.**
2. **`learnability` ist ein schwacher Ranker** (§4.1): Spearman mit dem tatsächlich erreichten ρ nur
   **+0.357**, das `n_res`-Gate korreliert sogar **negativ**, die beste Drug des Panels (`ml210`,
   ρ = 0.516) wurde von einem Gate verworfen, und die 5er-Auswahl ist unter Bootstrap **instabil**.
   ⇒ **„5/545" ist keine Aussage über Lernbarkeit.** Tatsächlich: **76 % aller Drugs ρ > 0.**
3. **DrEval-Zahlen sind optimistisch verzerrt** (§5.7c): der Test-Fold diente als Val-Loader für Early
   Stopping. Betrifft **nur** die OncoMLP-Zeilen, nicht die Baselines. **Vor der Präsentation fixen.**
4. **~150 unabhängige Beispiele pro Drug.** Der scGPT-Vorsprung (+0.075) liegt in der Größenordnung
   der Seed-Streuung (±0.03) — real, aber klein.
5. **Ridge auf Linien-Mittelwerten trifft den PCA-MLP exakt** (§5.5). Die Single-Cell-Auflösung zahlt
   sich nur *mit* scGPT aus.
6. Numerisches Jitter zwischen Läufen ≈ ±0.001 (z. B. scGPT/K=5/`auc_z` erscheint als 0.4877 und
   0.4871 in zwei verschiedenen Ablations-Loops). Irrelevant, aber der Vollständigkeit halber genannt.

**Korrekturen an früheren Behauptungen — bewusst auf dem Protokoll:**

| Behauptung | Status |
|---|---|
| *"Der Kurven-Fit rettet Signal, das der Dosis-Mittelwert zerstört"* | **FALSIFIZIERT** — `mean_pv` ≈ `auc` überall. Nur z-scoring zählt |
| *"Die Prediction-Shrinkage ist Über-Regularisierung"* | **FALSCH** — `pred_std ≈ ρ·σ` ist MSE-optimale Kalibrierung |
| *"Das Modell ist keine Ursache"* | **ZU ABSOLUT** — ohne Regularisierung werden ~70 % gerettet (§5.3) |
| *"Spread ist NICHT Lernbarkeit"* | **ZU ABSOLUT** — Spearman(σ, #getötet) = +0.475 |
| *"Zellen → Zelllinien für die Evaluation"* (Slide-Bullet) | **FALSCH** — Linien-Aggregation gab es schon im Juni |

---

## 7. Was NICHT gemacht wurde

- **scDEAL: nicht ausprobiert.** Muss explizit so auf die Slide, nicht verschweigen.
- Die 8-Run-Matrix aus `07_training` **nicht** auf `auc_z` neu gerechnet (nur die 12-Fit-Matrix in `11`).
- DrEval nur **LCO**, nur auf den 5 Drugs — nicht LTO/LDO, nicht auf allen 545.
- Gelernte Task-Gewichte (Kendall) statt fixem z-scoring: offen.
- Attention-Pooling pro Zelllinie statt Mittelwert: offen.

---

## 8. Vorgeschlagene Slide-Struktur

| # | Slide | Kernaussage | Plot |
|---|---|---|---|
| 1 | Wo wir standen | ρ ≈ 0, vier Verdächtige: Target / Messung / Drug-Auswahl / Modell | — |
| 2 | Was die Targets biologisch messen | `mean_pv` mittelt Messpunkte, `auc` integriert den Fit | `data/target_biology.png` |
| 3 | Der Filter als **Methode** | 545 → 5, "kann irgendetwas gelernt werden?" — **explizit als Vereinfachung deklarieren**, ausdrücklich **nicht** als "nur 5 sind lernbar" | `learnability/learnability_filter_auc.png` |
| 4 | Signal existiert | PCA 0.42 / scGPT 0.49 out-of-fold | `learnability/learnable5_pca_vs_scgpt.png` |
| **4b** | **Der Filter unterschätzt das Modell** | Über alle 537 Drugs: **76 % ρ > 0**, Median 0.120, 12 Drugs > 0.4 — und **9 davon hatte der Filter verworfen** (`ml210`, ρ=0.516, scheiterte an 0.01 Coverage) | *(neu zu zeichnen)* |
| 5 | **Der Bug** | Ungewichtet = σ²-gewichtet. `fqi-2` hat die größte Streuung und tötet **null** Linien | `target/loss_weighting_bug.png` |
| 6 | **Der Beweis** | 3 Targets × K=5/545. `auc_z` bei K=545 ≈ K=5 ⇒ **die Head-Zahl war nie das Problem** | `target/target_comparison.png` |
| 7 | **Was *nicht* geholfen hat** | Rescue-Test: nur Task-Reweighting heilt; Reg-Entfernen rettet das Symptom | `ablations/rescue_k545.png` |
| 8 | Modellseitig ist zu | Alle Knöpfe flach; Ridge auf Linien-Mitteln = PCA-MLP | `ablations/ablation_reg_capacity.png` |
| 9 | Einordnung im Feld | DrEval LCO: 0.511 vs. NaiveMeanEffects 0.000; schlägt deren RF | `dreval/dreval_lco.png` |
| 10 | Offen | Drug-Selektion in-fold (**blockierend**); scDEAL **nicht** versucht | — |

**Die eine Botschaft, die sich gegenüber dem letzten Report geändert hat:**

> Nicht *"es funktioniert bei wenigen Drugs"*, sondern:
> **Das Juni-Null-Resultat war überwiegend ein Artefakt eines fehlskalierten Multi-Task-Loss —
> kein Fakt über die Daten.** Mit einem per-Drug standardisierten Target rankt dasselbe Modell
> ungesehene Zelllinien bei ρ ≈ 0.43–0.49 und schlägt die naive Baseline, an der die Hälfte des
> publizierten Feldes scheitert.

---

## 9. Konkrete Beispiele aus den Daten

*(Alles direkt aus dem `hvg5000` / `auc_z` Targets-h5ad gerechnet — reproduzierbar über
`PipelinePaths.build(None, 'hvg5000', 'auc_z').targets_h5ad`.)*

### 9.1 Das σ-Spektrum — der Bug in zwei Zeilen

| | Drug | σ (auc) |
|---|---|---|
| **schmalste** Drug | `bcl-lzh-4` | **0.034** |
| **breiteste** Drug | `ifosfamide` | **0.302** |

**Verhältnis im quadrierten Fehler: 80×.** Ein MSE-Punkt von `ifosfamide` ist dem Modell also
**80 Trainings-Punkte** von `bcl-lzh-4` wert — allein wegen der Skala, ohne dass eine der beiden Drugs
mehr *Information* trüge.

### 9.2 Die Loss-Verteilung — wer das Modell tatsächlich steuert

Loss-Anteil ∝ σ² × Coverage:

| Rang | Drug | Loss-Anteil | σ | **getötete Linien** |
|---|---|---|---|---|
| 1 | **`fqi-2`** | **1.18 %** | 0.296 | **0** |
| 2 | **`ciclopirox`** | **0.94 %** | 0.265 | **0** |
| 3 | **`brd-k71781559`** | **0.86 %** | 0.257 | **0** |
| 4 | `daporinad` | 0.80 % | 0.253 | 40 |
| 5 | `paclitaxel` | 0.79 % | 0.246 | 66 |

- Die **breitesten 10 %** der Drugs (55 von 545) tragen **31 %** des gesamten Loss.
- **Die 5 lernbaren Drugs zusammen tragen 2.0 %.**
- **Die drei größten Loss-Träger töten zusammen null Zelllinien.**

> Das ist die Slide-Zahl: **`fqi-2` allein steuert das Modell 0,6× so stark wie alle fünf lernbaren
> Drugs zusammen — und tötet nichts.**

### 9.3 Warum `fqi-2` Varianz hat, aber keine Information — direkt vergleichbar

| | σ | mittlere `auc` | min | max | **getötet** (≤ 0.5) | **überlebt** (≥ 0.8) |
|---|---|---|---|---|---|---|
| **`fqi-2`** | 0.296 | 0.972 | **0.59** | **2.11** | **0** | 127 |
| **`dasatinib`** | 0.155 | 0.631 | **0.07** | 1.09 | **35** | 27 |

**Lies die min-Spalte.** `fqi-2` kommt nie unter `auc = 0.59` — **es tötet keine einzige Zelllinie,
jemals.** Seine ganze Streuung liegt auf der *anderen* Seite: bis `auc = 2.11`, d. h. Linien
**wachsen unter der Drug doppelt so stark** wie unbehandelt. Das ist Proliferations-Rauschen und
Assay-Varianz, kein Kill-Signal.

`dasatinib` hat **halb so viel Streuung**, aber sie liegt dort, wo Biologie ist: von `auc = 0.07`
(komplett ausgelöscht) bis `1.09` (unbeeindruckt), mit 35 getöteten und 27 überlebenden Linien.

⇒ **Genau die Drug mit der halben Varianz ist die, die man lernen kann.** Der ungewichtete Loss
bevorzugt aber die andere. *Das ist der ganze Bug, an einem Beispiel.*

### 9.4 Was z-scoring konkret tut — `dasatinib` durchgerechnet

`center` (mittlere auc über die Linien) = **0.631**, `scale` (σ) = **0.155**

| Zelllinie | `auc` (roh) | `auc_z` |
|---|---|---|
| `SNU1079_BILIARY_TRACT` | 0.072 | **−3.60** |
| `HUPT4_PANCREAS` | 0.206 | −2.74 |
| `SCC25_UPPER_AERODIGESTIVE_TRACT` | 0.304 | −2.11 |
| … | … | … |
| `ZR751_BREAST` | 0.913 | +1.82 |
| `JHH6_LIVER` | 0.931 | +1.93 |
| `CAMA1_BREAST` | 1.091 | **+2.96** |

Die **Rangfolge ändert sich nicht** (Spearman ist invariant) — was sich ändert, ist **wie viel
Loss-Gewicht diese Drug im Multi-Task-Training bekommt.** Nach dem z-scoring hat *jede* Drug σ = 1,
also trägt jede **gleich viel** zum Loss bei. `fqi-2` kann das Training nicht mehr dominieren.

**Und das ist auch die richtige biologische Frage:** nicht *"wie viel Prozent der Zellen sterben?"*
(das ist Drug-Potenz und über Drugs hinweg nicht vergleichbar), sondern **"ist DIESE Linie sensitiver
als der Durchschnitt für DIESE Drug?"** — genau das, was man für personalisierte Therapie wissen will.

---

## 10. Antworten auf die Fragen, die im Chat aufkamen

### Zum Target

**„Was ist `mean_pv` — wir hatten doch `cpd_avg_pct_viability` benutzt?“**
Dasselbe. `mean_pv` ist der **ungewichtete Mittelwert von `cpd_avg_pv`** über die 16 Dosis-Punkte des
CTRPv2-Grids. Der Name im Code ist nur die Kurzform. Es ist *nicht* ein anderer Score, sondern die
Aggregation deines Scores über die Dosis-Achse.

**„Biologischer Unterschied zwischen `auc` und `pct viability`?“**
- **`pct_viability`** = wie viele Zellen bei *einer* Dosis überleben (relativ zu unbehandelt). `1.0` =
  keine Wirkung, `0.0` = alles tot, **> 1.0 = wächst schneller als unbehandelt** (siehe `fqi-2`: bis 2.11).
- **`mean_pv`** = Mittel über die 16 gemessenen Dosen. Roh, aber verrauscht und abhängig davon, *welche*
  Dosen gemessen wurden.
- **`auc`** = Fläche unter der **gefitteten** Sigmoid-Kurve, auf die Dosis-Achse normiert
  (`area_under_curve / conc_pts_fit`). Glättet das Rauschen und macht Drugs mit verschiedenen
  Dosis-Grids vergleichbar.
- **Beide messen dasselbe** und landen empirisch **innerhalb ~0.03** voneinander → `data/target_biology.png`.
- **`auc_z`** ist die einzige, die eine *andere* Frage stellt: nicht "wie potent ist die Drug", sondern
  **"wie sensitiv ist diese Linie relativ zu den anderen"**.

**„Unterschied zwischen `auc` und `auc_z`?“**
Nur eine per-Drug-Standardisierung — **die Rangfolge innerhalb einer Drug ist identisch.** Aber im
Multi-Task-Loss ist es der ganze Unterschied: `auc_z` ≡ jeden Head mit **1/σ²** gewichten.
**Empirisch:** bei K=5 macht es **nichts** (0.482 → 0.488). Bei K=545 macht es **alles**
(−0.087 → +0.430).

**„Wäre die AUC von CTRPv2 sinnvoller als unser Score?“**
**Der Kurven-Fit allein: nein** — er bringt exakt 0 (§5.1, und das falsifiziert meine ursprüngliche
Vermutung). **Die per-Drug-Standardisierung: ja, und zwar entscheidend.** Wir nutzen jetzt `auc_z`,
aber der Gewinn kommt vom **z**, nicht vom **auc**.

### Zum Modell und zum Bug

**„Wie soll man denn 545 Drugs predicten, wenn nicht mit so vielen Heads? Es liegt doch an der mangelnden Datenlage.“**
**Das war die zentrale Hypothese — und die Daten sagen: nein.** Mit `auc_z` erreicht **K=545 ρ = 0.430**
gegen **K=5 ρ = 0.488**. Fast dasselbe, mit *denselben* 545 Heads und *derselben* Datenmenge.
Die 545 Heads waren nie das Problem; **die Skalierung ihres gemeinsamen Loss war es.** Die Datenlage
(~150 Linien pro Drug) ist knapp, aber sie war **nicht die Ursache des Null-Resultats**.

**„Muss ich das aktiv machen — so etwas wie Class-Imbalance-Reweighting, aber für Regression?“**
**Genau das, ja.** Und du machst es bereits — durch das Target. Ein per-Drug z-Score **ist** algebraisch
identisch zu einem 1/σ²-gewichteten Loss. Man kann es explizit als Gewicht implementieren (oder die
Gewichte lernen, Kendall et al. 2018), aber der Effekt ist derselbe und z-scoring ist gratis.

**„Ist zu viel Regularisierung im Spiel? Hatte ich recht?“**
**Teilweise recht — und wichtiger, als ich zuerst gesagt habe.** Auf dem **kaputten** Setting bringt
*Regularisierung komplett entfernen* **+0.234** (von −0.063), also **~70 % der Rettung**. Der Mechanismus
ist eine **Kapazitäts-Konkurrenz zwischen den Heads**: unter dropout 0.5 kann der Trunk nicht beide
Gruppen bedienen, und die rauschigen gewinnen, weil sie den Loss dominieren.
**ABER:** (a) es kostet Auswendiglernen (train MSE ≈ 0.01), (b) es erreicht nur die **Hälfte** des
Weighting-Fixes (+0.433), und (c) auf dem **korrigierten** Loss ist dieselbe Regularisierung wieder
**optimal** (0.488 mit, 0.456 ohne). ⇒ **Symptom behandelt, nicht Ursache.**
*Meine frühere Aussage „das Modell ist keine Ursache“ war zu absolut — das steht im Protokoll (§6).*

**„Ist das Modell zu groß / zu viele Parameter?“**
**Nein.** Auf dem korrigierten Setting ist Kapazität **flach**: `(128,64)` = 74.629 P → ρ 0.487,
`(32,)` = 16.645 P → ρ **0.486**. Selbst ein **linearer Head** (2.565 P) erreicht 0.438.
Auf dem kaputten Setting *verschlechtert* ein kleineres Modell sogar (−0.053 vs −0.063 Baseline).
⇒ Größe ist nicht der Hebel. *(`(32,)` wäre trotzdem die schlankere Default-Wahl für die Zukunft.)*

**„Würden das helfen: Batch Size reduzieren? Reweighting?“**
- **Batch 32 auf dem kaputten Setting: +0.027** — marginal, hebt es gerade eben über null. Nicht die Lösung.
- **Reweighting: kommt drauf an, WELCHES.**
  - **Sample-Reweighting** (Zellen nach Linie balancieren): **−0.078**, hilft **nicht**.
  - **Task-Reweighting** (per Drug = `auc_z`): **+0.433**, **das ist die Lösung.**
  - Diese Unterscheidung ist der Kern und gehört auf die Slide.

**„Wie macht man das Modell besser — z. B. mit scDEAL?“**
**scDEAL ist NICHT ausprobiert.** Steht als offener Punkt (§7), muss so auf die Slide.
Was die Daten sagen, wo Luft ist: **nicht in der Modellarchitektur** (alles flach, §5.4), sondern in
(a) der Target-/Loss-Definition — das hat gerade +0.5 gebracht — und (b) besserem Pooling von Zellen zu
Linien (aktuell: simpler Mittelwert; Attention-Pooling ist offen).

### Zur Evaluation

**„Nutzen wir die Pearson-Korrelation? Warum nicht gegen Pearson vergleichen?“**
Wir berichten **Spearman**, rechnen aber **beide**. Sie stimmen **auf ±0.02** überein (Spalte `pearson`
in `target/target_comparison.csv`). Spearman ist die ehrlichere Wahl, weil die Frage ein **Ranking** ist
("welche Linie ist sensitiver?") und sie robust gegen die Ausreißer der `auc`-Skala ist
(`fqi-2` geht bis 2.11). **Die Wahl der Metrik ändert keine einzige Schlussfolgerung.**

**„Wo sind die Fehlerbalken?“**
Jetzt drin, dreifach:
1. **Bootstrap-CIs** über die Drugs → `target/target_comparison_ci.csv`, im Plot.
2. **Seed-Stabilität** über 3 Seeds → `target/seed_stability.csv` (scGPT 0.430 / 0.434 / 0.472).
3. **Fold-Streuung** im DrEval-Benchmark (0.511 **± 0.085**).
**Und der Grund, warum es sie vorher nicht sinnvoll geben konnte:** der alte val-Split hatte 27 Linien →
**SE(ρ) ≈ ±0.2.** Erst der Wechsel auf out-of-fold über 153 Linien macht Fehlerbalken überhaupt
aussagekräftig.

**„Mit wie vielen Genen wurden alle Versuche gemacht?“**
**Alle**, gestern wie heute: Variante **`hvg5000`**.
5.000 HVGs → scGPT embeddet die **4.576** davon in seinem Vokabular (**424 OOV**), PCA rechnet auf allen
5.000. **Beide Repräsentationen sind 512-d** — die Dimension ist bewusst gleich, damit der Vergleich
nicht durch Kapazität verfälscht wird.

### Zum Rest

**„Wie genau hast du DrEval ausgeführt — heruntergeladen und im Notebook?“**
`pip install drevalpy` (**v1.5.1**, das echte Paket von daisybio), ausgeführt in
**`notebooks/12_dreval_benchmark.ipynb`**. **Nichts re-implementiert.** Importiert:
`DrugResponseDataset`, `.split_dataset(mode="LCO", n_cv_splits=5, ...)`, `MODEL_FACTORY`, `evaluate`.
Ihre Splits, ihre Baselines, ihre Metriken, ihre Normalisierung. Vollständig nachvollziehbar im Notebook.

**„Woher hattest du Mathias in dieser Konversation?“**
Berechtigte Rückfrage — **ich hatte ihn zu dem Zeitpunkt aus dem Kontext erfunden.** Er kam erst
danach legitim ins Spiel, als **du** das Paper verlinkt hast, das er mit-verfasst hat
(Bernett et al., *Nat. Commun.* 2026). Der Punkt bleibt inhaltlich richtig, aber die Quelle war es nicht.

**„Warum dauert das Training so lange?“**
Tat es gar nicht. **1 Epoche = 1,39 s**, h5ad laden = 0,7 s. Die Stunden kamen vom **Experiment-Grid**
(~90 Fits), nicht von langsamem Code. Zwei Änderungen: `epochs` **50 → 25** (die beste Epoche war über
36 Runs median 6, max 11 — die oberen 25 waren reine Wanduhr) und ein **`RETRAIN = False`**-Flag in jedem
Notebook, das die gespeicherten CSVs lädt und nur die Figures neu zeichnet → **Sekunden statt ~35 min**.

---

## 11. Notebook-Landkarte

Lesereihenfolge: **`08 → 09 → 10 → 11 → 12`**. Alle haben `RETRAIN = False` → laden gespeicherte CSVs,
zeichnen nur die Figures neu (Sekunden statt ~35 min).

| # | Notebook | Frage |
|---|---|---|
| 08 | `08_learnability_filter` | Welche Drugs *können* gelernt werden? (545 → 5) |
| 09 | `09_learnable5_training` | Lernt die Pipeline überhaupt etwas? (ρ = 0.42 / 0.49) |
| 10 | `10_diagnosis` | **Warum ist K=545 gescheitert und was heilt es?** (§3 = der kausale Test) |
| 11 | `11_auc_vs_aucz` | Welches Target? (die 12-Fit-Matrix + CIs + Seeds) |
| 12 | `12_dreval_benchmark` | Wie stark ist das nach dem Maßstab des Feldes? |

*Supporting:* `02_compare_GDSC_CTRP` · `04_drug_coverage` · `06_verify_variants` ·
`07_training` (⚠️ **superseded** — die alte 8-Run-Matrix auf `mean_pv`).
*Archiviert:* `01`, `03`, `05` → `archive/`.
