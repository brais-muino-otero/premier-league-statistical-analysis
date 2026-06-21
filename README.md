# Premier League Statistical Analysis (2022-23)

> A from-scratch statistical analysis of the Premier League 2022-23 season, synchronized across **R**, **R Markdown**, and **Jupyter (Python)**, plus an interactive **Shiny** app with a 6-language selector.

**Live demo:** [braismuinootero.shinyapps.io/app_premier_league_prediction_22_23](https://braismuinootero.shinyapps.io/app_premier_league_prediction_22_23/)

**Companion project:** [datamill-analytics](https://github.com/brais-muino-otero/datamill-analytics) — a production-style Dash/Python web app built on the same dataset, focused on software architecture (Docker, i18n, testing, leakage-safe ML pipeline) rather than statistical breadth. This repository is the academic, statistics-first counterpart: classical hypothesis testing, GLMs, classification, unsupervised learning and linear programming, each justified methodologically.

---

## What this project covers

A single guiding question — **what determines a Premier League match result, and can it be predicted?** — threaded through ten sections:

1. **Exploratory data analysis** — distributions and correlations across goals, shots, corners, fouls, cards.
2. **Statistical inference** — hypothesis testing on the home-advantage effect (goals and shots on target, home vs away).
3. **Dixon-Coles model** — maximum-likelihood estimation of attack/defense strength per team (`optim`, BFGS, in both R and Python independently); home advantage estimated at a **×1.34** multiplier on expected home goals.
4. **Poisson regression** — goals scored as a function of match statistics (`glm`, family `poisson`).
5. **Logistic regression** — home-win classification (`glm`, family `binomial`); away corners (AC) found significant and interpreted as a "scoreboard effect" (teams trailing push more corners in search of an equalizer), consistent with the negative home-corners coefficient in the Poisson model.
6. **Match outcome classification** — LDA, QDA, and k-NN (k selected via cross-validation on shared folds) predicting Home/Draw/Away, evaluated with confusion matrices, ROC curves, and 10-fold CV. Final test results: **Logistic / LDA 77.6%**, **QDA 72.4%**, **k-NN 64.5%** accuracy (k = 21), McFadden's pseudo-R² = 0.181 for the logistic model.
7. **PCA + K-means** — principal component analysis (`prcomp` / scikit-learn with n−1 scaling) and K-means clustering (centroid-based cluster relabeling: highest `GolesFavor` → elite, highest `GolesContra` → struggling) to discover team/match profiles, validated with the silhouette score (0.307).
8. **Linear programming** — a constrained optimization problem solved with `lpSolveAPI` (R) / SciPy `linprog` (Python), including shadow-price (dual value) interpretation.
9. **Monte Carlo season simulation** — 10,000-season replay using the Dixon-Coles parameters to estimate title / relegation probabilities.
10. **Conclusions** — synthesis of findings across all methods.

### A deliberate data-leakage fix worth mentioning

An earlier iteration of the classifiers included `EfficiencyHome`/`EfficiencyAway` (goals ÷ shots on target) as predictors. Since that ratio encodes the match result itself, it would have produced inflated, unrealistic accuracy. Both variables were excluded from every classifier, settling on `HomeWin ~ HST + AST + HC + AC + HF + AF` — a formula built only from match statistics that don't trivially leak the outcome. A `ShotsTargetDiff` variable was also excluded due to perfect collinearity with HST/AST, which would have caused QDA to fail with a rank-deficiency error.

## Reproducibility across R and Python

R is the canonical implementation. Four CSV bridges let the Jupyter notebook reuse R's exact random splits and simulation output instead of re-deriving them independently, so results match between languages wherever randomness is involved. The Dixon-Coles model itself is **not** bridged: it is estimated independently by maximum likelihood in both R (`optim`) and Python (`scipy.optimize.minimize`), and both converge to the same deterministic solution (home-advantage ≈ 0.294, ×1.34).

| File | Purpose |
|---|---|
| `montecarlo_probs.csv` | Monte Carlo simulation output (title/top-4/relegation probabilities per team) |
| `idx_train.csv` | Train/test split indices, shared so both languages evaluate on the same rows |
| `folds_knn.csv` | Cross-validation fold assignments used to select k for k-NN |
| `folds_cv.csv` | Cross-validation fold assignments used for the final model comparison |

Running `analisis_premier_22_23.R` first regenerates all four bridge files; the notebook reads them if present and falls back to an independent (and clearly flagged) computation otherwise, so it never fails silently.

## Interactive Shiny app

**Try it live:** [braismuinootero.shinyapps.io/app_premier_league_prediction_22_23](https://braismuinootero.shinyapps.io/app_premier_league_prediction_22_23/)

`app.R` ports the analysis into a 4-section interactive dashboard:

- **Dixon-Coles simulator** — pick a home/away team, see their expected goals (λ) and the resulting 1/X/2 outcome probabilities.
- **Live prediction panel** — GLM-based home-win probability and expected goals from custom match-stat sliders.
- **Tactical explorer** — interactive PCA + K-means cluster visualization (Plotly), hovering over each team shows its profile (Elite / Mid-table / Strugglers).
- **Technical optimizer** — the linear programming model (training-hours allocation), with shadow-price output, recalculated live as the sliders move.

All heavy models (`optim`, `glm`, `prcomp`, `kmeans`, `lpSolveAPI`) are fit **once**, in the global R scope, at app startup — so the UI stays instantly responsive across sessions instead of refitting on every interaction.

The app ships with a **6-language selector** (Español, English, Português, Italiano, Français, Deutsch): every label, plot title, axis, cluster name, and validation message is translated through a self-contained dictionary (no external i18n dependency), with input values preserved across language switches via `isolate()`.

### Running the app

```r
# From the project root, in R or RStudio:
shiny::runApp("app.R")
```
Or click **Run App** in RStudio.

> `epl_results_2022-23.csv` must be present in the working directory. The four CSV bridges are optional — if missing, the app still runs (Monte Carlo probabilities are simply not shown), but it's recommended to run `analisis_premier_22_23.R` first to generate them.

### Required R packages

```r
install.packages(c("shiny", "bslib", "plotly", "lpSolveAPI", "cluster",
                   "corrplot", "lmtest", "pscl", "MASS", "class", "pROC"))
```

## Project structure

```
premier-league-statistical-analysis/
├── app.R                          # Interactive Shiny dashboard (6 languages)
├── analisis_premier_22_23.R       # Consolidated R script (sections 1-10)
├── analisis_premier_22_23.Rmd     # R Markdown report (knit to HTML/Word/PDF)
├── analisis_premier_22_23.pdf     # Rendered final report (knit from the .Rmd)
├── analisis_premier_22_23.ipynb   # Python/Jupyter notebook (reads the R-exported CSV bridges)
├── epl_results_2022-23.csv        # Source dataset (380 matches, 23 variables)
├── montecarlo_probs.csv           # R → Python bridge: Monte Carlo probabilities
├── idx_train.csv                  # R → Python bridge: train/test split indices
├── folds_knn.csv                  # R → Python bridge: k-NN cross-validation folds
└── folds_cv.csv                   # R → Python bridge: final model-comparison folds
```

## Tech stack

**R**: `optim`, `glm`, `prcomp`, `kmeans`, `lpSolveAPI`, `shiny`, `bslib`, `plotly`, `rmarkdown`.
**Python**: `pandas`, `numpy`, `scikit-learn`, `scipy`, `statsmodels` (for the synchronized notebook).

## Author

**Brais Muiño Otero** — Data Science student (UDC), targeting ML/applied-science roles in the UK/Ireland.

🔗 **GitHub:** [github.com/brais-muino-otero](https://github.com/brais-muino-otero)  ·  **Live demo:** [Shiny app](https://braismuinootero.shinyapps.io/app_premier_league_prediction_22_23/)
