# Premier League Statistical Analysis (2022-23)

> A from-scratch statistical analysis of the Premier League 2022-23 season in R, synchronized across **R**, **R Markdown**, and **Jupyter (Python)**, plus an interactive **Shiny** app with a 6-language selector.

**Live demo:** [braismuinootero.shinyapps.io/analisis_premier_22_23](https://braismuinootero.shinyapps.io/analisis_premier_22_23/)

**Companion project:** [datamill-analytics](https://github.com/brais-muino-otero/datamill-analytics) — a production-style Dash/Python web app built on the same dataset, focused on software architecture (Docker, i18n, testing, leakage-safe ML pipeline) rather than statistical breadth. This repository is the academic, statistics-first counterpart: classical hypothesis testing, GLMs, classification, unsupervised learning and linear programming, each justified methodologically.

---

## What this project covers

A single guiding question — **what determines a Premier League match result, and can it be predicted?** — threaded through ten sections:

1. **Exploratory data analysis** — distributions and correlations across goals, shots, corners, fouls, cards.
2. **Statistical inference** — hypothesis testing on the home-advantage effect.
3. **Dixon-Coles model** — maximum-likelihood estimation of attack/defense strength per team (`optim`, BFGS); home advantage estimated at a **×1.34** multiplier on expected home goals.
4. **Poisson regression** — goals scored as a function of match statistics (`glm`, family `poisson`).
5. **Logistic regression** — over/under 2.5 goals classification (`glm`, family `binomial`); away corners (AC) found significant (p = 0.009) and interpreted as a "scoreboard effect" (teams trailing push more corners in search of an equalizer), consistent with the negative home-corners coefficient in the Poisson model.
6. **Match outcome classification** — LDA, QDA, and k-NN (k selected via cross-validation) predicting Home/Draw/Away, evaluated with confusion matrices, ROC curves, and 10-fold CV. Final test results: **LDA 77.6%**, **QDA 72.4%**, **k-NN 64.5%** accuracy, McFadden's pseudo-R² = 0.181 for the logistic baseline.
7. **PCA + K-means** — principal component analysis (`prcomp`, n−1 scaling) and K-means clustering (centroid-based cluster relabeling: highest `GolesFavor` → elite, highest `GolesContra` → struggling) to discover team/match profiles, validated with the silhouette score.
8. **Linear programming** — a constrained optimization problem solved with `lpSolveAPI` (R) / SciPy (Python bridge), including shadow-price (dual value) interpretation.
9. **Monte Carlo season simulation** — season replay using the Dixon-Coles parameters to estimate title / relegation probabilities.
10. **Conclusions** — synthesis of findings across all methods.

### A deliberate data-leakage fix worth mentioning

An earlier iteration of the classifiers included `EfficiencyHome`/`EfficiencyAway` (goals ÷ shots on target) as predictors. Since that ratio encodes the match result itself, it produced inflated accuracy (85-91%) and an unrealistic McFadden's R² of 0.70. Both variables were removed from every classifier, settling on `HomeWin ~ HST + AST + HC + AC + HF + AF` — a formula built only from match statistics that don't trivially leak the outcome. The honest numbers above (77.6% / 72.4% / 64.5%) reflect that correction. A `ShotsTargetDiff` variable was also excluded due to perfect collinearity with HST/AST, which had caused QDA to fail with a rank-deficiency error.

## Reproducibility across R and Python

Four CSV bridges let the Jupyter notebook reuse R's exact fitted parameters instead of re-deriving them independently, so deterministic results (coefficients, PCA loadings, the Dixon-Coles fit, the LP solution) match exactly between languages:

| File | Purpose |
|---|---|
| `parametros_dixon_coles.csv` | Fitted attack/defense strengths per team |
| `ventaja_local.csv` | Home-advantage multiplier |
| `montecarlo_probs.csv` | Monte Carlo simulation output (title/relegation probabilities) |
| `idx_train.csv` | Train/test split indices, shared so both languages evaluate on the same rows |
| `folds_cv.csv` / `folds_knn.csv` | Cross-validation fold assignments |

Stochastic results (CV fold splits, optimal k for k-NN: R selects k=21, Python's scikit-learn selects k=15) differ slightly by implementation but converge to the same test accuracy, and are reported as such rather than forced to match artificially.

## Interactive Shiny app

**Try it live:** [braismuinootero.shinyapps.io/analisis_premier_22_23](https://braismuinootero.shinyapps.io/analisis_premier_22_23/)

`app.R` ports the analysis into a 4-section interactive dashboard:

- **Dixon-Coles simulator** — explore attack/defense ratings and simulate matchups.
- **Live prediction panel** — GLM-based outcome prediction from custom inputs.
- **Tactical explorer** — PCA + K-means cluster visualization.
- **Technical optimizer** — the linear programming model, with sensitivity/shadow-price output.

All heavy models (`optim`, `glm`, `prcomp`, `kmeans`, `lpSolveAPI`) are fit **once**, in the global R scope, at app startup — so the UI stays instantly responsive across sessions instead of refitting on every interaction.

The app ships with a **6-language selector** (Español, English, Português, Italiano, Français, Deutsch): every label, plot title, axis, cluster name, and validation message is translated through a self-contained dictionary (no external i18n dependency), with input values preserved across language switches via `isolate()`.

### Running the app

```r
# From the project root, in R or RStudio:
shiny::runApp("app.R")
```
Or click **Run App** in RStudio.

> The four CSV bridges and `epl_results_2022-23.csv` must be present in the working directory — running `analisis_premier_22_23.R` first regenerates them if missing.

## Project structure

```
premier-league-statistical-analysis/
├── app.R                          # Interactive Shiny dashboard (6 languages)
├── analisis_premier_22_23.R       # Consolidated R script (sections 1-10)
├── analisis_premier_22_23.Rmd      # R Markdown report
├── analisis_premier_22_23.ipynb    # Python/Jupyter notebook (reads the R-exported CSVs)
├── analisis_premier_22_23.pdf      # Rendered final report
├── epl_results_2022-23.csv         # Source dataset (380 matches, 23 variables)
├── parametros_dixon_coles.csv      # R → Python bridge
├── ventaja_local.csv                # R → Python bridge
├── montecarlo_probs.csv             # R → Python bridge
├── idx_train.csv                    # R → Python bridge
├── folds_cv.csv                      # CV fold assignments
└── folds_knn.csv                     # k-NN fold assignments
```

## Tech stack

**R**: `optim`, `glm`, `prcomp`, `kmeans`, `lpSolveAPI`, `shiny`, `rmarkdown`.
**Python**: `pandas`, `numpy`, `scikit-learn`, `scipy` (for the synchronized notebook).

## Author

**Brais Muiño Otero** — Data Science student (UDC), targeting ML/applied-science roles in the UK/Ireland.

🔗 **GitHub:** [github.com/brais-muino-otero](https://github.com/brais-muino-otero)  ·  **Live demo:** [Shiny app](https://braismuinootero.shinyapps.io/analisis_premier_22_23/)
