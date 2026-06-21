###############################################################################
##  Análisis estadístico y predictivo de la Premier League 2022-2023
##  Aplicación interactiva en Shiny (R) · multilingüe · Autor: Brais Muiño Otero
##  Lanzar con:  shiny::runApp("app.R")   (o el botón "Run App" de RStudio)
##  (Guardar este archivo con codificación UTF-8.)
###############################################################################
##  Idiomas: Español, Português, English, Italiano, Français, Deutsch.
##  La LÓGICA (Dixon-Coles, GLM, PCA+K-means, programación lineal) es idéntica
##  a la del proyecto y los modelos pesados se ajustan UNA vez al arrancar
##  (ámbito global). El idioma solo cambia los textos de la interfaz (reactivo).
###############################################################################

.req <- c("shiny", "bslib", "plotly", "lpSolveAPI", "cluster")
.missing <- .req[!vapply(.req, requireNamespace, logical(1), quietly = TRUE)]
if (length(.missing))
  stop("Faltan paquetes. Ejecuta:\n  install.packages(c(",
       paste(sprintf('\"%s\"', .missing), collapse = ", "), "))")

library(shiny)
library(bslib)
library(plotly)
library(lpSolveAPI)
library(cluster)

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# =============================================================================
#  CARGA Y AJUSTE DE MODELOS (ÁMBITO GLOBAL -> UNA VEZ AL ARRANCAR)
# =============================================================================
DATASET <- "epl_results_2022-23.csv"
if (!file.exists(DATASET))
  stop("No se encuentra '", DATASET, "'. Colócalo junto a app.R.")

datos <- read.csv(DATASET, stringsAsFactors = FALSE)
datos$HomeWin <- as.integer(datos$FTR == "H")
message("Ajustando modelos (Dixon-Coles, GLM, PCA, K-means)... una sola vez.")

## ---- Dixon-Coles (optim BFGS) ----
equipos <- sort(unique(datos$HomeTeam)); n_eq <- length(equipos)
idx <- setNames(seq_along(equipos), equipos)
hi <- idx[datos$HomeTeam]; ai <- idx[datos$AwayTeam]; fthg <- datos$FTHG; ftag <- datos$FTAG
neg_log_ver <- function(p) {
  at <- p[1:n_eq]; de <- p[(n_eq + 1):(2 * n_eq)]; v <- p[2 * n_eq + 1]
  lh <- exp(v + at[hi] - de[ai]); la <- exp(at[ai] - de[hi])
  -(sum(dpois(fthg, lh, log = TRUE)) + sum(dpois(ftag, la, log = TRUE)))
}
.fit <- optim(c(rep(0, 2 * n_eq), 0.1), neg_log_ver, method = "BFGS", control = list(maxit = 10000))
ataque  <- setNames(.fit$par[1:n_eq], equipos)
defensa <- setNames(.fit$par[(n_eq + 1):(2 * n_eq)], equipos)
ventaja <- as.numeric(.fit$par[2 * n_eq + 1])

## ---- GLMs (temporada completa) ----
m_pois <- glm(FTHG    ~ HST + AST + HC + AC + HF + AF, data = datos, family = poisson)
m_log  <- glm(HomeWin ~ HST + AST + HC + AC + HF + AF, data = datos, family = binomial)

## ---- PCA (escalado n-1) + K-means ----
loc <- aggregate(cbind(GF = FTHG, GC = FTAG, Tiros = HS, TP = HST,
                       Corners = HC, Faltas = HF) ~ HomeTeam, datos, mean)
vis <- aggregate(cbind(GF = FTAG, GC = FTHG, Tiros = AS, TP = AST,
                       Corners = AC, Faltas = AF) ~ AwayTeam, datos, mean)
rownames(loc) <- loc$HomeTeam; loc$HomeTeam <- NULL
rownames(vis) <- vis$AwayTeam; vis$AwayTeam <- NULL
perfil <- as.matrix((loc[equipos, ] + vis[equipos, ]) / 2)

pca <- prcomp(perfil, scale. = TRUE)
scores <- pca$x[, 1:2]; load <- pca$rotation
if (load["GF", "PC1"] < 0)     { scores[, 1] <- -scores[, 1]; load[, "PC1"] <- -load[, "PC1"] }
if (load["Faltas", "PC2"] < 0) { scores[, 2] <- -scores[, 2]; load[, "PC2"] <- -load[, "PC2"] }
varexp <- (pca$sdev^2 / sum(pca$sdev^2)) * 100

set.seed(123)
km <- kmeans(scale(perfil), centers = 3, nstart = 25); cen <- km$centers
elite <- which.max(cen[, "GF"]); strug <- which.max(cen[, "GC"]); mid <- setdiff(1:3, c(elite, strug))
remap <- integer(3); remap[mid] <- 1L; remap[strug] <- 2L; remap[elite] <- 3L
cl <- remap[km$cluster]
sil <- mean(silhouette(km$cluster, dist(scale(perfil)))[, 3])

df_pca <- data.frame(
  Equipo = equipos, PC1 = scores[, 1], PC2 = scores[, 2],
  Cl = c("mid", "strug", "elite")[cl],          # código estable (no traducido)
  GF = round(perfil[, "GF"], 2), GC = round(perfil[, "GC"], 2),
  stringsAsFactors = FALSE
)

mc <- if (file.exists("montecarlo_probs.csv")) read.csv("montecarlo_probs.csv") else NULL

# =============================================================================
#  DICCIONARIO DE TRADUCCIONES  (es / pt / en / it / fr / de)
# =============================================================================
TR <- list(
  es = list(
    subtitle="Análisis estadístico y predictivo · Temporada 2022-2023", sec_label="Secciones",
    nav_dc="⚽ Simulador Dixon-Coles", nav_live="🔮 Predicción en Vivo",
    nav_pca="🗺️ Explorador Táctico", nav_opt="⏱️ Optimizador Técnico",
    home_team="🏠 Equipo local", away_team="✈️ Equipo visitante",
    dc_info="El modelo <b>Dixon-Coles</b> asigna a cada equipo fuerza de <b>ataque</b> y <b>debilidad defensiva</b>, con un factor de <b>ventaja local</b> (×1.34). Con los goles esperados (λ) de ambos se obtienen las probabilidades 1/X/2.",
    dc_lambda_local="λ %s (local)", dc_lambda_visit="λ %s (visitante)", dc_score="Marcador más probable",
    dc_lambda_title="Goles esperados (λ)", dc_probs_title="Probabilidad de resultado (1 / X / 2)", prob_y="Probabilidad (%)",
    res_local="Victoria local", res_draw="Empate", res_visit="Victoria visitante",
    mc_header="Probabilidades de toda la liga (Monte Carlo, 10.000 temporadas)", mc_title="Probabilidad de ser campeón (%)",
    live_info="Introduce las estadísticas del partido y obtén en tiempo real la <b>probabilidad de victoria local</b> (logística, modelo 7.3) y los <b>goles esperados del local</b> (Poisson, modelo 7.2).",
    sl_sot="Tiros a puerta", sl_corners="Córners", sl_fouls="Faltas",
    gauge_title="Probabilidad de victoria local", live_goals="Goles esperados (local)",
    live_read="Lectura del partido", live_shotsdiff="Diferencia de tiros a puerta",
    v_fav_local="Favorito claro el local", v_balanced="Partido equilibrado", v_fav_visit="Favorito el visitante / empate",
    or_header="Coeficientes del modelo logístico (Odds Ratios)", or_var="Variable", or_coef="Coeficiente (log-odds)", or_or="Odds Ratio exp(b)",
    pca_info="Cada equipo se proyecta sobre las dos primeras componentes principales. <b>PC1 (%.1f%%)</b> = dominio ofensivo (derecha = más ataque); <b>PC2 (%.1f%%)</b> = fricción/faltas. El K-means forma tres perfiles (silueta media = <b>%.3f</b>). Pasa el ratón por cada punto.",
    pca_show="Mostrar perfiles", elite="Élite", mid="Zona Media", strug="En apuros",
    pca_x="PC1 · Dominio ofensivo (%.1f%%)", pca_y="PC2 · Fricción / faltas (%.1f%%)", pca_legend="Perfil táctico",
    pca_groups_header="Composición de los clústeres", perfil_lbl="Perfil", gf_lbl="Goles a favor", gc_lbl="Goles en contra",
    pca_warn="Selecciona al menos un perfil.",
    opt_info="Reparte las horas de entrenamiento entre <b>ataque</b> (x₁) y <b>defensa</b> (x₂) maximizando <b>Z = 0.05·x₁ + 0.03·x₂</b>, con un máximo de horas y mínimos por área. Resuelto con lpSolve.",
    sl_hours="⏱️ Horas totales disponibles", sl_minatk="🎯 Mínimo de horas de ataque", sl_mindef="🛡️ Mínimo de horas de defensa",
    opt_x1="🎯 Horas de ataque (x₁)", opt_x2="🛡️ Horas de defensa (x₂)", opt_z="📈 Rendimiento máximo (Z)",
    opt_bar_title="Reparto óptimo de horas", opt_bar_y="Horas", opt_pie_title="Distribución del esfuerzo", pie_atk="Ataque", pie_def="Defensa",
    opt_note="Con <b>%d h</b> disponibles, el óptimo dedica <b>%.0f h al ataque</b> y <b>%.0f h a la defensa</b>, alcanzando <b>Z = %.2f</b>. Como el ataque rinde más (0.05 vs 0.03), conviene cubrir el mínimo defensivo y volcar el resto en ataque. Precio sombra de las horas: <b>%.2f</b> (cada hora extra añade %.2f).",
    opt_infeasible="Problema no factible: los mínimos (%d + %d = %d h) superan las %d h totales.",
    dc_distinct_warn="Selecciona dos equipos distintos."
  ),
  pt = list(
    subtitle="Análise estatística e preditiva · Época 2022-2023", sec_label="Secções",
    nav_dc="⚽ Simulador Dixon-Coles", nav_live="🔮 Previsão ao Vivo",
    nav_pca="🗺️ Explorador Tático", nav_opt="⏱️ Otimizador Técnico",
    home_team="🏠 Equipa da casa", away_team="✈️ Equipa visitante",
    dc_info="O modelo <b>Dixon-Coles</b> atribui a cada equipa força de <b>ataque</b> e <b>fragilidade defensiva</b>, com um fator de <b>vantagem caseira</b> (×1.34). Com os golos esperados (λ) de ambas obtêm-se as probabilidades 1/X/2.",
    dc_lambda_local="λ %s (casa)", dc_lambda_visit="λ %s (visitante)", dc_score="Resultado mais provável",
    dc_lambda_title="Golos esperados (λ)", dc_probs_title="Probabilidade de resultado (1 / X / 2)", prob_y="Probabilidade (%)",
    res_local="Vitória da casa", res_draw="Empate", res_visit="Vitória visitante",
    mc_header="Probabilidades de toda a liga (Monte Carlo, 10.000 épocas)", mc_title="Probabilidade de ser campeão (%)",
    live_info="Introduz as estatísticas do jogo e obtém em tempo real a <b>probabilidade de vitória da casa</b> (logística, modelo 7.3) e os <b>golos esperados da casa</b> (Poisson, modelo 7.2).",
    sl_sot="Remates à baliza", sl_corners="Cantos", sl_fouls="Faltas",
    gauge_title="Probabilidade de vitória da casa", live_goals="Golos esperados (casa)",
    live_read="Leitura do jogo", live_shotsdiff="Diferença de remates à baliza",
    v_fav_local="Casa claramente favorita", v_balanced="Jogo equilibrado", v_fav_visit="Visitante favorito / empate",
    or_header="Coeficientes do modelo logístico (Odds Ratios)", or_var="Variável", or_coef="Coeficiente (log-odds)", or_or="Odds Ratio exp(b)",
    pca_info="Cada equipa é projetada sobre as duas primeiras componentes principais. <b>PC1 (%.1f%%)</b> = domínio ofensivo (direita = mais ataque); <b>PC2 (%.1f%%)</b> = fricção/faltas. O K-means forma três perfis (silhueta média = <b>%.3f</b>). Passa o rato por cada ponto.",
    pca_show="Mostrar perfis", elite="Elite", mid="Zona Intermédia", strug="Em apuros",
    pca_x="PC1 · Domínio ofensivo (%.1f%%)", pca_y="PC2 · Fricção / faltas (%.1f%%)", pca_legend="Perfil tático",
    pca_groups_header="Composição dos clusters", perfil_lbl="Perfil", gf_lbl="Golos marcados", gc_lbl="Golos sofridos",
    pca_warn="Seleciona pelo menos um perfil.",
    opt_info="Distribui as horas de treino entre <b>ataque</b> (x₁) e <b>defesa</b> (x₂) maximizando <b>Z = 0.05·x₁ + 0.03·x₂</b>, com um máximo de horas e mínimos por área. Resolvido com lpSolve.",
    sl_hours="⏱️ Horas totais disponíveis", sl_minatk="🎯 Mínimo de horas de ataque", sl_mindef="🛡️ Mínimo de horas de defesa",
    opt_x1="🎯 Horas de ataque (x₁)", opt_x2="🛡️ Horas de defesa (x₂)", opt_z="📈 Rendimento máximo (Z)",
    opt_bar_title="Distribuição ótima de horas", opt_bar_y="Horas", opt_pie_title="Distribuição do esforço", pie_atk="Ataque", pie_def="Defesa",
    opt_note="Com <b>%d h</b> disponíveis, o ótimo dedica <b>%.0f h ao ataque</b> e <b>%.0f h à defesa</b>, alcançando <b>Z = %.2f</b>. Como o ataque rende mais (0.05 vs 0.03), convém cobrir o mínimo defensivo e colocar o resto no ataque. Preço-sombra das horas: <b>%.2f</b> (cada hora extra acrescenta %.2f).",
    opt_infeasible="Problema inviável: os mínimos (%d + %d = %d h) excedem as %d h totais.",
    dc_distinct_warn="Seleciona duas equipas diferentes."
  ),
  en = list(
    subtitle="Statistical and predictive analysis · 2022-2023 season", sec_label="Sections",
    nav_dc="⚽ Dixon-Coles Simulator", nav_live="🔮 Live Prediction",
    nav_pca="🗺️ Tactical Explorer", nav_opt="⏱️ Technical Optimizer",
    home_team="🏠 Home team", away_team="✈️ Away team",
    dc_info="The <b>Dixon-Coles</b> model assigns each team an <b>attack</b> strength and a <b>defensive weakness</b>, plus a <b>home-advantage</b> factor (×1.34). The expected goals (λ) of both yield the 1/X/2 probabilities.",
    dc_lambda_local="λ %s (home)", dc_lambda_visit="λ %s (away)", dc_score="Most likely score",
    dc_lambda_title="Expected goals (λ)", dc_probs_title="Outcome probability (1 / X / 2)", prob_y="Probability (%)",
    res_local="Home win", res_draw="Draw", res_visit="Away win",
    mc_header="League-wide probabilities (Monte Carlo, 10,000 seasons)", mc_title="Probability of winning the title (%)",
    live_info="Enter the match statistics and get, in real time, the <b>home-win probability</b> (logistic regression, model 7.3) and the <b>expected home goals</b> (Poisson regression, model 7.2).",
    sl_sot="Shots on target", sl_corners="Corners", sl_fouls="Fouls",
    gauge_title="Home-win probability", live_goals="Expected goals (home)",
    live_read="Match read", live_shotsdiff="Shots-on-target difference",
    v_fav_local="Home clear favorite", v_balanced="Balanced match", v_fav_visit="Away favorite / draw",
    or_header="Logistic model coefficients (Odds Ratios)", or_var="Variable", or_coef="Coefficient (log-odds)", or_or="Odds Ratio exp(b)",
    pca_info="Each team is projected onto the first two principal components. <b>PC1 (%.1f%%)</b> = offensive dominance (right = more attack); <b>PC2 (%.1f%%)</b> = friction/fouls. K-means forms three profiles (mean silhouette = <b>%.3f</b>). Hover over each point.",
    pca_show="Show profiles", elite="Elite", mid="Mid-table", strug="Strugglers",
    pca_x="PC1 · Offensive dominance (%.1f%%)", pca_y="PC2 · Friction / fouls (%.1f%%)", pca_legend="Tactical profile",
    pca_groups_header="Cluster composition", perfil_lbl="Profile", gf_lbl="Goals for", gc_lbl="Goals against",
    pca_warn="Select at least one profile.",
    opt_info="Split training hours between <b>attack</b> (x₁) and <b>defense</b> (x₂) maximizing <b>Z = 0.05·x₁ + 0.03·x₂</b>, subject to a maximum of hours and per-area minimums. Solved with lpSolve.",
    sl_hours="⏱️ Total hours available", sl_minatk="🎯 Minimum attack hours", sl_mindef="🛡️ Minimum defense hours",
    opt_x1="🎯 Attack hours (x₁)", opt_x2="🛡️ Defense hours (x₂)", opt_z="📈 Maximum performance (Z)",
    opt_bar_title="Optimal hour allocation", opt_bar_y="Hours", opt_pie_title="Effort distribution", pie_atk="Attack", pie_def="Defense",
    opt_note="With <b>%d h</b> available, the optimum assigns <b>%.0f h to attack</b> and <b>%.0f h to defense</b>, reaching <b>Z = %.2f</b>. Since attack yields more (0.05 vs 0.03), it is best to cover the defensive minimum and pour the rest into attack. Shadow price of hours: <b>%.2f</b> (each extra hour adds %.2f).",
    opt_infeasible="Infeasible problem: the minimums (%d + %d = %d h) exceed the %d total hours.",
    dc_distinct_warn="Select two different teams."
  ),
  it = list(
    subtitle="Analisi statistica e predittiva · Stagione 2022-2023", sec_label="Sezioni",
    nav_dc="⚽ Simulatore Dixon-Coles", nav_live="🔮 Previsione dal Vivo",
    nav_pca="🗺️ Esploratore Tattico", nav_opt="⏱️ Ottimizzatore Tecnico",
    home_team="🏠 Squadra di casa", away_team="✈️ Squadra ospite",
    dc_info="Il modello <b>Dixon-Coles</b> assegna a ogni squadra una forza di <b>attacco</b> e una <b>debolezza difensiva</b>, con un fattore di <b>vantaggio casalingo</b> (×1.34). Dai gol attesi (λ) di entrambe si ottengono le probabilità 1/X/2.",
    dc_lambda_local="λ %s (casa)", dc_lambda_visit="λ %s (ospite)", dc_score="Risultato più probabile",
    dc_lambda_title="Gol attesi (λ)", dc_probs_title="Probabilità di esito (1 / X / 2)", prob_y="Probabilità (%)",
    res_local="Vittoria casa", res_draw="Pareggio", res_visit="Vittoria ospite",
    mc_header="Probabilità dell'intero campionato (Monte Carlo, 10.000 stagioni)", mc_title="Probabilità di vincere il titolo (%)",
    live_info="Inserisci le statistiche della partita e ottieni in tempo reale la <b>probabilità di vittoria casalinga</b> (regressione logistica, modello 7.3) e i <b>gol attesi della squadra di casa</b> (regressione di Poisson, modello 7.2).",
    sl_sot="Tiri in porta", sl_corners="Calci d'angolo", sl_fouls="Falli",
    gauge_title="Probabilità di vittoria casalinga", live_goals="Gol attesi (casa)",
    live_read="Lettura della partita", live_shotsdiff="Differenza di tiri in porta",
    v_fav_local="Casa nettamente favorita", v_balanced="Partita equilibrata", v_fav_visit="Ospite favorito / pareggio",
    or_header="Coefficienti del modello logistico (Odds Ratio)", or_var="Variabile", or_coef="Coefficiente (log-odds)", or_or="Odds Ratio exp(b)",
    pca_info="Ogni squadra è proiettata sulle prime due componenti principali. <b>PC1 (%.1f%%)</b> = dominio offensivo (destra = più attacco); <b>PC2 (%.1f%%)</b> = attrito/falli. Il K-means forma tre profili (silhouette media = <b>%.3f</b>). Passa il mouse su ogni punto.",
    pca_show="Mostra profili", elite="Élite", mid="Zona Media", strug="In difficoltà",
    pca_x="PC1 · Dominio offensivo (%.1f%%)", pca_y="PC2 · Attrito / falli (%.1f%%)", pca_legend="Profilo tattico",
    pca_groups_header="Composizione dei cluster", perfil_lbl="Profilo", gf_lbl="Gol fatti", gc_lbl="Gol subiti",
    pca_warn="Seleziona almeno un profilo.",
    opt_info="Distribuisci le ore di allenamento tra <b>attacco</b> (x₁) e <b>difesa</b> (x₂) massimizzando <b>Z = 0.05·x₁ + 0.03·x₂</b>, con un massimo di ore e minimi per area. Risolto con lpSolve.",
    sl_hours="⏱️ Ore totali disponibili", sl_minatk="🎯 Ore minime di attacco", sl_mindef="🛡️ Ore minime di difesa",
    opt_x1="🎯 Ore di attacco (x₁)", opt_x2="🛡️ Ore di difesa (x₂)", opt_z="📈 Rendimento massimo (Z)",
    opt_bar_title="Ripartizione ottimale delle ore", opt_bar_y="Ore", opt_pie_title="Distribuzione dello sforzo", pie_atk="Attacco", pie_def="Difesa",
    opt_note="Con <b>%d h</b> disponibili, l'ottimo dedica <b>%.0f h all'attacco</b> e <b>%.0f h alla difesa</b>, raggiungendo <b>Z = %.2f</b>. Poiché l'attacco rende di più (0.05 vs 0.03), conviene coprire il minimo difensivo e destinare il resto all'attacco. Prezzo ombra delle ore: <b>%.2f</b> (ogni ora extra aggiunge %.2f).",
    opt_infeasible="Problema non ammissibile: i minimi (%d + %d = %d h) superano le %d h totali.",
    dc_distinct_warn="Seleziona due squadre diverse."
  ),
  fr = list(
    subtitle="Analyse statistique et prédictive · Saison 2022-2023", sec_label="Sections",
    nav_dc="⚽ Simulateur Dixon-Coles", nav_live="🔮 Prédiction en Direct",
    nav_pca="🗺️ Explorateur Tactique", nav_opt="⏱️ Optimiseur Technique",
    home_team="🏠 Équipe à domicile", away_team="✈️ Équipe à l'extérieur",
    dc_info="Le modèle <b>Dixon-Coles</b> attribue à chaque équipe une force d'<b>attaque</b> et une <b>faiblesse défensive</b>, avec un facteur d'<b>avantage à domicile</b> (×1.34). Les buts attendus (λ) des deux donnent les probabilités 1/N/2.",
    dc_lambda_local="λ %s (domicile)", dc_lambda_visit="λ %s (extérieur)", dc_score="Score le plus probable",
    dc_lambda_title="Buts attendus (λ)", dc_probs_title="Probabilité de résultat (1 / N / 2)", prob_y="Probabilité (%)",
    res_local="Victoire domicile", res_draw="Match nul", res_visit="Victoire extérieur",
    mc_header="Probabilités sur toute la ligue (Monte-Carlo, 10 000 saisons)", mc_title="Probabilité d'être champion (%)",
    live_info="Saisissez les statistiques du match et obtenez en temps réel la <b>probabilité de victoire à domicile</b> (régression logistique, modèle 7.3) et les <b>buts attendus à domicile</b> (régression de Poisson, modèle 7.2).",
    sl_sot="Tirs cadrés", sl_corners="Corners", sl_fouls="Fautes",
    gauge_title="Probabilité de victoire à domicile", live_goals="Buts attendus (domicile)",
    live_read="Lecture du match", live_shotsdiff="Différence de tirs cadrés",
    v_fav_local="Domicile nettement favori", v_balanced="Match équilibré", v_fav_visit="Extérieur favori / nul",
    or_header="Coefficients du modèle logistique (Odds Ratios)", or_var="Variable", or_coef="Coefficient (log-odds)", or_or="Odds Ratio exp(b)",
    pca_info="Chaque équipe est projetée sur les deux premières composantes principales. <b>PC1 (%.1f%%)</b> = domination offensive (droite = plus d'attaque) ; <b>PC2 (%.1f%%)</b> = friction/fautes. Le K-means forme trois profils (silhouette moyenne = <b>%.3f</b>). Survolez chaque point.",
    pca_show="Afficher les profils", elite="Élite", mid="Milieu de tableau", strug="En difficulté",
    pca_x="PC1 · Domination offensive (%.1f%%)", pca_y="PC2 · Friction / fautes (%.1f%%)", pca_legend="Profil tactique",
    pca_groups_header="Composition des clusters", perfil_lbl="Profil", gf_lbl="Buts pour", gc_lbl="Buts contre",
    pca_warn="Sélectionnez au moins un profil.",
    opt_info="Répartissez les heures d'entraînement entre <b>attaque</b> (x₁) et <b>défense</b> (x₂) en maximisant <b>Z = 0.05·x₁ + 0.03·x₂</b>, sous un maximum d'heures et des minimums par domaine. Résolu avec lpSolve.",
    sl_hours="⏱️ Heures totales disponibles", sl_minatk="🎯 Heures minimales d'attaque", sl_mindef="🛡️ Heures minimales de défense",
    opt_x1="🎯 Heures d'attaque (x₁)", opt_x2="🛡️ Heures de défense (x₂)", opt_z="📈 Performance maximale (Z)",
    opt_bar_title="Répartition optimale des heures", opt_bar_y="Heures", opt_pie_title="Répartition de l'effort", pie_atk="Attaque", pie_def="Défense",
    opt_note="Avec <b>%d h</b> disponibles, l'optimum consacre <b>%.0f h à l'attaque</b> et <b>%.0f h à la défense</b>, atteignant <b>Z = %.2f</b>. L'attaque rapportant davantage (0.05 vs 0.03), il vaut mieux couvrir le minimum défensif et mettre le reste en attaque. Prix dual des heures : <b>%.2f</b> (chaque heure supplémentaire ajoute %.2f).",
    opt_infeasible="Problème non réalisable : les minimums (%d + %d = %d h) dépassent les %d h totales.",
    dc_distinct_warn="Sélectionnez deux équipes différentes."
  ),
  de = list(
    subtitle="Statistische und prädiktive Analyse · Saison 2022-2023", sec_label="Abschnitte",
    nav_dc="⚽ Dixon-Coles-Simulator", nav_live="🔮 Live-Vorhersage",
    nav_pca="🗺️ Taktischer Explorer", nav_opt="⏱️ Technischer Optimierer",
    home_team="🏠 Heimmannschaft", away_team="✈️ Gastmannschaft",
    dc_info="Das <b>Dixon-Coles</b>-Modell weist jedem Team eine <b>Angriffsstärke</b> und eine <b>Abwehrschwäche</b> zu, plus einen <b>Heimvorteil</b>-Faktor (×1.34). Aus den erwarteten Toren (λ) beider ergeben sich die 1/X/2-Wahrscheinlichkeiten.",
    dc_lambda_local="λ %s (Heim)", dc_lambda_visit="λ %s (Gast)", dc_score="Wahrscheinlichstes Ergebnis",
    dc_lambda_title="Erwartete Tore (λ)", dc_probs_title="Ergebniswahrscheinlichkeit (1 / X / 2)", prob_y="Wahrscheinlichkeit (%)",
    res_local="Heimsieg", res_draw="Unentschieden", res_visit="Auswärtssieg",
    mc_header="Liga-weite Wahrscheinlichkeiten (Monte-Carlo, 10.000 Saisons)", mc_title="Meisterschaftswahrscheinlichkeit (%)",
    live_info="Gib die Spielstatistiken ein und erhalte in Echtzeit die <b>Heimsieg-Wahrscheinlichkeit</b> (logistische Regression, Modell 7.3) und die <b>erwarteten Heimtore</b> (Poisson-Regression, Modell 7.2).",
    sl_sot="Torschüsse", sl_corners="Eckbälle", sl_fouls="Fouls",
    gauge_title="Heimsieg-Wahrscheinlichkeit", live_goals="Erwartete Tore (Heim)",
    live_read="Spieleinschätzung", live_shotsdiff="Differenz der Torschüsse",
    v_fav_local="Heim klarer Favorit", v_balanced="Ausgeglichenes Spiel", v_fav_visit="Gast favorisiert / Unentschieden",
    or_header="Koeffizienten des logistischen Modells (Odds Ratios)", or_var="Variable", or_coef="Koeffizient (Log-Odds)", or_or="Odds Ratio exp(b)",
    pca_info="Jedes Team wird auf die ersten beiden Hauptkomponenten projiziert. <b>PC1 (%.1f%%)</b> = Offensivdominanz (rechts = mehr Angriff); <b>PC2 (%.1f%%)</b> = Reibung/Fouls. K-means bildet drei Profile (mittlere Silhouette = <b>%.3f</b>). Fahre über jeden Punkt.",
    pca_show="Profile anzeigen", elite="Spitze", mid="Mittelfeld", strug="Abstiegskampf",
    pca_x="PC1 · Offensivdominanz (%.1f%%)", pca_y="PC2 · Reibung / Fouls (%.1f%%)", pca_legend="Taktisches Profil",
    pca_groups_header="Cluster-Zusammensetzung", perfil_lbl="Profil", gf_lbl="Tore für", gc_lbl="Gegentore",
    pca_warn="Wähle mindestens ein Profil.",
    opt_info="Verteile die Trainingsstunden zwischen <b>Angriff</b> (x₁) und <b>Abwehr</b> (x₂) und maximiere <b>Z = 0.05·x₁ + 0.03·x₂</b>, mit einer Stundenobergrenze und Mindestwerten pro Bereich. Gelöst mit lpSolve.",
    sl_hours="⏱️ Verfügbare Gesamtstunden", sl_minatk="🎯 Mindeststunden Angriff", sl_mindef="🛡️ Mindeststunden Abwehr",
    opt_x1="🎯 Angriffsstunden (x₁)", opt_x2="🛡️ Abwehrstunden (x₂)", opt_z="📈 Maximale Leistung (Z)",
    opt_bar_title="Optimale Stundenverteilung", opt_bar_y="Stunden", opt_pie_title="Aufwandsverteilung", pie_atk="Angriff", pie_def="Abwehr",
    opt_note="Mit <b>%d h</b> verfügbar weist das Optimum <b>%.0f h dem Angriff</b> und <b>%.0f h der Abwehr</b> zu und erreicht <b>Z = %.2f</b>. Da der Angriff mehr bringt (0,05 vs 0,03), sollte man das Abwehr-Minimum abdecken und den Rest in den Angriff stecken. Schattenpreis der Stunden: <b>%.2f</b> (jede zusätzliche Stunde bringt %.2f).",
    opt_infeasible="Unzulässiges Problem: die Mindestwerte (%d + %d = %d h) überschreiten die %d Gesamtstunden.",
    dc_distinct_warn="Wähle zwei verschiedene Teams."
  )
)
LANGS <- c("Español" = "es", "Português" = "pt", "English" = "en",
           "Italiano" = "it", "Français" = "fr", "Deutsch" = "de")

# =============================================================================
#  FUNCIONES AUXILIARES (lógica idéntica al proyecto)
# =============================================================================
matriz_resultado <- function(lh, la, maxg = 10) {
  ph <- dpois(0:maxg, lh); pa <- dpois(0:maxg, la); M <- outer(ph, pa)
  ij <- which(M == max(M), arr.ind = TRUE)[1, ]
  list(local = sum(M[lower.tri(M)]), empate = sum(diag(M)),
       visit = sum(M[upper.tri(M)]), score = c(ij[1] - 1, ij[2] - 1))
}
resolver_lp <- function(horas, min_ata, min_def, r_ata = 0.05, r_def = 0.03) {
  lp <- make.lp(0, 2); invisible(lp.control(lp, sense = "max"))
  set.objfn(lp, c(r_ata, r_def))
  add.constraint(lp, c(1, 1), "<=", horas)
  add.constraint(lp, c(1, 0), ">=", min_ata)
  add.constraint(lp, c(0, 1), ">=", min_def)
  if (solve(lp) != 0) return(NULL)
  precio <- tryCatch(get.sensitivity.rhs(lp)$duals[1], error = function(e) NA_real_)
  list(x1 = get.variables(lp)[1], x2 = get.variables(lp)[2],
       z = get.objective(lp), precio = precio)
}
metric_card <- function(label, value) {
  div(class = "metric-card",
      div(class = "metric-label", label),
      div(class = "metric-value", value))
}

css <- "
.metric-card{background:#fff;border:1px solid #e6e6e6;border-radius:12px;
  padding:16px;text-align:center;box-shadow:0 1px 3px rgba(0,0,0,.05);margin-bottom:12px;}
.metric-label{color:#6c757d;font-size:.9rem;margin-bottom:4px;}
.metric-value{color:#1d3557;font-size:1.7rem;font-weight:700;line-height:1.2;}
.info-box{background:#e7f1ff;border-left:4px solid #1d3557;border-radius:6px;
  padding:12px 14px;margin-bottom:14px;color:#0b2545;font-size:.9rem;}
"

# =============================================================================
#  INTERFAZ DE USUARIO
# =============================================================================
ui <- fluidPage(
  theme = bs_theme(version = 5, bootswatch = "flatly", primary = "#1d3557"),
  tags$head(
    tags$style(HTML(css)),
    # Metaetiquetas mínimas para que WhatsApp/Telegram/Slack muestren una
    # tarjeta limpia (título + descripción) en vez de leer el HTML a lo loco.
    tags$meta(property = "og:title",       content = "⚽ Premier League 2022-2023 · Analítica"),
    tags$meta(property = "og:description", content = "Simulador Dixon-Coles, predicción en vivo, explorador táctico (PCA + K-means) y optimizador técnico. App interactiva en Shiny."),
    tags$meta(name = "description",        content = "Análisis estadístico y predictivo de la Premier League 2022-2023: Dixon-Coles, regresión de Poisson y logística, PCA, K-means y programación lineal.")
  ),
  titlePanel(
    title = div(
      span("⚽ Premier League 2022-2023", style = "font-weight:700;"),
      uiOutput("subtitle", inline = TRUE)
    ),
    windowTitle = "Premier League 2022-2023 · Analítica"  # texto plano para <title> y la pestaña del navegador
  ),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      selectInput("lang", "🌐 Idioma / Language", choices = LANGS, selected = "es"),
      radioButtons("sec", TR$es$sec_label,
                   choiceNames  = list(TR$es$nav_dc, TR$es$nav_live, TR$es$nav_pca, TR$es$nav_opt),
                   choiceValues = c("dc", "live", "pca", "opt"), selected = "dc"),
      tags$hr(),
      uiOutput("controls")
    ),
    mainPanel(width = 9, uiOutput("main"))
  )
)

# =============================================================================
#  SERVIDOR
# =============================================================================
server <- function(input, output, session) {

  L <- reactive(TR[[input$lang %||% "es"]])

  output$subtitle <- renderUI(span(style = "color:#6c757d;font-size:1rem;margin-left:10px;",
                                   paste("·", L()$subtitle)))

  # Actualizar etiquetas del menú lateral al cambiar de idioma (sin perder selección)
  observeEvent(input$lang, {
    tr <- L()
    updateRadioButtons(session, "sec", label = tr$sec_label,
                       choiceNames  = list(tr$nav_dc, tr$nav_live, tr$nav_pca, tr$nav_opt),
                       choiceValues = c("dc", "live", "pca", "opt"),
                       selected = input$sec %||% "dc")
  }, ignoreInit = TRUE)

  # --------- Controles por sección (se traducen; conservan valor al cambiar idioma) ---------
  output$controls <- renderUI({
    tr <- L(); sec <- input$sec %||% "dc"
    if (sec == "dc") {
      tagList(
        div(class = "info-box", HTML(tr$dc_info)),
        selectInput("local", tr$home_team, choices = equipos, selected = isolate(input$local) %||% "Man City"),
        selectInput("visit", tr$away_team, choices = equipos, selected = isolate(input$visit) %||% "Liverpool")
      )
    } else if (sec == "live") {
      tagList(
        div(class = "info-box", HTML(tr$live_info)),
        tags$b(tr$home_team),
        sliderInput("hst", paste0(tr$sl_sot, " (HST)"), 0, 20, isolate(input$hst) %||% 6),
        sliderInput("hc",  paste0(tr$sl_corners, " (HC)"), 0, 15, isolate(input$hc) %||% 5),
        sliderInput("hf",  paste0(tr$sl_fouls, " (HF)"), 0, 25, isolate(input$hf) %||% 11),
        tags$hr(), tags$b(tr$away_team),
        sliderInput("ast", paste0(tr$sl_sot, " (AST)"), 0, 20, isolate(input$ast) %||% 4),
        sliderInput("ac",  paste0(tr$sl_corners, " (AC)"), 0, 15, isolate(input$ac) %||% 5),
        sliderInput("af",  paste0(tr$sl_fouls, " (AF)"), 0, 25, isolate(input$af) %||% 11)
      )
    } else if (sec == "pca") {
      tagList(
        div(class = "info-box", HTML(sprintf(tr$pca_info, varexp[1], varexp[2], sil))),
        checkboxGroupInput("perfiles", tr$pca_show,
                           choiceNames  = list(tr$elite, tr$mid, tr$strug),
                           choiceValues = c("elite", "mid", "strug"),
                           selected = isolate(input$perfiles) %||% c("elite", "mid", "strug"))
      )
    } else {
      tagList(
        div(class = "info-box", HTML(tr$opt_info)),
        sliderInput("horas",   tr$sl_hours,  30, 150, isolate(input$horas) %||% 100, step = 5),
        sliderInput("min_ata", tr$sl_minatk,  0,  80, isolate(input$min_ata) %||% 20, step = 5),
        sliderInput("min_def", tr$sl_mindef,  0,  80, isolate(input$min_def) %||% 30, step = 5)
      )
    }
  })

  # --------- Cuerpo principal por sección ---------
  output$main <- renderUI({
    tr <- L(); sec <- input$sec %||% "dc"
    if (sec == "dc") {
      tagList(h3(tr$nav_dc), uiOutput("dc_metrics"),
              fluidRow(column(6, plotlyOutput("dc_lambda", height = "380px")),
                       column(6, plotlyOutput("dc_probs",  height = "380px"))),
              uiOutput("dc_mc_ui"))
    } else if (sec == "live") {
      tagList(h3(tr$nav_live),
              fluidRow(column(7, plotlyOutput("live_gauge", height = "360px")),
                       column(5, uiOutput("live_metrics"))),
              tags$h5(tr$or_header), tableOutput("live_or"))
    } else if (sec == "pca") {
      tagList(h3(tr$nav_pca), plotlyOutput("pca_plot", height = "560px"),
              tags$hr(), tags$h5(tr$pca_groups_header), uiOutput("pca_groups"))
    } else {
      tagList(h3(tr$nav_opt), uiOutput("opt_metrics"),
              fluidRow(column(6, plotlyOutput("opt_bar", height = "380px")),
                       column(6, plotlyOutput("opt_pie", height = "380px"))),
              uiOutput("opt_note"))
    }
  })

  # ---------------- 1) Dixon-Coles ----------------
  dc <- reactive({
    req(input$local, input$visit)
    validate(need(input$local != input$visit, L()$dc_distinct_warn))
    lh <- exp(ventaja + ataque[[input$local]] - defensa[[input$visit]])
    la <- exp(ataque[[input$visit]] - defensa[[input$local]])
    r <- matriz_resultado(lh, la)
    list(lh = lh, la = la, p_local = r$local, p_empate = r$empate, p_visit = r$visit, score = r$score)
  })

  output$dc_metrics <- renderUI({
    d <- dc(); tr <- L()
    fluidRow(
      column(4, metric_card(sprintf(tr$dc_lambda_local, input$local), sprintf("%.2f", d$lh))),
      column(4, metric_card(sprintf(tr$dc_lambda_visit, input$visit), sprintf("%.2f", d$la))),
      column(4, metric_card(tr$dc_score, sprintf("%d \u2013 %d", d$score[1], d$score[2])))
    )
  })
  output$dc_lambda <- renderPlotly({
    d <- dc(); tr <- L()
    plot_ly(x = c(input$local, input$visit), y = c(d$lh, d$la), type = "bar",
            marker = list(color = c("#1d3557", "#e63946")),
            text = sprintf("%.2f", c(d$lh, d$la)), textposition = "outside") %>%
      layout(title = tr$dc_lambda_title, showlegend = FALSE,
             yaxis = list(title = "\u03bb", range = c(0, max(d$lh, d$la) * 1.25)),
             xaxis = list(title = ""), plot_bgcolor = "white", paper_bgcolor = "white")
  })
  output$dc_probs <- renderPlotly({
    d <- dc(); tr <- L()
    plot_ly(x = c(tr$res_local, tr$res_draw, tr$res_visit),
            y = c(d$p_local, d$p_empate, d$p_visit) * 100, type = "bar",
            marker = list(color = c("#2ca02c", "#7f7f7f", "#d62728")),
            text = sprintf("%.1f%%", c(d$p_local, d$p_empate, d$p_visit) * 100),
            textposition = "outside") %>%
      layout(title = tr$dc_probs_title, showlegend = FALSE,
             yaxis = list(title = tr$prob_y, range = c(0, 100)),
             xaxis = list(title = ""), plot_bgcolor = "white", paper_bgcolor = "white")
  })
  output$dc_mc_ui <- renderUI({
    if (is.null(mc)) return(NULL)
    tagList(tags$hr(), tags$h5(L()$mc_header), plotlyOutput("dc_mc_plot", height = "380px"))
  })
  output$dc_mc_plot <- renderPlotly({
    req(!is.null(mc)); m <- mc[order(-mc$Campeon), ]; m <- m[m$Campeon > 0, ]
    plot_ly(x = m$Equipo, y = m$Campeon, type = "bar", marker = list(color = "#1d3557")) %>%
      layout(title = L()$mc_title, yaxis = list(title = "%"),
             xaxis = list(title = "", categoryorder = "total descending"),
             plot_bgcolor = "white", paper_bgcolor = "white")
  })

  # ---------------- 2) Predicción en vivo ----------------
  live <- reactive({
    req(input$hst, input$ast, input$hc, input$ac, input$hf, input$af)
    fila <- data.frame(HST = input$hst, AST = input$ast, HC = input$hc,
                       AC = input$ac, HF = input$hf, AF = input$af)
    list(prob  = as.numeric(predict(m_log,  fila, type = "response")),
         goles = as.numeric(predict(m_pois, fila, type = "response")))
  })
  output$live_gauge <- renderPlotly({
    p <- live()$prob * 100; tr <- L()
    plot_ly(type = "indicator", mode = "gauge+number", value = p, number = list(suffix = " %"),
            title = list(text = tr$gauge_title),
            gauge = list(axis = list(range = list(0, 100)), bar = list(color = "#1d3557"),
                         steps = list(list(range = c(0, 40), color = "#fde2e2"),
                                      list(range = c(40, 60), color = "#fff3cd"),
                                      list(range = c(60, 100), color = "#d8f3dc")),
                         threshold = list(line = list(color = "black", width = 3),
                                          thickness = 0.75, value = 50))) %>%
      layout(height = 360, margin = list(t = 60, b = 10))
  })
  output$live_metrics <- renderUI({
    l <- live(); tr <- L()
    v <- if (l$prob >= 0.6) tr$v_fav_local else if (l$prob >= 0.4) tr$v_balanced else tr$v_fav_visit
    tagList(metric_card(tr$live_goals, sprintf("%.2f", l$goles)),
            metric_card(tr$live_read, v),
            metric_card(tr$live_shotsdiff, sprintf("%+d", input$hst - input$ast)))
  })
  output$live_or <- renderTable({
    tr <- L(); co <- coef(m_log)
    df <- data.frame(names(co), round(co, 3), round(exp(co), 3))
    names(df) <- c(tr$or_var, tr$or_coef, tr$or_or); df
  }, striped = TRUE, bordered = TRUE, width = "100%")

  # ---------------- 3) Explorador táctico ----------------
  output$pca_plot <- renderPlotly({
    tr <- L(); df <- df_pca[df_pca$Cl %in% (input$perfiles %||% character(0)), , drop = FALSE]
    validate(need(nrow(df) > 0, tr$pca_warn))
    nom <- c(elite = tr$elite, mid = tr$mid, strug = tr$strug)
    df$Perfil <- factor(nom[df$Cl], levels = c(tr$elite, tr$mid, tr$strug))
    cmap <- setNames(c("#2ca02c", "#1f77b4", "#d62728"), c(tr$elite, tr$mid, tr$strug))
    plot_ly(df, x = ~PC1, y = ~PC2, color = ~Perfil, colors = cmap,
            type = "scatter", mode = "markers+text",
            text = ~Equipo, textposition = "top center", textfont = list(size = 10),
            hoverinfo = "text",
            hovertext = ~paste0("<b>", Equipo, "</b><br>", tr$perfil_lbl, ": ", Perfil,
                                "<br>", tr$gf_lbl, ": ", GF, " | ", tr$gc_lbl, ": ", GC),
            marker = list(size = 12, line = list(width = 1, color = "white"))) %>%
      layout(xaxis = list(title = sprintf(tr$pca_x, varexp[1]), zeroline = TRUE),
             yaxis = list(title = sprintf(tr$pca_y, varexp[2]), zeroline = TRUE),
             legend = list(title = list(text = tr$pca_legend)),
             plot_bgcolor = "white", paper_bgcolor = "white")
  })
  output$pca_groups <- renderUI({
    tr <- L(); nom <- c(elite = tr$elite, mid = tr$mid, strug = tr$strug)
    cols <- lapply(c("elite", "mid", "strug"), function(code) {
      eqs <- df_pca$Equipo[df_pca$Cl == code]
      column(4, tags$b(sprintf("%s (%d)", nom[[code]], length(eqs))),
             tags$p(paste(eqs, collapse = ", ")))
    })
    do.call(fluidRow, cols)
  })

  # ---------------- 4) Optimizador técnico ----------------
  opt <- reactive({
    req(input$horas, input$min_ata, input$min_def)
    validate(need(input$min_ata + input$min_def <= input$horas,
                  sprintf(L()$opt_infeasible, input$min_ata, input$min_def,
                          input$min_ata + input$min_def, input$horas)))
    r <- resolver_lp(input$horas, input$min_ata, input$min_def)
    validate(need(!is.null(r), "—")); r
  })
  output$opt_metrics <- renderUI({
    o <- opt(); tr <- L()
    fluidRow(column(4, metric_card(tr$opt_x1, sprintf("%.0f h", o$x1))),
             column(4, metric_card(tr$opt_x2, sprintf("%.0f h", o$x2))),
             column(4, metric_card(tr$opt_z, sprintf("%.2f", o$z))))
  })
  output$opt_bar <- renderPlotly({
    o <- opt(); tr <- L()
    plot_ly(x = c(paste0(tr$pie_atk, " (x\u2081)"), paste0(tr$pie_def, " (x\u2082)")),
            y = c(o$x1, o$x2), type = "bar", marker = list(color = c("#e63946", "#1d3557")),
            text = sprintf("%.0f h", c(o$x1, o$x2)), textposition = "outside") %>%
      layout(title = tr$opt_bar_title, showlegend = FALSE,
             yaxis = list(title = tr$opt_bar_y, range = c(0, input$horas * 1.1)),
             xaxis = list(title = ""), plot_bgcolor = "white", paper_bgcolor = "white")
  })
  output$opt_pie <- renderPlotly({
    o <- opt(); tr <- L()
    plot_ly(labels = c(tr$pie_atk, tr$pie_def), values = c(o$x1, o$x2), type = "pie",
            hole = 0.55, marker = list(colors = c("#e63946", "#1d3557"))) %>%
      layout(title = tr$opt_pie_title)
  })
  output$opt_note <- renderUI({
    o <- opt(); tr <- L()
    div(class = "info-box", style = "background:#d8f3dc;border-left-color:#2ca02c;",
        HTML(sprintf(tr$opt_note, input$horas, o$x1, o$x2, o$z, o$precio, o$precio)))
  })
}

# =============================================================================
shinyApp(ui = ui, server = server)
