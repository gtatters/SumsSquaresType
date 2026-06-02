library(shiny)

# https://hbctraining.github.io/Training-modules/RShiny/lessons/shinylive.html
# Run the shinylive::export line to populate the docs folder 
# so that shinylive works from github
#shinylive::export(appdir = "../SumsSquaresType/", destdir = "docs")
#httpuv::runStaticServer("docs/", port = 8008)

# ── Helpers ──────────────────────────────────────────────────────────────────

get_type1 <- function(dat) {
  fit <- lm(y ~ A * B, data = dat)
  a   <- anova(fit)
  data.frame(
    Term  = rownames(a),
    SS    = round(a$`Sum Sq`, 3),
    df    = a$Df,
    MS    = round(a$`Mean Sq`, 3),
    F_val = round(a$`F value`, 3),
    p_val = round(a$`Pr(>F)`, 4),
    stringsAsFactors = FALSE
  )
}

get_type2 <- function(dat) {
  fitFull_noA <- lm(y ~ B,     data = dat)
  fitFull_noB <- lm(y ~ A,     data = dat)
  fitAB       <- lm(y ~ A + B, data = dat)
  fitFull     <- lm(y ~ A * B, data = dat)
  
  ss_A   <- anova(fitFull_noA, fitAB)$`Sum of Sq`[2]
  ss_B   <- anova(fitFull_noB, fitAB)$`Sum of Sq`[2]
  ss_AB  <- anova(fitAB, fitFull)$`Sum of Sq`[2]
  ss_res <- sum(residuals(fitFull)^2)
  
  df_A   <- length(levels(dat$A)) - 1
  df_B   <- length(levels(dat$B)) - 1
  df_AB  <- df_A * df_B
  df_res <- nrow(dat) - length(levels(dat$A)) * length(levels(dat$B))
  
  ms_A <- ss_A/df_A; ms_B <- ss_B/df_B
  ms_AB <- ss_AB/df_AB; ms_res <- ss_res/df_res
  
  f_A  <- ms_A/ms_res;  f_B  <- ms_B/ms_res;  f_AB <- ms_AB/ms_res
  p_A  <- pf(f_A,  df_A,  df_res, lower.tail = FALSE)
  p_B  <- pf(f_B,  df_B,  df_res, lower.tail = FALSE)
  p_AB <- pf(f_AB, df_AB, df_res, lower.tail = FALSE)
  
  data.frame(
    Term  = c("A", "B", "A:B", "Residuals"),
    SS    = round(c(ss_A, ss_B, ss_AB, ss_res), 3),
    df    = c(df_A, df_B, df_AB, df_res),
    MS    = round(c(ms_A, ms_B, ms_AB, ms_res), 3),
    F_val = round(c(f_A, f_B, f_AB, NA), 3),
    p_val = round(c(p_A, p_B, p_AB, NA), 4),
    stringsAsFactors = FALSE
  )
}

get_type3 <- function(dat) {
  dat$A <- C(dat$A, contr.sum)
  dat$B <- C(dat$B, contr.sum)
  
  fitFull <- lm(y ~ A * B, data = dat)
  d1      <- drop1(fitFull, ~ A + B + A:B, test = "F")
  
  ss_A   <- d1["A",   "Sum of Sq"]
  ss_B   <- d1["B",   "Sum of Sq"]
  ss_AB  <- d1["A:B", "Sum of Sq"]
  ss_res <- sum(residuals(fitFull)^2)
  
  df_A   <- d1["A",   "Df"]
  df_B   <- d1["B",   "Df"]
  df_AB  <- d1["A:B", "Df"]
  df_res <- fitFull$df.residual
  
  ms_A   <- ss_A/df_A;  ms_B  <- ss_B/df_B
  ms_AB  <- ss_AB/df_AB; ms_res <- ss_res/df_res
  
  f_A  <- ms_A/ms_res;  f_B  <- ms_B/ms_res;  f_AB <- ms_AB/ms_res
  p_A  <- pf(f_A,  df_A,  df_res, lower.tail = FALSE)
  p_B  <- pf(f_B,  df_B,  df_res, lower.tail = FALSE)
  p_AB <- pf(f_AB, df_AB, df_res, lower.tail = FALSE)
  
  data.frame(
    Term  = c("A", "B", "A:B", "Residuals"),
    SS    = round(c(ss_A, ss_B, ss_AB, ss_res), 3),
    df    = c(df_A, df_B, df_AB, df_res),
    MS    = round(c(ms_A, ms_B, ms_AB, ms_res), 3),
    F_val = round(c(f_A, f_B, f_AB, NA), 3),
    p_val = round(c(p_A, p_B, p_AB, NA), 4),
    stringsAsFactors = FALSE
  )
}

make_data <- function(n11, n12, n21, n22,
                      effA, effB, interaction, noise, seed) {
  set.seed(seed)
  
  mu11 <- 50 + effA/2 + effB/2 + interaction/2
  mu12 <- 50 + effA/2 - effB/2 - interaction/2
  mu21 <- 50 - effA/2 + effB/2 - interaction/2
  mu22 <- 50 - effA/2 - effB/2 + interaction/2
  
  y <- c(rnorm(n11, mu11, noise),
         rnorm(n12, mu12, noise),
         rnorm(n21, mu21, noise),
         rnorm(n22, mu22, noise))
  A <- factor(c(rep("A1", n11 + n12), rep("A2", n21 + n22)))
  B <- factor(c(rep("B1", n11), rep("B2", n12),
                rep("B1", n21), rep("B2", n22)))
  data.frame(y = y, A = A, B = B)
}

sig_stars <- function(p) {
  ifelse(is.na(p), "",
         ifelse(p < .001, "***",
                ifelse(p < .01, "**",
                       ifelse(p < .05, "*",
                              ifelse(p < .1, ".", "")))))
}

cohen_label <- function(d) {
  d <- abs(d)
  if (d < 0.2)      "negligible"
  else if (d < 0.5) "small"
  else if (d < 0.8) "medium"
  else               "large"
}

cohen_colour <- function(d) {
  d <- abs(d)
  if (d < 0.2)      "#888888"
  else if (d < 0.5) "#D95F02"
  else if (d < 0.8) "#1B9E77"
  else               "#569BBD"
}

# What each term tests, by SS type
ss_tests <- list(
  "I" = list(
    "A"         = "y ~ A &rarr; y ~ 1",
    "B"         = "y ~ A + B &rarr; y ~ A",
    "A:B"       = "y ~ A * B &rarr; y ~ A + B",
    "Residuals" = "Full model error"
  ),
  "II" = list(
    "A"         = "y ~ A + B &rarr; y ~ B",
    "B"         = "y ~ A + B &rarr; y ~ A",
    "A:B"       = "y ~ A * B &rarr; y ~ A + B",
    "Residuals" = "Full model error"
  ),
  "III" = list(
    "A"         = "y ~ A * B &rarr; y ~ B + A:B",
    "B"         = "y ~ A * B &rarr; y ~ A + A:B",
    "A:B"       = "y ~ A * B &rarr; y ~ A + B",
    "Residuals" = "Full model error"
  )
)

# ── Interaction plot helper (reused for both orientations) ────────────────────
draw_interaction <- function(d, x_factor, group_factor,
                             show_points, show_ci, interaction_val) {
  
  # x_factor and group_factor are strings: "A" or "B"
  x_var   <- d[[x_factor]]
  grp_var <- d[[group_factor]]
  
  levX   <- levels(x_var)
  levGrp <- levels(grp_var)
  
  agg  <- aggregate(d$y, by = list(X = x_var, G = grp_var), FUN = mean)
  aggs <- aggregate(d$y, by = list(X = x_var, G = grp_var), FUN = sd)
  agg$sd <- aggs$x
  names(agg)[3] <- "y"
  
  ymin <- min(d$y) - 2
  ymax <- max(d$y) + 2
  
  group_cols <- setNames(c("#1B9E77", "#D95F02"), levGrp)
  pchs       <- setNames(c(16, 17),               levGrp)
  
  par(mar = c(4, 5, 3, 2), bg = "white")
  plot(NA, xlim = c(0.7, length(levX) + 0.3), ylim = c(ymin, ymax),
       xaxt = "n",
       xlab = paste("Factor", x_factor),
       ylab = "Mean of Y",
       main = paste0("Interaction Plot: Factor ", x_factor, " on x-axis"),
       cex.lab = 1.1, cex.main = 1.1, col.main = "#569BBD")
  axis(1, at = seq_along(levX), labels = levX, cex.axis = 1.1)
  grid(col = "#e0e0e0", lty = 1, lwd = 0.8)
  box(col = "#cccccc")
  
  for (j in seq_along(levGrp)) {
    gval <- levGrp[j]
    col  <- group_cols[gval]
    pch  <- pchs[gval]
    sub  <- agg[agg$G == gval, ]
    sub  <- sub[order(sub$X), ]
    xpos <- seq_along(levX)
    
    lines(xpos, sub$y,  col = col, lwd = 2.5, lty = j)
    points(xpos, sub$y, col = col, pch = pch, cex = 1.7)
    
    if (show_ci) {
      for (k in seq_along(levX)) {
        w <- 0.04
        segments(xpos[k], sub$y[k]-sub$sd[k],
                 xpos[k], sub$y[k]+sub$sd[k], col = col, lwd = 1.5)
        segments(xpos[k]-w, sub$y[k]-sub$sd[k],
                 xpos[k]+w, sub$y[k]-sub$sd[k], col = col, lwd = 1.5)
        segments(xpos[k]-w, sub$y[k]+sub$sd[k],
                 xpos[k]+w, sub$y[k]+sub$sd[k], col = col, lwd = 1.5)
      }
    }
  }
  
  if (show_points) {
    for (j in seq_along(levGrp)) {
      gval <- levGrp[j]
      col  <- group_cols[gval]
      for (k in seq_along(levX)) {
        pts  <- d$y[x_var == levX[k] & grp_var == gval]
        xjit <- jitter(rep(k, length(pts)), amount = 0.06)
        points(xjit, pts,
               col = adjustcolor(col, alpha.f = 0.35),
               pch = 16, cex = 0.9)
      }
    }
  }
  
  legend("topright",
         legend = paste0(group_factor, " = ", levGrp),
         col = group_cols, lty = 1:2, pch = pchs,
         lwd = 2, pt.cex = 1.3, cex = 0.9,
         bg = "white", box.col = "#cccccc")
  
  mtext(if (interaction_val == 0)
    "Interaction = 0 \u2192 lines should be parallel"
    else
      "Interaction \u2260 0 \u2192 lines are non-parallel by design",
    side = 3, line = 0.2, cex = 0.85, col = "#555555", font = 3)
}

# ── Venn diagram ─────────────────────────────────────────────────────────────
draw_venn <- function(ss_A, ss_B, ss_AB, ss_res, ss_type) {
  col_A   <- "#569BBD"
  col_B   <- "#1B9E77"
  col_AB  <- "#D95F02"
  col_res <- "#AAAAAA"
  
  total <- ss_A + ss_B + ss_AB + ss_res
  if (total <= 0) { plot.new(); text(0.5, 0.5, "No data"); return(invisible()) }
  
  r_scale <- 0.22
  rA  <- max(r_scale * sqrt(ss_A  / total), 0.015)
  rB  <- max(r_scale * sqrt(ss_B  / total), 0.015)
  rAB <- max(r_scale * sqrt(ss_AB / total), 0.015)
  
  cAx <- 0.35; cAy <- 0.62
  cBx <- 0.65; cBy <- 0.62
  cIx <- 0.50; cIy <- 0.38
  
  par(mar = c(1, 1, 3, 1), bg = "white")
  plot(NA, xlim = c(0, 1), ylim = c(0, 1),
       asp = 1, axes = FALSE, xlab = "", ylab = "",
       main = paste("Variance Partitioning \u2014 Type", ss_type, "SS"),
       cex.main = 1.15, col.main = "#569BBD")
  
  draw_circle <- function(cx, cy, r, col_fill, col_border = "white", alpha = 0.45) {
    theta <- seq(0, 2*pi, length.out = 200)
    polygon(cx + r * cos(theta), cy + r * sin(theta),
            col = adjustcolor(col_fill, alpha.f = alpha),
            border = col_border, lwd = 2)
  }
  
  r_res <- max(r_scale * sqrt(ss_res / total), 0.02)
  draw_circle(0.50, 0.55, min(r_res * 1.8, 0.28),
              col_fill = col_res, alpha = 0.15)
  
  draw_circle(cAx, cAy, rA,  col_fill = col_A)
  draw_circle(cBx, cBy, rB,  col_fill = col_B)
  draw_circle(cIx, cIy, rAB, col_fill = col_AB)
  
  text(cAx, cAy, paste0("A\n", round(ss_A, 1)),         cex = 0.85, font = 2, col = "white")
  text(cBx, cBy, paste0("B\n", round(ss_B, 1)),         cex = 0.85, font = 2, col = "white")
  text(cIx, cIy, paste0("A\u00d7B\n", round(ss_AB, 1)), cex = 0.85, font = 2, col = "white")
  text(0.50, 0.10, paste0("Residual: ", round(ss_res, 1)), cex = 0.85, col = col_res, font = 2)
  
  legend("bottomright",
         legend = c("Factor A", "Factor B", "A\u00d7B interaction", "Residual"),
         fill   = adjustcolor(c(col_A, col_B, col_AB, col_res), alpha.f = 0.6),
         border = "white", cex = 0.8, bty = "n")
  
  mtext("Circle size \u221d SS (schematic)", side = 1, line = -0.5,
        cex = 0.75, col = "#888888", font = 3)
}

# ── SS bar chart ──────────────────────────────────────────────────────────────
draw_ss_bars <- function(a1, a2, a3) {
  terms <- c("A", "B", "A:B")
  col1  <- "#569BBD"; col2 <- "#1B9E77"; col3 <- "#D95F02"
  
  mat <- matrix(NA, nrow = 3, ncol = 3,
                dimnames = list(c("Type I", "Type II", "Type III"), terms))
  for (trm in terms) {
    mat["Type I",   trm] <- a1$SS[a1$Term == trm]
    mat["Type II",  trm] <- a2$SS[a2$Term == trm]
    mat["Type III", trm] <- a3$SS[a3$Term == trm]
  }
  
  ymax <- max(mat, na.rm = TRUE) * 1.15
  if (ymax <= 0) ymax <- 1
  
  par(mar = c(5, 5, 3, 8), bg = "white")
  bp <- barplot(mat, beside = TRUE,
                col = c(col1, col2, col3), border = "white",
                ylim = c(0, ymax), ylab = "Sums of Squares",
                main = "SS by Term and Type",
                cex.main = 1.15, col.main = "#569BBD",
                cex.names = 1.1, cex.axis = 1.0, las = 1)
  
  for (j in seq_along(terms))
    for (i in 1:3) {
      val <- mat[i, j]
      if (!is.na(val) && val > 0)
        text(bp[i, j], val + ymax * 0.02,
             labels = round(val, 1), cex = 0.75, col = "#333333")
    }
  
  legend("topleft", inset = c(1.01, 0),
         legend = c("Type I", "Type II", "Type III"),
         fill = c(col1, col2, col3), border = "white",
         cex = 0.9, bty = "n", xpd = TRUE)
  
  abline(h = 0, col = "#cccccc")
}

# ── Marginal means plot ───────────────────────────────────────────────────────
draw_marginal <- function(means, sds, ns, all_y, all_groups,
                          factor_name, col, show_points, show_ci) {
  levs <- names(means)
  ymin <- min(all_y) - 2
  ymax <- max(all_y) + 2
  
  par(mar = c(4, 4, 3, 1), bg = "white")
  plot(NA, xlim = c(0.5, length(levs) + 0.5), ylim = c(ymin, ymax),
       xaxt = "n", xlab = factor_name, ylab = "Marginal Mean of Y",
       main = paste("Main Effect of", factor_name),
       sub  = "Error bars = \u00b11 SD of all observations at this level",
       cex.main = 1.1, col.main = "#569BBD",
       cex.sub  = 0.8, col.sub  = "#888888",
       cex.lab  = 1.0, cex.axis = 1.0)
  axis(1, at = seq_along(levs), labels = levs, cex.axis = 1.1)
  grid(col = "#e0e0e0", lty = 1, lwd = 0.8)
  box(col = "#cccccc")
  
  if (show_points) {
    for (i in seq_along(levs)) {
      pts  <- all_y[all_groups == levs[i]]
      xjit <- jitter(rep(i, length(pts)), amount = 0.06)
      points(xjit, pts,
             col = adjustcolor(col, alpha.f = 0.35),
             pch = 16, cex = 0.9)
    }
  }
  
  lines(seq_along(levs), means, col = col, lwd = 2, lty = 2)
  
  for (i in seq_along(levs)) {
    if (show_ci) {
      w <- 0.06
      segments(i, means[i] - sds[i], i, means[i] + sds[i], col = col, lwd = 1.5)
      segments(i-w, means[i]-sds[i], i+w, means[i]-sds[i], col = col, lwd = 1.5)
      segments(i-w, means[i]+sds[i], i+w, means[i]+sds[i], col = col, lwd = 1.5)
    }
    points(i, means[i], pch = 16, col = col, cex = 1.8)
    text(i, ymin + (ymax - ymin) * 0.04,
         paste0("n=", ns[i]), cex = 0.8, col = "#666666")
  }
}

# ── UI ───────────────────────────────────────────────────────────────────────
ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      body { background: #ffffff; font-size: 15px; color: #222; }

      .nav-tabs > li > a { color: #569BBD; font-size: 15px; }
      .nav-tabs > li.active > a,
      .nav-tabs > li.active > a:focus,
      .nav-tabs > li.active > a:hover {
        background-color: #569BBD !important;
        color: white !important;
        border-color: #569BBD;
      }

      .info-box {
        background: #f4f8fb;
        border-left: 4px solid #569BBD;
        border-radius: 3px;
        padding: 9px 13px;
        margin-bottom: 11px;
        font-size: 0.92em;
        line-height: 1.5;
      }
      .warn-box {
        background: #fff8e1;
        border-left: 4px solid #f0a500;
        border-radius: 3px;
        padding: 9px 13px;
        margin-bottom: 11px;
        font-size: 0.92em;
        line-height: 1.5;
      }
      .cohen-box {
        background: #f4f8fb;
        border-radius: 3px;
        padding: 8px 12px;
        margin-top: 4px;
        margin-bottom: 10px;
        font-size: 0.90em;
        line-height: 1.8;
      }

      .ss-table {
        width: 100%;
        border-collapse: collapse;
        font-size: 0.92em;
        margin-top: 4px;
      }
      .ss-table th {
        background-color: #569BBD;
        color: white;
        padding: 8px 10px;
        text-align: left;
        font-weight: 600;
      }
      .ss-table td { padding: 7px 10px; border: 1px solid #ddd; }
      .ss-table td.test-col { font-size: 0.85em; color: #555; font-style: italic; }
      .ss-table tr:nth-child(even) td { background-color: #f4f8fb; }
      .ss-table tr:last-child td     { font-weight: bold; }
      .ss-table tr.diff-row td       { background-color: #fff3cd !important; }

      .badge-sig {
        background: #2e7d32; color: white;
        border-radius: 3px; padding: 1px 5px; font-size: 0.82em;
      }
      .section-label {
        font-size: 0.78em; text-transform: uppercase;
        letter-spacing: .07em; color: #888;
        margin-bottom: 3px; margin-top: 8px;
      }
      .type-desc { font-size: 0.88em; color: #444; margin-top: 5px; line-height: 1.5; }

      h3 { color: #569BBD; font-size: 1.15em; margin-top: 0; }
    "))
  ),
  
  titlePanel("Sums of Squares: Type I, II, and III"),
  
  tabsetPanel(id = "tabs",
              
              # ── TAB 1: Build ──────────────────────────────────────────────────────
              tabPanel("\u2460 Build Data",
                       br(),
                       fluidRow(
                         column(6,
                                wellPanel(
                                  h3("Cell Sizes"),
                                  div(class = "info-box",
                                      "Set the number of observations in each of the four cells
               of the 2\u00d72 design independently. Equal cell sizes mean
               all three SS types give the same answer \u2014 unequal sizes
               are what make the choice matter."
                                  ),
                                  sliderInput("n11", "Cell A1 and B1:", min = 2, max = 30, value = 10, step = 1),
                                  sliderInput("n12", "Cell A1 and B2:", min = 2, max = 30, value = 10, step = 1),
                                  sliderInput("n21", "Cell A2 and B1:", min = 2, max = 30, value = 10, step = 1),
                                  sliderInput("n22", "Cell A2 and B2:", min = 2, max = 30, value = 10, step = 1),
                                  hr(),
                                  div(class = "section-label", "Seed (re-roll random data)"),
                                  sliderInput("seed", NULL, min = 1, max = 100, value = 42, step = 1)
                                ),
                                wellPanel(
                                  h3("Cell Counts"),
                                  tableOutput("cellCounts"),
                                  div(class = "info-box",
                                      "The table updates as you move the cell size sliders above.
               Try making the cells unequal to see how imbalance affects the SS types."
                                  )
                                )
                         ),
                         column(6,
                                wellPanel(
                                  h3("Effect Sizes"),
                                  div(class = "info-box",
                                      "These control how large each effect is in the same units as Y.
               Zero means no effect. Negative values reverse the direction.
               Noise is the within-group SD. Cohen's d = Effect \u00f7 Noise."
                                  ),
                                  sliderInput("effA",        "Effect of Factor A:",       min = -40, max = 40, value = 15, step = 1),
                                  sliderInput("effB",        "Effect of Factor B:",       min = -40, max = 40, value = 5,  step = 1),
                                  sliderInput("interaction", "Interaction (A \u00d7 B):", min = -40, max = 40, value = 0,  step = 1),
                                  sliderInput("noise",       "Noise (within-group SD):",  min = 1,   max = 20, value = 10, step = 1),
                                  uiOutput("cohenDisplay"),
                                  div(class = "warn-box",
                                      "\u26a0\ufe0f  When Interaction \u2260 0, Type II SS can be misleading
               because it ignores the interaction when testing main effects."
                                  )
                                )
                         )
                       )
              ),
              
              # ── TAB 2: Data Plot ──────────────────────────────────────────────────
              tabPanel("\u2461 Data Plot",
                       br(),
                       fluidRow(
                         column(3,
                                wellPanel(
                                  h3("Plot Options"),
                                  checkboxInput("showPoints", "Show individual data points", value = TRUE),
                                  checkboxInput("showCI",     "Show \u00b11 SD error bars",  value = TRUE),
                                  hr(),
                                  downloadButton("downloadData", "Download current data (.csv)"),
                                  br(), br(),
                                  div(class = "info-box",
                                      "The two interaction plots show identical data — only the
               x-axis factor differs. Notice how the same interaction
               pattern looks different depending on which factor is plotted
               on the x-axis. The choice is arbitrary.",
                                      br(), br(),
                                      "The options above apply to all plots on this tab."
                                  )
                                ),
                                wellPanel(
                                  h3("Cell Means"),
                                  tableOutput("cellMeansTable")
                                )
                         ),
                         column(9,
                                fluidRow(
                                  column(6,
                                         wellPanel(
                                           plotOutput("interactionPlotA", height = "300px")
                                         )
                                  ),
                                  column(6,
                                         wellPanel(
                                           plotOutput("interactionPlotB", height = "300px")
                                         )
                                  )
                                ),
                                uiOutput("mainEffectWarning"),
                                fluidRow(
                                  column(6,
                                         wellPanel(plotOutput("mainEffectA", height = "260px"))
                                  ),
                                  column(6,
                                         wellPanel(plotOutput("mainEffectB", height = "260px"))
                                  )
                                )
                         )
                       )
              ),
              
              # ── TAB 3: ANOVA Table ────────────────────────────────────────────────
              tabPanel("\u2462 ANOVA Table",
                       br(),
                       fluidRow(
                         column(4,
                                wellPanel(
                                  h3("Choose SS Type"),
                                  radioButtons("ssType", NULL,
                                               choices = c("Type I \u2014 Sequential"    = "I",
                                                           "Type II \u2014 Hierarchical" = "II",
                                                           "Type III \u2014 Marginal"    = "III"),
                                               selected = "I"
                                  ),
                                  uiOutput("ssDescription"),
                                  hr(),
                                  div(class = "info-box",
                                      strong("Key insight:"), " With equal cell sizes, all three types
               give identical SS. Imbalance is what makes the choice matter."
                                  )
                                )
                         ),
                         column(8,
                                wellPanel(
                                  h3(textOutput("tableTitle")),
                                  uiOutput("anovaTableUI")
                                ),
                                wellPanel(
                                  h3("SS Comparison Across All Three Types"),
                                  uiOutput("comparisonUI")
                                )
                         )
                       )
              ),
              
              # ── TAB 4: Visualise SS ───────────────────────────────────────────────
              tabPanel("\u2463 Visualise SS",
                       br(),
                       fluidRow(
                         column(12,
                                div(class = "info-box",
                                    strong("How to read these plots:"),
                                    " The bar chart (right) shows the SS for each term side by side
             for all three SS types \u2014 bars that differ in height are where
             the types disagree. The Venn diagram (left) gives a schematic
             sense of how the total variance is carved up; circle size is
             roughly proportional to SS. When cell sizes are equal all bars
             will be the same height and all circles the same size across types.",
                                    br(), br(),
                                    "Use the radio buttons to change which SS type the Venn diagram
             reflects \u2014 the bar chart always shows all three at once."
                                )
                         )
                       ),
                       fluidRow(
                         column(4,
                                wellPanel(
                                  h3("Venn Diagram"),
                                  radioButtons("vennType", "Show SS type:",
                                               choices  = c("Type I" = "I", "Type II" = "II", "Type III" = "III"),
                                               selected = "I", inline = TRUE
                                  ),
                                  plotOutput("vennPlot", height = "380px")
                                )
                         ),
                         column(8,
                                wellPanel(
                                  h3("SS by Term \u2014 All Three Types"),
                                  plotOutput("ssBarPlot", height = "380px")
                                )
                         )
                       )
              )
  )
)

# ── Server ───────────────────────────────────────────────────────────────────
server <- function(input, output, session) {
  
  dat <- reactive({
    make_data(input$n11, input$n12, input$n21, input$n22,
              input$effA, input$effB, input$interaction,
              input$noise, input$seed)
  })
  
  output$downloadData <- downloadHandler(
    filename = function() paste0("ss_data_seed", input$seed, ".csv"),
    content  = function(file) write.csv(dat(), file, row.names = FALSE)
  )
  
  output$cohenDisplay <- renderUI({
    noise <- input$noise
    if (is.null(noise) || noise == 0) return(NULL)
    dA  <- input$effA        / noise
    dB  <- input$effB        / noise
    dAB <- input$interaction / noise
    
    make_row <- function(label, d) {
      tags$div(
        tags$span(style = "color:#555; width:110px; display:inline-block;", label),
        tags$span(style = paste0("color:", cohen_colour(d), "; font-weight:bold;"),
                  paste0("d = ", round(d, 2), " (", cohen_label(d), ")"))
      )
    }
    div(class = "cohen-box",
        tags$strong("Approximate Cohen's d:"), br(),
        make_row("Factor A:",    dA),
        make_row("Factor B:",    dB),
        make_row("Interaction:", dAB)
    )
  })
  
  ss1 <- reactive({ tryCatch(get_type1(dat()), error = function(e) NULL) })
  ss2 <- reactive({ tryCatch(get_type2(dat()), error = function(e) NULL) })
  ss3 <- reactive({ tryCatch(get_type3(dat()), error = function(e) NULL) })
  
  current_ss <- reactive({
    switch(input$ssType, "I" = ss1(), "II" = ss2(), "III" = ss3())
  })
  
  output$cellCounts <- renderTable({
    d   <- dat()
    tab <- table(d$A, d$B)
    df  <- as.data.frame.matrix(tab)
    cbind(data.frame(` ` = rownames(tab), check.names = FALSE), df)
  }, striped = TRUE, bordered = TRUE, hover = TRUE)
  
  output$tableTitle <- renderText({
    paste("ANOVA Table \u2014 Type", input$ssType, "Sums of Squares")
  })
  
  output$ssDescription <- renderUI({
    switch(input$ssType,
           "I" = tagList(
             div(class = "type-desc",
                 strong("Type I is sequential."),
                 " SS are calculated in the order terms enter the model: A first,
            then B adjusted for A, then A\u00d7B adjusted for both. Swapping
            the order of A and B changes the results when groups are unequal."
             ),
             br(),
             div(class = "type-desc",
                 strong("When to use this type:"),
                 " Type I is rarely the first choice for unbalanced designs because
            the results depend on the order terms are entered \u2014 swapping
            A and B changes the SS for both. However, it is what base R's
            anova() function produces by default, and it is perfectly
            appropriate for balanced designs where all cell sizes are equal.
            It is also useful when you have a genuine scientific reason to
            believe one factor should be accounted for before another."
             ),
             br(),
             div(class = "type-desc",
                 em("To see order-dependence in action: download the data from the
              Data Plot tab and run anova(lm(y ~ A * B, data=d)) then
              anova(lm(y ~ B * A, data=d)) in R with unequal cell sizes.
              The SS for A and B will swap; Type II and III will not.")
             )
           ),
           "II" = tagList(
             div(class = "type-desc",
                 strong("Type II is hierarchical."),
                 " Each main effect is adjusted for the other, but not for the
            interaction. Only appropriate when the interaction is truly zero."
             ),
             br(),
             div(class = "type-desc",
                 strong("When to use this type:"),
                 " Type II is the default in several statistical packages and is
            common in fields where balanced designs are typical and interactions
            are not always expected. It has slightly more statistical power than
            Type III for detecting main effects when the interaction is truly
            zero, because it does not adjust for the interaction term. The key
            assumption is that the interaction is negligible \u2014 if a real
            interaction exists, Type II results for main effects can be misleading."
             )
           ),
           "III" = tagList(
             div(class = "type-desc",
                 strong("Type III is marginal."),
                 " Every term is tested controlling for all others, including the
            interaction. This is the default in SPSS and most software."
             ),
             br(),
             div(class = "type-desc",
                 strong("When to use this type:"),
                 " Type III is the default in SPSS and SAS, making it the most
            commonly encountered type in software output. It tests each term
            as if it were the last one entered, controlling for everything
            else including the interaction. This makes it symmetric and
            consistent regardless of design balance. The philosophical
            criticism is that testing a main effect while controlling for an
            interaction it is part of is conceptually awkward \u2014 if
            A\u00d7B is real, the effect of A averaged across B may not be
            meaningful on its own."
             )
           )
    )
  })
  
  output$anovaTableUI <- renderUI({
    ss    <- current_ss()
    a1    <- ss1(); a2 <- ss2(); a3 <- ss3()
    type  <- input$ssType
    tests <- ss_tests[[type]]
    
    if (is.null(ss))
      return(div("Cannot compute \u2014 ensure each cell has at least 2 observations."))
    
    diff_rows <- character(0)
    if (!is.null(a1) && !is.null(a2) && !is.null(a3)) {
      for (trm in c("A", "B", "A:B")) {
        v1 <- a1$SS[a1$Term == trm]
        v2 <- a2$SS[a2$Term == trm]
        v3 <- a3$SS[a3$Term == trm]
        if (length(v1) && length(v2) && length(v3) &&
            !(abs(v1-v2) < 0.01 && abs(v1-v3) < 0.01))
          diff_rows <- c(diff_rows, trm)
      }
    }
    
    rows_html <- ""
    for (i in seq_len(nrow(ss))) {
      row  <- ss[i, ]
      cls  <- if (row$Term %in% diff_rows) " class='diff-row'" else ""
      fval <- if (is.na(row$F_val)) "&mdash;" else row$F_val
      pval <- if (is.na(row$p_val)) "&mdash;" else {
        star <- sig_stars(row$p_val)
        pstr <- formatC(row$p_val, digits = 4, format = "f")
        if (row$p_val < .05)
          paste0("<span class='badge-sig'>", pstr, " ", star, "</span>")
        else
          paste0(pstr, " ", star)
      }
      test_label <- if (!is.null(tests[[row$Term]])) tests[[row$Term]] else ""
      rows_html <- paste0(rows_html,
                          "<tr", cls, ">",
                          "<td>", row$Term, "</td>",
                          "<td class='test-col'>", test_label, "</td>",
                          "<td>", row$SS,   "</td>",
                          "<td>", row$df,   "</td>",
                          "<td>", row$MS,   "</td>",
                          "<td>", fval,     "</td>",
                          "<td>", pval,     "</td>",
                          "</tr>")
    }
    
    tbl <- paste0(
      "<table class='ss-table'>",
      "<thead><tr>",
      "<th>Term</th>",
      "<th>What is being tested?</th>",
      "<th>SS</th><th>df</th><th>MS</th><th>F</th><th>p</th>",
      "</tr></thead><tbody>", rows_html, "</tbody></table>"
    )
    
    note <- if (length(diff_rows) > 0)
      div(class = "warn-box", style = "margin-top:10px;",
          paste0("\u26a0\ufe0f  Highlighted rows (",
                 paste(diff_rows, collapse = ", "),
                 ") differ across SS types because cell sizes are unequal."))
    else
      div(class = "info-box", style = "margin-top:10px;",
          "\u2713  All three types agree here \u2014 cell sizes are equal.
         Try making them unequal in the Build tab.")
    
    tagList(HTML(tbl), note)
  })
  
  output$comparisonUI <- renderUI({
    a1 <- ss1(); a2 <- ss2(); a3 <- ss3()
    if (is.null(a1) || is.null(a2) || is.null(a3))
      return(div("Comparison not available."))
    
    rows_html <- ""
    for (trm in c("A", "B", "A:B")) {
      v1 <- round(a1$SS[a1$Term == trm], 2)
      v2 <- round(a2$SS[a2$Term == trm], 2)
      v3 <- round(a3$SS[a3$Term == trm], 2)
      if (!length(v1)) next
      cls <- if (!(abs(v1-v2) < 0.01 && abs(v1-v3) < 0.01)) " class='diff-row'" else ""
      rows_html <- paste0(rows_html,
                          "<tr", cls, ">",
                          "<td>", trm, "</td>",
                          "<td>", v1,  "</td>",
                          "<td>", v2,  "</td>",
                          "<td>", v3,  "</td>",
                          "</tr>")
    }
    
    HTML(paste0(
      "<table class='ss-table'>",
      "<thead><tr>",
      "<th>Term</th><th>Type I SS</th><th>Type II SS</th><th>Type III SS</th>",
      "</tr></thead><tbody>", rows_html, "</tbody></table>"
    ))
  })
  
  output$cellMeansTable <- renderTable({
    d   <- dat()
    agg <- aggregate(y ~ A + B, data = d,
                     FUN = function(x) paste0(round(mean(x), 1),
                                              " (n=", length(x), ")"))
    colnames(agg)[3] <- "Mean (n)"
    agg
  }, striped = TRUE, bordered = TRUE)
  
  output$mainEffectWarning <- renderUI({
    if (input$interaction != 0)
      div(class = "warn-box",
          "\u26a0\ufe0f  An interaction is present. The marginal means below
         average across levels of the other factor, which can be misleading
         when the lines above are non-parallel. Interpret main effects
         with caution."
      )
    else
      div(class = "info-box",
          "\u2139\ufe0f  No interaction is set. The marginal means below
         reliably summarise the effect of each factor independently."
      )
  })
  
  # Interaction plot — Factor A on x-axis
  output$interactionPlotA <- renderPlot({
    draw_interaction(dat(), "A", "B",
                     input$showPoints, input$showCI, input$interaction)
  })
  
  # Interaction plot — Factor B on x-axis
  output$interactionPlotB <- renderPlot({
    draw_interaction(dat(), "B", "A",
                     input$showPoints, input$showCI, input$interaction)
  })
  
  output$mainEffectA <- renderPlot({
    d <- dat()
    draw_marginal(
      means       = tapply(d$y, d$A, mean),
      sds         = tapply(d$y, d$A, sd),
      ns          = tapply(d$y, d$A, length),
      all_y       = d$y,
      all_groups  = d$A,
      factor_name = "Factor A",
      col         = "#569BBD",
      show_points = input$showPoints,
      show_ci     = input$showCI
    )
  })
  
  output$mainEffectB <- renderPlot({
    d <- dat()
    draw_marginal(
      means       = tapply(d$y, d$B, mean),
      sds         = tapply(d$y, d$B, sd),
      ns          = tapply(d$y, d$B, length),
      all_y       = d$y,
      all_groups  = d$B,
      factor_name = "Factor B",
      col         = "#1B9E77",
      show_points = input$showPoints,
      show_ci     = input$showCI
    )
  })
  
  output$vennPlot <- renderPlot({
    ss <- switch(input$vennType, "I" = ss1(), "II" = ss2(), "III" = ss3())
    if (is.null(ss)) { plot.new(); return(invisible()) }
    draw_venn(ss$SS[ss$Term == "A"],   ss$SS[ss$Term == "B"],
              ss$SS[ss$Term == "A:B"], ss$SS[ss$Term == "Residuals"],
              input$vennType)
  })
  
  output$ssBarPlot <- renderPlot({
    a1 <- ss1(); a2 <- ss2(); a3 <- ss3()
    if (is.null(a1) || is.null(a2) || is.null(a3)) { plot.new(); return(invisible()) }
    draw_ss_bars(a1, a2, a3)
  })
}

shinyApp(ui = ui, server = server)