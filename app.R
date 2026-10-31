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

# ── SS as model comparisons (Tab 4) ──────────────────────────────────────────
# Every SS term is the squared distance between the fitted values of two
# models: a smaller model without the term and a larger model with it.
# Within a cell every observation shares the same fitted values, so
# SS_term = sum over cells of n_cell * (larger fit - smaller fit)^2.

ss_models <- list(
  "I"   = list("A"   = c("y ~ 1",       "y ~ A"),
               "B"   = c("y ~ A",       "y ~ A + B"),
               "A:B" = c("y ~ A + B",   "y ~ A * B")),
  "II"  = list("A"   = c("y ~ B",       "y ~ A + B"),
               "B"   = c("y ~ A",       "y ~ A + B"),
               "A:B" = c("y ~ A + B",   "y ~ A * B")),
  "III" = list("A"   = c("y ~ B + A:B", "y ~ A * B"),
               "B"   = c("y ~ A + A:B", "y ~ A * B"),
               "A:B" = c("y ~ A + B",   "y ~ A * B"))
)

term_cols <- c("A" = "#569BBD", "B" = "#1B9E77", "A:B" = "#D95F02",
               "Error" = "#009499")

# Fitted values of the smaller and larger model behind one SS term
fit_pair <- function(dat, type, term) {
  f <- function(fm) fitted(lm(fm, data = dat))
  if (term == "A:B") return(list(small = f(y ~ A + B), big = f(y ~ A * B)))
  if (type == "I") {
    if (term == "A") return(list(small = f(y ~ 1), big = f(y ~ A)))
    return(list(small = f(y ~ A), big = f(y ~ A + B)))
  }
  if (type == "II") {
    if (term == "A") return(list(small = f(y ~ B), big = f(y ~ A + B)))
    return(list(small = f(y ~ A), big = f(y ~ A + B)))
  }
  # Type III: sum-to-zero coding, then drop the main-effect column only
  d3 <- dat
  d3$A <- C(d3$A, contr.sum)
  d3$B <- C(d3$B, contr.sum)
  full <- lm(y ~ A * B, data = d3)
  X    <- model.matrix(full)
  keep <- setdiff(colnames(X), if (term == "A") "A1" else "B1")
  list(small = fitted(lm(dat$y ~ X[, keep] - 1)), big = fitted(full))
}

ss_cell_x <- c("A1.B1" = 0.70, "A1.B2" = 1.40, "A2.B1" = 2.25, "A2.B2" = 2.95)
ss_pch    <- c("B1" = 19, "B2" = 17)   # circle = B1, triangle = B2
ss_cell_pch <- function(cl) ss_pch[sub(".*\\.", "", cl)]

ss_axes <- function(ylim) {
  plot(NA, xlim = c(0.35, 3.35), ylim = ylim, xaxt = "n", bty = "l",
       xlab = "", ylab = "y", las = 1, cex.axis = 1.0, cex.lab = 1.1)
  # faint band behind each cell, and a divider between A1 and A2
  usr <- par("usr")
  rect(ss_cell_x - 0.28, usr[3], ss_cell_x + 0.31, usr[4], col = "#f4f4f4", border = NA)
  abline(v = mean(ss_cell_x[2:3]), col = "gray80", lty = 3)
  box(bty = "l")
  axis(1, at = ss_cell_x, labels = c("B1", "B2", "B1", "B2"),
       tick = FALSE, line = -0.6, cex.axis = 1.0)
  axis(1, at = c(mean(ss_cell_x[1:2]), mean(ss_cell_x[3:4])), labels = c("A1", "A2"),
       tick = FALSE, line = 0.6, font = 2, cex.axis = 1.1)
}

# Fixed horizontal jitter for each observation (does not touch the data RNG)
ss_jitter <- function(n) ((seq_len(n) - 0.5) / n - 0.5) * 0.3

draw_ss_term <- function(dat, type, term, ylim) {
  fp   <- fit_pair(dat, type, term)
  cell <- interaction(dat$A, dat$B, sep = ".")
  ss   <- sum((fp$big - fp$small)^2)
  col  <- term_cols[term]
  
  par(mar = c(3.2, 4.2, 5.2, 0.8))
  ss_axes(ylim)
  for (cl in names(ss_cell_x)) {
    i  <- which(cell == cl)
    x0 <- ss_cell_x[cl]
    points(x0 + ss_jitter(length(i)), dat$y[i],
           pch = ss_cell_pch(cl), col = adjustcolor("gray30", 0.55), cex = 0.95)
    ys <- fp$small[i][1]
    yb <- fp$big[i][1]
    segments(x0 - 0.2, ys, x0 + 0.2, ys, col = "gray35", lwd = 2.5, lty = 2)
    segments(x0 - 0.2, yb, x0 + 0.2, yb, col = "black", lwd = 3)
    segments(x0 + 0.23, ys, x0 + 0.23, yb, col = col, lwd = 7, lend = 1)
    text(x0 + 0.23, max(ys, yb), paste0("×", length(i)),
         pos = 3, offset = 0.3, col = col, cex = 0.9, font = 2)
  }
  title(main = bquote(SS[.(term)] == .(formatC(ss, format = "f", digits = 1))),
        col.main = col, line = 3.3, cex.main = 1.35)
  mtext(paste0(ss_models[[type]][[term]][1], "  →  ",
               ss_models[[type]][[term]][2]),
        side = 3, line = 1.8, cex = 1.0, col = "gray25")
  mtext(expression("dashed = " * hat(Y) * " without term;   solid = " * hat(Y) * " with term"),
        side = 3, line = 0.2, cex = 0.85, col = "gray45")
}

draw_ss_error <- function(dat, ylim) {
  full <- fitted(lm(y ~ A * B, data = dat))
  cell <- interaction(dat$A, dat$B, sep = ".")
  ss   <- sum((dat$y - full)^2)
  col  <- term_cols["Error"]
  
  par(mar = c(3.2, 4.2, 5.2, 0.8))
  ss_axes(ylim)
  for (cl in names(ss_cell_x)) {
    i  <- which(cell == cl)
    x0 <- ss_cell_x[cl]
    xs <- x0 + ss_jitter(length(i))
    segments(xs, dat$y[i], xs, full[i], col = adjustcolor(col, 0.8), lwd = 1.3)
    points(xs, dat$y[i], pch = ss_cell_pch(cl), col = adjustcolor("gray30", 0.6), cex = 0.95)
    segments(x0 - 0.2, full[i][1], x0 + 0.2, full[i][1], col = "black", lwd = 3)
  }
  title(main = bquote(SS[error] == .(formatC(ss, format = "f", digits = 1))),
        col.main = col, line = 3.3, cex.main = 1.35)
  mtext("each point → its cell mean", side = 3, line = 1.8,
        cex = 1.0, col = "gray25")
  mtext(expression("black = " * hat(Y) * " of full model;  identical for all Types"),
        side = 3, line = 0.2, cex = 0.85, col = "gray45")
}

draw_ss_panels <- function(dat, type) {
  ylim <- range(dat$y) + c(-2, 0.12 * diff(range(dat$y)))
  par(mfrow = c(1, 4), bg = "white", cex = 1)  # stop mfrow shrinking text
  for (tm in c("A", "B", "A:B")) {
    draw_ss_term(dat, type, tm, ylim)
    if (tm == "A")
      legend("topright", pch = ss_pch, legend = names(ss_pch), col = "gray30",
             bty = "n", cex = 0.9, horiz = TRUE, inset = c(0, -0.02))
  }
  draw_ss_error(dat, ylim)
}

# Stacked bar: SS total vs the selected Type's terms (plus any unassigned gap)
draw_partition <- function(ss_tab, ss_total, type) {
  vals <- c(ss_tab$SS[ss_tab$Term == "A"],   ss_tab$SS[ss_tab$Term == "B"],
            ss_tab$SS[ss_tab$Term == "A:B"], ss_tab$SS[ss_tab$Term == "Residuals"])
  names(vals) <- c("A", "B", "A:B", "Error")
  cols <- term_cols[names(vals)]
  gap  <- ss_total - sum(vals)
  ytop <- max(ss_total, sum(vals)) * 1.3
  
  par(mar = c(3, 5, 3.5, 9), bg = "white")
  plot(NA, xlim = c(0.3, 2.7), ylim = c(0, ytop), xaxt = "n", bty = "l",
       xlab = "", ylab = "Sum of squares", las = 1)
  axis(1, at = 1:2, labels = c("SS total", paste("Type", type, "terms")),
       tick = FALSE, cex.axis = 1.05, font = 2)
  title(main = paste("Type", type, "sums of squares"),
        cex.main = 1.15, col.main = "#569BBD")
  
  rect(0.6, 0, 1.4, ss_total, col = "gray80", border = "white")
  text(1, ss_total, formatC(ss_total, format = "f", digits = 1), pos = 3, font = 2)
  
  bottoms <- c(0, cumsum(vals)[-length(vals)])
  for (k in seq_along(vals)) {
    rect(1.6, bottoms[k], 2.4, bottoms[k] + vals[k], col = cols[k], border = "white")
    if (vals[k] > ytop * 0.04)
      text(2, bottoms[k] + vals[k] / 2,
           formatC(vals[k], format = "f", digits = 1), col = "white", font = 2, cex = 0.9)
  }
  if (gap > ss_total * 0.001) {
    rect(1.6, sum(vals), 2.4, ss_total, col = NA, border = "gray40", lty = 2, lwd = 1.5)
    gap_lab <- if (type == "II") "shared by A and B:\ncredited to neither\n" else
      "not assigned\nto any term\n"
    text(2, ss_total, paste0(gap_lab, formatC(gap, format = "f", digits = 1)),
         pos = 3, cex = 0.8, col = "gray25")
  } else if (gap < -ss_total * 0.001) {
    segments(1.5, ss_total, 2.5, ss_total, lty = 2, col = "gray40", lwd = 1.5)
    text(2, ss_total, "SS total", pos = 1, cex = 0.8, col = "gray25")
  }
  legend("topleft", inset = c(1.01, 0), xpd = TRUE, bty = "n", cex = 0.9,
         fill = c(cols, "gray80"), border = "white",
         legend = c("Factor A", "Factor B", "A×B", "Error", "SS total"))
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

# ── Safari-safe plotting ─────────────────────────────────────────────────────
# Draws a base-R plot as a PNG image. Stand-in for renderPlot(), which
# currently fails in Safari under Shinylive. Use with imageOutput().
render_base <- function(id, session, draw) {
  renderImage({
    cd <- session$clientData
    w  <- cd[[paste0("output_", id, "_width")]]
    h  <- cd[[paste0("output_", id, "_height")]]
    req(w > 0, h > 0)
    pr <- cd$pixelratio
    if (is.null(pr)) pr <- 1
    f <- tempfile(fileext = ".png")
    png(f, width = w * pr, height = h * pr, res = 72 * pr)
    tryCatch(draw(), finally = dev.off())
    list(src = f, contentType = "image/png", width = w, height = h)
  }, deleteFile = TRUE)
}

# ── UI ───────────────────────────────────────────────────────────────────────
ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      body { background: #ffffff; font-size: 15px; color: #222; }

      .nav-tabs > li.active > a,
      .nav-tabs > li.active > a:hover,
      .nav-tabs > li.active > a:focus {
        color: #009499 !important;
        font-size: 12px; font-weight: bold;
        background-color: #f0fafa !important;
      }
      .nav-tabs > li > a {
        color: #666;
        font-size: 12px;
        padding: 5px 9px;
      }
      .nav-tabs > li > a:hover {
        color: #009499;
        background-color: #f0fafa;
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

      .ss-strip { display: flex; gap: 10px; margin-top: 4px; }
      .ss-strip .cell {
        flex: 1; background: #f7f7f7; border-radius: 3px; padding: 6px 10px;
      }
      .ss-strip .lab { font-size: 0.82em; color: #555; }
      .ss-strip .num { font-size: 24px; font-weight: 700; line-height: 1.2; }
      .ss-strip .note { font-size: 0.78em; color: #777; }

      h3 { color: #569BBD; font-size: 1.15em; margin-top: 0; }
      #cellCounts { margin-bottom: 0; }
      #cellCounts table { margin-bottom: 0; }
    "))
  ),
  
  titlePanel("Sums of Squares: Type I, II, and III"),
  
  tabsetPanel(id = "tabs",
              
              # ── TAB 1: Build ──────────────────────────────────────────────────────
              tabPanel("1. Build Data",
                       br(),
                       fluidRow(
                         column(6,
                                wellPanel(
                                  h3("Cell Counts"),
                                  tags$table(style = "width:100%; border:none;",
                                             tags$tr(
                                               tags$td(style = "vertical-align:top; padding-right:16px; width:1%;",
                                                       tableOutput("cellCounts")
                                               ),
                                               tags$td(style = "vertical-align:middle;",
                                                       div(class = "info-box",
                                                           "The table updates as you move the cell size sliders below.
               Try making the cells unequal to see how imbalance affects the SS types."
                                                       )
                                               )
                                             )
                                  )
                                ),
                                wellPanel(
                                  h3("Cell Sizes"),
                                  div(class = "info-box",
                                      "Set the number of observations in each of the four cells
               of the 2\u00d72 design independently. Equal cell sizes mean
               all three SS types give the same answer \u2014 unequal sizes
               are what make the choice matter."
                                  ),
                                  sliderInput("n11", "Cell A1 and B1:", min = 2, max = 30, value = 12, step = 1),
                                  sliderInput("n12", "Cell A1 and B2:", min = 2, max = 30, value = 4, step = 1),
                                  sliderInput("n21", "Cell A2 and B1:", min = 2, max = 30, value = 6, step = 1),
                                  sliderInput("n22", "Cell A2 and B2:", min = 2, max = 30, value = 14, step = 1),
                                  hr(),
                                  div(class = "section-label", "Seed (re-roll random data)"),
                                  sliderInput("seed", NULL, min = 1, max = 100, value = 42, step = 1)
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
                                  sliderInput("interaction", "Interaction (A \u00d7 B):", min = -40, max = 40, value = 10, step = 1),
                                  sliderInput("noise",       "Noise (within-group SD):",  min = 1,   max = 20, value = 5, step = 1),
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
              tabPanel("2. Data Plot",
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
                                           imageOutput("interactionPlotA", height = "300px")
                                         )
                                  ),
                                  column(6,
                                         wellPanel(
                                           imageOutput("interactionPlotB", height = "300px")
                                         )
                                  )
                                ),
                                uiOutput("mainEffectWarning"),
                                fluidRow(
                                  column(6,
                                         wellPanel(imageOutput("mainEffectA", height = "260px"))
                                  ),
                                  column(6,
                                         wellPanel(imageOutput("mainEffectB", height = "260px"))
                                  )
                                )
                         )
                       )
              ),
              
              # ── TAB 3: ANOVA Table ────────────────────────────────────────────────
              tabPanel("3. ANOVA Table",
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
              tabPanel("4. Visualise SS",
                       br(),
                       fluidRow(
                         column(4,
                                wellPanel(
                                  h3("Choose SS Type"),
                                  radioButtons("vizType", NULL,
                                               choices  = c("Type I — Sequential"    = "I",
                                                            "Type II — Hierarchical" = "II",
                                                            "Type III — Marginal"    = "III"),
                                               selected = "I",
                                  ),
                                  uiOutput("vizDescription")
                                )
                         ),
                         column(8,
                                div(class = "info-box",
                                    strong("How to read these plots:"),
                                    " Every sum of squares compares two models: a smaller one ",
                                    "without the term and a larger one with it. In each cell, the ",
                                    strong("dashed grey bar"), " is the predicted value (\u0176) from the smaller model and the ",
                                    strong("solid black bar"), " is the predicted value (\u0176) from the larger model. ",
                                    "The ", strong("coloured bar"), " is the gap between them. Every observation in that ",
                                    "cell shares the same gap, so it counts n times (×n). The SS is the sum of the ",
                                    "squared gaps for every observation. With unequal cell sizes, notice which cells ",
                                    "dominate each SS, and how the dashed bars move when you change Type."
                                ),
                                uiOutput("ssStrip")
                         )
                       ),
                       fluidRow(
                         column(12,
                                wellPanel(imageOutput("ssPanels", height = "420px"))
                         )
                       ),
                       fluidRow(
                         column(5,
                                wellPanel(
                                  h3("Where Does SS Total Go?"),
                                  imageOutput("partitionPlot", height = "360px")
                                )
                         ),
                         column(7,
                                wellPanel(
                                  h3("SS by Term — All Three Types"),
                                  imageOutput("ssBarPlot", height = "360px")
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
  output$interactionPlotA <- render_base("interactionPlotA", session, function() {
    draw_interaction(dat(), "A", "B",
                     input$showPoints, input$showCI, input$interaction)
  })
  
  # Interaction plot — Factor B on x-axis
  output$interactionPlotB <- render_base("interactionPlotB", session, function() {
    draw_interaction(dat(), "B", "A",
                     input$showPoints, input$showCI, input$interaction)
  })
  
  output$mainEffectA <- render_base("mainEffectA", session, function() {
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
  
  output$mainEffectB <- render_base("mainEffectB", session, function() {
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
  
  viz_ss <- reactive({
    switch(input$vizType, "I" = ss1(), "II" = ss2(), "III" = ss3())
  })
  
  ss_total <- reactive({
    d <- dat()
    sum((d$y - mean(d$y))^2)
  })
  
  output$ssPanels <- render_base("ssPanels", session, function() {
    d <- dat()
    if (is.null(viz_ss())) { plot.new(); return(invisible()) }
    draw_ss_panels(d, input$vizType)
  })
  
  output$partitionPlot <- render_base("partitionPlot", session, function() {
    ss <- viz_ss()
    if (is.null(ss)) { plot.new(); return(invisible()) }
    draw_partition(ss, ss_total(), input$vizType)
  })
  
  output$ssStrip <- renderUI({
    ss <- viz_ss()
    if (is.null(ss)) return(NULL)
    get_ss <- function(trm) ss$SS[ss$Term == trm]
    tot  <- ss_total()
    summ <- sum(ss$SS)
    cell <- function(lab, val, col, note = NULL) {
      div(class = "cell",
          div(class = "lab", HTML(lab)),
          div(class = "num", style = paste0("color:", col, ";"),
              formatC(val, format = "f", digits = 1)),
          if (!is.null(note)) div(class = "note", note))
    }
    gap_note <- if (abs(tot - summ) < 0.05) "= SS total" else
      paste0("SS total = ", formatC(tot, format = "f", digits = 1))
    div(class = "ss-strip",
        cell("SS<sub>A</sub>",   get_ss("A"),         term_cols["A"]),
        cell("SS<sub>B</sub>",   get_ss("B"),         term_cols["B"]),
        cell("SS<sub>A:B</sub>", get_ss("A:B"),       term_cols["A:B"]),
        cell("SS<sub>error</sub>", get_ss("Residuals"), term_cols["Error"]),
        cell("Sum of all four",  summ, "#222", gap_note))
  })
  
  output$vizDescription <- renderUI({
    switch(input$vizType,
           "I" = div(class = "type-desc",
                     strong("Type I (sequential):"),
                     " A is compared with a model containing nothing but the grand mean, so A
                     gets credit for everything it can explain. B is then judged only against
                     what A left over. The four SS add up exactly to SS total."),
           "II" = div(class = "type-desc",
                      strong("Type II (hierarchical):"),
                      " each main effect is judged after the other main effect is already in the
                      model (but without the interaction). With unequal cells, variation that A
                      and B could both explain is credited to neither, so the four SS no longer
                      add up to SS total."),
           "III" = div(class = "type-desc",
                       strong("Type III (marginal):"),
                       " each main effect is judged after everything else, including the
                       interaction. The smaller model for A forces the A1 and A2 averages of the
                       cell means to be equal, giving every cell equal weight whatever its n.
                       That is why small cells can get surprisingly large gaps. In the bar chart
                       below, the dashed box is variation not assigned to any term. Unlike Type II,
                       it is not purely the overlap between A and B, because the Type III
                       comparison models are less tidy; in some designs the terms can even add up
                       to more than SS total.")
    )
  })
  
  output$ssBarPlot <- render_base("ssBarPlot", session, function() {
    a1 <- ss1(); a2 <- ss2(); a3 <- ss3()
    if (is.null(a1) || is.null(a2) || is.null(a3)) { plot.new(); return(invisible()) }
    draw_ss_bars(a1, a2, a3)
  })
}

shinyApp(ui = ui, server = server)