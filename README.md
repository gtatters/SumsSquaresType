# Two-Way ANOVA Interaction Explorer

An interactive Shiny app for BIOL 3P96 (Biostatistics) at Brock University.

## What this app does

Two-way ANOVA asks whether two categorical factors (A and B) affect a
response variable, and whether they interact. An interaction means the
effect of one factor depends on the level of the other — the lines in an
interaction plot are non-parallel.

This app lets you simulate data with any combination of main effects and
an interaction, then explore the results through several lenses:
  
  - **Data Plot** — raw data with group means
- **Interaction Plots** — the classic lines plot, with Factor A or B on the x-axis
- **Marginal Means** — main-effect summaries (with a warning when an interaction
                                              is present)
- **Variance Partitioning** — Venn diagram showing how SS is divided among
Factor A, Factor B, the interaction, and residuals
- **Type I / II / III SS** — side-by-side comparison of the three sums-of-squares
approaches, with explanations of when each is appropriate

## How to use

1. Use the sliders to set the effect of Factor A, Factor B, and the interaction.
2. Adjust group sizes (balanced or unbalanced) to see how SS types diverge.
3. Switch tabs to explore different representations of the same data.
4. Read the explanations below each plot for guidance on interpretation.

## Learning goals

- Understand what an interaction is and how to spot it in a plot
- See why marginal means can be misleading when an interaction is present
- Understand the difference between Type I, II, and III sums of squares
- Connect the ANOVA table to the visual partitioning of variance

## Course context

Developed for BIOL 3P96 — Biostatistics, Brock University.
Built with R and Shiny (base R graphics only).