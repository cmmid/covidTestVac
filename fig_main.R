
suppressPackageStartupMessages({
  require(data.table)
  require(ggplot2)
  require(metR)
  require(patchwork)
})

.args <- if (interactive()) c(
  "main.png"
) else commandArgs(trailingOnly = TRUE)

PPDmul <- function(TNR, TPR, seropos) TNR + seropos*(1+TPR-TNR)

#' if we think about contour lines correspond to PPDmul == 1,
#' these are the (exclusive) minimum TPRs required to achieve
#' any expansion in coverage
#' TPR = (1 - TNR)*(1/seropos - 1)
#' @examples
#' TPRcontour <- function(TNR, seropos) (1 - TNR)*(1/seropos - 1)
#' surface.dt <- data.table(expand.grid(seropos = seq(0,1,by=0.01), TNR = seq(0,1,by=0.01)))
#' surface.dt[, TPR := TPRcontour(TNR, seropos)]
#' ggplot(surface.dt[TPR <= 1]) + geom_contour_filled(aes(seropos, TNR, z=TPR))

#' alternatively, if we consider some baseline coverage & some fixed (high) TNR
#' how does TPR expand that coverage?
#' newcoverage = basecoverage*(TNR + seropos*(1+TPR-TNR)) =>
#' (newcoverage/basecoverage - TNR)/seropos - (1-TNR) = TPR

covaxcoverage <- .2
refspec <- .99
TPRcovcontour <- function(
  seropos, expcoverage, basecoverage = covaxcoverage, TNR = refspec
) ((expcoverage/basecoverage) - TNR)/seropos - (1-TNR)

seroposscale <- seq(0.2,0.8,by=0.0001)

ref <- data.table(expand.grid(
  seropos = seroposscale,
  expcoverage = seq(0, .5,by=0.0001)
))
ref[, TPR := TPRcovcontour(seropos, expcoverage) ]

TPRres <- function(x) fifelse(x < 0.5, (x/0.5)*0.2, 0.2+(x-0.5)/0.5*0.8)
inv.rescale <- function(x) fifelse(x < 0.2, x/0.2*0.5, 0.5+(x-0.2)/0.8*0.5)

ref[, TPRrescale := TPRres(TPR) ]

TPRlines <- data.table(expand.grid(
  seropos = seq(0.2, 0)
))

plottercore <- function(
  dt, aesydef
) ggplot(dt[between(TPR,0,1)]) +
  aes(x = seropos, fill = TPRrescale) +
  aesydef + geom_raster() +
  geom_contour(
    aes(z=TPR),
    breaks = c(0.1, 0.3, 0.5, 0.7, 0.9, 0.95, 0.99),
    show.legend = FALSE,
    color = "black"
  ) +
  geom_text_contour(
    aes(z=TPR, label = sprintf("%i%%", round(stat(level)*100))),
    color = "white", stroke.color = "black",
    breaks = c(0.1, 0.3, 0.5, 0.7, 0.9, 0.95, 0.99),
    skip = 0,
    show.legend = FALSE,
    stroke = 0.2, check_overlap = TRUE
  ) + 
  scale_x_continuous(
    "% Seropositive",
    minor_breaks = NULL,
    labels = scales::percent_format(accuracy = 1)
  ) +
  scale_fill_continuous(
    "Test Sensitivity",
    breaks = TPRres(c(0.1, 0.3, 0.5, 0.7, 0.9, 0.95, 0.99)),
    labels = function(x) scales::percent(inv.rescale(x), accuracy = 1)
  )

pA <- plottercore(ref, aes(y=expcoverage)) +
  coord_cartesian(
    xlim=c(0.2, 0.8),
    ylim=c(0, covaxcoverage*2),
    expand = FALSE
  ) +
  scale_y_continuous(
    "Total Immunized %, Assuming COVAX-like Baseline Coverage",
    breaks = seq(0.2, 0.4, by=0.05),
    minor_breaks = NULL,
    labels = scales::percent_format(accuracy = 1)
  )

#' (T/V + (1-TNR))/seropos - (1-TNR) <= TPR  

refB <- data.table(expand.grid(
  seropos = seroposscale,
  TVratio = seq(0, 1, by=0.001)
))

refB[, TPR := (TVratio + (1-refspec))/seropos - (1-refspec) ]
refB[, TPRrescale := TPRres(TPR) ]

pB <- plottercore(refB, aes(y=TVratio)) +
  coord_cartesian(
    xlim=c(0.2, 0.8),
    ylim=c(0, 1),
    expand = FALSE
  ) +
  scale_y_continuous(
    "Testing Cost % of Vaccination Dose Cost",
    breaks = seq(0, 1, by=0.25),
    minor_breaks = NULL,
    labels = scales::percent_format(accuracy = 1)
  )

res.p <- pA / pB &
  theme_minimal() &
  plot_layout(guides = "collect") &
  plot_annotation(tag_levels = "A")

# one color bar fine
# embiggen labels
# align contour labels - manually place at 65 (all except 95, 99 - put those at 70, 75)

ggsave(tail(.args, 1), res.p, height = 15, width = 8, dpi = 600, bg = "white")
