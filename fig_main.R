
suppressPackageStartupMessages({
  require(data.table)
  require(ggplot2)
  require(ggrepel)
  require(patchwork)
})

#' want two panels, which provide two complementary "maps" of how to think
#' about introducing an Ab-testing scheme
#'  - how much would that scheme expand some baseline coverage level?
#'  - when does it make more sense to buy additional coverage with testing vs vaccination?
#'  
#' for these "maps", the reader would work from measurable values, navigate the "maps"
#' and find the answer to these questions
#' 
#' we're choosing the following the navigation scheme:
#'  - scan across to likely seroprevalence in a population
#'  - go up (perpendicular to x-axis) to test sensitivity contour line
#'  - go left (parallel to x-axis) to look up value on y-axis
#'  
#' As discussed in the manuscript, we present "maps" for
#' particular specificity & reference baseline coverage assumptions,
#' since the maps are to illustrate the result (i.e. the relationship equations)
#' not themselves the result
#' 
#' if you want different "maps", use different arguments when invoking this
#' script or change the defaults below / in the Makefile. Also, we choose
#' particular questions, but others are answerable with the relations after
#' a bit of algebra (e.g., given a test price, what is required sensitivity
#' for it to be a preferable purchase)

.defaults <- c(specificity = 0.99, coverage = 0.20)
.args <- if (interactive()) c(
  .defaults, "main.png"
) else commandArgs(trailingOnly = TRUE)

refspec <- as.numeric(.args[1])
covaxcoverage <- as.numeric(.args[2])

PPDmul <- function(TNR, TPR, seropos) TNR + seropos*(1+TPR-TNR)

targetTPRcontours = c(seq(1,9,by=2)/10, 0.95, 0.99)

#' relationship is linear in seropositive fraction,
#' so can just get line endpoints
seroposscale <- c(0.2, 0.8)

refCoverageMap <- data.table(expand.grid(
  TPR = targetTPRcontours,
  seropos = seroposscale
))

refCoverageMap[, TVratio := PPDmul(refspec, TPR, seropos)-1 ]
refCoverageMap[, coverage := covaxcoverage*(TVratio+1) ]

plottercore <- function(dt, aesydef, yl) ggplot(dt) +
  aes(x = seropos, color = "black") +
  aesydef + 
  geom_line(aes(group=TPR)) +
  geom_text_repel(
    aes(label = scales::percent(TPR, accuracy = 1)),
    data = function(ldt) ldt[seropos == max(seropos)],
    direction = "y",
    hjust = "left",
    size = 2,
    nudge_x = 0.01,
    box.padding = 0.01,
    xlim = c(NA, Inf),
    show.legend = FALSE
  ) +
  coord_cartesian(
    ylim = yl, xlim = range(seroposscale),
    clip = "off"
  ) +
  scale_x_continuous(
    "% Seropositive",
    minor_breaks = NULL,
    labels = scales::percent_format(accuracy = 1)
  ) + scale_color_manual(
    name = NULL,
    values = c(black = "black"),
    labels = sprintf("Test sensitivity, X%%\nspecificity %i%%", round(refspec*100))
  )

pA <- function(
  dt,
  yl.cov=c(0, covaxcoverage*2),
  ylres.cov=diff(yl.cov)/4,
  ...
) plottercore(dt, aes(y=coverage), yl.cov) +
  scale_y_continuous(
    "Total Immunized %\nAssuming COVAX-like Baseline Coverage",
    breaks = seq(yl.cov[1], yl.cov[2], by=ylres.cov),
    minor_breaks = NULL,
    labels = scales::percent_format(accuracy = 1)
  )

pB <- function(
  dt,
  yl.cost = c(0, 1),
  ylres.cost = diff(yl.cost)/4,
  ...
) plottercore(dt[TPR != 0.01], aes(y=TVratio), yl.cost) +
  scale_y_continuous(
    "Testing Cost % of\nVaccination Dose Cost",
    breaks = seq(yl.cost[1], yl.cost[2], by=ylres.cost),
    minor_breaks = NULL,
    labels = scales::percent_format(accuracy = 1)
  )

fullplotter <- function(...) {
  upper <- pA(...)
  lower <- pB(...)
  (upper / lower + plot_layout(guides = "collect") + plot_annotation(tag_levels = "A")) &
  theme_minimal(base_size = 10) & theme(
    plot.margin = margin(r = 1, unit = "line"),
    legend.position = "bottom"
  ) 
}

res.p <- fullplotter(dt = refCoverageMap)

ggsave(tail(.args, 1), res.p, height = 7.5, width = 4, dpi = 600, bg = "white")
