# Post-communist collapses: year-specific dummies for the former Soviet
# republics (with Albania) and for the former Yugoslav republics.
#
# Each group has its own column: for a member country in its collapse window
# (*_COLLAPSE_WINDOWS), the year of the collapse ("1" for the first year, "2"
# for the second, ...), and "none" otherwise. Entered as
#   i(fsu_collapse_year, ref = "none") + i(yug_collapse_year, ref = "none")
# every collapse year gets its own shift, common to the group's countries and on
# top of the year fixed effects. Collapse years are counted from each country's
# own start, so a collapse that began earlier (Albania, 1990) shares the shifts
# of the same stage of the others' collapse (1991 onwards) rather than of the
# same calendar year. For the former Soviet and Yugoslav republics, which all
# start in 1991, collapse year k is calendar year 1990 + k.
#
# Former Soviet republics (the Baltic states included): 1991-1995, the years in
# which their average growth of GDP per capita was negative in PWT and WB (-9,
# -26, -12, -14 and -4 percent); every republic contracted in 1991-93, and their
# individual troughs fall between 1993 and 1999. Albania: 1990-1992, its
# contraction (-11, -33 and -7 percent in PWT, against the republics' average
# -9, -26 and -12 in their first three collapse years) before growth resumed
# in 1993.
#
# Former Yugoslav republics: 1991-1993, the contraction after the break-up
# (Croatia -23, -12, -7; Serbia -13, -33, -37; Montenegro -12, -27, -46 percent
# in PWT), with growth resuming in 1994. WB has no Serbia before 1996 and no
# Montenegro before 1998; Kosovo has no economic data in either panel.
FORMER_SOVIET <- c("ARM", "AZE", "BLR", "EST", "GEO", "KAZ", "KGZ", "LTU", "LVA",
                   "MDA", "RUS", "TJK", "TKM", "UKR", "UZB")
FORMER_YUGOSLAV <- c("SVN", "HRV", "BIH", "SRB", "MNE", "MKD", "XKO")
SOVIET_COLLAPSE_WINDOWS <- c(
  stats::setNames(rep(list(1991:1995), length(FORMER_SOVIET)), FORMER_SOVIET),
  list(ALB = 1990:1992)
)
YUGOSLAV_COLLAPSE_WINDOWS <- stats::setNames(rep(list(1991:1993), length(FORMER_YUGOSLAV)),
                                             FORMER_YUGOSLAV)

# Adds `column`: the collapse year (1, 2, ...) for a country in its collapse
# window, "none" otherwise.
add_collapse_dummies <- function(panel, windows, column) {
  stage <- mapply(function(country, year) {
    position <- match(year, windows[[country]])
    if (is.na(position)) "none" else as.character(position)
  }, panel$GID_0, panel$year, USE.NAMES = FALSE)
  panel[[column]] <- stage
  panel
}

add_soviet_collapse <- function(panel, windows = SOVIET_COLLAPSE_WINDOWS) {
  add_collapse_dummies(panel, windows, "fsu_collapse_year")
}

add_yugoslav_collapse <- function(panel, windows = YUGOSLAV_COLLAPSE_WINDOWS) {
  add_collapse_dummies(panel, windows, "yug_collapse_year")
}
