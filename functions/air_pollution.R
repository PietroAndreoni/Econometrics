# National sulphur dioxide emissions, used as a proxy for changes in air
# pollution.
#
# data/so-emissions-by-world-region-in-million-tonnes.csv (CEDS via Our World in
# Data): Entity, Code, Year, "Sulfur dioxide". Despite the file name the values
# are tonnes of SO2 (USA 1980: 25.9 million). Rows without an ISO3 code and the
# OWID_* aggregates are dropped; OWID codes Kosovo as OWID_KOS, so it is lost.
#
# dlSO2 is the log change in emissions between consecutive years, split into
# its increases dlSO2_pos = max(dlSO2, 0) and decreases dlSO2_neg =
# min(dlSO2, 0), so dlSO2 = dlSO2_pos + dlSO2_neg; <var>_l<k> is the k-th lag of
# each. Each is missing when a year it needs is absent or its emissions are not
# positive; nothing is filled.
load_so2_emissions <- function(lags = 1:3) {
  raw <- utils::read.csv(
    data_file("so-emissions-by-world-region-in-million-tonnes.csv"),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  so2 <- tibble::tibble(
    GID_0 = raw$Code,
    year = as.integer(raw$Year),
    SO2 = as.numeric(raw[["Sulfur dioxide"]])
  ) %>%
    filter(nzchar(GID_0), !startsWith(GID_0, "OWID_")) %>%
    arrange(GID_0, year) %>%
    group_by(GID_0) %>%
    mutate(
      lSO2 = if_else(SO2 > 0, log(SO2), NA_real_),
      dlSO2 = lSO2 - lag_by_year(lSO2, year, 1),
      dlSO2_pos = pmax(dlSO2, 0),
      dlSO2_neg = pmin(dlSO2, 0)
    )
  for (v in c("dlSO2", "dlSO2_pos", "dlSO2_neg")) {
    for (k in lags) {
      so2 <- mutate(so2, "{v}_l{k}" := lag_by_year(.data[[v]], year, k))
    }
  }
  ungroup(so2)
}
