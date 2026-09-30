## =====================================================================
## Step05_AMU quantification.R
## Oral AMU quantification: UDD, mg/PCU, TI_UDD, TI_DDD
##
## Input : data/VN_AMU_Quantification.csv  - one row per farm x per treatment day x antimicrobial record per AMD
##         data/DDDvet.csv                  - EMA DDDvet (mg/kg/day) defined daily dose for broilers (DDDvet) 
## Output:                                    - mg/PCU by purpose
##         output/Table_UDD_by purpose.docx  - UDD by purpose (therapeutic / preventive)
##         output/AMU quantification.docx   - DDD, UDD, TI_UDD, TI_DDD per antimicrobial
##
## Key definitions used in this script
##   UDD         = Used Daily Dose (mg/kg/day), UDD per farm and per treatment day per AMD has been calculated and saved in the dataset. We aim to analyze their distribution.
##   amount_mg   = Dose_in_mg_kg_day x alive_start_of_day x predicted_weight_kg
##                 (= the amount of AMDs (active moiety) used on that treatment day)

##   mg/PCU      = total mg per farm / PCU_kg
##   The numerator of mg/PCU function -> Sum (amount_mg) = Sum (Dose_in_mg_kg_day x alive_start_of_day x predicted_weight_kg)
##                                      (= the amount of AMDs (active moiety) used on all treatment days)
##   The denominator of mg/PCU function -> PCU_kg = a standardized theoretical weight of 1 kg for a broiler chicken at treatment x N_broilers (Population Correction Unit)

##   TI          = total mg / (UDD or DDDvet x days_at_risk x kg_at_risk) x 100
##                 (= number of UDD or DDDvet per 100 bird-days at risk)
##   kg_at_risk (chicken biomass at risk) = mean treatment body weight (kg) x N_broilers
## =====================================================================

############## 1. Data import ##############
# Strings that should be read as missing values (NA), including Excel error codes
na_codes <- c("", "NA", "#VALUE!", "#N/A", "#DIV/0!", "NaN")

# Main AMU dataset: all oral treatment records
data1 <- read.csv("data/VN_AMU_Quantification.csv",
                  na.strings = na_codes,
                  stringsAsFactors = FALSE)   # keep text as character, not factor

# Reference DDDvet values (one row per active ingredient)
data_ddd <- read.csv("data/DDDvet.csv")


############## 2. Data cleaning ##############
to_num <- function(x, comma = c("thousands", "decimal", "none"),
                   var = deparse(substitute(x))) {
  comma <- match.arg(comma)                          # validate the 'comma' option
  if (is.numeric(x)) return(x)                       # already numeric: do nothing
  x_chr <- gsub("\\s", "", as.character(x))          # drop all spaces
  if (comma == "thousands") x_chr <- gsub(",", "", x_chr)
  if (comma == "decimal")   x_chr <- gsub(",", ".", x_chr)
  out <- suppressWarnings(as.numeric(x_chr))         # convert; we report problems ourselves
  bad <- !is.na(x_chr) & x_chr != "" & is.na(out)    # text that failed to convert
  if (any(bad)) warning(var, ": ", sum(bad), " non-numeric value(s) set to NA, e.g. ",
                        paste(head(unique(x_chr[bad]), 5), collapse = ", "))
  out
}

data1 <- data1 %>%
  mutate(
    # Drop leading/trailing spaces in drug names so joins match; empty -> NA
    Administration_antimicrobials = na_if(trimws(Administration_antimicrobials), ""),
    predicted_weight_kg = to_num(predicted_weight_kg, "decimal"),    # e.g. "0,85" -> 0.85
    alive_start_of_day  = to_num(alive_start_of_day,  "thousands"),  # e.g. "5,000" -> 5000
    Dose_in_mg_kg_day   = to_num(Dose_in_mg_kg_day,   "none"),
    Days_cycle          = to_num(Days_cycle,          "none"),
    # Purpose of use, coded du_why_use: 1 = therapeutic, 2 = preventive
    # (any other code, or a missing one, becomes NA)
    purpose = factor(trimws(as.character(du_why_use)), levels = c("1", "2"),
                     labels = c("Therapeutic", "Preventive"))
  )

# Clean the DDDvet reference table: trim names, make DDDvet numeric,
# drop empty names and exact duplicate rows
data_ddd <- data_ddd %>%
  mutate(Administration_antimicrobials = trimws(Administration_antimicrobials),
         DDDvet = to_num(DDDvet, "none")) %>%
  distinct(Administration_antimicrobials, DDDvet)

# Each antimicrobial must have exactly one DDDvet. 
dup_ddd <- unique(data_ddd$Administration_antimicrobials[
  duplicated(data_ddd$Administration_antimicrobials)])

############## 3. AMEG category ##############
# Lookup table: antimicrobial -> EMA AMEG category (B / C / D).
# The row order here sets the display order in all output tables.
amr_cat <- tibble::tribble(
  ~Administration_antimicrobials, ~Category,
  "Colistin","B","Enrofloxacin","B","Norfloxacin","B","Ceftiofur","B",
  "Florfenicol","C","Tilmicosin","C","Tylosin","C","Neomycin","C",
  "Gentamicin","C","Cefalexin","C","Apramycin","C","Lincomycin","C",
  "Azithromycin","C","Erythromycin","C",
  "Ampicillin","D","Amoxicillin","D","Doxycycline","D","Sulfamonomethoxine","D",
  "Sulfamonomethoxine_TMP","D","Sulfaclozine","D","Sulfadimidine","D","Trimethoprim","D",
  "Sulfadiazine","D","Sulfamethoxazole","D","Oxytetracycline","D","Sulfadimerazine","D",
  "Spectinomycin","D","Sulphaquinoxaline","D","Sulphaguanidine","D", "Bacitracin","D"
)
# Character vector of drug names in display order (used as factor levels)
amd_levels <- amr_cat$Administration_antimicrobials


############## 4. Treatment records ##############
# Keep only true treatment records (drug name AND dose present), then compute
# the amount of active ingredient given to the flock on that day (mg)
data_treat <- data1 %>%
  filter(!is.na(Administration_antimicrobials), !is.na(Dose_in_mg_kg_day)) %>%
  mutate(amount_mg = Dose_in_mg_kg_day * alive_start_of_day * predicted_weight_kg)


############## 5. Weight calculation ##############
# Mean treatment weight (kg), weighted by the number of birds alive on each treatment day: sum(weight x birds) / sum(birds).
weight_treat_kg <- data_treat %>%
  filter(!is.na(predicted_weight_kg), !is.na(alive_start_of_day)) %>%
  summarise(w = sum(predicted_weight_kg * alive_start_of_day) / sum(alive_start_of_day)) %>%
  pull(w)                                           # extract the single value as a number


############## 6. Farm-level denominators for both mg/PCU and TI caculation##############
# One row per farm (ALL farms in data1, including farms without any treatment):
#   N_broilers   : flock size = maximum number of birds alive during the cycle
#   days_at_risk : cycle length (Days_cycle) of that farm
#   PCU_kg       : population correction unit = 1 kg x N_broilers
#   kg_at_risk   : biomass at risk = mean treatment weight x N_broilers

weight_pcu_kg <- 1 # Standard weight for the PCU (kg per bird)

farm_denom <- data1 %>%
  group_by(farm_id) %>%
  summarise(
    # if every bird count is missing, return NA instead of -Inf from max()
    N_broilers   = if (all(is.na(alive_start_of_day))) NA_real_
    else max(alive_start_of_day, na.rm = TRUE),
    days_at_risk = first(na.omit(Days_cycle)),     # first non-missing cycle length
    .groups = "drop"
  ) %>%
  mutate(PCU_kg     = weight_pcu_kg   * N_broilers,
         kg_at_risk = weight_treat_kg * N_broilers)


############## 7. Formatting functions ##############
# Helpers that turn a numeric vector into one formatted text cell for tables.
# fmt_iqr(): "median (Q1–Q3)", d = decimals
fmt_iqr <- function(x, d = 2) {
  x <- x[is.finite(x)]; if (!length(x)) return("-")
  q <- unname(quantile(x, c(.25, .75)))            # 1st and 3rd quartiles
  sprintf("%.*f (%.*f\u2013%.*f)", d, median(x), d, q[1], d, q[2])   # \u2013 = en dash
}

# fmt_minmax(): "min-max"
fmt_minmax <- function(x, d = 2) {
  x <- x[is.finite(x)]; if (!length(x)) return("-")
  sprintf("%.*f-%.*f", d, min(x), d, max(x))
}

# fmt_meansd(): "mean ± SD"; only the mean if there is a single value (SD undefined)
fmt_meansd <- function(x, d = 1) {
  x <- x[is.finite(x)]; if (!length(x)) return("-")
  if (length(x) < 2) return(sprintf("%.*f", d, mean(x)))
  sprintf("%.*f \u00b1 %.*f", d, mean(x), d, sd(x))                 # \u00b1 = ±
}

# fmt_range(): "min – max" (with spaces and an en dash)
fmt_range <- function(x, d = 1) {
  x <- x[is.finite(x)]; if (!length(x)) return("-")
  sprintf("%.*f \u2013 %.*f", d, min(x), d, max(x))
}

# fmt_ms(): shortcut for "mean ± SD" with 1 decimal
fmt_ms <- function(x) fmt_meansd(x, d = 1)

# fmt_mr(): "median (min-max)"
fmt_mr <- function(x, d = 1) {
  x <- x[is.finite(x)]; if (!length(x)) return("-")
  sprintf("%.*f (%.*f-%.*f)", d, median(x), d, min(x), d, max(x))
}


############## 8. UDD per antimicrobial ##############
# Distribution of the used daily dose (mg/kg/day) across all treatment records of each antimicrobial. 
#UDD_mean and UDD_median are used later for TI_UDD.
udd_summary <- data_treat %>%
  group_by(Administration_antimicrobials) %>%
  summarise(
    UDD_mean   = mean(Dose_in_mg_kg_day),
    UDD_median = median(Dose_in_mg_kg_day),
    UDD_min    = min(Dose_in_mg_kg_day),
    UDD_max    = max(Dose_in_mg_kg_day),
    n          = n(),                              # number of treatment records
    .groups = "drop"
  )

# Overall UDD per antimicrobial, rounded for display
udd_table_overall <- udd_summary %>%
  transmute(
    Administration_antimicrobials,
    `UDD (mean)`            = round(UDD_mean, 2),
    `UDD (median, min-max)` = paste0(round(UDD_median, 2),
                                     " (", round(UDD_min, 2), "-", round(UDD_max, 2), ")")
  )

# UDD per antimicrobial split by purpose (long format: one row per drug x purpose)
udd_long <- data_treat %>%
  filter(!is.na(purpose)) %>%
  group_by(Administration_antimicrobials, purpose) %>%
  summarise(med_iqr = fmt_iqr(Dose_in_mg_kg_day, 2),
            min_max = fmt_minmax(Dose_in_mg_kg_day, 2),
            .groups = "drop")

# Wide format: one row per drug; therapeutic and preventive columns side by side
udd_table_purpose <- udd_long %>%
  pivot_wider(id_cols = Administration_antimicrobials,
              names_from = purpose,
              values_from = c(med_iqr, min_max),
              names_glue = "{purpose}_{.value}",   # e.g. "Therapeutic_med_iqr"
              names_expand = TRUE) %>%             # keep columns even if a purpose has no data
  transmute(                                       # rename to short column names
    amd      = Administration_antimicrobials,
    ther_med = Therapeutic_med_iqr, ther_mm = Therapeutic_min_max,
    prev_med = Preventive_med_iqr,  prev_mm = Preventive_min_max
  ) %>%
  mutate(across(-amd, ~replace_na(.x, "-"))) %>%   # drug never used for a purpose -> "-"
  left_join(amr_cat, by = c("amd" = "Administration_antimicrobials")) %>%  # add AMEG category
  mutate(Category = replace_na(Category, "-")) %>% # drug not in amr_cat -> "-"
  arrange(factor(amd, levels = amd_levels)) %>%    # order B -> C -> D as in amr_cat
  relocate(Category, .after = amd)                 # Category as 2nd column

grp_rows <- which(udd_table_purpose$Category != dplyr::lead(udd_table_purpose$Category))

# Word table with a two-level header:
ft_udd <- flextable(udd_table_purpose) %>%
  set_header_labels(amd = "", Category = "",                      # lower header row
                    ther_med = "UDD median (IQR)", ther_mm = "UDD (min-max)",
                    prev_med = "UDD median (IQR)", prev_mm = "UDD (min-max)") %>%
  add_header_row(values = c("AMD (active ingredient)", "Category",  # upper header row
                            "Therapeutic purpose", "Preventive purpose"),
                 colwidths = c(1, 1, 2, 2)) %>%
  merge_at(i = 1:2, j = 1, part = "header") %>%    # merge both header rows in column 1
  merge_at(i = 1:2, j = 2, part = "header") %>%    # ... and in column 2
  theme_booktabs() %>%                             # clean academic style
  align(align = "center", part = "all") %>%
  align(j = 1, align = "left", part = "all") %>%   # drug names left-aligned
  valign(valign = "center", part = "header")
if (length(grp_rows) > 0)                          # draw B/C/D separator lines if any
  ft_udd <- hline(ft_udd, i = grp_rows, border = fp_border(width = 1.2))
ft_udd <- autofit(ft_udd)                          # fit column widths to content

save_as_docx(ft_udd, path = "output/Table_UDD_by purpose.docx")


############## 9. mg/PCU by purpose ##############
# Farm x purpose: total mg / PCU.
mg_pcu_purpose <- data_treat %>%
  filter(!is.na(purpose)) %>%                      # records with a known purpose only
  group_by(farm_id, purpose) %>%
  summarise(total_mg = sum(amount_mg, na.rm = TRUE), .groups = "drop") %>%
  complete(farm_id = farm_denom$farm_id, purpose, fill = list(total_mg = 0)) %>%
  left_join(farm_denom %>% select(farm_id, PCU_kg), by = "farm_id") %>%
  mutate(mg_PCU = total_mg / PCU_kg, purpose = as.character(purpose))

# Summary across farms: median (IQR), mean ± SD and range of mg/PCU per purpose
mg_pcu_summary <- mg_pcu_purpose %>%
  mutate(purpose = factor(purpose, c("Preventive", "Therapeutic"))) %>%  # display order
  group_by(`Purpose of usage` = purpose) %>%
  summarise(`Median (IQR)`      = fmt_iqr(mg_PCU, 1),
            `Mean ± SD`         = fmt_meansd(mg_PCU, 1),
            `Range (min - max)` = fmt_range(mg_PCU, 1),
            .groups = "drop") %>%
  arrange(`Purpose of usage`)
print(mg_pcu_summary)


############## 10. TI_UDD per antimicrobial ##############
# Farm x antimicrobial: total amount used (mg), joined with the farm's denominators (flock size, cycle length, kg at risk)
farm_amd <- data_treat %>%
  group_by(farm_id, Administration_antimicrobials) %>%
  summarise(total_amount_mg = sum(amount_mg, na.rm = TRUE), .groups = "drop") %>%
  left_join(farm_denom %>% select(farm_id, N_broilers, days_at_risk, kg_at_risk),
            by = "farm_id")

# Treatment incidence based on UDD (days treated per 100 days at risk), using the mean and the median UDD of each antimicrobial
ti_udd <- farm_amd %>%
  left_join(udd_summary %>% select(Administration_antimicrobials, UDD_mean, UDD_median),
            by = "Administration_antimicrobials") %>%
  mutate(
    TI_UDD_mean   = total_amount_mg / (UDD_mean   * days_at_risk * kg_at_risk) * 100,
    TI_UDD_median = total_amount_mg / (UDD_median * days_at_risk * kg_at_risk) * 100
  )


############## 11. TI_DDD per antimicrobial ##############
# Same formula as TI_UDD, but with the EMA DDDvet as the standard daily dose.
# Drugs without a DDDvet get TI_DDD = NA (shown as "-" in the table).
ti_ddd <- farm_amd %>%
  left_join(data_ddd, by = "Administration_antimicrobials") %>%
  mutate(TI_DDD = total_amount_mg / (DDDvet * days_at_risk * kg_at_risk) * 100)


############## 12. Summary table: DDD, UDD, TI_UDD, TI_DDD per antimicrobial ##############
# UDD display text: "median (min-max)"
udd_disp <- udd_summary %>%
  transmute(Administration_antimicrobials,
            UDD = sprintf("%.1f (%.1f-%.1f)", UDD_median, UDD_min, UDD_max))

# TI_UDD across the farms that used each drug: mean ± SD and median (min-max)
udd_tab <- ti_udd %>%
  group_by(Administration_antimicrobials) %>%
  summarise(TI_UDD_ms = fmt_ms(TI_UDD_mean),
            TI_UDD_mr = fmt_mr(TI_UDD_mean),
            .groups = "drop")

# TI_DDD across the farms that used each drug: mean ± SD and median (min-max)
ddd_tab <- ti_ddd %>%
  group_by(Administration_antimicrobials) %>%
  summarise(TI_DDD_ms = fmt_ms(TI_DDD),
            TI_DDD_mr = fmt_mr(TI_DDD),
            .groups = "drop")

# DDD column = DDDvet reference value
add_tab <- data_ddd %>% rename(DDD = DDDvet)

# Start from amr_cat so every listed drug appears once, in B -> C -> D order,
# then attach all indicators by drug name
final_tab <- amr_cat %>%
  left_join(add_tab,  by = "Administration_antimicrobials") %>%
  left_join(udd_disp, by = "Administration_antimicrobials") %>%
  left_join(udd_tab,  by = "Administration_antimicrobials") %>%
  left_join(ddd_tab,  by = "Administration_antimicrobials") %>%
  mutate(Administration_antimicrobials =
           factor(Administration_antimicrobials, levels = amd_levels)) %>%
  arrange(Category, Administration_antimicrobials)

# Safety check: each drug must appear only once (no duplication from joins)
stopifnot(!anyDuplicated(final_tab$Administration_antimicrobials))

# Compact version: one row per drug, median (min-max) only
final_one_row <- final_tab %>%
  filter(!is.na(UDD)) %>%                          # drop drugs without UDD
  select(Administration_antimicrobials, Category, DDD, UDD,
         TI_UDD = TI_UDD_mr, TI_DDD = TI_DDD_mr) %>%
  mutate(DDD = as.character(DDD)) %>%              # numeric -> text, so missing DDD can show "-"
  mutate(across(where(is.character), ~replace_na(.x, "-")))  # remaining NA -> "-"

# Export to Word (Times New Roman, 9 pt, booktabs style)
ft <- flextable(final_one_row) %>%
  set_header_labels(
    Administration_antimicrobials = "Antimicrobial",
    Category = "Category",
    DDD      = "DDD",
    UDD      = "UDD\n(median, min-max)",
    TI_UDD   = "TI_UDD\n(median, min-max)",
    TI_DDD   = "TI_DDD\n(median, min-max)"
  ) %>%
  theme_booktabs() %>%
  align(align = "center", part = "all") %>%
  align(j = 1, align = "left", part = "all") %>%
  fontsize(size = 9, part = "all") %>%
  font(fontname = "Times New Roman", part = "all") %>%
  autofit()

save_as_docx(ft, path = "output/AMU quantification_table2.docx")