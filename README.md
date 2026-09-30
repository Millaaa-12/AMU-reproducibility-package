### Reproducibility package: Antimicrobial usage in slow-growing broiler chicken production cycles in northern Vietnam in 2022 to 2023
This repository contains the R code and supporting files for reproducing the analyses for:
**Antimicrobial usage in slow-growing broiler chicken production cycles in northern Vietnam in 2022 to 2023**

### Creator
[Mila (Chen) Xin] - ORCID [0009-0004-0257-7732] - [City University of Hong Kong]

### License
Code: BSD-3-Clause (see `LICENSE`, https://opensource.org/licenses/BSD-3-Clause).
Without a license, re-users would have no permission to use, modify or share this code (https://choosealicense.com/no-permission/).

### Dates
- Code last updated: [2026-09-29]
- Analysis last run (fresh R installation, one sitting): [2026-09-29]

### How to run
1. Install R [4.4.3] and JAGS [4.3.2] (https://mcmc-jags.sourceforge.io/).
2. Set the working directory to this folder.
3. `source("Master_Run_All.R")` - runs Step00 -> Step05 in one session and writes
   `log/Master_log.txt` (all printed results, sessionInfo(), package versions, run time).

| Order | File | Content |
|---|---|---|
| 0 | Step00_data_process.R | builds logistic/Poisson datasets |
| 1 | Step01_GAMM.R | Generalized Additive Mixed Models (GAMMs) |
| 2 | Step02_Markov_preventive.R | Bayesian two-state Markov model (SMP), preventive AMU |
| 3 | Step03_Markov_therapeutic.R | Bayesian two-state Markov model (SMP), therapeutic AMU |
| 4 | Step04_All_AMU_simulation.R | Joint posterior-predictive simulation of P-AMU, T-AMU and any AMU (PPC) and risk ratio of initiating preventive vs therapeutic AMU by broiler age |
| 5| Step05_AMU_quantification.R | AMU quantification, including UDD, mg/PCU, TI_DDD and TI_UDD

### File encoding
All scripts are UTF-8.

### Data and codebook
| File | Unit of observation | Desciption|
|---|---|---|
| `data/VN_Daily_farm_records.csv` | farm-day | Daily farm records |
| `VN_AMU_Quantification.csv` | farm-treatment-day administered drug | Antimicrobial-use quantification records |

`VN_AMU_Quantification.csv` may contain multiple records for one farm-day when multiple antimicrobial products or active substances are administered.

### Variables in daily farm records
| Variable | Description | Values / units |
|---|---|---|
| farm_id | farm identifier | 1-60 |
| Admin_days_cycle | day of age of the broiler | days |
| Drug_usage_final | any antimicrobial used that day | 0 = no, 1 = yes |
| preventive / therapeutic | antimicrobial used for prevention / terapeutic purpose | 0/1 |
| dead| the number of birds dead that day | count |
| alive_start_of_day | the number of alive birds at start of day | count |
| Days_cycle | length of the production cycle | days |

| *Derived* all_AMU | factor of Drug_usage_final (Step01) | 0/1 |
| *Derived* mortality_risk | dead / alive_start_of_day (Step01) | proportion |
| *Derived* event_id, duration | treatment-event number and length (Step00) | integer / days |


### Variables in antimicrobial-use quantification records
# Core identifiers and production-cycle variables

| Variable | Description | Values / units |
|---|---|---|
| `farm_id` | Unique farm identifier | 1–60 |
| `Days_cycle` | Total length of the production cycle | days |
| `Date_day_zero` | Calendar date corresponding to day 0 of the production cycle | date, `YYYY-MM-DD` |
| `Date_last_day` | Calendar date of the last production-cycle day | date, `YYYY-MM-DD` |
| `date` | Calendar date of the record or antimicrobial administration | date, `YYYY-MM-DD` |
| `day` | Broiler age on the record date, calculated from `Date_day_zero` | days |
| `predicted_weight_kg` | Predicted mean live body weight of a broiler on that day | kg/broiler |
| `entire_batch` | whether drugs has been given to birds in the full production batch or not |  0 = no, 1 = yes |

# Antimicrobial-use variables

| Variable | Description | Values / units |
|---|---|---|
| `Administration_antimicrobials` | The administrated antimicrobial | text, antimicrobial |
| `du_why_use` | Reported reason for antimicrobial use | categorical; e.g., preventive, therapeutic =1, preventive =2 |
| `du_ifthera_disease` | Reported disease or clinical indication when antimicrobial use is therapeutic | text or coded disease category; `NA` if not therapeutic |
| `route` | Route of antimicrobial administration | text; mixed in feed, mixed in drinking water, other |
| `brand` | Commercial product/brand name | text |
| `substance` | Antimicrobial active substance | standardized substance name |
| `concentration` | Product-label concentration as originally recorded | text; e.g., `%`, mg/g, or mg/mL |
| `concentration_mg_1gram_or_1ml` | Amount of active substance per gram or millilitre of product | mg/g or mg/mL |
| `concentration_colistin_IU_gram` | Colistin activity concentration, where applicable | IU/g; `NA` for non-colistin products |
| `du_amt_g_ml` | Unit used to record the administered product quantity | `g` or `mL` |
| `du_amt` | Total amount of antimicrobial product administered | g or mL, according to `du_amt_g_ml` |

# Feed- and water-administration variables

| Variable | Description | Values / units |
|---|---|---|
| `route_amtfeed_kg` | Amount of feed used to deliver the antimicrobial | kg |
| `route_amtwater_litres` | Amount of drinking water used to deliver the antimicrobial | litres |
| `the_concentration_in_diluted_amtwater_mgL` | Concentration of active substance in diluted medicated drinking water | mg/L |
| `waterintake_24h_litter` | Total flock drinking-water intake over 24 hours | litres/24 h |
| `waterintake_per_broiler_24h_ml` | Actual drinking-water intake per broiler over 24 hours | mL/broiler/24 h |
| `water_prepare_water_intake` | The comporision between the amount 24-hour drinking-water intake and the amount of water that used to deliver the antimicrobial | real number |

# Active-moiety and dose variables

| Variable | Description | Values / units |
|---|---|---|
| `Active_moiety_per_mg_of_derivative` | Conversion factor from antimicrobial derivative to active moiety | mg active moiety/mg derivative, or IU/mg where applicable |
| `Note_for_active_moiety_per_mg_of_derivative` | Notes or source supporting the active-moiety conversion factor | text |
| `total_active_flock_dose_mg_or_IU` | Total amount of active antimicrobial moiety administered to the flock on that treatment day | mg or IU/flock/day |
| `Dose_in_mg_or_IU_per_chicken` | Active antimicrobial dose administered per treated chicken | mg or IU/chicken/day |
| `Dose_in_mg_kg_day` | Used Daily Dose (UDD): active antimicrobial dose per kilogram body weight per day | mg/kg/day |

# Derived analytical variables and definitions

| Variable / indicator | Definition | Unit |
|---|---|---|
| `UDD` | Used Daily Dose; equivalent to `Dose_in_mg_kg_day` for each farm, treatment day, and antimicrobial. | mg/kg/day |
| `amount_mg` | `Dose_in_mg_kg_day × alive_start_of_day × predicted_weight_kg`; amount of active antimicrobial moiety used on a treatment day. | mg/day |
| `PCU_kg` | Population Correction Unit: `1 kg × N_broilers`. | kg |
| `mg_PCU` | Total active antimicrobial use divided by PCU: `Σ(amount_mg) / PCU_kg`. | mg/PCU |
| `kg_at_risk` | Chicken biomass at risk: mean body weight at treatment × number of broilers. | kg |
| `TI` | Treatment incidence: `total mg / (UDD or DDDvet × days_at_risk × kg_at_risk) × 100`. Interpreted as the number of UDDs or DDDvets per 100 bird-days at risk. | UDDs or DDDvets/100 bird-days at risk |

### Coding and quality rules

- Use `NA` for missing or not-applicable values; do not use `0` when the value is unknown.
- Record dates consistently as `YYYY-MM-DD`.
- Use standardized active-substance names in `substance`, while retaining the original product name in `brand`.
