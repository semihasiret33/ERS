# CLAUDE.md - AI Assistant Guide for ERS Repository

## 📋 Repository Overview

**Repository Name:** ERS (Extreme Response Style)
**Full Title:** Aşırı Yanıt Tarzı (ERS) Düzeltme Yöntemlerinin PISA Verileri Üzerindeki Etkisi: Simülasyon ve Uygulama
**Purpose:** Araştırma projesi - ERS düzeltme yöntemlerinin (Z-skor, kovaryans, ML-MNRM, IRTree) PISA gibi büyük ölçekli değerlendirmelerdeki ülke ortalamaları ve sıralamaları üzerindeki etkilerini simülasyon ve gerçek veri analizi ile karşılaştırmak
**Status:** Proje kurulum aşaması
**Language:** R (primary), possibly Python for simulations

This document serves as a comprehensive guide for AI assistants (like Claude) working on this research project. It contains essential information about the project structure, statistical methods, conventions, and workflows to follow.

---

## 🏗️ Codebase Structure

```
ERS/
├── .git/                           # Git version control
├── CLAUDE.md                       # This file - AI assistant guide
├── README.md                       # Project overview (Turkish & English)
├── data/                           # Data files (not committed if large)
│   ├── raw/                        # Raw PISA data files
│   ├── processed/                  # Processed/cleaned data
│   └── simulated/                  # Simulated data from studies
├── R/                              # R source code
│   ├── simulation/                 # Simulation study scripts
│   │   ├── 01_data_generation.R    # ML-MNRM data generation
│   │   ├── 02_correction_methods.R # ERS correction implementations
│   │   └── 03_evaluation.R         # Bias and RMSE calculations
│   ├── analysis/                   # PISA real data analysis
│   │   ├── 01_data_preparation.R   # PISA data loading and cleaning
│   │   └── 02_ers_measurement_and_models.R  # ERS indices + Models 0-4
│   ├── functions/                  # Reusable utility functions
│   │   ├── ers_indices.R           # ERS calculation functions
│   │   ├── correction_methods.R    # Correction method implementations
│   │   ├── evaluation_metrics.R    # Bias, RMSE, correlation functions
│   │   └── irtree_functions.R      # IRTree model custom functions
│   └── visualization/              # Plotting and visualization scripts
│       ├── simulation_plots.R      # Simulation results plots
│       └── country_comparisons.R   # Country ranking visualizations
├── output/                         # Analysis outputs
│   ├── figures/                    # Generated plots and charts
│   ├── tables/                     # Result tables (CSV, LaTeX)
│   └── reports/                    # Generated reports
├── docs/                           # Documentation
│   ├── methodology.md              # Detailed methodology notes
│   ├── references.bib              # Bibliography
│   └── notes/                      # Research notes and decisions
├── tests/                          # Unit tests for R functions
└── renv/                           # R package management (if using renv)
```

### Directory Organization Principles

- **R/** - All R source code, organized by study phase (simulation vs. analysis)
- **data/** - Data files organized by type (raw, processed, simulated)
- **output/** - Generated outputs (figures, tables, reports)
- **docs/** - Documentation and research notes
- **tests/** - Unit tests for custom R functions

---

## 🔧 Technology Stack

### Core Technologies

- **Primary Language:** R (version 4.0+)
- **Statistical Packages:**
  - `mirt` - Multidimensional Item Response Theory models
  - `TAM` - Test Analysis Modules (for ML-MNRM)
  - `lavaan` - Latent Variable Analysis (for SEM/covariate models)
  - `tidyverse` - Data manipulation and visualization
  - `psych` - Psychological scales and ERS indices

### Data Analysis
- **PISA Data:** OECD PISA 2018/2022 student questionnaire data
- **Data Format:** SPSS (.sav) or CSV, converted to R data frames
- **Sample Size:** Multi-country datasets (5+ countries)

### Simulation Tools
- **ML-MNRM Implementation:** TAM package or custom JAGS/Stan models
- **Replication:** 1000 iterations per condition
- **Random Seed Management:** For reproducibility

### Dependency Management
- **renv** (recommended) - R environment management for reproducibility
- Alternative: Document required packages in `DESCRIPTION` or `packages.R`

### Version Control
- **Git** - Primary version control
- **Large File Storage:** Git LFS for large PISA datasets (if needed)

---

## 🚀 Development Workflows

### Git Workflow

**Branch Naming Convention:**
- Feature branches: `claude/claude-md-[session-id]` (for Claude AI sessions)
- Feature branches: `feature/[feature-name]`
- Bug fixes: `fix/[bug-description]`
- Hotfixes: `hotfix/[issue]`

**Commit Message Format:**
```
[Type]: Brief description

Detailed explanation if needed
- Bullet points for multiple changes
- Keep it clear and concise
```

**Types:**
- `feat`: New feature
- `fix`: Bug fix
- `docs`: Documentation changes
- `refactor`: Code refactoring
- `test`: Adding or updating tests
- `chore`: Maintenance tasks
- `style`: Code style changes (formatting, etc.)

### Pull Request Process

1. Create feature branch from main
2. Make changes with clear, atomic commits
3. Push to remote: `git push -u origin <branch-name>`
4. Create PR with descriptive title and summary
5. Wait for review and approval
6. Merge when approved

---

## 🤖 Key Conventions for AI Assistants

### Code Quality Standards

1. **Write Clean, Readable Code**
   - Use meaningful variable and function names
   - Keep functions small and focused (single responsibility)
   - Comment complex logic, but prefer self-documenting code
   - Follow DRY (Don't Repeat Yourself) principle

2. **Error Handling**
   - Always handle errors appropriately
   - Use try-catch blocks where necessary
   - Provide meaningful error messages
   - Log errors for debugging

3. **Testing**
   - Write tests for new features
   - Ensure existing tests pass before committing
   - Aim for meaningful test coverage (quality over quantity)

4. **Security**
   - Never commit secrets, API keys, or credentials
   - Use environment variables for sensitive data
   - Validate and sanitize all inputs
   - Follow security best practices (OWASP guidelines)

### File Operations

**Prefer Editing Over Creating:**
- Always check if a file exists before creating a new one
- Use Edit tool for modifications to existing files
- Only create new files when absolutely necessary

**Reading Before Writing:**
- Always read a file before modifying it
- Understand the existing code structure
- Maintain consistency with existing patterns

### Communication

**Be Concise:**
- Keep explanations clear and brief
- Focus on what changed and why
- Use bullet points for clarity

**Code References:**
- Reference specific files and line numbers
- Format: `file_path:line_number`
- Example: `src/services/auth.js:45`

---

## 📝 Common Development Tasks

### Setting Up R Environment

```r
# Install required packages
install.packages(c("tidyverse", "mirt", "TAM", "lavaan", "psych", "here"))

# If using renv for reproducibility
renv::init()
renv::snapshot()  # After installing all packages
renv::restore()   # To restore environment
```

### Running Simulation Study

```r
# Navigate to project root
setwd("/path/to/ERS")

# Run simulation scripts in order
source("R/simulation/01_data_generation.R")
source("R/simulation/02_correction_methods.R")
source("R/simulation/03_evaluation.R")

# Or run all at once
source("R/run_simulation_study.R")
```

### Running PISA Analysis

```r
# Run PISA analysis scripts in order
source("R/analysis/01_data_preparation.R")
source("R/analysis/02_ers_measurement_and_models.R")  # All 5 models (0-4)

# Generate visualizations
source("R/visualization/country_comparisons.R")
```

### Generating Outputs

```r
# Generate all figures
source("R/visualization/simulation_plots.R")
source("R/visualization/country_comparisons.R")

# Export tables
source("R/export_tables.R")
```

### Running Tests

```r
# Run unit tests (if using testthat)
testthat::test_dir("tests/")

# Or run specific test file
testthat::test_file("tests/test_ers_indices.R")
```

### Starting a New Analysis Component

```bash
# Create and switch to feature branch
git checkout -b feature/new-analysis-component

# Create new R script
touch R/analysis/08_new_component.R

# Make your changes, test, commit
git add R/analysis/08_new_component.R
git commit -m "feat: Add new analysis component"

# Push to remote
git push -u origin feature/new-analysis-component
```

---

## 🎯 Project-Specific Guidelines

### R Coding Conventions

**Naming Conventions:**
- **Variables:** `snake_case` (e.g., `ers_index`, `country_means`)
- **Functions:** `snake_case` with verb prefixes (e.g., `calculate_ers_index()`, `apply_zscore_correction()`)
- **Constants:** `UPPER_SNAKE_CASE` (e.g., `N_REPLICATIONS = 1000`)
- **Files:** `snake_case.R` with number prefixes for ordered scripts (e.g., `01_data_generation.R`)
- **Data frames:** `snake_case` (e.g., `pisa_data`, `simulation_results`)

**R Code Style (tidyverse style guide):**
- **Indentation:** 2 spaces (no tabs)
- **Line length:** Maximum 80 characters
- **Assignment:** Use `<-` not `=` for assignment
- **Quotes:** Use double quotes `"` for strings
- **Pipes:** Use `%>%` (magrittr) or `|>` (base R 4.1+) for data pipelines
- **Spacing:** Space after commas, around operators

**Example:**
```r
# Good
calculate_ers_index <- function(data, items) {
  ers_scores <- data %>%
    select(all_of(items)) %>%
    rowwise() %>%
    mutate(
      ers_index = sd(c_across(everything()), na.rm = TRUE)
    ) %>%
    pull(ers_index)

  return(ers_scores)
}

# Bad
calculateERSIndex<-function(data,items){ers_scores=data%>%select(all_of(items))%>%rowwise()%>%mutate(ers_index=sd(c_across(everything()),na.rm=TRUE))%>%pull(ers_index);return(ers_scores)}
```

### Statistical Analysis Guidelines

1. **Reproducibility:**
   - Always set random seeds: `set.seed(12345)`
   - Document package versions
   - Use `sessionInfo()` in output

2. **ERS Measurement:**
   - Use Greenleaf ERS Index as standard
   - Document alternative indices if tested
   - Handle missing data appropriately

3. **Model Naming:**
   - Model 0: `model0_baseline` or `raw_scores`
   - Model 1: `model1_zscore` or `zscore_corrected`
   - Model 2: `model2_covariate` or `ers_controlled`
   - Model 3: `model3_mlmnrm` or `irt_corrected`
   - Model 4: `model4_irtree` or `irtree_corrected`

4. **Evaluation Metrics:**
   - Calculate bias: `bias = mean(estimate - true_value)`
   - Calculate RMSE: `rmse = sqrt(mean((estimate - true_value)^2))`
   - Use rank-order correlation for country rankings

### Code Organization

**Script Structure:**
```r
# ============================================================================
# Script: 01_data_generation.R
# Purpose: Generate simulated data using ML-MNRM
# Author: [Name]
# Date: 2025-11-28
# ============================================================================

# 1. Setup ----
library(tidyverse)
library(TAM)

set.seed(12345)

# 2. Parameters ----
N_COUNTRIES <- 2
N_STUDENTS_PER_COUNTRY <- 500
N_ITEMS <- 10

# 3. Functions ----
generate_mlmnrm_data <- function(...) {
  # Function implementation
}

# 4. Main Analysis ----
# Analysis code here

# 5. Save Results ----
saveRDS(results, "data/simulated/simulation_data.rds")
```

### Documentation Standards

1. **Function Documentation:**
   - Use roxygen2 style comments
   - Document parameters, return values, examples

```r
#' Calculate Greenleaf ERS Index
#'
#' Computes the Greenleaf Extreme Response Style index for each respondent
#' based on their standard deviation across Likert items.
#'
#' @param data Data frame containing response data
#' @param items Character vector of item column names
#' @return Numeric vector of ERS index values
#' @references Greenleaf (1992)
#' @examples
#' ers_scores <- calculate_greenleaf_ers(pisa_data, c("item1", "item2"))
calculate_greenleaf_ers <- function(data, items) {
  # Implementation
}
```

2. **Inline Comments:**
   - Explain "why", not "what"
   - Use `# ----` for section breaks
   - Comment complex statistical procedures

---

## 🔍 Important Notes for AI Assistants

### Before Making Changes

1. **Read the relevant files first** - Never propose changes to code you haven't read
2. **Understand the statistical context** - Know the methodology and how changes affect statistical validity
3. **Check for existing patterns** - Follow established R coding conventions
4. **Search for related code** - Ensure consistency across simulation and analysis scripts

### When Implementing Statistical Analysis

1. **Avoid over-engineering** - Keep statistical methods simple and interpretable
2. **Don't add unrequested features** - Stick to the four models (0-3) unless explicitly asked
3. **Preserve reproducibility** - Never remove or change random seeds without updating documentation
4. **Statistical validity first** - Ensure corrections are methodologically sound

### R-Specific Considerations

1. **Package Dependencies:**
   - Check if packages are already loaded before adding new `library()` calls
   - Document any new package dependencies
   - Consider computational efficiency for large PISA datasets

2. **Data Handling:**
   - PISA data is large - use efficient data.table or tidyverse operations
   - Always check for missing values in survey data
   - Respect PISA sampling weights if required

3. **Model Estimation:**
   - ML-MNRM models can be computationally intensive
   - Consider running complex models on subsets first
   - Save intermediate results to avoid re-running lengthy computations

4. **Output Management:**
   - Always save results to appropriate directories (`output/`)
   - Use descriptive file names with dates: `simulation_results_2025-11-28.rds`
   - Export both R objects (.rds) and readable formats (.csv, .txt)

### Git Operations Best Practices

1. **Branch naming:** Must start with `claude/` for AI sessions
2. **Always use:** `git push -u origin <branch-name>`
3. **Retry on network errors:** Up to 4 times with exponential backoff (2s, 4s, 8s, 16s)
4. **Never skip hooks:** Unless explicitly requested
5. **Verify before force operations:** Never force push to main/master

### Commit Strategy

1. **Clear, descriptive messages** - Explain the "why" not just the "what"
2. **Atomic commits** - One logical change per commit
3. **Follow the commit format** - Use conventional commit types
4. **Review before committing** - Check `git status` and `git diff`

---

## 📚 Resources and Documentation

### Key References

**Primary Literature:**
1. Ulitzsch et al. (2023) - ERS in PISA 2015 questionnaire data
2. Lu & Bolt (2015) - ML-MNRM for ERS in PISA attitude-achievement paradox
3. Schoenmakers et al. (2023) - Model choice for ERS correction
4. Greenleaf (1992) - ERS index methodology

**PISA Resources:**
- [PISA Data Download](https://www.oecd.org/pisa/data/)
- [PISA 2018 Technical Report](https://www.oecd.org/pisa/data/pisa2018technicalreport/)
- [PISA Questionnaire Framework](https://www.oecd.org/pisa/pisaproducts/)

**R Package Documentation:**
- [mirt package](https://cran.r-project.org/web/packages/mirt/)
- [TAM package](https://cran.r-project.org/web/packages/TAM/)
- [lavaan package](https://lavaan.ugent.be/)

### Internal Documentation

- **Methodology:** `docs/methodology.md` - Detailed statistical methods
- **Bibliography:** `docs/references.bib` - Full reference list
- **Research Notes:** `docs/notes/` - Decision logs and analysis notes
- **Code Documentation:** Inline roxygen2 comments in R functions

---

## 🔄 Updating This Document

This document should be updated whenever:
- New technologies or frameworks are added
- Development workflows change
- New conventions are established
- Project structure evolves significantly

**Last Updated:** 2025-11-28
**Version:** 2.0.0 (Updated for ERS Research Project)

---

## 🎓 Quick Reference for Common Scenarios

### Scenario 1: Adapting R code from another study
1. Read the provided R code thoroughly
2. Identify which components map to this study (simulation vs. PISA analysis)
3. Adapt variable names to match our conventions (snake_case)
4. Update file paths to match our directory structure
5. Preserve statistical methods but modernize R syntax (tidyverse)
6. Document changes and rationale
7. Test adapted code with sample data
8. Commit with descriptive message

### Scenario 2: Implementing a new ERS correction method
1. Review literature reference for the method
2. Check if similar methods exist in `R/functions/correction_methods.R`
3. Implement as a standalone function with roxygen2 documentation
4. Add to both simulation and PISA analysis pipelines
5. Create unit tests in `tests/`
6. Update comparison scripts to include new method
7. Commit and push

### Scenario 3: Debugging simulation or model estimation
1. Check random seed settings for reproducibility
2. Verify input data structure and dimensions
3. Test with smaller sample size first
4. Check for convergence issues in IRT models
5. Review error messages for package-specific issues
6. Document solution in `docs/notes/`
7. Fix and commit

### Scenario 4: Generating results for manuscript
1. Ensure all analysis scripts have run successfully
2. Run visualization scripts to generate figures
3. Export tables in both CSV and LaTeX formats
4. Check output/ directory for all required files
5. Generate summary report with `sessionInfo()`
6. Commit outputs (or document how to reproduce)
7. Create tagged release for manuscript version

### Scenario 5: User provides R code to adapt
1. Save provided code to appropriate location (e.g., `data/external_code/`)
2. Read and analyze the code structure
3. Identify reusable functions vs. study-specific code
4. Map to our project structure (simulation/ vs. analysis/)
5. Refactor to match our conventions
6. Test thoroughly
7. Document adaptations in commit messages
8. Remove original external code after successful adaptation

---

**Remember:** This is a living document. As the ERS project grows and evolves, keep this guide updated to reflect the current state of the codebase and best practices.
