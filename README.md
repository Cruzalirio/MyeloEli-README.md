# Survival Analysis Using Parametric and Semiparametric Models

This project focuses on modeling survival data using both parametric and semiparametric techniques. The analysis employs regression models to assess the influence of various clinical variables on the time to the event of interest (e.g., survival or progression).

## Key Steps

### 1. **Importing Libraries and Loading Data**
   - **Libraries Used**: `survival`, `spBayesSurv`, `flexsurv`, and others for survival analysis and Bayesian modeling.
   - **Data**: Loaded from an Excel file, with specific patient health-related variables selected for analysis.

### 2. **Parametric Models (Weibull Model)**
   - **Weibull Survival Model**: Fitted using the `survreg` function in R, with `OS_Time` as the dependent variable (time to event) and `OS_Censor` as the censoring indicator. Predictor variables include `Pet_Ct`, `AGE`, `B2M`, `ALB`, etc.
   - **Interpretation**: Coefficients indicate the influence of each variable on survival. For example, a negative coefficient for `Pet_Ct` suggests that higher values are associated with shorter time to event.

### 3. **Bayesian Models**
   - **Bayesian Semiparametric Model**: Implemented using the `survregbayes` function, leveraging Markov Chain Monte Carlo (MCMC) methods for parameter inference. This approach handles parameter uncertainty and provides robust inferences.
   - **Variable Selection**: The model automatically selects relevant variables based on their influence on survival time, using a probability-based approach.

### 4. **Evaluation of Predictive Capacity**
   - **Survival Curves**: Generated using `ggsurvplot` to visualize survival probabilities as a function of categorical variables (e.g., `Pet_Ct`). Confidence intervals are included to show estimate uncertainty.

### 5. **Flexsurv Model**
   - **Flexible Survival Model**: Fitted using `flexsurvreg` to explore alternative survival time distributions and compare results with the traditional Weibull model.

## Summary
This project demonstrates a comprehensive approach to survival analysis, combining parametric (Weibull) and semiparametric (Bayesian) models to evaluate the impact of clinical variables on survival time. The analysis includes model fitting, variable selection, and predictive evaluation, with visualizations to aid interpretation.