# MOPF
modelling psychometric functions

## R Code Status

### Logistic Regression

| File | Notes |
|------|-------|
|[one threshold](r/analysis/logreg_one.R) | should be ok |
|[one threshold non-linear](r/analysis/logreg_one_nl.R) | <ul><li>necessary for truncated prior</li><li>check non-linear model</li><li>param recov better than linear</li><li>[prior changes nothing](r/analysis/compare_nl.R)</li></ul> |
|[one threshold ranef](r/analysis/logreg_one_ranef.R) | <ul><li>should be ok</li><li>bad recovery of correlation</li></ul> |
|[pre/post threshold](r/analysis/logreg_prepost.R) | <ul><li>should be ok</li></ul> |
|[pre/post threshold ranef](r/analysis/logreg_prepost_ranef) | bad param recovery |

### Lapse Model Wichmann

| File | Notes |
|------|-------|
|[one threshold](r/analysis/wichmann_one.R) | should be ok |
|[one threshold ranef](r/analysis/wichmann_one_ranef.R) | <ul><li> conditional effects look very wrong with (recomended) beta priors for laps and guess</li><li> look ok with (improper) studentt priors</li><li>param recov still bad</li></ul> |
|[pre/post threshold](r/analysis/wichmann_prepost.R) | <ul><li>predictors on guess and lapse seem to improperly move psychometric function</li><li>conditional effects between 0 and 2</li><li>bad param recov</li></ul> |

### Data Simulation

| File | Notes |
|------|-------|
|[data_simulation](r/simulation/data_simulation.R) | <ul><li>check covariance matrix</li></ul> 
