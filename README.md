# predictability
A simple bayesian analysis of the compartimental SCIR model introduced in Reference: **Predictability: Can the turning point and end of an expanding epidemic be precisely forecast while the epidemic is still spreading?** using epidemiological data for Spain.
## Codes:
- All codes require the library rjags, invoked using: 
```{r}
require(rjags)
```
- Main codes:
  + jags-SCIR-realtime.R: Main code. Bayesian implementation of the SCIR model using JAGS and real-time data available while the epidemic was ongoing. Also contains auxiliary functions to plot the solution. 
  + jags-SCIR-synthetic.R: Same implementation but using synthetic data generated with the code deSolve-SCIR.R (see below).
  + jags-SCIR-forensic-analysis.R: Analysis of the end of the epidemic after the peak was reached using updated data by the end of july.
  + jags-gompertz.R: Bayesian implementation of the Gompertz model using updated data by the end of july.
 
- To generate the synthetic data we have used the library deSolve
```{r}
require(deSolve)
```
- Auxiliary codes:
  + plot.SCIR.linear.output.R: Auxiliary function to plot the output of jags-SCIR-realtime.R but using a linear vertical scale.
  + deSolve-SIQR.R: Integraion of the Ordinary Differential Equation for the SCIR model.


## Data:
All the datasets have been obtained from publicly available repositories. 
- Data sets used in real-time forecasting (late march), created from Spanish Ministry of Health and curated by DATADISTA: 
https://github.com/datadista/datasets/tree/master/COVID%2019/old_series
  + confirmed-march31.csv: Total number of confirmed cases (original data available until the 31st of march)
  + deaths-march31.csv: Total number of confirmed deaths (original data available until the 31st of march)
  + recovered-march31.csv: Total number of confirmed recovered (original data available until the 31st of march)


- Also splitted by spanish autonomous regions:
  + aut-regions-confirmed-march31.csv: Total number of confirmed cases (original data available until the 31st of march)
  + aut-regions-deaths-march31.csv: Total number of confirmed deaths (original data available until the 31st of march)
  + aut-regions-recovered-march31.csv: Total number of confirmed recovered (original data available until the 31st of march)
  + names.csv


- Updated Datasets (late july), created from Spanish Ministry of Health and curated by El-Pais-Data:
https://www.epdata.es/datos/coronavirus-china-datos-graficos/498
  + covid-19-es.csv: Different columns contain information about confirmed, recovered and deaths.

## References
- V1: https://arxiv.org/abs/2004.08842
- V2: 
