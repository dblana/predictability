# predictability
## Codes:
- All codes require the library rjags, invoked using: 
```{r}
dim(iris)
```
  + jags-SCIR-realtime.R
  + jags-SCIR-synthetic.R
  + jags-SCIR-forensic-analysis.R
  + jags-gompertz.R
  + plot.SCIR.linear.output.R
  + deSolve-SIQR.R


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
