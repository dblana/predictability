# predictability
A simple bayesian analysis of the compartimental SCIR model introduced in Reference: **Predictability: Can the turning point and end of an expanding epidemic be precisely forecast while the epidemic is still spreading?** using epidemiological data for Spain.
## Codes:
- All codes require the library rjags, invoked using: 
```{r}
require(rjags)
```
- Main codes:
  + [jags-SCIR-realtime.R](https://github.com/mariocastro73/predictability/blob/master/jags-SCIR-realtime.R): Main code. Bayesian implementation of the SCIR model using JAGS and real-time data available while the epidemic was ongoing. Also contains auxiliary functions to plot the solution. 
  + [jags-SCIR-synthetic.R](https://github.com/mariocastro73/predictability/blob/master/jags-SCIR-synthetic.R): Same implementation but using synthetic data generated with the code deSolve-SCIR.R (see below).
  + [jags-SCIR-post-peak-analysis.R](https://github.com/mariocastro73/predictability/blob/master/jags-SCIR-post-peak-analysis.R): Analysis of the end of the epidemic after the peak was reached using updated data by the end of July.
  + [jags-gompertz.R](https://github.com/mariocastro73/predictability/blob/master/jags-gompertz.R): Bayesian implementation of the Gompertz model using updated data by the end of July.
 
- To generate the synthetic data we have used the library deSolve
```{r}
require(deSolve)
```
- Auxiliary codes:
  + [plot.SCIR.linear.output.R](https://github.com/mariocastro73/predictability/blob/master/plot.SCIR.linear.output.R): Auxiliary function to plot the output of jags-SCIR-realtime.R but using a linear vertical scale.
  + [deSolve-SCIR.R](https://github.com/mariocastro73/predictability/blob/master/deSolve-SCIR.R): Integraion of the Ordinary Differential Equation for the SCIR model.


## Data:
All the datasets (see [datasets folder](https://github.com/mariocastro73/predictability/tree/master/datasets)) have been obtained from publicly available repositories. 
- Data sets used in real-time forecasting (late March), officially reported by Spanish Ministry of Health and curated by DATADISTA: 
https://github.com/datadista/datasets/tree/master/COVID%2019/old_series
  + confirmed-march31.csv: Total number of confirmed cases (original data available until March 31st)
  + deaths-march31.csv: Total number of confirmed deaths (original data available until March 31st)
  + recovered-march31.csv: Total number of confirmed recovered (original data available until March 31st)


- Also splitted by Spanish Autonomous Regions:
  + aut-regions-confirmed-march31.csv: Total number of confirmed cases (original data available until March 31st)
  + aut-regions-deaths-march31.csv: Total number of confirmed deaths (original data available until March 31st)
  + aut-regions-recovered-march31.csv: Total number of confirmed recovered (original data available until March 31st)
  + names.csv


- Updated Datasets (post-peak data), officially reported by Spanish Ministry of Health and curated by El-Pais-Data by late July, 2020:
https://www.epdata.es/datos/coronavirus-china-datos-graficos/498
  + covid-19-es.csv: Different columns contain information about confirmed, recovered and deaths.

**Note:** The data series for recovered was interrupted on 17th May. 

## Sample outputs
In the [output folder](https://github.com/mariocastro73/predictability/tree/master/output), we include some samples produced by the code and used to generate the figures in the paper.

## References
- V1: https://arxiv.org/abs/2004.08842v1
- V2: https://arxiv.org/abs/2004.08842 
