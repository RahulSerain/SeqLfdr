Sequential Large-Scale Multiple Testing
This repository contains R implementations and simulation programs for
sequential large-scale multiple testing procedures. The examples cover
several hypothesis-testing settings and illustrate the procedures
through Monte Carlo simulation.
Repository Structure
``` text
.
├── Functions/
│   ├── EstNull.func.R.txt
│   ├── Seq.adap.txt
│   ├── adaptZ.func.R.txt
│   ├── adpt.cutz.R.txt
│   ├── epsest.func.R.txt
│   ├── jin.cai.pi0.R.txt
│   ├── lfdr.gen.R.txt
│   └── lin.itp.R.txt
│
├── Example1_NormalMean/
│   ├── ex1_normal_mean_simulation.R
│   └── norm.hypo.R
├── Example2_Exponential/
│   ├── ex2_exponential_scale_simulation.R
│   └── exp.hypo.R
├── Example3_NormalVariance/
│   ├── Norm.sigma.hypo.R
│   └── ex3_normal_variance_simulation.R
├── Example4_BernoulliProp/
│   ├── Bern.Hypo.alt.R
│   └── ex4_bernoulli_simulation.R
├── Example5_One_Sample_t_test/
│   ├── ex5_one_sample_t_simulation.R
│   └── t test.R
├── Example6_CauchyLocation/
│   ├── cauchy.hypo.R
│   └── ex6_cauchy_location_simulation.R
└── Example7_Two_Sample_t_test/
    ├── ex7_two_sample_t_simulation.R
    └── t.2.sample.R
```
Common Functions
The `Functions/` directory contains functions shared by the simulation
examples:
`EstNull.func.R.txt` --- estimation of the null component.
`epsest.func.R.txt` --- estimation of the tuning parameter.
`lin.itp.R.txt` --- linear interpolation.
`adpt.cutz.R.txt` --- adaptive cutoff calculation.
`adaptZ.func.R.txt` --- adaptive procedure based on transformed test
statistics.
`Seq.adap.txt` --- sequential adaptive procedure.
`lfdr.gen.R.txt` --- local false discovery rate estimation.
`jin.cai.pi0.R.txt` --- estimation of the proportion of true null
hypotheses.
Simulation Examples
Example 1: Normal Mean
`Example1_NormalMean/` contains the hypothesis/data-generation functions
in `norm.hypo.R` and the simulation driver
`ex1_normal_mean_simulation.R`.
Example 2: Exponential Scale
`Example2_Exponential/` contains the model-specific functions in
`exp.hypo.R` and the simulation driver
`ex2_exponential_scale_simulation.R`.
Example 3: Normal Variance
`Example3_NormalVariance/` contains the model-specific functions in
`Norm.sigma.hypo.R` and the simulation driver
`ex3_normal_variance_simulation.R`.
Example 4: Bernoulli Proportion
`Example4_BernoulliProp/` contains the hypothesis/data-generation
functions in `Bern.Hypo.alt.R` and the simulation driver
`ex4_bernoulli_simulation.R`.
Example 5: One-Sample Student's t-Test
`Example5_One_Sample_t_test/` contains the testing functions in
`t test.R` and the simulation driver `ex5_one_sample_t_simulation.R`.
Example 6: Cauchy Location
`Example6_CauchyLocation/` contains the model-specific functions in
`cauchy.hypo.R` and the simulation driver
`ex6_cauchy_location_simulation.R`.
Example 7: Two-Sample Student's t-Test
`Example7_Two_Sample_t_test/` contains the testing functions in
`t.2.sample.R` and the simulation driver
`ex7_two_sample_t_simulation.R`.
Running a Simulation
Each example has a dedicated simulation driver. Open the required
example directory, check the simulation parameters, and run the
corresponding simulation script in R.
The simulation programs use the `parallel` package for Monte Carlo
replications. A typical setup is:
``` r
library(parallel)

ncore = detectCores() - 2

cl = makeCluster(ncore)

clusterSetRNGStream(cl, 12345)
```
The required common functions are sourced on the cluster workers by the
example-specific driver.
Monte Carlo Simulation
The simulation drivers evaluate the procedures over specified
configurations of:
number of hypotheses,
probabilities of the different hypothesis states,
model parameters,
error-rate parameters, and
initial sample size where applicable.
Multiple Monte Carlo replications are performed for each configuration.
The resulting quantities are summarized using Monte Carlo estimates and,
where implemented, Monte Carlo standard errors.
Reproducibility
Parallel random-number streams are initialized using:
``` r
clusterSetRNGStream(cl, 12345)
```
The seed may be changed when a different Monte Carlo realization is
required.
The number of parallel workers can be adjusted through:
``` r
ncore = detectCores() - 2
```
according to the available computing resources.
File Naming
The files in `Functions/` retain their `.txt` extensions where
applicable. These files contain R code and are loaded using `source()`.
The example directories separate model-specific functions from
simulation drivers, while the `Functions/` directory contains common
components used across examples.
Citation
If you use this repository in research, please cite the associated
methodological paper or research work for which these implementations
and simulation examples were developed.
