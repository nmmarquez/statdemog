# Homework 4  
## Neal Marquez  



2) Use the Lee-Carter method to obtain probabilistic forecasts of mortality in 
Mexico for 2015-2020.  Obtain probabilistic forecasts of the mortality index,
kt.  Hence find a confidence interval for the age-specific mortality rate, nmx,
for women aged 75-80 in 2015-2020.

The Lee Crater parameters were calculated in two seperate ways, once by least
squares and another by SVD. The models proced very similiar mean forecasts 
however gretaly differed in their uncertainty as seen by the standard deviation
in the random walk of the $k_t$ terms which were `.0522` & `.3018` for the least
squares and SVD methods respectively. The plots for each forecast method are 
shown against the  actual mexico female mortality data in log space below.

![](/home/nmarquez/Documents/Classes/statdemog/week5/allleecarterls.jpg "")
![](/home/nmarquez/Documents/Classes/statdemog/week5/allleecartersvd.jpg "")

The uncetainty from the random walk model for each age group is shown below 
for the least squares model fit and the svd model fit respectively where ymin
is the 2.5% confidence interval and ymax is the 97.5% confidence interval.

### Least Squares Model fit Estimates and Uncertainty

|     lnmxt| age| year|model |      ymin|       ymax|
|---------:|---:|----:|:-----|---------:|----------:|
| -4.154011|   0| 2015|ls    | -4.256414| -4.0516086|
| -7.366427|   1| 2015|ls    | -7.468829| -7.2640241|
| -8.436859|   5| 2015|ls    | -8.539261| -8.3344560|
| -8.574348|  10| 2015|ls    | -8.676751| -8.4719459|
| -8.281179|  15| 2015|ls    | -8.383581| -8.1787760|
| -7.970710|  20| 2015|ls    | -8.073113| -7.8683078|
| -7.716431|  25| 2015|ls    | -7.818834| -7.6140286|
| -7.415483|  30| 2015|ls    | -7.517886| -7.3130805|
| -7.022140|  35| 2015|ls    | -7.124542| -6.9197373|
| -6.559744|  40| 2015|ls    | -6.662147| -6.4573413|
| -6.080090|  45| 2015|ls    | -6.182493| -5.9776874|
| -5.600567|  50| 2015|ls    | -5.702969| -5.4981641|
| -5.132307|  55| 2015|ls    | -5.234710| -5.0299047|
| -4.671966|  60| 2015|ls    | -4.774369| -4.5695635|
| -4.217219|  65| 2015|ls    | -4.319621| -4.1148162|
| -3.760514|  70| 2015|ls    | -3.862916| -3.6581110|
| -3.306348|  75| 2015|ls    | -3.408750| -3.2039452|
| -2.778054|  80| 2015|ls    | -2.880457| -2.6756517|
| -2.324761|  85| 2015|ls    | -2.427164| -2.2223587|
| -1.912305|  90| 2015|ls    | -2.014708| -1.8099025|
| -1.537097|  95| 2015|ls    | -1.639500| -1.4346948|
| -1.090144| 100| 2015|ls    | -1.192547| -0.9877415|

### SVD Model fit Estimates and Uncertainty

|     lnmxt| age| year|model |      ymin|       ymax|
|---------:|---:|----:|:-----|---------:|----------:|
| -4.121788|   0| 2015|svd   | -4.713411| -3.5301656|
| -7.314561|   1| 2015|svd   | -7.906183| -6.7229382|
| -8.393671|   5| 2015|svd   | -8.985294| -7.8020485|
| -8.538055|  10| 2015|svd   | -9.129678| -7.9464326|
| -8.243947|  15| 2015|svd   | -8.835569| -7.6523240|
| -7.932990|  20| 2015|svd   | -8.524613| -7.3413676|
| -7.679733|  25| 2015|svd   | -8.271355| -7.0881099|
| -7.380616|  30| 2015|svd   | -7.972239| -6.7889937|
| -6.990557|  35| 2015|svd   | -7.582179| -6.3989341|
| -6.532980|  40| 2015|svd   | -7.124602| -5.9413571|
| -6.057571|  45| 2015|svd   | -6.649193| -5.4659481|
| -5.581362|  50| 2015|svd   | -6.172985| -4.9897394|
| -5.115113|  55| 2015|svd   | -5.706735| -4.5234902|
| -4.656113|  60| 2015|svd   | -5.247736| -4.0644902|
| -4.202307|  65| 2015|svd   | -4.793930| -3.6106845|
| -3.746838|  70| 2015|svd   | -4.338461| -3.1552158|
| -3.294020|  75| 2015|svd   | -3.885642| -2.7023970|
| -2.762654|  80| 2015|svd   | -3.354277| -2.1710316|
| -2.309824|  85| 2015|svd   | -2.901447| -1.7182013|
| -1.899113|  90| 2015|svd   | -2.490735| -1.3074900|
| -1.526774|  95| 2015|svd   | -2.118397| -0.9351517|
| -1.086289| 100| 2015|svd   | -1.677912| -0.4946665|

3) The purpose of this question is to give you some practice in fitting the
simplest form of Bayesian hierarchical model, namely the random intercept model,
or Bayesian random effects one-way analysis of variance model.
This question uses data from the 1975 U.S. Sustaining Effects Study of
elementary education, available as the `egsingle` dataset in the `mlmRev` R
package. This gives data on 1,721 students in 60 schools. We will take
math (Mathematics achievement score) as the outcome variable, and schoolid
(the code for the school the student attends) as the grouping variable. 
We will use only the data for year 0.5 (for which there are data on 1672
students). Your task is to estimate the Bayesian random effects 
one-way analysis of variance model for these data and interpret the results.

a. Write out the Bayesian random effects one-way analysis of variance model for
these data. What are the unknown parameters to be estimated?

This hierarchical bayesian model will model the mean of all students across
all schools, the variance across all students, as well as the variation
accounted for by the effect of the schools themselves. These parameters are 
$\mu_{\alpha}$, $\sigma_{\epsilon}$, and $\sigma_{\alpha}$ respectively. The
model will also produce a fit for $\alpha_j$ which are a school specific 
effects that deviate away from $\mu_{\alpha}$ with a standard deviation of 
$\sigma_{\alpha}$. The data that this model is fitting is math test scores
$y_{i,j}$ where $i$ is the student and $j$ is the school. The functional
form of the model is as follows

$$
y_{i,j} = \alpha_{j} + \epsilon_{i,j} \\
\epsilon_{i,j} \stackrel{iid}{\sim} \mathcal{N}(0, \sigma^{2}_{\epsilon}) \\
\alpha_{j} \stackrel{iid}{\sim} \mathcal{N}(\mu_{\alpha},\sigma^{2}_{\alpha}) \\
$$

b. Specify a reasonable prior distribution for the parameters. Explain your
reasoning.  

We will need to specify a prior distribution for the model such that it covers 
a suitable paramter space and can reasonably estimated in the posterior. I will 
use prior distributions with high variance in order to reflect that I have 
little information on how I believe these parameters are dispersed. Given that
our response variable has a variance of `7.2142` this should be large enough to
reach a good posterior estimate that converges.

$$
\sigma_{y} \sim Uniform(0,100)  \\
\sigma_{\alpha} \sim Uniform(0,100)  \\
\mu_{\alpha} \sim \mathcal{N}(0,100^2)  \\
$$


c. Estimate the model in a Bayesian way via Markov chain Monte Carlo.  

The model will be estimated using PyMC 3 which uses a Hamiltonian Markov Chain 
varient of the Markov Chain Monte Carlo Algorithm with a No U-Turn Sampling 
method. This sampling method is ideal for estimating continuous parameters where
information on the first order gradient can be used to inform step size. The
parameter start points are selected by estimating the maximum a posteriori (MAP)
via automatic differentiation variational inference (ADVI).  

d. Assess the convergence of your algorithm & whether it has run for enough
iterations.

The model appers to reach converegence as the samples seem to fluctuate at
random around a stable value relative to its intil start value. The images below
show all 5000 iterations as well as only the last 4500 iterations, where 500 
were removed for burn in, to demonstrate this convergence graphically for the 
parameters we wish to estimate. The iteration number is on the x-axis while the
sampling value for that iteration is on the y-axis.

# ![Meow](/home/nmarquez/Documents/Classes/statdemog/week5/alltrace.png "")

# ![](/home/nmarquez/Documents/Classes/statdemog/week5/burntrace.png "")

e. Summarize the posterior distribution you obtain in graphical and table form.  

The posterior distributions are plotted above for the full set of iterations as 
well as the set acounting for 500 iterations of burn in. The distributions by in
large have a centered mean and appear to be well estimated. Below is a table for
the three paramters most of interest and descriptive statitsics of the posterior
distribution including the mean, standard deviation, min, max, and the 2.5, 50, 
& 97.5 quantiles.

| stat | mu_alpha  | sigma_alpha | sigma_epsilon |
|----- | --------- | ----------- | ------------- |
|mean  | -0.312602 | 0.514002    | 1.098387      |
|std   |  0.073985 | 0.061300    | 0.019521      |
|min   | -0.569113 | 0.333738    | 1.036662      |
|2.5%  | -0.461543 | 0.404085    | 1.061083      |
|50%   | -0.311958 | 0.510248    | 1.097888      |
|97.5% | -0.166686 | 0.643797    | 1.136900      |
|max   | -0.023175 | 0.757744    | 1.166740      |
 
These values are in line with the estimates that are given by estimating a 
linear mixed effects model via restricted maximum likelihood (REML).

```
 mu_alpha     sigma_alpha   sigma_epsilon 
-0.3111669     0.5019289     1.0972490 
```
