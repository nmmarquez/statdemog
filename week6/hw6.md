# Homework Week 6  
### Neal Marquez  

1. Redo Homework 3, Question 3(b) once again. This time, solve it both
analytically (expressing the predictive distribution as a normal distribution,
conditional on the estimated parameters), and by simulation. Show your working
clearly.  

The Lee Carter model was fit two times, once using a Least Squares approach and another time using a SVD approach. For each of the models values for $a_x$,$b_x$, $\mu_k$, and $k_t$ were obtained and used to determine the mean forecast as well as $\sigma_k$, the error for the autoregressive $k$ terms and $\sigma_{\epsilon}$ the overall error in the model. Predictions were made for Mexican females aged 75-80 for the 2015-2020 time period. Uncertainty was generated by simulating data using 10000 simulations for the two error terms. The model was then forecasted to 2015 as follows

$log(m_{75,2015}) = a_{75} + b_{75} * k_{2015} + \epsilon_{75,2015}$  
$k_{2015} = k_{2010} + \mu_k + \epsilon_{2015}$  

The uncertainty for 2015 estimates of mortality can also be calculated analytically. The first source of error comes from the $k_t$ which has the following distribution.

$k_{2015} \sim \mathcal{N}(k_{2010} + \mu_k, \sigma_k)$

the $k_{2015}$ is then multiplied by the scalar value for b_{75} which alters the variance for the distribution such that the resulting variance is $b_{75}^2 \sigma_k^2$. This resulting distribution is normal and can be combined with the variance from the overall error term which is also normal and the variances can be added such that the final variance of the estimates are now $b_{75}^2 \sigma_k^2 + \sigma^2_{\epsilon}$. This total variance can be converted to a standard deviation and 1.96 times the standard deviation can then be added and subtracted from the mean values in order to get the 2.5% and the 97.5% quantiles.

The values generated by the two methods are shown below first for the least squares approach and then for the svd approach where ymin and ymax are the intervals that were simulated and ymin2 and ymax2 are the intervals that were created analytically. Mxt is the median mortality rate for women in mexico age 75 in 2015.

### least squares method
|     mxt|      ymin|      ymax|     ymin2|     ymax2|
|---------:|---------:|---------:|---------:|---------:|
| 0.0366498| 0.0321386| 0.0416743| 0.0322563| 0.0416417|

### svd method
|     mxt|      ymin|      ymax|     ymin2|     ymax2|
|---------:|---------:|---------:|---------:|---------:|
| 0.0371044| 0.0328375| 0.0419152| 0.0327884| 0.0419885|

Below are the plots for the least squares model and the svd model respectively.

![](/home/nmarquez/Documents/Classes/statdemog/week6/lsleecarter.jpg "")
![](/home/nmarquez/Documents/Classes/statdemog/week6/svdleecarter.jpg "")

2. Obtain the values of TFR for Honduras for 1950-2015 from the 2015 World Population Prospects.

a. Fit a version of the double logistic decline model by nonlinear least squares to
these data, assuming that the error variance remains constant over time (this
will just give one set of double logistic parameter values and an estimated error
variance). Note: Honduras is still in Phase II up to 2015.

The double logistic model uses two logistic equations in order to capture the change between time periods of measured total fertility rate as a country moves from phase two to three. The model is fit using the `nlminb` optimizer in `R` in order to fit the 5 parameters of the model as well as the error variance term. Because Honduras is still in Phase 2 the model was unable to converge on a final answer as the parameters $\Delta_{1,2,3,4}$ dictate the time at which phase 3 starts and without data the model fails to choose a single set of parameters to converge on. Nevertheless the fit is appropriate for the data.

![](/home/nmarquez/Documents/Classes/statdemog/week6/tfrdiff.jpg "")