# Homework Week 1  

1. Suppose the force of mortality at age x for a cohort is $\mu(x) = x^2$.
Denote by X the age at death of a person randomly chosen from this cohort (a
random variable).  

    a. Find the Cumulitive hazard function.  

The cumulitive hazard function is $\Lambda(x) = \int_{0}^{x} \mu(u)du$. Given 
that our function is $\mu(x) = x^2$ we can take the integral and find that 
$\Lambda(x) =  \frac{x^3}{3} - \frac{0}{3} = \frac{x^3}{3}$.  

    b. Find the survival function.  

The survival function is defined as $S(x) = exp(- \Lambda (x))$ which makes 
it a function of the cumulitive hazard function. For our hazard function the 
survival function is $S(x) = exp(\frac{x^3}{-3})$.  

    c. Find and graph the probability density function of X.

The pdf can be found be taking the additive inverse of the derivative of the 
Survival function such the $f(x) = - \frac{d}{dx}S(x)$ which gives us 
$f(x) = exp(\frac{x^3}{-3}) ~ x^2$.  

![](/home/nmarquez/Documents/Classes/statdemog/week1/pdfX.jpg "")  

    d. What is the name of this distribution?  

Chi squared distribution.  

    e. Find the life expectancy at birth of a member of this cohort.  
    
Life expectancy at birth is given as $e_0 = E(X) = \int_{0}^{\infty} S(u)du$. 
this simplifies to 
$e_0 = \frac{\Gamma(\frac{1}{3})}{3^{(\frac{2}{3})}} \approx 1.2879$  

    f. Find the life expectancy at age 2 of a member of this cohort.  
    
Life Expectancy at any age can be calculated using 
$e_x = \frac{\int_{x}^{\infty} S(u)du}{S(x)}$. For the case of life expectancy 
at age 2 we get $e_2 = \frac{\int_{2}^{\infty} S(u)du}{S(2)}$ which can be
reduced to $e_2 = \frac{\Gamma(8/3, 1/3)}{3^{2/3} \times S(2)}$ and finally 
$e_2 \approx .2087$.

    g. Find 2q2 for this cohort.  
    
Given that $_nq_x = S(x) - S(x+n)$ We can directly compute this as 
$exp(2^3/-3) - exp(4^3/-3) \approx .0695$.  

2. Estimating life tables.  

    a. From the UN’s 2015 World Population Prospects, extract the estimates of 
    the age-specific mortality rates nmx for females in Mexico in 2005-2010. 
    Plot them against age on the raw and logarithmic scales, and comment on any 
    unusual features.  
    
![](/home/nmarquez/Documents/Classes/statdemog/week1/gomp.jpg "")  
![](/home/nmarquez/Documents/Classes/statdemog/week1/loggomp.jpg "")

    b. Using these, derive a life table for this population and time period. 
    Show your work.

Work is shown in accompanying code.  

| age|   nMx|   nqx|   npx|        lx|       ndx|       nLx|         Tx|     ex|
|---:|-----:|-----:|-----:|---------:|---------:|---------:|----------:|------:|
|   0| 0.018| 0.018| 0.982| 100000.00|  1765.137|  98764.40| 7827612.83| 78.276|
|   1| 0.001| 0.004| 0.996|  98234.86|   423.821| 392091.81| 7728848.43| 78.677|
|   5| 0.000| 0.002| 0.998|  97811.04|   161.092| 488652.48| 7336756.62| 75.009|
|  10| 0.000| 0.001| 0.999|  97649.95|   126.007| 487934.73| 6848104.14| 70.129|
|  15| 0.000| 0.002| 0.998|  97523.94|   166.918| 487202.42| 6360169.41| 65.216|
|  20| 0.000| 0.002| 0.998|  97357.02|   228.551| 486213.74| 5872966.99| 60.324|
|  25| 0.001| 0.003| 0.997|  97128.47|   293.964| 484907.45| 5386753.25| 55.460|
|  30| 0.001| 0.004| 0.996|  96834.51|   389.128| 483199.72| 4901845.80| 50.621|
|  35| 0.001| 0.006| 0.994|  96445.38|   549.950| 480852.03| 4418646.07| 45.815|
|  40| 0.002| 0.009| 0.991|  95895.43|   820.308| 477426.38| 3937794.04| 41.063|
|  45| 0.003| 0.013| 0.987|  95075.12|  1252.741| 472243.76| 3460367.66| 36.396|
|  50| 0.004| 0.021| 0.979|  93822.38|  1929.156| 464289.02| 2988123.90| 31.849|
|  55| 0.007| 0.032| 0.968|  91893.23|  2962.095| 452060.89| 2523834.88| 27.465|
|  60| 0.010| 0.050| 0.950|  88931.13|  4482.106| 433450.39| 2071774.00| 23.296|
|  65| 0.016| 0.078| 0.922|  84449.02|  6619.208| 405697.10| 1638323.61| 19.400|
|  70| 0.026| 0.121| 0.879|  77829.82|  9394.206| 365663.57| 1232626.51| 15.837|
|  75| 0.040| 0.183| 0.817|  68435.61| 12522.000| 310873.05|  866962.94| 12.668|
|  80| 0.063| 0.271| 0.729|  55913.61| 15142.485| 241711.84|  556089.89|  9.946|
|  85| 0.098| 0.387| 0.613|  40771.13| 15782.650| 164399.00|  314378.05|  7.711|
|  90| 0.150| 0.527| 0.473|  24988.47| 13177.349|  91999.00|  149979.05|  6.002|
|  95| 0.227| 0.679| 0.321|  11811.13|  8017.496|  39011.89|   57980.04|  4.909|
| 100| 0.400| 1.000| 0.000|   3793.63|  3793.630|  18968.15|   18968.15|  5.000|

    c. Find the life expectancy at birth and at age 10 for this population.  
    
Life expectancy at birth is $78.276$ while at age ten it is $70.129$.