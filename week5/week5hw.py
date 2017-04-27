import pymc3 as pm
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

DF = pd.read_csv("/home/nmarquez/Documents/Classes/statdemog/week5/school.csv")
DF = DF.query("year == .5")[["schoolid", "math"]]
DF["schoolfactor"] = pd.factorize(DF["schoolid"])[0]
J = len(np.unique(DF.schoolfactor))
N = DF.shape[0]


with pm.Model() as BYHM:
    mu_alpha = pm.Normal("mu_alpha", mu=0, sd=100.)
    sigma_alpha = pm.Uniform("sigma_alpha", lower=0., upper=100.)
    sigma_epsilon = pm.Uniform("sigma_epsilon", lower=0., upper=100.)

    alpha_j = pm.Normal("alpha_j", mu_alpha, sd=sigma_alpha, shape=J)

    mathhat = alpha_j[DF.schoolfactor.values]

    ylike = pm.Normal("ylike", mu=mathhat, sd=sigma_epsilon,
                      observed=DF.math.values)

with BYHM:
    step = pm.NUTS()
    trace = pm.sample(5000, step=step)

pm.traceplot(trace[0:])
plt.suptitle("All 5000 traces", fontsize=24, y=.98)
plt.subplots_adjust(top=0.92)
plt.savefig('/home/nmarquez/Documents/Classes/statdemog/week5/alltrace.png')
with BYHM:
    # Use ADVI for initialization
    mu, sds, elbo = pm.variational.advi(n=100000)
    step = pm.NUTS(scaling=BYHM.dict_to_array(sds)**2, is_cov=True)
    ht = pm.sample(5000, step, start=mu)

burn = 500
pm.traceplot(ht[burn:])
plt.suptitle("Trace with 500 burn-in", fontsize=24, y=.98)
plt.subplots_adjust(top=0.92)
plt.savefig('/home/nmarquez/Documents/Classes/statdemog/week5/burntrace.png')

params = ["mu_alpha", "sigma_alpha", "sigma_epsilon"]
perc = [.025, .5, .975]
desc = pd.DataFrame({k: pd.Series(ht[k][burn:]).describe(percentiles=perc)[1:]
                     for k in params})
desc
