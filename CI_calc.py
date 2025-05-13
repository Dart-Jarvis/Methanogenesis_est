import numpy as np
from scipy import stats

# Your data
data = [
9.249,7.561,7.039,8.657,9.358,9.854,7.161,9.987,8.768,9.330,
7.940,8.637,9.246,7.648,9.986,7.019,9.493,9.043,9.989,7.931,
8.686,8.315,8.739,9.235,9.227,9.513,10.307,9.276,10.012,10.988,
9.994
]

# Sample size, mean, and standard error
n = len(data)
mean = np.mean(data)
sem = stats.sem(data)  # Standard error of the mean

# 95% CI (two-tailed t-distribution)
confidence = 0.95
h = sem * stats.t.ppf((1 + confidence) / 2., n-1)

# Confidence interval
ci_lower = mean - h
ci_upper = mean + h

print(f"Mean = {mean:.3f}")
print(f"95% CI = ({ci_lower:.3f}, {ci_upper:.3f})")
print(h)
