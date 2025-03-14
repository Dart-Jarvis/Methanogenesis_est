import math
import numpy as np

d13C1,dD1,D13CD1,DDD1=-74.604, 793.236, -5.266, -15.509
x13C1=0.0112372*(d13C1/1000+1)
x13C1=x13C1/(1+x13C1)
xD1=0.00015576*(dD1/1000+1)
xD1=xD1/(1+xD1)
sto_13CD=4*x13C1*xD1*(1-xD1)**3
sto_DD=6*(1-x13C1)*(xD1**2)*((1-xD1)**2)
sto_CH4=(1-x13C1)*(1-xD1)**4
r13CD_sto=sto_13CD/sto_CH4
rDD_sto=sto_DD/sto_CH4
r13CD=(D13CD1/1000+1)*r13CD_sto
rDD=(DDD1/1000+1)*rDD_sto
d13CD=(r13CD/6.34600e-06-1)*1000
dDD=(rDD/1.28866e-07-1)*1000
print(d13CD,dDD)

print(np.sqrt((0.028)**2+(0.0027)**2+(0.022)**2+(0.019)**2+(0.037)**2))