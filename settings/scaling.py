import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('results/scaling.txt',delimiter=',')

cores = data[:,0]
time = data[:,1]
error = np.sqrt(data[:,2])/1000 #arreglar esto

time_1_core = time[0]

Speedup = time_1_core/time
Ideal_s = np.linspace(cores[0],cores[3],len(cores))
Efficiency = Speedup/Ideal_s
Ideal_e = np.linspace(1,1,len(cores))

fig, (ax0, ax1) = plt.subplots(nrows=1, ncols=2, figsize=(12, 6))
ax0.set_title('Speedup')
ax0.errorbar(cores,Speedup,yerr=error,label='Speedup')
ax0.plot(cores,Ideal_s,color='red',linestyle='-.',label='Ideal speedup')
ax0.set_xlabel('# Cores')
ax0.set_xticks(cores)
ax0.legend()

ax1.set_title('Parallel efficiency')
ax1.errorbar(cores,Efficiency,yerr=error,label='Parallel efficiency')
ax1.plot(cores,Ideal_e,color='red',linestyle='-.',label='Ideal parallel efficiency')
ax1.set_xlabel('# Cores')
ax1.set_xticks(cores)
ax1.legend()

plt.show()