from scipy.optimize import minimize
import numpy as np
import matplotlib.pylab as plt

def getData(freqs, phases, time):
    totalData = np.zeros(200000)
    for signalInc in range(0, freqs.size):
        y = np.sin(2 * np.pi * freqs[signalInc] * time + phases[signalInc])
        totalData = np.add(totalData, y)
    return totalData, max(totalData)

def maxSignal(phi):
    time = np.linspace(0, 100 / 100000000., 200000)
    signalFreqs = np.array([65, 75, 85, 95]) * 1000000
    data, maxData = getData(signalFreqs, [0, phi[0], phi[1], phi[2]], time)
    return maxData
#result = minimize(maxSignal, np.array([1, 0, 0]))
time = np.linspace(0, 100 / 100000000., 200000)
signalFreqs = np.array([65, 75, 85, 95]) * 1000000
bestPhases = [0, 3.80277662, -0.47467704, -0.27651992]
bestPhases2 = [0, 3.80277662, 5.808508267, 6.006665]

# 4 at zeros: 3.72mW
# 4 at found best phases: 5.196

data, maxD = getData(signalFreqs, bestPhases2, time)
plt.plot(time, data)
#plt.show()
signalFreqs = np.array([65, 75, 85, 95]) * 1000000
time = np.linspace(0, 100 / 100000000., 200000)
data = [None]*7
maxSignal = [None]*7

data[0], maxSignal[0] = getData(signalFreqs, [0, 0, 0, 0], time)
data[1], maxSignal[1] = getData(signalFreqs, [1, 0, 0, 0], time)
data[2], maxSignal[2] = getData(signalFreqs, [2, 0, 0, 0], time)
data[3], maxSignal[3] = getData(signalFreqs, [3, 0, 0, 0], time)
data[4], maxSignal[4] = getData(signalFreqs, [0, 1, 2, 3], time)
data[5], maxSignal[5] = getData(signalFreqs, [3.1415, 0, 0, 0], time)
data[6], maxSignal[6] = getData(signalFreqs, [0, 1, 1, 1], time)
#data[7], maxSignal[7] = getData(signalFreqs, [np.pi/2, np.pi/2, np.pi/2, np.pi/2], time)

#print(maxSignal[7])
#print(maxSignal[0])

plt.plot(time, data[0])
plt.plot(time, data[1])
plt.show()

powers = [3.669, 4.543, 5.31, 5.29, 4.984, 5.28, 4.307]
print(len(maxSignal), len(powers))
plt.plot(maxSignal, powers, linestyle="none", marker='o')
plt.xlabel('(Relative) Maximum Voltage in outputted signal')
plt.ylabel('Power observed in rail (mW)')
plt.show()

fig1 = plt.figure(1)
"""
plt.plot(t, totalData, color='b', label="2,0,0,0")

phases = np.array([0, 0, 0, 0])
totalData = np.zeros(1000)
for signalInc in range(0, freqs.size):
    y = np.sin(2 * np.pi * freqs[signalInc] * t + phases[signalInc])
    totalData = np.add(totalData, y)
plt.plot(t, totalData, color='r', label='0,0,0,0')

phases = np.array([1, 0, 0, 0])
totalData = np.zeros(1000)
for signalInc in range(0, freqs.size):
    y = np.sin(2 * np.pi * freqs[signalInc] * t + phases[signalInc])
    totalData = np.add(totalData, y)
plt.plot(t, totalData, color='g', label='1,0,0,0')

plt.show()
"""