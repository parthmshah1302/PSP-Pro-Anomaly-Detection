import matplotlib.pyplot as plt
import numpy as np
delta=0.5
delta2=1
delta3=2
delta4=4
lamda = np.linspace(0, 30, 1000)
# kavya demoo 


power=1/6
driftv=pow((delta*lamda)/2,power)
driftv2=pow((delta2*lamda)/2,power)
driftv3=pow((delta3*lamda)/2,power)
driftv4=pow((delta4*lamda)/2,power)

fig, ax = plt.subplots()
plt.xlabel("Drift velociy (v)")
plt.ylabel("Lambda (λ)")
ax.plot(driftv, lamda)
ax.plot(driftv2, lamda)
ax.plot(driftv3, lamda)
ax.plot(driftv4, lamda)
plt.show()