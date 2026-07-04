import matplotlib.pyplot as plt
import numpy as np


test = 0 # Should be same as in T_E.c

potentials_file = f"./T_E_test_{test}/T_E_test_{test}_potentials.txt"
if(test == 0):
    potentials_file = "./poisson_sio2.txt"

x1 = []
y1 = []
with open(potentials_file, "r") as file:
    for line in file:
        if not line.strip():
            continue
        columns = line.split()
        x1.append(float(columns[0]))
        # Convert Joules to eV
        y1.append(float(columns[2]) / 1.60218e-19)

N = int(max(x1)) 


data2 = np.loadtxt('transmission_data.txt')


fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

ax1.step(x1, y1, where='mid', color='darkblue', linewidth=2, label='Potentials')
ax1.bar(x1, y1, width=1.0, alpha=0.3, color='royalblue', edgecolor='darkblue', align='center', label='Slices')

ax1.axvline(x=0.5, color='black', linestyle='--', linewidth=1.5)

ax1.axvline(x=N - 0.5, color='black', linestyle='--', linewidth=1.5)

y_pos = max(y1) * 0.9 if y1 else 0.9  
ax1.text(0.0, y_pos/2, 'Channel', color='black', fontsize=10, fontweight='bold', ha='center')
ax1.text(N / 2, y_pos/2, 'Barrier/Oxide', color='black', fontsize=10, fontweight='bold', ha='center')
ax1.text(N, y_pos/2, 'Gate', color='black', fontsize=10, fontweight='bold', ha='center')

ax1.set_xlabel('Slice Number i')
ax1.set_ylabel('Potential $V_i$ (in eV)')
ax1.set_title('Potential Profile (Slices)')
ax1.grid(True, linestyle='--', alpha=0.5)
ax1.legend()

ax2.plot(data2[:, 0], data2[:, 1], color='crimson', linewidth=2, label='Transmission Probability T(E)')
ax2.set_xlabel('Electron Energy (eV)')
ax2.set_ylabel('Transmission Probability')
ax2.set_title('Transmission Probability vs Electron Energy')
ax2.grid(True, linestyle='--', alpha=0.5)
ax2.legend()

plt.tight_layout()
plt.show()
