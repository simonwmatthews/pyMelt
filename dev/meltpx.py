# %%
import pyMelt as m
import numpy as np
import matplotlib.pyplot as plt

# %%
meltpx = m.lithologies.lambart.MeltPX(
            alk=3.09, Mgn=0.6, Al2O3=9.66, CaO=13.09, 
            high_F_prod=0.003
            )

# %%
P = np.linspace(0,6,100)
Tsol = np.zeros(np.shape(P))
Tliq = np.zeros(np.shape(P))

for i in range(len(P)):
    Tsol[i] = meltpx.TSolidus(P[i])
    Tliq[i] = meltpx.TLiquidus(P[i])

fig, ax = plt.subplots()
ax.plot(Tsol, P)
ax.plot(Tliq, P)
ax.invert_yaxis()
plt.show()


# %%
P = 4.0
T = np.linspace(1100, 1800, 500)
F = np.zeros(np.shape(T))

for i in range(len(T)):
    F[i] = meltpx.F(P, T[i])

fig, ax = plt.subplots()
ax.plot(T, F)
ax.axvline(meltpx.TSolidus(P), c='k')
ax.axvline(meltpx.TLiquidus(P), c='k')
plt.show()


# %%
lz = m.lithologies.katz.lherzolite()
mantle = m.mantle([lz, meltpx], [10.0, 1.0], ['lz', 'px'])
col = mantle.adiabaticMelt(1350.0, Pend=0.1, Pstart=6.0)

# %%
col.plot()


# %%
meltpx.F(3.0, 1400.0)
# %%
A = meltpx._A(3.0)
B = meltpx._B(3.0)
# %%
(- B + np.sqrt(B**2 - 20 * A))/10.0

# %%

(- B - np.sqrt(B**2 - 20 * A))/10.0
# %%
