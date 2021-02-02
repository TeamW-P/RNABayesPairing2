from matplotlib import pyplot as plt
import matplotlib
import seaborn as sns
import math
sns.set(style="white")

#x = [3,9,15]
#y1 = [11.29,39.63, 50.45]
#y3 = [13.79,43.55,60.01]
#y6 = [14.07, 45.68, 65.25]
#y9 = [15.63,46.55,68.35]
#y12 = [15.77,47.40,68.66]
#y15 = [17.03,48.44,71.28]

x =  [1, 3, 6, 9, 12, 15]
ya = [11.29, 13.79, 14.07, 15.63, 15.77, 17.03]
yb = [39.63, 43.55, 45.68, 46.55, 47.40, 48.44]
yc = [50.45, 60.01, 65.25, 68.35, 68.66, 71.28]
matplotlib.rcParams.update({"xtick.labelsize":15,"ytick.labelsize":15})


f, ax = plt.subplots(figsize=(15,7))
#f, ax = plt.subplots()

plt.subplot(121)
plt.title("BayesPairing2 execution time scaling",fontsize=19)
plt.ylabel("time(s)",fontsize=19)
plt.xlabel("Number of modules searched",fontsize=19)
plt.tight_layout()

#palette = sns.cubehelix_palette(10, start=2.6, rot=-0.2)[4:]
palette = sns.color_palette("muted")

plt.plot(x,ya, color=palette[0],linewidth=3, label="3 seqs", alpha=0.99)
plt.plot(x,yb, color=palette[1],linewidth=3, label="9 seqs", alpha =0.99)
plt.plot(x,yc, color=palette[2],linewidth=3, label="15 seqs", alpha = 0.99)
#plt.plot(x,y9, color=palette[3],linewidth=1.9, label="9 modules", alpha = 0.99)
#plt.plot(x,y12, color=palette[4],linewidth=1.9, label="12 modules", alpha = 0.99)
#plt.plot(x,y15, color=palette[5],linewidth=1.9, label="15 modules", alpha = 0.99)
plt.xticks([1, 3, 6, 9, 12, 15])

plt.legend(fontsize=19)
#plt.show()

#sns.set()
#here x is the number of modules
x = [1,3,6]
bp2_ss = [0.631,1.42,2.25]
bp1 = [21.413,49.653,135.937]
bp2 = [11.29,13.79,14.073]


palette = sns.color_palette("colorblind")
plt.subplot(122)
plt.plot(x,bp2_ss, color=palette[0],linewidth=3, label="BP2 with structure")
plt.plot(x,bp2, color=palette[1],linewidth=3, label="BP2")
plt.plot(x,bp1, color=palette[2],linewidth=3, label="BP1")
plt.xticks([1,3,6])
#ax.tick_params(axis='x', labelsize=15)
#ax.tick_params(axis='y', labelsize=15)

plt.title("BayesPairing execution time",fontsize=19)
plt.ylabel("time(s)",fontsize=19)
plt.xlabel("Number of modules searched in 15 seq. (len 200)",fontsize=19)
plt.yscale("log")
plt.legend(fontsize=19)
plt.tight_layout()
plt.savefig("speed.pdf",format="pdf")
plt.show()