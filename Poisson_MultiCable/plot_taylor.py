from IPython import embed
import matplotlib.pyplot as plt
from os import system

def plot_rates_1(lines, i, color,marker):
    epsilons = [float(k) for k in lines[i+1].split()]
    res_0 = [float(k) for k in lines[i+2].split()]
    res_1 = [float(k) for k in lines[i+3].split()]
    rate_0 = [float(k) for k in lines[i+4].split()]
    rate_1 = [float(k) for k in lines[i+5].split()]
    # plt.semilogx(epsilons[1:], rate_0,color=color, linestyle="dashed")
    plt.semilogx(epsilons[1:], rate_1,color=color, linestyle="solid",
                 marker=marker,label="# Cells: %s" %lines[i].split()[-1])

def plot_rates_0(lines, i, color):
    epsilons = [float(k) for k in lines[i+1].split()]
    res_0 = [float(k) for k in lines[i+2].split()]
    res_1 = [float(k) for k in lines[i+3].split()]
    rate_0 = [float(k) for k in lines[i+4].split()]
    rate_1 = [float(k) for k in lines[i+5].split()]
    plt.semilogx(epsilons[1:], rate_0,color=color, linestyle="solid",
                 marker="o",label="# Cells: %s" %lines[i].split()[-1])




filename = "output/taylor_data.txt"

i = 0
for res in [2,3,4,5,7]:
    if i == 0:
        suffix = ">"
    else:
        suffix = ">>"
    system("python3 refine_mesh.py {0:1d}".format(res))
    system("python3 taylor_test.py" + suffix + filename)
    i+=1

taylor_data = open(filename, "r")
lines = taylor_data.readlines()
indicies = [i for i, s in enumerate(lines) if '#Cells' in s]
colors = ["r","b", "g","k","c","m"]
markers = ["o", "s", "h", "D","*","X"]
counter=0
plt.figure()
for i in indicies:
    plot_rates_1(lines,i,colors[counter],markers[counter])
    counter+=1
ax = plt.gca()
ax.legend(fontsize="x-large")
plt.grid(True)
plt.ylim((1.6,2.1))
ax.tick_params(labelsize="x-large")
plt.xlabel(r"Perturbation length ($\epsilon$)",fontsize="x-large")
plt.ylabel(r"Convergence rate",fontsize="x-large")
plt.tight_layout()
plt.savefig("output/taylor_test.png",dpi=250)
import os
os.system("convert output/taylor_test.png -trim output/taylor_test.png")
