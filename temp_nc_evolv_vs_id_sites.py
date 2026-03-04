
import matplotlib.pyplot as plt

fig, ax = plt.subplots()

id_sites_av = [6.6, 7.2, 5.7, 4.8, 4.11, 7.67, 8.62, 4.79, 7.04, 5.3]
evo_av = [22, 23, 24, 27, 30, 12, 13, 22, 14, 12]

for i, (x, y) in enumerate(zip(id_sites_av, evo_av), start=1):
    ax.scatter(x, y, color=f"C{i-1}", label=i)

ax.legend()
ax.set_xlabel("avg. num of identical bp pair sites in NC neighborhood")
ax.set_ylabel("avg. nc neigh. diversity")
plt.savefig("neigh_div_over_id_sites_averages.pdf", format="pdf")