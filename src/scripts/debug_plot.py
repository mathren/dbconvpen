import numpy as np
import paths
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import sys
from MESAreader import read_MESA_header, getSrcCol


def plot_cell_T_rho(pfile, fname=None):
    src, col = getSrcCol(pfile)
    # for c in col: print(c)
    z = src[:, col.index("zone")]
    logT = src[:, col.index("logT")]
    logRho = src[:, col.index("logRho")]
    logr = src[:, col.index("logR")]
    he4 = src[:, col.index("he4")]
    eps_he = src[:, col.index("tri_alfa")]
    fig = plt.figure(figsize=(10,10))
    gs = gridspec.GridSpec(100, 100)
    ax = fig.add_subplot(gs[:, :])
    ax.scatter(z, logT, label="logT")
    ax.scatter(z, logRho, label="logRho")
    ax.plot(z, logr, lw=1, ls='--', c="#808080", zorder=0)
    bx = ax.twinx()
    bx.scatter(z, eps_he, c='k')
    ax.legend()
    ax.axvline(2652, 0,1, ls='--', lw=3)
    if fname:
        plt.savefig(paths.figures/str(fname))
    else:
        plt.show()
    plt.close()

if __name__ == "__main__":
    plot_cell_T_rho(sys.argv[1], fname="test.png")
