import numpy as np
import paths
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import sys
from MESAreader import read_MESA_header, GetSrcCol


def plot_cell_T_rho(pfile, fname=None):
    src, col = GetSrcCol(pfile)
    z = src[:, col.index("zone")]
    logT = src[:, col.index("logT")]
    logRho = src[:, col.index("logRho")]
    fig = plt.figure(figsize=(10,10))
    gs = gridspec.GridSpec(100, 100)
    ax = fig.add_subplot(gs[:, :])
    ax.plot(z, logT, label="logT")
    ax.plot(z, logRho, label="logRho")
    ax.legend()
    if fname:
        plt.savefig(paths.figures/str(fname))
    else:
        plt.show()
    plt.close()

if __name__ == "__main__":
    plot_cell_T_rho(sys.argv[1])
