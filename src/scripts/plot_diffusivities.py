import numpy as np
import paths
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import sys
from MESAreader import read_MESA_header, GetSrcCol



def get_age_from_profile(pfile):
    src, col = read_MESA_header(pfile)
    return float(src[col.index("star_age")])


def plot_Dmix(pfile, fname=None):
    fig = plt.figure(figsize=(10,10))
    gs = gridspec.GridSpec(100, 100)
    ax = fig.add_subplot(gs[:, :])

    age = get_age_from_profile(pfile)
    title = f"{age:.2f}"+r"\,\mathrm{Myr}$"

    src, col = getSrcCol(pfile)
    m = src[:, col.index("mass")]
    log_D_mix = src[:, col.index('log_D_mix')]
    log_D_mix_non_rotation = src[:, col.index('log_D_mix_non_rotation')]
    log_D_conv = src[:, col.index('log_D_conv')]
    log_D_semi = src[:, col.index('log_D_semi')]
    log_D_thrm = src[:, col.index('log_D_thrm')]

    ax.plot(m, log_D_mix, lw=10, label="$\mathrm{Total}$", c="c")
    ax.plot(m, log_D_mix_non_rotation, label="$\mathrm{Excluding\ rotation}$", c='b')
    ax.plot(m, log_D_conv, label="$\mathrm{Convection}$", c="r")
    ax.plot(m, log_D_semi, label="$\mathrm{Semiconvection}$", c="m")
    ax.plot(m, log_D_thrm, label="$\mathrm{Thermohaline}$", c="#CC99BB")
    # ax.legend()
    # ax.set_title(title, fontsize=30)
    ax.text(0.5,0.92, title, fontsize=30, va="center", ha="center", transform=ax.transAxes)
    if fname:
        plt.savefig(paths.figures / fname)
    else:
        plt.show()
    plt.close()


if __name__ == "__main__":
    try:
        plot_Dmix(sys.argv[1], sys.argv[2])
    except:
        plot_Dmix(sys.argv[1], fname=None)
