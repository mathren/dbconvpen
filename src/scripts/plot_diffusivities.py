import numpy as np
import paths
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import sys


def getSrcCol(f, clean=True, convert=True, bin_fname=None):
    """
    Read MESA output (history or profiles) and optionally
    save a copy in binary format for faster access later.

    Parameters:
    ----------
    clean:   `bool` if True use log_scrubber to remove retry steps before converting to binary
              and parsing the output.
    convert: `bool`, if True the file f is saved to
              binary format for faster reading in the future.
    bin_fname:   `str` name of the binary file if needed

    Returns:
    --------
    src: `np.array`, shape (number of timesteps,
          number of columns), the data
    col: `list`, column names from the header
    """
    # print(f)
    # TODO: maybe one day I'll update this to be a pandas dataframe
    # should work both for history and profiles
    # read header
    with open(f, "r") as P:
        for i, line in enumerate(P):
            if i == 5:
                col = line.split()
                break
    # check if binary exists
    if bin_fname == None:
        bin_fname = str(str(f)[:-4]) + ".npy"
    if os.path.isfile(bin_fname):
        # read the column and binary
        src = reader(f, len(col), 6, bin_fname)
    else:  # binary file does not exist
        print("... Binary file does not yet exist")
        if ("history.data" in str(f)) and clean:
            scrub(f)
        if convert:
            src = reader(f, len(col), 6, bin_fname)
        else:
            src = np.genfromtxt(f, skip_header=6)
    return src, col


def read_MESA_header(fname):
    """reader for the header of MESA profile or history files.
    Parameters:
    ----------
    fname : `string`, path to the file to open
    Returns:
    --------
    src : `list` values, some cannot be converted to float
    col : `list`, columns names
    """
    with open(fname, "r") as f:
        for i, line in enumerate(f):
            if i == 1:
                col = line.split()
            if i == 2:
                src = line.split()
                break
    return src, col


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
