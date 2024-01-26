import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from plot_tools import *
from scipy.interpolate import interp1d
from ase.io import read
'''
def add_lines_to_core(core_ax, s, E,
                      orientation="bottom",
                      color="C0",
                      x_shift = 0.05):

    fig = core_ax.get_figure()

    # get bbox in in figure coordinates
    bbox = core_ax.get_position()
    # transform to display coordinates
    bbox = bbox.transformed(fig.transFigure)
    # transfrom to Data coordinates for ploting
    bbox = bbox.transformed(ax.transData.inverted())

    x0 = bbox.xmin
    y0 = bbox.ymin

    x1 = bbox.xmax
    y1 = bbox.ymax

    x0 += x_shift
    x1 -= x_shift

    if orientation == "left":

        plot_x1 = x0
        plot_y1 = y1

        plot_x2 = x0
        plot_y2 = y0

    elif orientation == "right":

        plot_x1 = x1
        plot_y1 = y1

        plot_x2 = x0
        plot_y2 = y1

    elif orientation == "top":

        plot_x1 = x0
        plot_y1 = y1

        plot_x2 = x1
        plot_y2 = y1

    elif orientation == "bottom":

        plot_x1 = x0
        plot_y1 = y0

        plot_x2 = x1
        plot_y2 = y0


    ax.plot((plot_x1, s), (plot_y1, E),
            lw=2.0, color=color, linestyle="dashed", clip_on=False)

    ax.plot((plot_x2, s), (plot_y2, E),
            lw=2.0, color=color, linestyle="dashed", clip_on=False)
'''
def add_core_plot(fig, i, rect,
                  pure_glide_averaged_images,
                  QM_MM_bulk,
                  alat=3.17,
                  core_positions=None,
                  legend=False,
                  core_size=75, spine_color="black",
                  spine_lw=1.0, core_lw=1.0,
                  legend_bbox=(0.375, 1.0, 0.5, 0.5),
                  xlim=(-1.9, 4.3)):

    if core_positions is not None:
        core_x, core_y = core_positions.T

        f_core = interp1d(core_x, core_y, kind='cubic')
        xnew_core = np.linspace(core_x.min(), core_x.max(), num=41, endpoint=True)

    core1 = fig.add_axes(rect)

    core1.set_xticks(())
    core1.set_yticks(())

    image = pure_glide_averaged_images[i]
    sd.plot_vitek(image, QM_MM_bulk, alat=alat, plot_axes=core1, xyscale=6)
    core1.set_xlim(xlim)
    core1.set_ylim((-1.3, 2.0))

    #core1.scatter(core_x, core_y, marker="^", facecolor='None', s=10, edgecolor="C3")

    #core1.scatter(core_x[i], core_y[i], marker = "^", s=core_size, c="C3", zorder=10)


    #core = core1.scatter([], [], s=100, marker="^",
                   #color="C3", label="Core position")
    #W = core1.scatter([], [], marker="o", c="white", s=100, edgecolors="k", label="W atoms")
    #core1.scatter(core_x, core_y, marker="^", facecolor='None', s=15, edgecolor="C3")

    #core1.plot(xnew_core, f_core(xnew_core),
    #           linestyle="dashed", color="C3", lw=core_lw,
    #           label="Core trajectory")

    core1.set_aspect("equal")

    if legend:
        core1.legend(bbox_to_anchor=legend_bbox,
                     fontsize=12, loc="lower center")

    for axis in ['top','bottom','left','right']:
        core1.spines[axis].set_linewidth(spine_lw)
        core1.spines[axis].set_color(spine_color)


    # core1.set_axis_off()

    return core1

if __name__ == '__main__':

    pure_glide_averaged_images = read("pure_glide_averaged_images.xyz", index=":")
    QM_MM_core_positions = np.loadtxt("QMMM_core_positions.txt")
    pure_glide_dataframe = pd.read_csv("pure_glide_averaged_energies.csv")

    QM_MM_bulk = read("QMMM_bulk.xyz")

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111)

    s_W, E_W, QMMM = plot_averaged_barrier(pure_glide_averaged_images, pure_glide_dataframe,
                                           label="QM/MM: VASP/EAM3", ms=7.5,
                                           ax=ax, s_norm=True, lw=2)

    MM, = ax.plot([], [], label="MM:", alpha=0.0)

    s, E = np.loadtxt("eam3_s_E.txt")
    Sfit, Efit = np.loadtxt("eam3_fit.txt")

    norm = s.max()
    s /= norm
    Sfit /= norm

    ax.plot(Sfit, Efit, 'C5-', lw=2)
    Ef = max(E) - E[0]
    eam3, = ax.plot(s, E, "sC5", label="MM: EAM3")

    lines = np.loadtxt("eam3_lines.txt")
    proper_lines = []

    #    for line in lines.reshape((9, 6)):
    #        proper_lines.append(line.reshape((2,3)))

    #    for x, y in proper_lines:
    #        x/=norm
    #        ax.plot(x, y, '-C5')

    GAP = np.genfromtxt("GAP_barrier.csv", delimiter=',')
    x = GAP.T[0]
    y = GAP.T[1]
    #    ax.scatter(x, y, z, facecolors='none', edgecolors='C2')
    f2 = interp1d(x, y, kind='cubic')
    xnew = np.linspace(x.min(), x.max(), num=101, endpoint=True)
    #GAP, = ax.plot(xnew, f2(xnew), color="C2", label="GAP: SBC")

    # Create a legend for the first line.

    proxy_eam3, = ax.plot([-10], [-10], "C5s-", label="Cluster MM: EAM3")
    proxy_QMMM = ax.errorbar([-10], [-10], fmt="C0o-", yerr=[0], label="Cluster QM/MM: VASP/EAM3")

    first_legend = ax.legend(handles=[proxy_QMMM, proxy_eam3],
                             loc="upper right")
                             #bbox_to_anchor=(0.6, 0.5, 0.5, 0.5))

    # Add the legend manually to the current Axes.
    plt.gca().add_artist(first_legend)

    DFT, = ax.plot([], [], label="Quadrupole DFT:", alpha=0.0)

    mvg_siesta = np.genfromtxt("MVG_siesta.csv", delimiter=',')
    x = mvg_siesta.T[0]
    y = mvg_siesta.T[1]/1.0e3
    siesta = ax.scatter(x, y, label="MVG SIESTA GGA", facecolors='none', edgecolors='C4')
    f2 = interp1d(x, y, kind='cubic')
    xnew = np.linspace(x.min(), x.max(), num=41, endpoint=True)
    ax.plot(xnew, f2(xnew), linestyle="dashed", color="C4", lw=2)

    DVC_pwscf = np.genfromtxt("DVC_PWSCF.csv", delimiter=',')
    x = DVC_pwscf.T[0]
    y = DVC_pwscf.T[1]/1.0e3
    siesta = ax.scatter(x, y, label="DVC PWSCF GGA", facecolors='none', edgecolors='C2')
    f2 = interp1d(x, y, kind='cubic')
    xnew = np.linspace(x.min(), x.max(), num=41, endpoint=True)
    ax.plot(xnew, f2(xnew), linestyle="dashed", color="C2", lw=2)

    """
    mvg_pwscf = np.genfromtxt("MVG_PWSGCF_GGA.csv", delimiter=',')
    x = mvg_pwscf.T[0]
    y = mvg_pwscf.T[1]/1.0e3
    PWSCF = ax.scatter(x, y, label="MVG PWSCF GGA", facecolors='none', edgecolors='C1')
    f2 = interp1d(x, y, kind='cubic')
    xnew = np.linspace(x.min(), x.max(), num=41, endpoint=True)
    ax.plot(xnew, f2(xnew), linestyle="dashed", color="C1", lw=2)


    GAP_DFT = np.genfromtxt("GAP_DFT.csv", delimiter=',')
    x = GAP_DFT.T[0]
    y = GAP_DFT.T[1]
    CASTEP = ax.scatter(x, y, label="SBC CASTEP", facecolors='none', edgecolors='C2')
    f2 = interp1d(x, y, kind='cubic')
    xnew = np.linspace(x.min(), x.max(), num=41, endpoint=True)
    ax.plot(xnew, f2(xnew), linestyle="dashed", color="C2")
    """
    s, E = np.loadtxt("s_E_VASP.txt")
    Sfit, Efit = np.loadtxt("fit_VASP.txt")

    norm = s.max()
    s /= norm
    Sfit /= norm

    ax.plot(Sfit, Efit/2., 'C3-', linestyle="dashed", lw=2)
    Ef = max(E) - E[0]
    VASP = ax.scatter(s, E/2., edgecolors='C3', facecolors="None", label="This work VASP")

    kcd_vasp = np.genfromtxt("KCD_VASP.csv", delimiter=',')
    x = kcd_vasp.T[0]
    y = kcd_vasp.T[1]
    # we use our reaction coordinate since they do not report it
    KCD_VASP = ax.scatter(s, y, label="MVG PWSCF GGA", facecolors='none', edgecolors='C1')
    f2 = interp1d(s, y, kind='cubic')
    xnew = np.linspace(x.min(), x.max(), num=41, endpoint=True)
    ax.plot(xnew, f2(xnew), linestyle="dashed", color="C1", lw=2)

    ax.set_xlim(-0.15, 1.15)
    ax.set_ylim(-0.015, 0.12)

    ax.grid(True, linestyle="dashed")


    ax.set_xlabel('Normalized NEB coordinate')
    ax.set_ylabel('Energy (eV)')

    proxy_siesta, = ax.plot([], [], "C4o--", markerfacecolor='none', label="MVG SIESTA GGA")
    proxy_PWSCF, = ax.plot([], [], "C1o--", markerfacecolor='none', label="MVG PWSCF GGA")
    proxy_KCD_VASP, = ax.plot([], [], "C1o--", markerfacecolor='none', label="KCD VASP GGA")
    proxy_PWSCF_D, = ax.plot([], [], "C2o--", markerfacecolor='none', label="DVC PWSCF GGA")
    # proxy_CASTEP, = ax.plot([], [], "C2o--", markerfacecolor='none', label="SBC CASTEP")
    proxy_VASP, = ax.plot([], [], "C3o--", markerfacecolor='none', label="This work VASP")

    ax.legend(loc="upper left", handles=[DFT,
                                         proxy_siesta,
                                         proxy_PWSCF_D,
                                        # proxy_CASTEP,
                                         proxy_KCD_VASP,
                                         proxy_VASP],
              fontsize=12)

    rect = (0.095, 0.465, 0.28, 0.14)
    i = 0
    first_pure_core = add_core_plot(fig, i, rect,
                                    pure_glide_averaged_images,
                                    QM_MM_bulk,
                                    QM_MM_core_positions,
                                    spine_color="C0",
                                    spine_lw=2.0, core_lw=2.0)

    add_lines_to_core(ax, first_pure_core, s_W[0], E_W[0],
                      orientation="bottom",
                      color="C0", x_shift=0.07)

    rect = (0.375, 0.12, 0.28, 0.14)
    i = 3
    second_pure_core = add_core_plot(fig, i, rect,
                                     pure_glide_averaged_images,
                                     QM_MM_bulk,
                                     QM_MM_core_positions,
                                     spine_color="C0",
                                     spine_lw=2.0, core_lw=2.0)

    add_lines_to_core(ax, second_pure_core , s_W[3], E_W[3],
                      orientation="top",
                      color="C0", x_shift=0.065)

    rect = (0.65, 0.465, 0.28, 0.14)
    i = -1
    last_pure_core = add_core_plot(fig, i, rect,
                                   pure_glide_averaged_images,
                                   QM_MM_bulk,
                                   QM_MM_core_positions,
                                   spine_color="C0",
                                   spine_lw=2.0, core_lw=2.0,
                                   legend=True,
                                   legend_bbox=(0.20, 0.95, 0.5, 0.5))

    add_lines_to_core(ax, last_pure_core , s_W[-1], E_W[-1],
                      orientation="bottom",
                      color="C0", x_shift=0.07)


    #plt.show()

    fig.savefig("Averaged_QMMM_vs_DFT_MD.pdf")
    fig.savefig("Averaged_QMMM_vs_DFT_MD.png", dpi=300)
