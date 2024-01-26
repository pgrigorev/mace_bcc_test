from scipy.interpolate import interp1d
from ase.utils.geometry import find_mic
from scipy.integrate import cumtrapz
import numpy as np
import matplotlib.pyplot as plt
from ase.neb import fit0
from ase.utils.forcecurve import fit_raw
import sys

sys.path.insert(-1, "../LIB/")
from nebForceIntegrator import NEB
import matscipy.dislocation as sd
from ase.calculators.singlepoint import SinglePointCalculator

def plot_H_at_core(H_image, bulk,
                   core_position=None,
                   ax=None,
                   marker="d",
                   imp_size=300, imp_color="C4",
                   core_size=150,
                   imp_symbol="H",
                   legend_bbox=(1.0, 0.95),
                   title="QM/MM VASP/EAM3",
                   xyscale=6,
                   insert=False):

    symbols = np.array(H_image.get_chemical_symbols())
    W_atoms = H_image[symbols == "W"]
    H_pos = H_image[symbols == imp_symbol].positions.T

    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 6))

    ax.scatter(H_pos[0], H_pos[1], s=imp_size, c=imp_color,
               marker=marker, edgecolor="k", zorder=15)

    core = ax.scatter([], [], s=100, marker="^",
                          color="C3", label="core position")

    W = ax.scatter([], [], marker="o", c="white", s=100, edgecolors="k", label="W atoms")

    if imp_symbol is not None:
        H = ax.scatter([], [], marker=marker, c=imp_color, s=125, edgecolor="k",
                       label=imp_symbol + ' atom')

        handles = [core, W, H]

    else:

        handles = [W, core]

    sd.plot_vitek(W_atoms, bulk, alat=3.17,
                  plot_axes=ax, xyscale=xyscale)

    if not insert:

        ax.legend(loc="upper right", handles=handles,
                  ncol=1, bbox_to_anchor=legend_bbox)
        ax.set_xlim((-2.5, 5.5))
        ax.set_ylim((-4, 4))
        ax.set_xlabel("<112> coordinate, $\AA$", fontsize=15)
        ax.set_ylabel("<110> coordinate, $\AA$", fontsize=15)
        ax.set_title(title, fontsize=15)

    ax.set_aspect('equal')

    if core_position is not None:
        core_x, core_y = core_position
        ax.scatter(core_x, core_y, s=core_size, marker="^", color="C3", zorder=10)

    if ax is None:
        plt.show()


def plot_band_forces(images, forces_label, label="averaged barrier",
                     s_norm=True,
                     ax=None, **kwargs):

    if ax is None:
        fig, ax = plt.subplots()

    for im in images:
        im.set_calculator(SinglePointCalculator(im, forces=im.get_array(forces_label)))

    neb = NEB(images, force_only=True)
    neb.get_forces()

    E = neb.get_potential_energies()

    R = [atoms.positions for atoms in images]
    F =  [image.get_forces() for image in images]
    A = images[0].cell
    pbc = images[0].pbc
    s, __, Sfit, Efit, lines = fit_raw(E, F, R, A, pbc)

    s = np.array(s)

    if s_norm:

        norm = s.max()
        s /= norm
        Sfit /= norm

    s_e_plot = ax.plot(Sfit, Efit, label=label, **kwargs)
    color = s_e_plot[0].get_color()
    ax.errorbar(s, E, fmt='o',
                yerr=[0.0, 0.0, 0.0, 0.0279*E[3], 0.0, 0.0, 0.039*E[3]],
                c=color,
                **kwargs)

    for x, y in lines:
        if s_norm:
            x /= norm
        ax.plot(x, y, linestyle="dashed", c=color, **kwargs)

    ax.set_xlabel("Normalised NEB coordinate")
    ax.set_ylabel("Energy (eV)")

    ax.grid(True, linestyle="dashed")

    return Sfit, Efit, s, E, F, R, A, pbc


def plot_averaged_barrier(images, final_dataframe, label="label", marker="o",
                          ax=None, s_norm=False, s_shift=0.0, **kwargs):

    add_legend=False

    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 5))
        add_legend = True

    #E = final_dataframe["E"].copy().values
    s = final_dataframe["s"].copy().values

    for im in images:
        im.set_calculator(SinglePointCalculator(im, forces=im.get_array("averaged_forces")))

    neb = NEB(images, force_only=True)
    neb.get_forces()

    E = neb.get_potential_energies()
    R = [atoms.positions for atoms in images]
    F = neb.forces
    A = images[0].cell
    pbc = images[0].pbc
    __, __, Sfit, Efit, __ = fit0(E, F, R, A, pbc)

    if s_norm:

        norm = s.max()
        s /= norm
        Sfit /= norm

    s += s_shift
    Sfit += s_shift

    s_e_plot = ax.plot(Sfit, Efit, label=label, **kwargs)
    color = s_e_plot[0].get_color()

    result = ax.errorbar(s, E, fmt=marker,
                         yerr=final_dataframe["std_E"].values,
                         xerr=final_dataframe["std_s"].values,
                         c=color,
                         label=label,
                         **kwargs)

    __, __, __, Efit_plus, __ = fit0(E + final_dataframe["std_E"].values,
                                     F, R, A, pbc)

    __, __, __, Efit_minus, __ = fit0(E - final_dataframe["std_E"].values,
                                      F, R, A, pbc)

    ax.fill_between(Sfit,
                    Efit_plus,
                    y2=Efit_minus, color=color, alpha=0.25)

    if add_legend:
        ax.legend(loc="best")

    if s_norm:
        ax.set_xlabel("Normalised NEB coordinate", fontsize=15)
    else:
        ax.set_xlabel("Absolute NEB coordinate", fontsize=15)

    ax.set_ylabel("Energy (eV)", fontsize=15, labelpad=0)

    ax.grid(True, linestyle="dashed")

    """
    E = neb.get_potential_energies()
    R = [atoms.positions for atoms in neb.images]
    F = neb.forces
    A = neb.images[0].cell
    pbc = neb.images[0].pbc
    s, E, Sfit, Efit, __ = fit0(E, F, R, A, pbc)

    s = np.array(s)
    if s_norm:

        norm = s.max()
        s /= norm
        Sfit /= norm

    s += s_shift

    ax.scatter(s, E, marker="x", c="k", zorder=10,
           label="neb.forces")
    """
    return s, E, result

def add_lines_to_core(ax, core_ax, s, E,
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
