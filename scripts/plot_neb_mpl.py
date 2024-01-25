import os

import numpy as np
import matplotlib.pyplot as plt

import matscipy.dislocation as sd

from ase.utils.forcecurve import fit_images
from scipy.interpolate import interp1d

from ase.io import read
from ase.neb import NEB

def plot_neb(images, ax=None, label="", forces_only=False, show=True,
             symmetrize=False, hollow=False, marker="o", linestyle="dashed",
             mlw=1.5,  lw=2.0, **kwargs):

    if ax is None:
        fig, ax = plt.subplots()

    force_fit = fit_images(images)


    if forces_only:
        s, E, _ = NEB(images).integrate_forces()
        if symmetrize:
            E = 0.5 * (E + E[::-1])
            E -= E.min()

        fitted_E = interp1d(s, E)


        ax.plot(s, E, linestyle=linestyle, marker=None, lw=lw, **kwargs)
        ax.plot(force_fit.path / force_fit.path[-1], 
                   fitted_E(force_fit.path / force_fit.path[-1]), linestyle="none",
                   mew=mlw,  marker=marker, **kwargs)
        color = ax.get_lines()[-1].get_color()
    else:
        ax.plot(force_fit.fit_path / force_fit.path[-1], force_fit.fit_energies, lw=lw,
        linestyle=linestyle, **kwargs)
        color = ax.get_lines()[-1].get_color()
        if hollow:
            if "zorder" in kwargs:
                zorder = kwargs.get("zorder")
                del kwargs["zorder"]
            ax.scatter(force_fit.path / force_fit.path[-1], force_fit.energies, mew=mlw, facecolors="white", edgecolor=color, s=ms, marker=marker, zorder=zorder+1, **kwargs)
        else:
            ax.plot(force_fit.path / force_fit.path[-1], force_fit.energies, mew=mlw, linestyle="none", marker=marker, **kwargs)


    if forces_only:
        handle, = ax.plot([], [], marker=marker, linestyle=linestyle,
        label=label, lw=0.75*lw, markeredgewidth=0.75*mlw, **kwargs)
    else:
        if hollow:
            handle, = ax.plot([], [], maker=marker,
                              color=color, linestyle=linestyle, label=label, lw=0.75*lw,  
                              markeredgewidth=0.75*mlw, markerfacecolor="white", markeredgecolor=color)
        else:
            handle, = ax.plot([], [], marker=marker, linestyle=linestyle, 
                              label=label, lw=0.75*lw, markeredgewidth=0.75*mlw, **kwargs)

    if show:
        ax.legend()
        ax.set_ylabel("Energy (eV)")
        ax.set_xlabel("Normalised NEB path")
    else:
        return handle

colors = {#"WHff": "Wang et. al.",
          "MCM2": "black",
          "Mason17": "#0036FA",
          "Hiremath22": "#FFA12B",
          "TabGAP_v2": "#FF0019",
          "Juslin05": "#00831D"}


labels = {#"WHff": "Wang et. al.",
          "MCM2": "EAM",
          "Mason17": "EAM/ZBL",
          "Hiremath22": "MEAM",
          "TabGAP_v2": "TabGAP",
          "Juslin05": "ABOP"}
ms=6.5
marker_styles = {#"WHff": "Wang et. al.",
                 "MCM2": {"marker": "s", "ms" : f"{ms}"},
                 "Mason17": {"marker": "^", "fillstyle": "none", "ms" : f"{ms}"},
                 "Hiremath22": {"marker": ">", "ms" : f"{ms}"},
                 "TabGAP_v2": {"marker": "o", "fillstyle": "none", "ms" : f"{ms}"},
                 "Juslin05": {"marker": "*", "fillstyle": "none", "markersize": f"{1.2 * ms}"}}


if __name__ == '__main__':

    qmml_images = {}

    fig, axes = plt.subplots(figsize=(7, 4), ncols=2, sharey=True, sharex=True)

    lw = 1.9

    for potential in snakemake.params.potentials:
        for ax, dislo_name in zip(axes, ["junction", "junction110"]):
            junction_images = read(f"../results/{dislo_name}_glide_neb_{potential}.xyz", index=":")
            plot_neb(junction_images, ax=ax, label=labels[potential], show=False, lw=lw, color=colors[potential], **marker_styles[potential])
            #plot_neb(refitted_LML_images[dislo_name], ax=ax, label="Retrained LML $\Theta_0+\delta\Theta$", show=False, marker="o", hollow=False, lw=lw, zorder=10)

    for ax, dislo_name in zip(axes, ["junction", "junction110"]):
        qmml_images = read(f"../data/{dislo_name}_glide_neb_QMML.xyz", index=":")
        plot_neb(qmml_images, ax=ax, forces_only=True, linestyle="solid", marker="x",
                 label="QM/ML", show=False, symmetrize=True, ms=7.5, mlw=2.0, c="grey", lw=2.0, zorder=11)
        ax.set_xlabel("Normalised NEB path", fontsize=13)
        ax.set_xlim((0.0, 1.0))

    axes[0].set_ylabel("Energy (eV)", fontsize=13)
    #axes[1].set_ylabel("Energy (eV)", fontsize=13)
    axes[0].text(0.1, 0.8, "a)", fontsize=15)
    axes[1].text(0.8, 0.8, "b)", fontsize=15)
    axes[0].set_title(r"Junction $ \leftangle $100$ \rightangle ${001} glide")
    axes[1].set_title(r"Junction $ \leftangle $100$ \rightangle ${011} glide")
    fig.tight_layout()


    # axes[1].legend(loc="lower left", ncol=4, bbox_to_anchor=(-1.4, 1.075, 0.5, 0.5), framealpha=1.0)
    # axes[1].legend(loc="upper center", fontsize=10)
    axes[1].legend(loc="upper center", fontsize=10, ncol=2,
    bbox_to_anchor=(-0.25, 0.475, 0.5, 0.5), framealpha=1.0)
    axes[0].set_ylim((-0.32, 1.1))

    fig.savefig(snakemake.output.junction_glide_plot)
    fig.savefig("../results/junction_glide.png")
    fig, ax = plt.subplots(figsize=(4, 4))

    for potential in snakemake.params.potentials:
        screw_images = read(f"../results/screw_glide_neb_{potential}.xyz", index=":")
        plot_neb(screw_images, ax=ax, label=labels[potential], show=False, lw=lw, color=colors[potential], **marker_styles[potential])

    qmml_images = read(f"../data/screw_glide_neb_QMML.xyz", index=":")
    plot_neb(qmml_images, ax=ax, forces_only=True, linestyle="solid", marker="x",
             label="QM/ML", show=False, symmetrize=True, ms=7.5, mlw=2.0, c="grey", lw=2.0, zorder=11)
    ax.set_xlabel("Normalised NEB path", fontsize=13)
    ax.set_ylabel("Energy (eV)", fontsize=13)
    ax.set_xlim((0.0, 1.0))
    ax.set_ylim((-0.21, 0.25))
    ax.legend(loc="upper center", fontsize=10, ncol=2)

    fig.tight_layout()
    fig.savefig(snakemake.output.screw_glide_plot)
    fig.savefig("../results/screw_glide.png")
