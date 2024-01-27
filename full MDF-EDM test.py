# -*- coding: utf-8 -*-
"""
Created on Wed Oct 11 15:09:48 2023

@author: galai
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from equilibrator_api import Q_, ComponentContribution
from matplotlib.backends.backend_pdf import PdfPages
from sbtab import SBtab
from tqdm import tqdm

from equilibrator_pathway import EnzymeCostModel, ThermodynamicModel


with PdfPages("ecm_results.pdf") as pdf:
    comp_contrib = ComponentContribution()

    # %% example for MDF analysis

    pp = ThermodynamicModel.from_sbtab(
        "ecoli_noor_2016_mdf.tsv", comp_contrib=comp_contrib
    )

    mdf_result = pp.mdf_analysis()

    fig1 = mdf_result.plot_concentrations()
    pdf.savefig(fig1)

    fig2 = mdf_result.plot_driving_forces()
    pdf.savefig(fig2)

    # %% example for ECM analysis

    model = EnzymeCostModel.from_sbtab(
        "ecoli_noor_2016_ecm.tsv", comp_contrib=comp_contrib
    )
    model.add_validation_data("ecoli_noor_2016_reference.tsv")

    print("Solving ECM problem")
    np.random.seed(1982)
    ecm_sol = model.optimize_ecm()
    res_sbtab = ecm_sol.to_sbtab()

    fig3, ax = plt.subplots(1, 1, figsize=(7, 5))
    ax.set_title("ECM solution")
    ecm_sol.plot_enzyme_demand_breakdown(ax, plot_measured=True)
    ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))
    fig3.tight_layout()
    pdf.savefig(fig3)

    fig4 = plt.figure(figsize=(16, 5))
    ax = fig4.add_subplot(1, 3, 1, xscale="log", yscale="log")
    ax.set_title("Metabolite Concentrations")
    ecm_sol.validate_metabolite_conc(ax)

    ax = fig4.add_subplot(1, 3, 2, xscale="log", yscale="log")
    ax.set_title("Enzyme Concentrations")
    ecm_sol.validate_enzyme_conc(ax)

    ax = fig4.add_subplot(1, 3, 3)
    ecm_sol.plot_volumes_pie(ax=ax)

    pdf.savefig(fig4)

output_fname = "ecm_results.tsv"
print(f"Writing the ECM results to file: {output_fname}")
res_sbtab.write(output_fname)
