# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     custom_cell_magics: kql
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.11.2
#   kernelspec:
#     display_name: Julia 1.10.5
#     language: julia
#     name: julia-1.10
# ---

# %% [markdown]
# # Generate KLB-1 phase diagram from Holland et al. (2018)
# ## Part 1 : Run MAGEMin over grid
#
# **This notebook must be run with julia**
#
# This notebook generates a phase diagram grid using MAGEMin, supported by the `magemin_functions.jl` script. Ideally this would be in the same notebook as the workflow for generating the phaseDiagram object, but python and julia will not play nicely on my machine.

# %%
include("../magemin_functions.jl")

# %%
Xoxides = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "Fe2O3"; "K2O"; "Na2O"; "TiO2"; "Cr2O3"; "H2O"];
klb1 = [44.84, 3.51, 3.07, 39.52, 8.20*(1-0.032), 8.20*0.032/(55.845+16)*(55.845*2+16*3)/2, 0.01, 0.3, 0.11, 0.01, 0.0]
sys_in = "wt"

T0 = 1100.0
T1 = 2100.0
nT = 200
P0 = 0.0
P1 = 80.0
nP = 200

# %% [markdown]
# Note that I have to specify the exact path to place the output file.

# %%
df = run_tp_grid(klb1, T0, T1, nT, P0, P1, nP, Xoxides, sys_in, "/Users/sm905/repos/pyMelt/phaseDiagrams/klb1_holland2018/table_klb1_holland2018.csv")


# %%
