# Construct the KG1 phase diagram from Holland et al. (2018)
include("../magemin_functions.jl")

Xoxides = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "Fe2O3"; "K2O"; "Na2O"; "TiO2"; "Cr2O3"; "H2O"];
kg1 = [46.97, 9.75, 7.35, 23.57, 9.77*(1-0.032), 9.77*0.032/(55.845+16)*(55.845*2+16*3)/2, 0.12, 1.52, 0.78, 0.17, 0.0]
sys_in = "wt"

T0 = 1100.0
T1 = 2100.0
nT = 200
P0 = 0.0
P1 = 80.0
nP = 200

df = run_tp_grid(kg1, T0, T1, nT, P0, P1, nP, Xoxides, sys_in, "/Users/sm905/repos/pyMelt/phaseDiagrams/kg1_holland2018/table_kg1_holland2018.csv")