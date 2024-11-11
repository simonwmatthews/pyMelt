# Construct the G2 phase diagram from Holland et al. (2018)
include("../magemin_functions.jl")

Xoxides = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "Fe2O3"; "K2O"; "Na2O"; "TiO2"; "Cr2O3"; "H2O"];
g2 = [50.05, 15.76, 11.74, 7.90, 9.35*(1-0.032), 9.35*0.032/(55.845+16)*(55.845*2+16*3)/2, 0.03, 3.04, 1.97, 0.3, 0.0]
sys_in = "wt"

T0 = 1000.0
T1 = 1800.0
nT = 250
P0 = 0.0
P1 = 50.0 
nP = 200

df = run_tp_grid(g2, T0, T1, nT, P0, P1, nP, Xoxides, sys_in, "/Users/sm905/repos/pyMelt/phaseDiagrams/g2_holland2018/table_g2_holland2018.csv")