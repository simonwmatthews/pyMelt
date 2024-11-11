# This script contains some functions that help in using MAGEMin to construct
# phase diagram grids for use with pyMelt.
#
# Simon Matthews (simonm@hi.is), September 2024

using MAGEMin_C
using DataFrames
using CSV
using Printf

function create_table()
    df = DataFrame(P=Float64[], T=Float64[]) 
    return df
end

function add_output_to_table(df, out, P, T, rowstart)
    row = rowstart
    rowtemp = Dict(zip(out.ph, out.ph_frac_wt))
    row = merge(row, rowtemp)
    row["P"] = P
    row["T"] = T
    row["density"] = out.rho
    row["enthalpy"] = out.enthalpy
    row["entropy"] = out.entropy
    row["thermal_expansivity"] = out.alpha
    row["heat_capacity"] = out.cp
    row["fO2"] = out.fO2
    row["fO2_QFM"] = out.dQFM

    ssveccount = 0
    for i in 1:size(out.ph, 1)
        if out.ph_type[i] == 1
            ssveccount = ssveccount + 1

            phname = out.ph[i]
            # OXIDE WTPT
            oxnames = out.oxides
            oxnames = [phname*"_wtpt_"*ox for ox in oxnames]
            oxwtpt = out.SS_vec[ssveccount].Comp_wt * 100
            rowtemp = Dict(zip(oxnames, oxwtpt))
            row = merge(row, rowtemp)

            # ENDMEMBERS
            # Get names of model endmembers
            emnames = out.SS_vec[ssveccount].emNames
            # Add the mineral name in front of endmember name
            emnames = [phname*"_"*en for en in emnames]
            # Get the endmember fractions
            emfracs = out.SS_vec[ssveccount].emFrac
            # Create dictionary and merge it with row
            rowtemp = Dict(zip(emnames, emfracs))
            row = merge(row, rowtemp)
        end
    end

    for kk in keys(row)
        if hasproperty(df, kk) == false
            df[!, kk] .= 0.0
        end
    end

    for kk in names(df)
        if haskey(row, kk) == false
            row[kk] = 0.0
        end
    end

    push!(df, row)
end

function run_tp_grid(comp, T0, T1, nT, P0, P1, nP, Xoxides, sys_in, fname)
    data = Initialize_MAGEMin("ig", verbose=false)
    df = create_table()
    for i in 1:nP
        @printf("Pressure %i of %i\n", i, nP)
        P = (i-1) / (nP - 1) * P1 + (1 - (i-1)/(nP-1)) * P0
        if P == 0
            P = 0.01
        end
        for j in 1:nT
            T = (j-1) / (nT - 1) * T1 + (1 - (j-1)/(nT-1)) * T0
            out = single_point_minimization(P, T, data, X=comp, Xoxides=Xoxides, sys_in=sys_in)
            add_output_to_table(df, out, P, T, Dict([]))
        end
        CSV.write(fname, df)
    end
    return df
end