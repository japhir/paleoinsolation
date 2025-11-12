# Plotting Scripts for PaleoInsolation
# This script generates the figures in the paper.
# Copyright (C) 2025 Ilja J. Kocken

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.


################################################################################
###                            packages and globals                          ###
################################################################################

using CSV
using DataFrames
# I have both, one for interactive exploration, 
using GLMakie
# the other for final figs.
using CairoMakie
using AlgebraOfGraphics
using PaleoInsolation # the package I wrote for this paper
using Arrow
using RCall # for astrochron and palinsol
# global settings
set_theme!(theme_latexfonts())
update_theme!(
    Theme(Axis = (xgridvisible = false, ygridvisible = false)),
    fontsize=12,
    joinstyle=:round,
)
# update_theme!(fontsize=12) # TODO: make legend patches smaller?
inch = 96
pt = 4/3
cm = inch / 2.54
# GLMakie.activate!() # for interactive exploration
CairoMakie.activate!() # to save PDFs for manuscript
global_S0 = 1367.1 # Menviel 2019, same as models
# summer insolation label
sil = rich("65°N peak summer\ninsolation (Wm", superscript("−2"),")")




################################################################################
###                               read data                                  ###
################################################################################

# snvec  c: computed PT-ZB18a_1-1.dat
PT_ZB18a = CSV.read("dat/PT-ZB18a_1-1.dat",
                    DataFrame;
                    delim = ' ',
                    stripwhitespace = true,
                    ignorerepeated = true,
                    header = [:time,
                              :eccentricity,
                              :obliquity,
                              :precession,
                              :lpx,
                              :climatic_precession]
                    )
PT_ZB18a.year .= PT_ZB18a.time * 1e3 .+ 2000
PT_ZB18a.insolation .= PaleoInsolation.insolation.(
    PT_ZB18a.eccentricity,
    PT_ZB18a.obliquity,
    mod.(PT_ZB18a.lpx .- pi, 2pi),
    longitude = pi/2, latitude = 65*pi/180,
    S0 = global_S0, H = nothing
)

# snvec  c: computed PT-ZB20a_1-1.dat
PT_ZB20a = CSV.read("dat/PT-ZB20a_1-1.dat",
                    DataFrame;
                    delim = ' ',
                    stripwhitespace = true,
                    ignorerepeated = true,
                    header = [:time,
                              :eccentricity,
                              :obliquity,
                              :precession,
                              :lpx,
                              :climatic_precession]
                    )

## get ZB18a(1,1) 65° N ins from Richard's webtool
# https://astrocyclo.soest.hawaii.edu/
zb_65 = CSV.read("dat/ZB18a_1-1_65_90.dat",
                 DataFrame,
                 comment = "%", header = [:time, :insolation],
                 delim = ' ', ignorerepeated = true)

ins_1361 = PaleoInsolation.insolation.(PT_ZB18a.eccentricity,
                         PT_ZB18a.obliquity,
                         # x.lpx; # this MUST be rewrapped!
                         mod.(PT_ZB18a.lpx .- pi, 2pi),
                         latitude = 65*pi/180,
                         longitude = pi/2,
                         S0 = 1361, H = nothing)
ins_1365 = PaleoInsolation.insolation.(PT_ZB18a.eccentricity,
                         PT_ZB18a.obliquity,
                         # x.lpx; # this MUST be rewrapped!
                         mod.(PT_ZB18a.lpx .- pi, 2pi),
                         latitude = 65*pi/180,
                         longitude = pi/2,
                         S0 = 1365, H = nothing)

## get pmip3
pmip3 = CSV.read("dat/Schmidt2011_pmip3_lm_orbital_parameters.txt",
                 DataFrame,
                 skipto = 23, delim = ' ', ignorerepeated = true, ignoreemptyrows = true,
                 header = [:year, :eccentricity, :obliquity, :perimin180])
pmip3.time = (pmip3.year .- 1950) .* 1e-3 # convert to time relative to 1950 in kyr
pmip3.lpx = deg2rad.(pmip3.perimin180 .+ 180)
pmip3.insolation = PaleoInsolation.insolation.(pmip3.eccentricity,
                                     deg2rad.(pmip3.obliquity),
                                     pmip3.lpx,
                                     # mod.(pmip3.lpx .- pi, 2pi),
                                     S0 = global_S0)

# I prepped these in R using palinsol::ber78 and palinsol::ber90
# therefore, the lpx is already - pi'd!
ber_df = CSV.read("dat/Ber78_from_palinsol.csv", DataFrame,
                  skipto = 2,
                  header = [:time, :obliquity, :eccentricity, :lpx, :epsp])
# ber_df.climatic_precession = ber_df.eccentricity .* sin.(ber_df.lpx)
# must be re-wrapped for this calc
ber_df.climatic_precession = ber_df.eccentricity .* sin.(mod.(ber_df.lpx .- pi, 2pi))
ber_df.insolation = PaleoInsolation.insolation.(ber_df.eccentricity,
                                      ber_df.obliquity,
                                      ber_df.lpx,
                                      longitude = pi/2, latitude = 65*pi/180,
                                      S0 = global_S0, H = nothing)
ber_df.year = ber_df.time * 1e3 .+ 1950
ber91_df = CSV.read("dat/Ber91_from_palinsol.csv", DataFrame,
                  skipto = 2,
                  header = [:time, :obliquity, :eccentricity, :lpx, :epsp])
# ber91_df.climatic_precession = ber91_df.eccentricity .* sin.(ber91_df.lpx)
ber91_df.climatic_precession = ber91_df.eccentricity .* sin.(mod.(ber91_df.lpx .- pi, 2pi))
ber91_df.insolation = PaleoInsolation.insolation.(ber91_df.eccentricity,
                                        ber91_df.obliquity,
                                        ber91_df.lpx,
                                        longitude = pi/2, latitude = 65*pi/180,
                                        S0 = global_S0, H = nothing)
ber91_df.year = ber91_df.time * 1e3 .+ 1950

# La04 from raw fortran included file orbital parameters
las_xss = CSV.read("La04/INSOLN.LA2004.BTL.100.ASC",
                    DataFrame,
                    delim = ' ', ignorerepeated = true,
                    header = [:time, :eccentricity, :obliquity, :lpx])
las_xss.year = las_xss.time * 1e3 .+ 2000
# but because we specify time in kyr, that 5% diff isn't going to make the difference
@. las_xss.climatic_precession = las_xss.eccentricity * sin(las_xss.lpx)
las_xss.insolation = PaleoInsolation.insolation.(las_xss.eccentricity,
                                       las_xss.obliquity,
                                       # las_xss.lpx; # this MUST be rewrapped!
                                       mod.(las_xss.lpx .- pi, 2pi),
                                       latitude = 65*pi/180,
                                       longitude = pi/2,
                                       S0 = global_S0, H = nothing)
# we also provide some other sources of La04
tend = -100
times_las = tend:1:0
times_kyr = 0.0:-1:tend
# this is quite slow
las_ins_65_90_R = rcopy(R"Map(\(x) palinsol::Insol(palinsol::la04(x), S0 = $global_S0), $times_las*1e3)")

# other sources of las
# the binary file with insola probably has least rounding
# get La04 insolation by running insola function
# from their own website
# https://vo.imcce.fr/insola/earth/online/earth/La2004
# get La04 using the web tool:
# https://vo.imcce.fr/insola/earth/online/earth/online/index.php
# from -1 Myr to 0 Myr every 1000 years S0 = 1368, mean daily insolation / true longitude
las_xss2 = CSV.read("dat/La04_web_65_90_true.dat",
                     DataFrame,
                     delim = ' ', ignorerepeated = true,
                     header = [:time, :insolation])
# generated using insola Fortran codes
las_xss3 = CSV.read("La04/insolation_mean90_65.dat",
                     DataFrame,
                     delim = ' ', ignorerepeated = true,
                     header = [:time, :insolation])
las_xss3.insolation = parse.(Float64, replace.(las_xss3.insolation, "D" => "E"))
# insola using true longitude
las_xss5 = CSV.read("La04/insolation_true90_65.dat",
                    DataFrame,
                    delim = ' ', ignorerepeated = true,
                    header = [:time, :insolation])
las_xss5[!, :insolation] = parse.(Float64, replace.(las_xss5[!, :insolation], "D" => "E"))
las_xss6 = CSV.read("La04/insola_mean_daily_true_lon_1367.1.txt",
                    DataFrame,
                    delim = ' ', ignorerepeated = true,
                    header = [:time, :insolation])
las_xss6[!, :insolation] = parse.(Float64, replace.(las_xss6[!, :insolation], "D" => "E"))


# La10b = CSV.read(download("http://vo.imcce.fr/insola/earth/online/earth/La2010/La2010b_ecc3L.dat"),
#                  DataFrame,
#                  delim = ' ',
#                  stripwhitespace = true,
#                  ignorerepeated = true,
#                  header = ["time", "eccentricity"]
#                  )
# cache
# Arrow.write("out/La10b.arrow", La10b)
La10b = DataFrame(Arrow.Table("out/La10b.arrow"))

# to contrast, we also calculate insolation using palinsol
# ber78_ins_palinsol = rcopy(R"""
# Map(\(x) palinsol::Insol(orbit = x, long = pi/2, lat = 65*pi/180, S0 = 1367.1, H = NULL),
#     list(varpi = $(ber_df.lpx), eps = $(ber_df.obliquity), ecc = $(ber_df.eccentricity))
# )
# """)
# this was too hard with the mapping and the interface between Julia and R...
ber_ins_R = rcopy(R"Map(\(x) palinsol::Insol(palinsol::ber78(x), S0 = $global_S0), $times_las*1e3)")
# note that these ages are offset a little because palinsol uses 1950


################################################################################
###                                 plots                                    ###
################################################################################

# Figure 1
zbl = "ZB18a(1,1)"
fsol = Figure()
ax_ins = Axis(fsol[1,1], 
              ylabel = sil,
              xlabel = "Time (Myr)")
# ins
lines!(ax_ins,
       (PT_ZB18a.year .- 2000) * 1e-6,
       PT_ZB18a.insolation,
       color = Makie.wong_colors()[1],
       linewidth = 3, label = zbl)
lines!(ax_ins,
       (ber_df.year .- 2000) * 1e-6,
       ber_df.insolation,
       label = "Ber78", color = Makie.wong_colors()[4], alpha = 0.8)
# lines!(ax_ins,
#        (ber91_df.year .- 2000) * 1e-6,
#        ber91_df.insolation,
#        label = "BL91", color = Makie.wong_colors()[6], alpha = 0.8)
lines!(ax_ins,
       (las_xss.year .- 2000) * 1e-6,
       las_xss.insolation,
       label = "La04", linewidth = 1.5, linestyle = :dash,
       color = Makie.wong_colors()[2]
       )

# no La10b for this one
Legend(fsol[0,1], ax_ins, orientation = :horizontal, halign = :right)
rowgap!(fsol.layout, 5)
xlims!(-2.5, -2)

save("imgs/compare_solutions_insolation.pdf",
     fsol, size=(5.5inch, 2inch), px_per_unit = 300/inch)

fsol


# The extended figure
# Figure A1: fsol
zbl = "ZB18a(1,1)"
fsol = Figure()
ax_ecc = Axis(fsol[1,1],
              ylabel = "Eccentricity (-)\nClimatic precession (-)", xlabel = "", xticklabelsvisible = false)
ax_obl = Axis(fsol[2,1],
              ylabel = "Obliquity (°)", xlabel = "", xticklabelsvisible = false)
ax_ins = Axis(fsol[3,1], 
              ylabel = sil,
              xlabel = "", xticklabelsvisible = false)
ax_insdiff = Axis(fsol[4,1], 
              ylabel = rich("Δ insolation\n(Wm", superscript("−2"),")"),
                  xlabel = "Time (Myr)",
                  yticks = -90:30:90)
linkxaxes!(ax_ecc, ax_obl, ax_ins, ax_insdiff)

# ecc
lines!(ax_ecc,
       # we put everything in Myr before 2000 explicitly
       (PT_ZB18a.year .- 2000) * 1e-6,
       PT_ZB18a.eccentricity,
       color = Makie.wong_colors()[1],
       linewidth = 3, label = zbl)
lines!(ax_ecc,
       (ber_df.year .- 2000) * 1e-6, 
       # ber_df.year * 1e-6 .- 5.0e-5, # convert to t0 = 2000
       # (ber_df.year .- 1950) * 1e-6,
       ber_df.eccentricity, label = "Ber78", color = Makie.wong_colors()[4], alpha = 0.8)
lines!(ax_ecc,
       (ber91_df.year .- 2000) * 1e-6,
       ber91_df.eccentricity, label = "BL91",
       color = Makie.wong_colors()[6], alpha = 0.8)
lines!(ax_ecc, (las_xss.year .- 2000) * 1e-6, las_xss.eccentricity,
       label = "La04", linewidth = 1.5, linestyle = :dash,
       color = Makie.wong_colors()[2]
       )
# cp
lines!(ax_ecc, (PT_ZB18a.year .- 2000) * 1e-6,
       PT_ZB18a.climatic_precession,
       color = Makie.wong_colors()[1],
       linewidth = 3, label = zbl)
lines!(ax_ecc,
       (ber_df.year .- 2000) * 1e-6, 
       ber_df.climatic_precession,
       label = "Ber78", color = Makie.wong_colors()[4], alpha = 0.8)
lines!(ax_ecc,
       (ber91_df.year .- 2000) * 1e-6,
       ber91_df.climatic_precession,
       label = "BL91", color = Makie.wong_colors()[6], alpha = 0.8)
lines!(ax_ecc, (las_xss.year .- 2000) * 1e-6,
       las_xss.climatic_precession,
       label = "La04", linewidth = 1.5, linestyle = :dash,
       color = Makie.wong_colors()[2]
       )

# obl
lines!(ax_obl,
       (PT_ZB18a.year .- 2000) * 1e-6,
       rad2deg.(PT_ZB18a.obliquity),
       color = Makie.wong_colors()[1],
       linewidth = 3, label = zbl)
lines!(ax_obl,
       (ber_df.year .- 2000) * 1e-6,
       rad2deg.(ber_df.obliquity),
       label = "Ber78", color = Makie.wong_colors()[4], alpha = 0.8)
lines!(ax_obl,
       (ber91_df.year .- 2000) * 1e-6,
       rad2deg.(ber91_df.obliquity), label = "BL91", color = Makie.wong_colors()[6], alpha = 0.8)
lines!(ax_obl,
       (las_xss.year .- 2000) * 1e-6,
       rad2deg.(las_xss.obliquity),
       label = "La04", linewidth = 1.5, linestyle = :dash,
       color = Makie.wong_colors()[2]
       )
# ins
lines!(ax_ins,
       (PT_ZB18a.year .- 2000) * 1e-6,
       PT_ZB18a.insolation,
       color = Makie.wong_colors()[1],
       linewidth = 3, label = zbl)
lines!(ax_ins,
       (ber_df.year .- 2000) * 1e-6,
       ber_df.insolation,
       label = "Ber78", color = Makie.wong_colors()[4], alpha = 0.8)
lines!(ax_ins,
       (ber91_df.year .- 2000) * 1e-6,
       ber91_df.insolation,
       label = "BL91", color = Makie.wong_colors()[6], alpha = 0.8)
lines!(ax_ins,
       (las_xss.year .- 2000) * 1e-6,
       las_xss.insolation,
       label = "La04", linewidth = 1.5, linestyle = :dash,
       color = Makie.wong_colors()[2]
       )

# linearly interpolate our data to Ber78 and Ber91 timescales
# in order to calculate differences
using DataInterpolations
PT_ZB18a_interp = LinearInterpolation(
    # they need the time in increasing order for interpolation
    reverse(PT_ZB18a.insolation), reverse(PT_ZB18a.year))

lines!(ax_insdiff,
       (ber_df.year .- 2000) * 1e-6,
       ber_df.insolation .- PT_ZB18a_interp.(ber_df.year),
       label = "Ber78", color = Makie.wong_colors()[4])
lines!(ax_insdiff,
       (ber91_df.year .- 2000) * 1e-6,
       ber91_df.insolation .- PT_ZB18a_interp.(ber91_df.year),
       label = "BL91", color = Makie.wong_colors()[6])
lines!(ax_insdiff,
       (las_xss.year .- 2000) * 1e-6,
       las_xss.insolation .- PT_ZB18a_interp(las_xss.year),
       label = "La04", linewidth = 1.5, linestyle = :dash,
       color = Makie.wong_colors()[2]
       )

# no La10b for this one
Legend(fsol[0,1], ax_ins, orientation = :horizontal, halign = :right)
rowgap!(fsol.layout, 5)
ylims!(-70, 70)

xlims!(-2.5, -2)

save("imgs/compare_solutions_insolation_detailed.pdf",
     fsol, size=(5.5inch, 5inch), px_per_unit = 300/inch)



# Figure 2: fecc
zbl = "ZB18a"
fecc = Figure()
ax_ecc1 = Axis(fecc[1,1])
ax_ecc2 = Axis(fecc[2,1])
ax_ecc3 = Axis(fecc[3,1])
ax_ecc4 = Axis(fecc[4,1], xlabel = "Time (Myr)")
Label(fecc[1:4,0], "Eccentricity (-)", rotation = pi/2)
for ax in [ax_ecc1, ax_ecc2, ax_ecc3, ax_ecc4]
    lines!(ax, PT_ZB18a.time .* 1e-3, PT_ZB18a.eccentricity,
           label = zbl, linewidth = 3,
           color = Makie.wong_colors()[1]
           )
    lines!(ax, La10b.time * 1e-3,
           La10b.eccentricity,
           color = Makie.wong_colors()[4],
           linestyle = :dash,
           label = "La10b")
    lines!(ax, las_xss.time * 1e-3, las_xss.eccentricity,
           label = "La04", #linestyle = :dash,
           # linewidth = 3,
           color = Makie.wong_colors()[2]
           )
end
Legend(fecc[0,1], ax_ecc1, orientation = :horizontal, halign = :right)
rowgap!(fecc.layout, 1, 5)
# linkxaxes!([ax_ecc1, ax_ecc2, ax_ecc3, ax_ecc4])
xlims!(ax_ecc1, -50, -45)
xlims!(ax_ecc2, -45, -40)
xlims!(ax_ecc3, -40, -35)
xlims!(ax_ecc4, -35, -30)

save("imgs/compare_eccentricity_-50_-30.pdf", fecc, size=(5.5inch, 5inch), px_per_unit = 300/inch)

fecc


# Figure 3: fcomp
update_theme!(
    joinstyle=:miter,
)
# contrast different insolation calculations for La04
# figure for paper?
fcomp = Figure()
ax = Axis(fcomp[1,1], xlabel = "Time (Myr)", ylabel = sil, title = "S₀ = $(global_S0)")
linesegments!(ax, times_las .- 0.050, ber_ins_R,
              # color = :darkred,
              color = :lightgreen,
              linecap = :round,
              label = "palinsol::ber78() |> palinsol::Insol()")
linesegments!(ax, (ber_df.year .- 2000) * 1e-3, ber_df.insolation,
              color = Makie.wong_colors()[4],
              linecap = :round,
              # color = :red,
              # linestyle = :dash, 
              label = "palinsol::ber78() |> PaleoInsolation.insolation()")
# lines!(ax, times, las_ins_65_90,
#        linewidth = 3,
#        color = :purple, label = "SNVec.insolation")
lines!(ax, las_xss6.time, las_xss6.insolation,
       color = Makie.wong_colors()[2], linewidth = 7, label = "insola La04 mean daily true longitude")
scatter!(ax, las_xss.time, las_xss.insolation,
              # color = :orange4,
         # color = (:orange4, 0),
         # strokecolor = :darkorange4, 
         color = :darkorange4, 
         markersize = 8,
         # strokewidth = 1,
         marker = :diamond,
              # linewidth = 0.5, 
              # linecap = :round,
              label = "insola La04 |> PaleoInsolation.insolation()")
# NOTE: palinsol standardizes on 1950
scatter!(ax, times_las .- 0.050, las_ins_65_90_R,
              # color = Makie.wong_colors()[2], 
              # color = :darkblue, 
         color = :darkorange3,
         # strokecolor = Makie.wong_colors()[2], 
         # color = Makie.wong_colors()[2], 
         markersize = 8,
         # strokewidth = 1,
         # marker = :utriangle,
         marker = :star4,
              # linewidth = 0.5, 
              # linecap = :round,
              # linestyle = :dash,
              label = "palinsol::la04() |> palinsol::Insol()")
scatterlines!(ax, PT_ZB18a.time, PT_ZB18a.insolation,
              color = Makie.wong_colors()[1],
              markercolor = :white,
              strokecolor = Makie.wong_colors()[1],
              markersize = 3,
              linewidth = 2,
              strokewidth = 1,
              label = "ZB18a(1,1) |> PaleoInsolation.insolation()")

# lines!(ax, las_xss.time, las_xss.insolation,
#        linewidth = 5,
#        color = (:purple, 0.4), label = "raw f90 ASCI + SNVec.insolation S₀ = $(global_S₀)")
# # almost the same, minor diffs
# lines!(ax, las_xss5.time, las_xss5.insolation,
       # color = :brown, linestyle = :dash, label = "insola true longitude S₀ = 1368.0")
# alsmost identical to above
# lines!(ax, @lift($xss.time), @lift($xss.insolation),
#        color = (:green, 0.6), linewidth = 5,
#        label = "ZB18a(1,1) SNVec.jl true longitude S₀ = $(global_S₀)")
axislegend(position=:lt)
# # different sub figure for different S0 value?
# ax2 = Axis(fcomp[2,1], xlabel = "Time (kyr)", ylabel = sil, title = "S₀ = 1365")
# # almost the same for young interval, out of phase for older interval
# lines!(ax2, zb_65.time, zb_65.insolation,
#        color = (:red, 0.6), linewidth = 5,
#        label = "ZB18a(1,1) astrocyclo webtool S₀ = 1365")
# # same pattern but different H0 value??
# # lines!(ax, PT_ZB18a.time, ins_1361, label = "ZB18a(1,1) SNVec.jl S₀ = 1361")
# lines!(ax2, PT_ZB18a.time, ins_1365, label = "ZB18a(1,1) PaleoInsolation.jl S₀ = 1365")
# lines!(ax, las_xss3.time, las_xss3.insolation,
#        color = :orange, linestyle = :dash, label = "insola mean longitude S₀ = 1368.0")
# # these have different phasing!!
# lines!(ax, pmip3.time, pmip3.insolation,
#        color = :black,
#        label = "PMIP3 Ber78 from 0 to 2100 AD")
# # higher res, falls slightly above
# DataInspector()
# xlims!(-127, 1)
xlims!(-71, 1)
ylims!(460, 550)

save("imgs/compare_insolation_calculation.pdf", fcomp,
     size=(6inch, 5inch), px_per_unit = 300/inch)

update_theme!(
    joinstyle=:round,
)





# Figure 4: modern insolation
xi = 1
x = @view PT_ZB18a[xi,:]

lats = deg2rad.(-90:90 |> collect)
lons = deg2rad.(0:359 |> collect)

ins = zeros(length(lons), length(lats))
for (j, lat) in enumerate(lats), (i, lon) in enumerate(lons)
    ins[i,j] = PaleoInsolation.insolation(
        x.eccentricity[xi],
        x.obliquity[xi],
        (x.lpx[xi] .- pi) .%  2pi,
        # longitude = day2l(x.eccentricity[xi], x.lpx[xi], lon), # if you want a "modern" calendar instead of angular calendar
        longitude = lon,
        latitude = lat,
        S0 = global_S0, H = nothing)
end

fmod = Figure()
ax_ins = Axis(fmod[1,1],
              xlabel = "True solar longitude (°)",
              # xlabel = "Mean solar longitude (°)",
              ylabel = "Earth’s latitude (°)")
cb1 = image!(ax_ins,
             extrema(lons) ./ pi .* 180,
             extrema(lats) ./ pi .* 180,
             ins,
             # ins2,
             interpolate = false,
             colormap = :solar, colorrange = (0, 600))
contour!(ax_ins,
         extrema(lons) ./ pi .* 180,
         extrema(lats) ./ pi .* 180,
         ins;
         labels = true, levels = 100:100:600, color = :white)
scatter!(ax_ins, [90], [65])
Colorbar(fmod[1,2], cb1, vertical = true, label = L"Insolation (Wm$^{-2}$)")

save("imgs/ZB18a_insolation_t0.png",fmod, size=(6inch, 3inch), px_per_unit = 300/inch)

# Figure 5: fmd maximum diffs
# these were generated in my other script in test_snvec.jl, where I have access
# to an interactive figure where I can tweak Ed and Td easily using my SNVec.jl
# WIP package. It loops over timesteps in 2 kyr increments over the past 50 Myr,
# computes insolation for the full grid for each solutions
# calculates differences in insolation
# takes extrema of these differences
# and puts these in the tables that I'm reading in here.
using Arrow
ber_diffs = DataFrame(Arrow.Table("out/ber_diffs2_50Myr_ins_2kyr.arrow"))
ber91_diffs = DataFrame(Arrow.Table("out/ber91_diffs_50Myr_ins_2kyr.arrow"))
las_diffs = DataFrame(Arrow.Table("out/las_diffs2_50Myr_ins_2kyr.arrow"))
las_diffs_phase = DataFrame(Arrow.Table("out/las_diffs_phase_50.arrow"))

# plot diffs
fmd = Figure()
ax_diff = Axis(fmd[1,1],
               xlabel = "Time (Myr)",
               ylabel = rich("ΔInsolation (Wm", superscript("−2"), ")"))
band!(ax_diff, ber_diffs.time * 1e-3, ber_diffs.lwr,  ber_diffs.upr, label = "Ber78", color = :red)
band!(ax_diff, ber91_diffs.time * 1e-3, ber91_diffs.lwr,  ber91_diffs.upr, label = "BL91", color = :darkred)
band!(ax_diff, las_diffs.time * 1e-3, las_diffs.lwr,  las_diffs.upr, label = "La04", color = :orange)
band!(ax_diff, las_diffs_phase.time * 1e-3, las_diffs_phase.lwr,  las_diffs_phase.upr, label = "La04 phase", color = :black)


# band!(ax_diff, lds.time * 1e-3, lds.lwr,  lds.upr, label = "La04", color = :orange)
# band!(ax_diff, las_diffs5.time * 1e-3, las_diffs5.lwr,  las_diffs5.upr, label = "La04", color = :orange)
# Legend(fmd[1,2], ax_diff, tellwidth = true, tellheight = false, unique = true, merge = true)
Legend(fmd[0,1], ax_diff, tellwidth = false, tellheight = true, unique = true, merge = true, halign = :right, orientation = :horizontal)
rowgap!(fmd.layout, 5)

# band!(ax_diff, las_diffs3.time * 1e-3, las_diffs3.lwr,  las_diffs3.upr, label = "La04 phase", color = :black)
# band!(ax_diff, las_diffs4.time * 1e-3, las_diffs4.lwr,  las_diffs4.upr, label = "La04 phase", color = :black)
# # DataInspector() # doesn't play nice with band
# Legend(fmd[1,2], ax_diff, tellwidth = true, tellheight = false, unique = true, merge = true)
save("imgs/extrema_insolation.png", fmd, size=(5.5inch, 2.5inch), px_per_unit = 300/inch)
# # deactivate_interaction!(ax_diff, :rectanglezoom)
# # make it possible to click and get coordinates
# # pt = select_point(ax_diff.scene)

save("imgs/extrema_insolation.pdf", fmd, size=(5.5inch, 2.5inch), px_per_unit = 300/inch)
save("imgs/extrema_insolation.png", fmd, size=(5.5inch, 2.5inch), px_per_unit = 300/inch)

# Figure 6: ins differences combined
# this uses my SNVec.jl interactive figure
fins
# to keep this document relatively clean, I'm not including it here
# but it's available at https://github.com/japhir/InsolationExplorer.jl
# after realizing I shouldn't interpolate
set_close_to!(sg.sliders[3], -934)
d1 = copy(ber_insdiff[])
set_close_to!(sg.sliders[3], -2108)
d1_2 = copy(ber_insdiff[])
set_close_to!(sg.sliders[3], -5038)
d2 = copy(ber_insdiff[])
set_close_to!(sg.sliders[3], -9892)
xlims!(ax_lpx, [-10200, -9500])
d3 = copy(las_insdiff[])
set_close_to!(sg.sliders[3], -14824)
xlims!(ax_lpx, [-15200, -14500])
d4 = copy(las_insdiff[])
set_close_to!(sg.sliders[2], 1.146564088532723) # 0--20 La04 65 ins 2 kyr
set_close_to!(sg.sliders[3], -9892)
xlims!(ax_lpx, [-10200, -9500])
d3f = copy(las_insdiff[])
extrema(d3f)
set_close_to!(sg.sliders[3], -14824)
xlims!(ax_lpx, [-15200, -14500])
d4f = copy(las_insdiff[])
extrema(d4f)

# after figuring out that with a slightly different td value La04 is fully aligned
set_close_to!(sg.sliders[3], -40000)
xlims!(ax_lpx, [-40350, -39650])
extrema(las_insdiff[])
set_close_to!(sg.sliders[3], -30000)
xlims!(ax_lpx, [-30350, -29650])
set_close_to!(sg.sliders[3], -3206)


# new insdiff figure
finsdiff = Figure()
ax1 = Axis(finsdiff[1,1], title = "Ber78 − ZB18a\n−934 kyr", aspect = DataAspect())
cb1 = image!(ax1,
             extrema(lons) ./ pi .* 180,
             extrema(lats) ./ pi .* 180,
             d1,
             interpolate = false,
             colormap = :diff, colorrange = (-100, 100))
contour!(ax1,
         extrema(lons) ./ pi .* 180,
         extrema(lats) ./ pi .* 180,
         d1;
         labels = true,
         levels = [-20, -15, -10, -5, 5, 10, 15, 20, 25],
         color = :black)
ax1_2 = Axis(finsdiff[2,1], title = "Ber78 − ZB18a\n−2108 kyr", aspect = DataAspect())
cb1_2 = image!(ax1_2,
             extrema(lons) ./ pi .* 180,
             extrema(lats) ./ pi .* 180,
             d1_2,
             interpolate = false,
             colormap = :diff, colorrange = (-100, 100))
contour!(ax1_2,
         extrema(lons) ./ pi .* 180,
         extrema(lats) ./ pi .* 180,
         d1_2;
         labels = true,
         levels = [-50, -40, -30, -20, -10, 10, 20, 30, 40, 50, 60],
         color = :black)
ax2 = Axis(finsdiff[3,1], title = "Ber78 − ZB18a\n−5038 kyr", aspect = DataAspect())
cb2 = image!(ax2,
             extrema(lons) ./ pi .* 180,
             extrema(lats) ./ pi .* 180,
             d2,
             interpolate = false,
             colormap = :diff, colorrange = (-100, 100))
contour!(ax2,
         extrema(lons) ./ pi .* 180,
         extrema(lats) ./ pi .* 180,
         d2;
         labels = true,
         levels = -600:20:600,
         color = :white)

Colorbar(finsdiff[1,2], cb1, vertical = true, label = L"$\Delta$Insolation (Wm$^{-2}$)")
Label(finsdiff[4,1], "True solar longitude (°)", tellwidth = false, tellheight = true)
Label(finsdiff[1:3,0], "Earth’s latitude (°)", tellheight = false, tellwidth = true,
      rotation = 0.5pi)
hidexdecorations!(finsdiff.content[1], ticks = false)
hidexdecorations!(finsdiff.content[2], ticks = false)
# hidexdecorations!(finsdiff.content[3], ticks = false)
# hidexdecorations!(finsdiff.content[4], ticks = false)

save("imgs/insolation_differences.png", finsdiff,
     size=(4inch, 6inch), px_per_unit = 300/inch)

# which I then combined with panels of eccentricity/cp/obliquity in Inkscape



# Figure A: f obliquity
f, ax = lines(1e-3 * PT_ZB18a.time,
              rad2deg.(PT_ZB18a.obliquity),
              label = "ZB18a(1,1)",
              axis = (; xlabel = "Time (Myr)",
                      ylabel = "Obliquity (°)"))
lines!(ax,
       1e-3 * PT_ZB20a.time,
       rad2deg.(PT_ZB20a.obliquity),
       # alpha = 0.8,
       label = "ZB20a(1,1)")
xlims!(ax, (-72,0))
vlines!(ax, [-71,-58], color = :black, linestyle = :dash)
Legend(f[0,1], ax, tellwidth = false, tellheight = true,
       halign = :right, orientation = :horizontal)
rowgap!(f.layout, 6)

save("imgs/obliquity.png", f,
     size=(4.5inch, 2inch), px_per_unit = 300/inch)
save("imgs/obliquity.pdf", f,
     size=(4.5inch, 2inch), px_per_unit = 300/inch)



# Figure A: fobl
# obliquity ranges
obl = DataFrame(
category_labels = [
    # fill("ZB18a(1,1)\npast 100 Myr", nrow(x1_1));
    # fill("ZB18a(1,1)\n−60 to −50 Myr", length(obl_Eocene));
    # fill("ZB18a(1,1)\n−20 to −10 Myr", length(obl_Miocene));
    "Vahlenkamp et al., 2018\nmiddle Eocene (43–46 Ma)",
    "Lunt et al., 2011\nKeery et al., 2018\nTierney et al., 2022\nPETM/Eocene",
    "Goldner et al., 2014\nMMCO (17–14.5 Ma)"
],
low = [
    # x1_1.obliquity * 180/pi;
    # obl_Eocene;
    # obl_Miocene;
    22.1, # vahlenkamp
    22, # Lunt
    22, # Goldner 2014a
    # keery.obliquity;
],
high = [
    24.5, # vahlenkamp
    24.5, # Lunt
    25, # Goldner 2014a
])
obl.mean = obl.low .+ (obl.high .- obl.low) /2
obl.timeperiod = [:midEocene, :PETMEocene, :MMCO]

# not really Eocene/Miocene but older and younger
obl_Eocene = PT_ZB20a[PT_ZB20a.time .>= -60e3 .&& PT_ZB20a.time .<= -50e3, :obliquity] * 180 / pi
obl_Miocene = PT_ZB18a[PT_ZB18a.time .>= -20e3 .&& PT_ZB18a.time .<= -10e3, :obliquity] * 180 / pi
obl_modern = PT_ZB20a[PT_ZB20a.time .>= -10e3 , :obliquity] * 180 / pi
extrema(obl_modern)

# override with actual ranges for direct comparison
# full time period
# obl_Eocene = PT_ZB20a[PT_ZB20a.time .>= -56e3 .&& PT_ZB20a.time .<= -33.9e3, :obliquity] * 180 / pi
# obl_Miocene = PT_ZB18a[PT_ZB18a.time .>= -23.03e3 .&& PT_ZB18a.time .<= -5.333e3, :obliquity] * 180 / pi
# time period for respective studies
obl_Eocene = PT_ZB20a[PT_ZB20a.time .>= -46e3 .&& PT_ZB20a.time .<= -43e3, :obliquity] * 180 / pi
obl_PETM = PT_ZB20a[PT_ZB20a.time .>= -60e3 .&& PT_ZB20a.time .<= -52e3, :obliquity] * 180 / pi
obl_Miocene = PT_ZB18a[PT_ZB18a.time .>= -17e3 .&& PT_ZB18a.time .<= -14.5e3, :obliquity] * 180 / pi

# eae1, eae2 = extrema(PT_ZB18a[PT_ZB18a.time .>= -58e3, :obliquity] * 180 / pi)
# tae1, tae2 = extrema(PT_ZB20a[PT_ZB20a.time .>= -71e3, :obliquity] * 180 / pi)
Eae1, Eae2 = extrema(obl_Eocene)
Pae1, Pae2 = extrema(obl_PETM)
Mae1, Mae2 = extrema(obl_Miocene)

cE = colorant"rgba(253,180,108,1)"
cP = colorant"rgba(253,100,50,1)"
cM = colorant"rgba(220,220,0,1)"

pl = data(obl) *
    mapping(:category_labels => "",
            :mean, :low, :high,
            color = :timeperiod => "",
            ) *
                visual(CrossBar; orientation = :horizontal,
                       show_midline = false, width = 0.2)
f = draw(pl,
         # order of colors is expermentally determined, kinda dumb
         scales(;Color = (; palette = [cM, cP, cE]));
         legend = (;show = false))
Label(f.figure[2,1], "Obliquity (°)",
      tellheight = true, tellwidth = false)
# rangebars!(f.figure[1,1], [5], eae1, eae2,
#            whiskerwidth = 10, direction = :x,
#            label = "ZB18a(1,1)\n−58 to 0 Myr",
#            color = :black)
# rangebars!(f.figure[1,1], [5.02], tae1, tae2,
#            whiskerwidth = 10, direction = :x,
#            label = "ZB20a(1,1)\n−71 to 0 Myr",
#            color = :black)
rangebars!(f.figure[1,1], [1.15], Mae1, Mae2,
           whiskerwidth = 10, direction = :x,
           label = "ZB18a(1,1)\n−17 to −14.5 Myr",
           color = cM)
rangebars!(f.figure[1,1], [2.15], Pae1, Pae2,
           whiskerwidth = 10, direction = :x,
           label = "ZB20a(1,1)\n−60 to −52 Myr",
           color = cP)
rangebars!(f.figure[1,1], [3.15], Eae1, Eae2,
           whiskerwidth = 10, direction = :x,
           label = "ZB20a(1,1)\n−46 to −43 Myr",
           color = cE)
text!(f.figure[1,1], 22.7, 1.2; text = "ZB18a(1,1)\n−17 to −14.5 Myr")
text!(f.figure[1,1], 22.7, 2.2; text = "ZB20a(1,1)\n−60 to −52 Myr")
text!(f.figure[1,1], 22.7, 3.2; text = "ZB20a(1,1)\n−46 to −43 Myr")
# axislegend(position = :rt)
# legend!(f.figure[1,2], f)
# Legend(f.figure[0,1], f.figure.content[1], "",
#        tellheight = true, tellwidth = false)
ylims!(0.3, 4.3)

save("imgs/obliquity_ranges.pdf", f,
     size=(4.5inch, 3inch), px_per_unit = 300/inch)


# Figure A: fsol extension
fsol = Figure()
zbl = "ZB18a(1,1)"
ax_ins1 = Axis(fsol[1,1])
ax_ins2 = Axis(fsol[2,1])
ax_ins3 = Axis(fsol[3,1])
ax_ins4 = Axis(fsol[4,1], xlabel = "Time (Myr)")
Label(fsol[1:4,0], sil, rotation = pi/2)
lns = Vector{Lines}(undef, 4)
for ax in [ax_ins1, ax_ins2, ax_ins3, ax_ins4]
    lns[1] = lines!(ax, ber91_df.time * 1e-3, ber91_df.insolation,
           label = "BL91",
           color = Makie.wong_colors()[6]
           )
    lns[2] = lines!(ax, ber_df.time * 1e-3, ber_df.insolation,
           label = "Ber79",
           # linestyle = :dash,
           color = Makie.wong_colors()[4]
           )
    lns[3] = lines!(ax, PT_ZB18a.time .* 1e-3, PT_ZB18a.insolation,
           label = zbl, linewidth = 3,
           color = Makie.wong_colors()[1]
           )
    lns[4] = lines!(ax, las_xss.time * 1e-3, las_xss.insolation,
           label = "La04",
           linestyle = :dash,
           color = Makie.wong_colors()[2]
           )
end
bl, be, zb, ls = lns
Legend(fsol[0,1],
       [
           be,
           bl,
           ls,
           zb
       ],
       ["Ber78",
        "BL91",
        "La04",
        "ZB18a(1,1)"],
       orientation = :horizontal,
       halign = :right)
rowgap!(fsol.layout, 1, 5)
# linkxaxes!([ax_ins1, ax_ins2, ax_ins3, ax_ins4])
xlims!(ax_ins1, -2, -1.5)
xlims!(ax_ins2, -1.5, -1)
xlims!(ax_ins3, -1, -0.5)
xlims!(ax_ins4, -0.5, 0)

save("imgs/compare_solutions_past2Myr.pdf", fsol,
     size=(5.5inch, 5.5inch), px_per_unit = 300/inch)
