## libraries
using CSV
using DataFrames
using GLMakie
# using CairoMakie
using Downloads
using RCall

# "package" with insolation and readbindata functions
using PaleoInsolation
# hacky way to include insolation function
# includet("src/insolation.jl")
# includet("src/readbindata.jl")

inch = 96
pt = 4/3
cm = inch / 2.54

update_theme!(
    GLMakie.Theme(Axis = (xgridvisible = false, ygridvisible = false)),
    fontsize=12,
    joinstyle=:round,
)

## load data
# HNBody c: input ZB18a-plan3.dat
ZB18a = CSV.read("dat/ZB18a-plan3.dat",
                 DataFrame,
                 delim = ' ',
                 stripwhitespace = true,
                 ignorerepeated = true,
                 skipto = 18,
                 header = [:time, #1
                           :semimajor, #7
                           :eccentricity, #8
                           :inclination, #9
                           :long_periapse, #12
                           :long_ascnode, #10
                           :arg_periapse, #11
                           :mean_anomaly, #15
                           ])

# HNBody c: input ZB20a-plan3.dat
ZB20a = CSV.read("dat/ZB20a-plan3.dat",
                 DataFrame,
                 delim = ' ',
                 stripwhitespace = true,
                 ignorerepeated = true,
                 skipto = 18,
                 header = [:time, #1
                           :semimajor, #7
                           :eccentricity, #8
                           :inclination, #9
                           :long_periapse, #12
                           :long_ascnode, #10
                           :arg_periapse, #11
                           :mean_anomaly, #15
                           ])

# snvec  c: computed PT-ZB18a_1-1.bin
ZB18a_bin = readbindata()

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
                              # :precession_unwrapped,
                              # :lph,
                              :lpx,
                              :climatic_precession]
                    )

# snvec  c: computed PT-ZB20a_1-1.bin
ZB20a_bin = readbindata("dat/PT-ZB20a_1-1.bin")

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


# astrocycle web tool pre-computed reference insolation
ZB18a_insolation_ref = CSV.read("dat/ZB18a_insolation_ref.dat",
                                DataFrame,
                                skipto = 2,
                                delim = ' ',
                                stripwhitespace = true,
                                ignorerepeated = true,
                                header = [:time, :insolation],
                                )

# paleoinsolation Fortran: 65deg N summer insolation dat/ZB18a_insolation.dat
ZB18a_insolation = CSV.read("dat/ZB18a_insolation.dat",
                            DataFrame;
                            delim = ' ',
                            stripwhitespace = true,
                            ignorerepeated = true,
                            header = [:time,
                                      :eccentricity,
                                      :obliquity,
                                      :precession,
                                      :lpx,
                                      :climatic_precession,
                                      :insolation]
                            )


# interpolated output from Fortran
# vcat(
    # CSV.read("interp_1620.dat",
    #                     DataFrame;
    #                     delim = ' ',
    #                     stripwhitespace = true,
    #                     ignorerepeated = true,
    #                     header = [:yearAD,
    #                               :eccentricity,
    #                               :obliquity,
    #                               :lpx]
    #                         ),
    # CSV.read("interp_-9849.dat",
    #                     DataFrame;
    #                     delim = ' ',
    #                     stripwhitespace = true,
    #                     ignorerepeated = true,
    #                     header = [:yearAD,
    #                               :eccentricity,
    #                               :obliquity,
    #                               :lpx]
    #          ),
# )
ZB18a_interp = CSV.read("interp_-40000.dat",
                        DataFrame;
                        delim = ' ',
                        stripwhitespace = true,
                        ignorerepeated = true,
                        header = [:yearAD,
                                  :eccentricity,
                                  :obliquity,
                                  :lpx,
                                  :insolation]
                            )
ZB18a_interp.time = 1e-3 * (ZB18a_interp.yearAD .- 2000)
ZB18a_interp.insol_65 = insolation.(
    ZB18a_interp.eccentricity,
    ZB18a_interp.obliquity,
    mod.(ZB18a_interp.lpx .- pi,2pi);
    S0 = 1365.0
)

ZB18a_interp_eccmax = CSV.read("interp_-250_-170.dat",
                        DataFrame;
                        delim = ' ',
                        stripwhitespace = true,
                        ignorerepeated = true,
                        header = [:yearAD,
                                  :eccentricity,
                                  :obliquity,
                                  :lpx,
                                  :insolation]
                            )
ZB18a_interp_eccmax.time = 1e-3 * (ZB18a_interp_eccmax.yearAD .- 2000)



### plotting

# eccentricity
p = Figure()
ax = Axis(p[1,1], xlabel = "Time (kyr)", ylabel = "Eccentricity (-)")
lines!(ax, ZB18a.time / 365.25e3, ZB18a.eccentricity,
       label = "HNBody output")
lines!(ax, ZB18a_bin.time, ZB18a_bin.eccentricity, label = "snvec.c output")
axislegend()
xlims!(-250,-170)
ylims!(0.035, 0.051)
save("imgs/ecc_jitters.png",p)
# lines!(ZB20a.time / 365.25e3, ZB20a.eccentricity)

# HNBody c output
lines(ZB18a.time[begin:100] / 365.25e3, ZB18a.eccentricity[begin:100])
lines!(ZB20a.time[begin:100] / 365.25e3, ZB20a.eccentricity[begin:100])
# add PT processed output
lines!(PT_ZB18a.time[begin:100], PT_ZB18a.eccentricity[begin:100])
lines!(PT_ZB20a.time[begin:100], PT_ZB20a.eccentricity[begin:100])
# add insolation from web-tool
lines!(ZB18a_insolation.time[begin:100], ZB18a_insolation.eccentricity[begin:100])
# add Fortran linterp
scatterlines!(ZB18a_interp.time, ZB18a_interp.eccentricity)
# xlims!(-0.391,-0.37)
# ylims!(0.01684,0.016854)

# inclination
# HNBody c output
lines(ZB18a.time[begin:100] / 365.25e3, ZB18a.inclination[begin:100])
lines!(ZB20a.time[begin:100] / 365.25e3, ZB20a.inclination[begin:100])
# seems smooth enough
# full time interval
lines(ZB18a.time / 365.25e3, ZB18a.inclination)
lines!(ZB20a.time / 365.25e3, ZB20a.inclination)

# longitude of perihelion
lines(ZB18a.time[begin:100] / 365.25e3, ZB18a.long_periapse[begin:100])

lines(ZB18a.time / 365.25e3, ZB18a.long_periapse)
lines!(ZB20a.time / 365.25e3, ZB20a.long_periapse)
# also very subtly jittery

# precession
# snvec  c output
f, ax = scatterlines(
    PT_ZB18a.time[begin:100],
    PT_ZB18a.climatic_precession[begin:100],
      label = "ZB18a(1,1)",
      axis = (;xlabel = "Time (kyr)",
              ylabel = "Climatic precession (-)"))
# lines!(ax, PT_ZB20a.time, PT_ZB20a.climatic_precession, label = "ZB20a(1,1)")
# scatter!(ax, PT_ZB20a.time, PT_ZB20a.climatic_precession, label = "ZB20a(1,1)")
scatterlines!(ax,
       ZB18a_interp.time,
       # ZB18a_interp.eccentricity .* sin.(ZB18a_interp.lpx),
       ZB18a_interp.eccentricity .* sin.(mod.(ZB18a_interp.lpx, 2pi)),
         label = "ZB18a(1,1) interp",
         color = :darkred)
# xlims!(ax, -0.391,-0.37)
# ylims!(ax, 0.01672, 0.01676)

save("imgs/interpolate.pdf", f,
     size=(4.5inch, 2inch), px_per_unit = 300/inch)
# SHIT wth is happening here?



xlims!(ax, -3,0)
ylims!(ax, 0.014, 0.0174)
axislegend(ax, merge=true, position = :rb)

save("imgs/need_to_interpolate?.png", f,
     size=(4.5inch, 2inch), px_per_unit = 300/inch)

# lines!(PT_ZB18a.time, PT_ZB18a.eccentricity)
# lines!(PT_ZB20a.time, PT_ZB20a.eccentricity)

# plot obliquity for paper
f, ax = scatterlines(PT_ZB18a.time[begin:100],
                     rad2deg.(PT_ZB18a.obliquity[begin:100]),
              label = "ZB18a(1,1)",
              axis = (; #xlabel = "Time (Myr)",
                      xlabel = "Time (kyr)",
                      ylabel = "Obliquity (°)"))
scatterlines!(ax, # 1e-3*
              PT_ZB20a.time[begin:100],
              rad2deg.(PT_ZB20a.obliquity[begin:100]),
       alpha = 0.8,
       label = "ZB20a(1,1)")

scatterlines!(ax,
         ZB18a_interp.time,
         rad2deg.(ZB18a_interp.obliquity),
         label = "ZB18a(1,1) interp")
# xlims!(ax, -0.391,-0.37)
# ylims!(ax, 23.4871, 23.4905)

# plot lpx for paper
# scatterlines(PT_ZB18a.time, PT_ZB18a.lph)
# scatterlines!(PT_ZB18a.time, PT_ZB18a.lphu)
scatterlines(PT_ZB18a.time, PT_ZB18a.precession)
# scatterlines(PT_ZB18a.time, PT_ZB18a.precession_unwrapped)
scatterlines(PT_ZB18a.time, PT_ZB18a.lpx)
scatterlines(PT_ZB18a.time,
             mod.(PT_ZB18a.lpx,2pi))
scatterlines(PT_ZB18a.time, PT_ZB18a.climatic_precession)
# scatterlines(PT_ZB20a.time, PT_ZB20a.climatic_precession)
# check
lines!(ZB18a_insolation.time, ZB18a_insolation.climatic_precession,
       color = :yellow)
# scatterlines(PT_ZB18a.time, PT_ZB18a.obliquity)

f, ax = scatterlines(
    PT_ZB18a.time,
    PT_ZB18a.lpx,
    label = "PT-ZB18a(1,1) raw")
scatterlines!(
    ax,
    PT_ZB18a.time,
    mod.(PT_ZB18a.lpx,2pi),
    label = "PT-ZB18a(1,1) mod")
scatterlines!(
    ax,
    ZB18a_interp.time,
    ZB18a_interp.lpx,
    label = "ZB18a(1,1) interp raw")
scatterlines!(
    ax,
    ZB18a_interp.time,
    mod.(ZB18a_interp.lpx,2pi),
    label = "ZB18a(1,1) interp mod")
# scatterlines(ZB18a.time,
#              ZB18a.long_periapse,
#              label = "ZB18a(1,1) interp")
xlims!(-40,0)
ylims!(-pi,2.1pi)
axislegend()

f, ax = scatterlines(PT_ZB18a.time[begin:100] # * 1e-3
              ,
              mod.(rad2deg.(PT_ZB18a.lpx[begin:100]), 360),
              label = "ZB18a(1,1)",
              axis = (; #xlabel = "Time (Myr)",
                      xlabel = "Time (kyr)",
                      ylabel = "LPX (°)"))
scatterlines!(ax, # 1e-3*
       PT_ZB20a.time[begin:100],
       mod.(rad2deg.(PT_ZB20a.lpx[begin:100]), 360),
       alpha = 0.8,
       label = "ZB20a(1,1)")

scatterlines!(ax,
         ZB18a_interp.time,
         mod.(rad2deg.(ZB18a_interp.lpx), 360),
         label = "ZB18a(1,1) interp")
axislegend()



# plot insolation for paper
f, ax = scatterlines(ZB18a_insolation_ref.time[begin:100],
              ZB18a_insolation_ref.insolation[begin:100],
              label = "ZB18a(1,1) web",
              axis = (; #xlabel = "Time (Myr)",
                      xlabel = "Time (kyr)",
                      ylabel = "65°N summer insolation (Wm¯²)"))
scatterlines!(ax,
              ZB18a_insolation.time,
              ZB18a_insolation.insolation,
              label = "ZB18a(1,1) Fortran")
scatterlines!(ax,
         ZB18a_interp.time,
         ZB18a_interp.insolation,
         color = :darkred,
         label = "ZB18a(1,1) Fortran interp")
scatterlines!(ax,
         ZB18a_interp_eccmax.time,
         ZB18a_interp_eccmax.insolation,
         color = :darkred,
         label = "ZB18a(1,1) Fortran interp")
# scatterlines!(ax,
#          ZB18a_interp.time,
#          ZB18a_interp.insol_65,
#          color = :darkred,
#          label = "ZB18a(1,1) interp")
axislegend(merge=true)
xlims!(-254,0)
save("imgs/interpolation_ok.png",f)
# OK so it does differ very slightly!
# let's try to interpolate the whole interval

# vlines!(ax, [-58, -71], color = Makie.wong_colors()[1:2])
# hlines!(ax, collect(rad2deg.(extrema(PT_ZB18a.obliquity))),
        # color = fill(Makie.wong_colors()[1],2))
# hlines!(ax, collect(rad2deg.(extrema(PT_ZB18a.obliquity[PT_ZB18a.time .>= -58e3]))),
        # color = fill(Makie.wong_colors()[1],2))

xlims!(ax, (-72,0))
Legend(f[0,1], ax, tellwidth = false, tellheight = true,
       halign = :right, orientation = :horizontal)
rowgap!(f.layout, 6)

save("imgs/obliquity.png", f,
     size=(4.5inch, 2inch), px_per_unit = 300/inch)

xlims!(ax, -11, -8)
ylims!(ax, 24.2,24.238)
save("imgs/need_to_interpolate_obliquity?.png", f,
     size=(4.5inch, 2inch), px_per_unit = 300/inch)


lines(PT_ZB18a.time, PT_ZB18a.lpx)
lines!(PT_ZB20a.time, PT_ZB20a.lpx)

lines(PT_ZB18a.time, mod.(PT_ZB18a.lpx .- pi,2pi))
lines!(PT_ZB20a.time, mod.(PT_ZB20a.lpx .- pi,2pi))

# paleoinsolation Fortran output
# note that fortran already wraps lpx again!
lines(ZB18a_insolation.time, ZB18a_insolation.lpx)
lines!(PT_ZB18a_interp.time, mod.(PT_ZB18a.lpx .- pi,2pi))

lines(ZB18a_insolation.time, ZB18a_insolation.insolation)
lines(ZB18a_interp.time, ZB18a_interp.insol_65)
lines!(ZB18a_insolation_ref.time, ZB18a_insolation_ref.insolation)


## tripple-check modern value

# Joussaume & Bracconot 1997
# based on either Berger 1978 or Laskar 1993
# 1950: ω - 180° = 102.04°

# Otto-Bliesner et al., 2017 1850
# based on Berger & Loutre 91
# 1850: ω - 180° = 100.3°

# Direct Berger & Loutre source:
# https://doi.pangaea.de/10.1594/PANGAEA.56040

# ZB18a(1,1)
# 2000: ω =
rad2deg(PT_ZB18a.lpx[1]) # 102.9°
# 2000: ω - 180° =
rad2deg(mod(PT_ZB18a.lpx[1]-pi,2pi)) # 282.917
rad2deg(ZB18a_insolation.lpx[1]) # this is wrapped and -180! 282.917

# vs palinsol in R
# Berger 1978
# 1950: ω - 180° =
rcopy(R"palinsol::ber78(0)[3] * 180 / pi")
# 282.039
# 1950: ω - 180° =
rcopy(R"palinsol::la04(0)[3] * 180 / pi")
# 282.0633
# 1950: ω - 180° =
rcopy(R"palinsol::ber90(0)[3] * 180 / pi")
# 281.3828
# we know palinsol preps La04 by already subtracting pi
# https://github.com/mcrucifix/palinsol/blob/master/data-raw/LA04.R#L50

# vs. La04 raw Fortran output
# http://vo.imcce.fr/insola/earth/online/earth/La2004/INSOLN.LA2004.BTL.100.ASC
La04 = CSV.read(Downloads.download("http://vo.imcce.fr/insola/earth/online/earth/La2004/INSOLN.LA2004.BTL.100.ASC"),
                DataFrame,
               delim = ' ',
               stripwhitespace = true,
               ignorerepeated = true,
               header = ["time", "eccentricity", "obliquity", "pibar"]
               )
# 2000 ω =
rad2deg(La04.pibar[1])
# 102.9179
# 2000 ω - 180 =
rad2deg(mod(La04.pibar[1] - pi, 2pi))
# 282.9179


# piControl = 1850
rcopy(R"palinsol::ber90(-100, degree = TRUE)['varpi']")
#> 281.381
R".Last.value - 180"
#> 101.3811
# midHolocene = 6 ka
rcopy(R"palinsol::ber90(-6000, degree = TRUE)['varpi]")
#> 281.2798
R".Last.value - 180"
#> 101.2798
# lig127k = 127 ka
rcopy(R"palinsol::ber90(-127000, degree = TRUE)['varpi']")
#> 279.20363
R".Last.value - 180"
#> 99.20363
