lout_general = Layout(xaxis=attr(mirror="all", ticks="inside", linecolor=:black, showgrid=false, zeroline=false), yaxis=attr(mirror="all", ticks="inside", linecolor=:black, showgrid=false, zeroline=false), font=attr(size=14, color="Black"), plot_bgcolor=:white)

lout_thz = Layout(xaxis=attr(range=[-2, 2], title="Idő (ps)"), yaxis=attr(range=[-2, 4] .- 0.0, title="Transzverzális koord. (mm)"))

lout_dist = Layout(xaxis=attr(range=[0, 8], title="Kristályhossz (mm)"), yaxis=attr(range=[-1, 5] .- 0.0, title="Transzverzális koord. (mm)"))

lout_dist_spect = Layout(xaxis=attr(range=[0, 8], title="Kristályhossz (mm)"), yaxis=attr(range=[-1, 5] .- 0.0, title="Transzverzális koord. (mm)"))

lout_thz_spectr = Layout(xaxis=attr(range=[0, 5], title="Frekvencia (THz)"), yaxis=attr(range=[-2, 4] .- 0.0, title="Transzverzális koord. (mm)"))

lout_thz1D = Layout(xaxis=attr(range=[-2, 2], title="Idő (ps)"), yaxis=attr(title="Elektromos térerősség (kV/cm)"), legend=attr(x=0.7, y=0.95))

lout_int1D = Layout(xaxis=attr(range=[-2,5],title="Transzverzális koordináta (mm)"), yaxis=attr(title="Energiasűrűség (rel. egys.)",range=[0,1.05]))

lout_thz1D_spectr = Layout(xaxis=attr(range=[0, 7], title="Frekvencia (THz)"), yaxis=attr(title="Spektrális amplitúdó (rel. egys.)", range=[0, 1]), legend=attr(x=0.7, y=0.95))

lout_imax = Layout(xaxis=attr(range=[0, 2], title="Kristályhossz (mm)"), yaxis=attr(title=L"Intenzitás (GW/cm2$)", range=[70, 105]))

lout_effic = Layout(yaxis=attr(title="Hatásfok (%)", range=[0, 4.5]), xaxis=attr(title="Kristályhossz (mm)", range=[0, 8]))

colorbar_field = attr(
    title="Elektromos tér (kV/cm)", # title here
    titleside="right",
    titlefont=attr(
        size=18,
        family="sans-serif",
        color=:black
    ),
    tickvals=range(-600, 600, 9)
)

colorbar_spectr = attr(
    title="Spektrális amplitúdó (rel. egys)", # title here
    titleside="right",
    titlefont=attr(
        size=18,
        family="sans-serif",
        color=:black
    ),
    tickvals=range(0, 1, 9)
)

colorbar_dist = attr(
    title="Intenzitás sűrűség (rel. egys)", # title here
    titleside="right",
    titlefont=attr(
        size=18,
        family="sans-serif",
        color=:black
    ),
    tickvals=range(0, 1, 9)
)
colorbar_dist_spect = attr(
    title="Központi frekvencia (THz)", # title here
    titleside="right",
    titlefont=attr(
        size=18,
        family="sans-serif",
        color=:black
    ),
    tickvals=range(0, 3, 5)
)
