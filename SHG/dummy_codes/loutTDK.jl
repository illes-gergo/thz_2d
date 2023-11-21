lout_general = Layout(xaxis=attr(mirror="all", ticks="inside", linecolor=:black, showgrid=false, zeroline=false), yaxis=attr(mirror="all", ticks="inside", linecolor=:black, showgrid=false, zeroline=false), font=attr(size=14, color="Black"), plot_bgcolor=:white)

lout_thz = Layout(xaxis=attr(range=[-3, 3], title="Time (ps)"), yaxis=attr(range=[-7, 10], title="Transversal pos. (mm)"))

lout_thz_spectr = Layout(xaxis=attr(range=[0, 1], title="Frequency (THz)"), yaxis=attr(range=[-7,10], title="Transversal pos. (mm)"))

lout_thz1D = Layout(xaxis=attr(range=[-5, 5], title="Time (ps)"), yaxis=attr(title="Electric field (kV/cm)"), legend=attr(x=0.8, y=1))

lout_thz1D_spectr = Layout(xaxis=attr(range=[0, 5], title="Frequency (THz)"), yaxis=attr(title="Spectral Amplitude (arb. u.)",range=[0,14e-6]), legend=attr(x=0.8, y=1))

lout_imax=Layout(xaxis=attr(range=[0,2],title="Kristályhossz (mm)"),yaxis=attr(title=L"Intenzitás (GW/cm2$)",range=[70,105]))

lout_effic = Layout(yaxis=attr(title="Hatásfok (%)", range=[0,2.5]), xaxis=attr(title="Kristályhossz (mm)",range=[0,8]))