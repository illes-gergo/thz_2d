l_generic = Layout(xaxis_linecolor=:black, yaxis_linecolor=:black,xaxis_mirror="all",yaxis_mirror="all",xaxis_ticks="inside",yaxis_ticks="inside", plot_bgcolor=:white, xaxis_showgrid=false,yaxis_showgrid=false, xaxis_zeroline=false, yaxis_zeroline=false,font=attr(size=14,color="Black"),showlegend=false)

l_efficsh = Layout(xaxis=attr(title="Kristályhossz (mm)"),yaxis=attr(title="Másodharmonikus keltés hatásfoka (%)"))

l_heatmap = Layout(coloraxis = attr(cmin=0,cmax=1))

l_compare2 = Layout(xaxis=attr(range=[-3e-12,3e-12],title="Idő (s)"),yaxis=attr(min=0))
