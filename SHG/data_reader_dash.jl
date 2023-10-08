using HDF5, Dash, PlotlyJS, DelimitedFiles, DashBootstrapComponents, LazyGrids, Base.Threads
include("data_reader_dash_funcs.jl")

calcdir = "/home/illesg/cst/2d-calculations/"


app = dash(external_stylesheets=[dbc_themes.GRID])


p1s = (scatter(x=sort(rand(10)),y=rand(10)))
p2s = (scatter(x=sort(rand(10)),y=rand(10)))
p3s = (scatter(x=sort(rand(10)),y=rand(10)))
p4s = (scatter(x=sort(rand(10)),y=rand(10)))

p1 = heatmap(x=(rand(10,10)),y=rand(10,10),z=rand(10,10))
p2 = heatmap(x=(rand(10,10)),y=rand(10,10),z=rand(10,10))
p3 = heatmap(x=(rand(10,10)),y=rand(10,10),z=rand(10,10))
p4 = heatmap(x=(rand(10,10)),y=rand(10,10),z=rand(10,10))



tupofcalcs = createtup(calcdir)

app.layout = html_div((
        html_h1("Vezérlés"),
        html_h2("Számolás kiválasztása"),
        dbc_row((
            dbc_col((
                    dcc_radioitems(
                    id="calcselector",
                    options=tupofcalcs,
                    value=tupofcalcs[1].value,
                    labelStyle=Dict("display" => "inline-block")
                )
                ), width=12)
           
        )),
	html_h2("Kristályhossz választása"),
	dbc_row((
	
	dbc_col((

		 dcc_slider(
			id="crystalslider",
			min=1,
			step=1,
			max=401,
		       marks=Dict([i=>"$(i)" for i in 1:50:401]),
		       value=101
			    )

		 ),width=10),
	dbc_col((
		 html_button("2D ábrák frissítése",id="crystalupdate")
		 ))

		 )),
         html_h1("Összefoglaló értékek"),
         dbc_row((
            html_div(
                "Maximális Intenzitást:"*"20GW"
            )
         )),
        html_h1("2D ábrák"),
        dbc_row((
            dbc_col((
                html_h2("THz-es alak"),
                dcc_graph(
                    id="thz_t",
                    #figure=Plot(heatmap(z=zeros(2,2)))
                    figure=plot(heatmap(x=(rand(10,10)), y=rand(10,10),z=zeros(10,10)))
                ),
                html_h2("THz-es spektrum"),
                dcc_graph(
                    id="thz_o",
                    figure=plot(heatmap(x=(rand(10,10)), y=rand(10,10),z=zeros(10,10)))
                ))),
            dbc_col((
                html_h2("Pumpa alak"),
                dcc_graph(
                    id="pump_t",
                    figure=plot(heatmap(x=(rand(10,10)), y=rand(10,10),z=zeros(10,10)))
                ),
                html_h2("Pumpa spektrum"),
                dcc_graph(
                    id="pump_o",
                    figure=plot(heatmap(x=(rand(10,10)), y=rand(10,10),z=zeros(10,10)))
                )))
        )),
        html_h1("Metszetek"),
	dbc_row((
		 dbc_col((
			dcc_slider(
			id="xslider",
			min=1,
			step=1,
			max=2048,
			marks=Dict([i=>"$(i)" for i in 1:floor(Int,2048/8):2048]),
		       value=1024
			    )
			  ))
		 )),
        dbc_row((
            dbc_col((
                html_h2("THz-es alak"),
                dcc_graph(
                    id="thz_t_single",
                    figure=plot(p1s)
                ),
                html_h2("THz-es spektrum"),
                dcc_graph(
                    id="thz_o_single",
                    figure=plot(p2s)
                ))),
            dbc_col((
                html_h2("Pumpa alak"),
                dcc_graph(
                    id="pump_t_single",
                    figure=plot(p3s)
                ),
                html_h2("Pumpa spektrum"),
                dcc_graph(
                    id="pump_o_single",
                    figure=plot(p4s)
                )))
        ))
    ), style=Dict("width" => "95vW"))



callback!(app, Output("thz_o", "figure"),Output("thz_t", "figure"), Output("pump_o", "figure"),Output("pump_t", "figure"),Input("crystalupdate", "n_clicks"), State("calcselector","value"),  State("crystalslider","value"),prevent_initial_call=true) do clicks, calc, cslider
    
    FID = h5open(calc,"r")
    t_val = read(FID["/t"])
    x_val = read(FID["/x"])
    o_val = read(FID["/omega"])
    omega0 = read(FID["/omega0"])
    tt,xx = ndgrid(t_val,x_val)
    oo,_ = ndgrid(o_val.-o_val[1],x_val)
    #println(keys(FID))
    p1 =@spawn plot(heatmap(x=oo,y=xx,z=abs.(read(FID["/$(cslider)/ATHz_xo"])),colorscale="Jet")) 
    p2 =@spawn plot(heatmap(x=tt,y=xx,z=real.(read(FID["/$(cslider)/ATHz_xt"])),colorscale="Jet")) 
    p3 =@spawn plot(heatmap(x=oo,y=xx,z=abs.(read(FID["/$(cslider)/Aop"])),colorscale="Jet")) 
    p4 =@spawn plot(heatmap(x=tt,y=xx,z=abs.(read(FID["/$(cslider)/Eop"])),colorscale="Jet")) 
    
    wait.([p1,p2,p3,p4])
    return (p1.result,p2.result,p3.result,p4.result)
end

callback!(app, Output("thz_o_single", "figure"),Output("thz_t_single", "figure"), Output("pump_o_single", "figure"),Output("pump_t_single", "figure"),Input("xslider", "value"),Input("crystalupdate", "n_clicks"), State("calcselector","value"),  State("crystalslider","value"),prevent_initial_call=true) do x_ind, nclicks, calc, cslider
    FID = h5open(calc,"r")
    t_val = read(FID["/t"])
    x_val = read(FID["/x"])
    o_val = read(FID["/omega"])
    omega0 = read(FID["/omega0"])
   # tt,xx = ndgrid(t_val,x_val)
   # oo,_ = ndgrid(o_val.-o_val[1],x_val)
    #println(keys(FID))

    p1s["data"]["x"]=t_val
    p1s["data"]["y"]=read(FID["/$(cslider)/ATHz_xt"])[:,x_ind]

    p2s["data"]["x"]=t_val
    p2s["data"]["y"]=read(FID["/$(cslider)/ATHz_xt"])[:,x_ind]

    p3s["data"]["x"]=t_val
    p3s["data"]["y"]=read(FID["/$(cslider)/ATHz_xt"])[:,x_ind]

    p4s["data"]["x"]=t_val
    p4s["data"]["y"]=read(FID["/$(cslider)/ATHz_xt"])[:,x_ind]

    return (plot(p1s),plot(p2s),plot(p3s),plot(p4s))
end


run_server(app, "0.0.0.0", 12345, debug=false)
