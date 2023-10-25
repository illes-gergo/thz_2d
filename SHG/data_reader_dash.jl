using HDF5, Dash, PlotlyJS, DelimitedFiles, DashBootstrapComponents, LazyGrids, Base.Threads
include("data_reader_dash_funcs.jl")

calcdir = "/home/illesg/cst/2d-calculations/"


app = dash(external_stylesheets=[dbc_themes.GRID])


p1s = Plot(scatter(x=sort(rand(10)), y=rand(10)))
p2s = Plot(scatter(x=sort(rand(10)), y=rand(10)))
p3s = Plot(scatter(x=sort(rand(10)), y=rand(10)))
p4s = Plot(scatter(x=sort(rand(10)), y=rand(10)))

p1 = Plot(heatmap(x=(rand(10, 10)), y=rand(10, 10), z=rand(10, 10), colorscale="Jet"))
p2 = Plot(heatmap(x=(rand(10, 10)), y=rand(10, 10), z=rand(10, 10), colorscale="Jet"))
p3 = Plot(heatmap(x=(rand(10, 10)), y=rand(10, 10), z=rand(10, 10), colorscale="Jet"))
p4 = Plot(heatmap(x=(rand(10, 10)), y=rand(10, 10), z=rand(10, 10), colorscale="Jet"))

p1e = Plot(scatter(x=1:10, y=rand(10)))
p2e = Plot(scatter(x=1:10, y=rand(10)))

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
        dbc_row((dbc_col((
                    dcc_slider(
                    id="crystalslider",
                    min=1,
                    step=1,
                    max=401,
                    marks=Dict([i => "$(i)" for i in 1:50:401]),
                    value=101
                )
                ), width=6),
            dbc_col((
                    dcc_input(
                    id="crystal-input",
                    min=1,
                    max=401,
                    step=1,
                    type="number",
                    debounce=true,
                    value=101
                )
                ), width=2),
            dbc_col((
                    html_button("2D ábrák frissítése", id="crystalupdate")
                ), width=2))),
        html_h1("Összefoglaló értékek"),
        dbc_row((
            html_div(
            "Maximális Intenzitást:" * "20GW"
        )
        )),
        dbc_row((
            dbc_col((
                dcc_graph(
                id="effic-graph",
                figure=p1e
            )
            )),
            dbc_col((
                dcc_graph(
                id="efficSH-graph",
                figure=p2e
            )
            ))
        )),
        html_h1("2D ábrák"),
        dbc_row((
            dbc_col((
                html_h2("THz-es alak"),
                dcc_graph(
                    id="thz_t",
                    #figure=Plot(heatmap(z=zeros(2,2)))
                    figure=(p1)
                ),
                html_h2("THz-es spektrum"),
                dcc_graph(
                    id="thz_o",
                    figure=(p2)
                ))),
            dbc_col((
                html_h2("Pumpa alak"),
                dcc_graph(
                    id="pump_t",
                    figure=(p3)
                ),
                html_h2("Pumpa spektrum"),
                dcc_graph(
                    id="pump_o",
                    figure=(p4)
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
                marks=Dict([i => "$(i)" for i in 1:floor(Int, 2048 / 8):2048]),
                value=1024
            )
            ), width=6)
        )),
        dbc_row((
            dbc_col((
                html_h2("THz-es alak"),
                dcc_graph(
                    id="thz_t_single",
                    figure=(p1s)
                ),
                html_h2("THz-es spektrum"),
                dcc_graph(
                    id="thz_o_single",
                    figure=(p2s)
                ))),
            dbc_col((
                html_h2("Pumpa alak"),
                dcc_graph(
                    id="pump_t_single",
                    figure=(p3s)
                ),
                html_h2("Pumpa spektrum"),
                dcc_graph(
                    id="pump_o_single",
                    figure=(p4s)
                )))
        ))
    ), style=Dict("width" => "95vW"))



callback!(app, Output("thz_o", "figure"), Output("thz_t", "figure"), Output("pump_o", "figure"), Output("pump_t", "figure"), Input("crystalupdate", "n_clicks"), State("calcselector", "value"), State("crystalslider", "value"), prevent_initial_call=true) do clicks, calc, cslider

    FID = h5open(calc, "r")
    t_val = read(FID["/t"])
    x_val = read(FID["/x"])
    o_val = read(FID["/omega"])
    omega0 = read(FID["/omega0"])
    tt, xx = ndgrid(t_val, x_val)
    oo, _ = ndgrid(o_val .- o_val[1], x_val)
    #println(keys(FID))
    t1 = @spawn begin
        z_raw = read(FID["/$(cslider)/ATHz_xo"])
        _, I = findmax(abs.(z_raw))
        z_shifted = circshift(z_raw, [0, 1024 - I[2]])
        p1.data[1][:x] = oo[1:floor(Int, end / 8), floor(Int, 3 * end / 8):floor(Int, 5 * end / 8)] ./ 2 ./ pi
        p1.data[1][:y] = xx[1:floor(Int, end / 8), floor(Int, 3 * end / 8):floor(Int, 5 * end / 8)]
        p1.data[1][:z] = abs.(z_shifted[1:floor(Int, end / 8), floor(Int, 3 * end / 8):floor(Int, 5 * end / 8)])
    end

    t2 = @spawn begin
        z_raw = read(FID["/$(cslider)/ATHz_xt"])
        _, I = findmax(abs.(z_raw))
        z_shifted = circshift(z_raw, [floor(Int, size(z_raw, 1) / 2), floor(Int, size(z_raw, 2) / 2)] .- [I[1], I[2]])
        z_processed = real.(z_shifted[floor(Int, 3 * end / 8):floor(Int, 5 * end / 8), floor(Int, 3 * end / 8):floor(Int, 5 * end / 8)])
        p2.data[1][:x] = tt[floor(Int, 3 * end / 8):floor(Int, 5 * end / 8), floor(Int, 3 * end / 8):floor(Int, 5 * end / 8)]
        p2.data[1][:y] = xx[floor(Int, 3 * end / 8):floor(Int, 5 * end / 8), floor(Int, 3 * end / 8):floor(Int, 5 * end / 8)]
        p2.data[1][:z] = real.(z_processed)
    end

    t3 = @spawn begin
        z_raw = read(FID["/$(cslider)/Aop"])
        _, I = findmax(abs.(z_raw))
        z_shifted = circshift(z_raw, [0, 1024 - I[2]])
        p3.data[1][:x] = oo[floor(Int, 3 * end / 8):floor(Int, 5 * end / 8), floor(Int, 3 * end / 8):floor(Int, 5 * end / 8)] ./ 2 ./ pi
        p3.data[1][:y] = xx[floor(Int, 3 * end / 8):floor(Int, 5 * end / 8), floor(Int, 3 * end / 8):floor(Int, 5 * end / 8)]
        p3.data[1][:z] = abs.(z_shifted[floor(Int, 3 * end / 8):floor(Int, 5 * end / 8), floor(Int, 3 * end / 8):floor(Int, 5 * end / 8)])
    end

    t4 = @spawn begin
        z_raw = read(FID["/$(cslider)/Eop"])
        _, I = findmax(abs.(z_raw))
        z_shifted = circshift(z_raw, [floor(Int, size(z_raw, 1) / 2), floor(Int, size(z_raw, 2) / 2)] .- [I[1], I[2]])
        z_processed = abs.(z_shifted[floor(Int, 3 * end / 8):floor(Int, 5 * end / 8), floor(Int, 3 * end / 8):floor(Int, 5 * end / 8)]) .^ 2
        p4.data[1][:x] = tt[floor(Int, 3 * end / 8):floor(Int, 5 * end / 8), floor(Int, 3 * end / 8):floor(Int, 5 * end / 8)]
        p4.data[1][:y] = xx[floor(Int, 3 * end / 8):floor(Int, 5 * end / 8), floor(Int, 3 * end / 8):floor(Int, 5 * end / 8)]
        p4.data[1][:z] = z_processed
    end
    wait.([t1, t2, t3, t4])
    return (p1, p2, p3, p4)
end

callback!(app, Output("thz_o_single", "figure"), Output("thz_t_single", "figure"), Output("pump_o_single", "figure"), Output("pump_t_single", "figure"), Input("xslider", "value"), Input("crystalupdate", "n_clicks"), State("calcselector", "value"), State("crystalslider", "value"), prevent_initial_call=true) do x_ind, nclicks, calc, cslider
    FID = h5open(calc, "r")
    t_val = read(FID["/t"])
    x_val = read(FID["/x"])
    o_val = read(FID["/omega"])
    omega0 = read(FID["/omega0"])
    # tt,xx = ndgrid(t_val,x_val)
    # oo,_ = ndgrid(o_val.-o_val[1],x_val)
    #println(keys(FID))

    p1s.data[1][:x] = (o_val .- o_val[1]) ./ 2 ./ pi
    p1s.data[1][:y] = abs.(read(FID["/$(cslider)/ATHz_xo"])[:, x_ind])

    p2s.data[1][:x] = t_val
    _, I = findmax(abs.(read(FID["/$(cslider)/ATHz_xt"])))
    p2s.data[1][:y] = circshift(real.(read(FID["/$(cslider)/ATHz_xt"])[:, x_ind]), 512 - I[1])

    p3s.data[1][:x] = o_val ./ 2 ./ pi
    p3s.data[1][:y] = abs.(read(FID["/$(cslider)/Aop"])[:, x_ind])

    p4s.data[1][:x] = t_val
    _, I = findmax(abs.(read(FID["/$(cslider)/Eop"])[:, x_ind]))
    p4s.data[1][:y] = circshift(abs.(read(FID["/$(cslider)/Eop"])[:, x_ind]) .^ 2, 512 - I[1])

    return (p1s, p2s, p3s, p4s)
end

callback!(app, Output("effic-graph", "figure"), Input("calcselector", "value")) do calc
    FID = h5open(calc, "r")
    z_val = read(FID["/z"])
    En0 = sum(read(FID["/En"])[1, :])
    EnTHz = sum(read(FID["/EnTHz"]), dims=2)
    effic = EnTHz ./ En0
    p1e.data[1][:x] = z_val
    p1e.data[1][:y] = effic
    return (p1e)
end


callback!(app, Output("crystalslider", "value"), Input("crystal-input", "value")) do clength
    return clength
end

callback!(app, Output("crystal-input", "value"), Input("crystalslider", "value")) do clength
    return clength
end

run_server(app, "0.0.0.0", 12345, debug=false)
