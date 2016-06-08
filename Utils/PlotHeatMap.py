import plotly.plotly as py
import plotly.graph_objs as go

# z = data, 
def PlotHeatMap(z, x, y, HeatMapTitle):

    annotations = []
    for n, row in enumerate(z):
        for m, val in enumerate(row):
            var = z[n][m]
            annotations.append(
                dict(
                    text=str(val),
                    x=x[m], y=y[n],
                    xref='x1', yref='y1',
                    font=dict(color='white'),
                    showarrow=False)
                )


    trace = go.Heatmap(x=x, y=y, z=z, showscale=True,
                       colorbar=dict(title='Avg. Time Per Step (sec)',titleside='right')
                        )
    layout = go.Layout(title='Global Font',font=dict(family='Arial, monospace', size=18))
    fig = go.Figure(data=[trace],layout=layout)
    fig['layout'].update(
        title=HeatMapTitle,
        annotations=annotations,
        xaxis=dict(title='Cell Size',ticks='', side='bottom'),
        # ticksuffix is a workaround to add a bit of padding
        yaxis=dict(title='Num PEs Per Node', ticks='', ticksuffix=''),
        width=700,
        height=700,
        autosize=False,
        titlefont=dict(family='Arial, monospace', size=22)
    )
    url = py.plot(fig, filename=HeatMapTitle, height=750)



# x = ['A', 'B', 'C', 'D', 'E']
# y = ['W', 'X', 'Y', 'Z']

# #       x0    x1    x2    x3    x4
# z = [[2.00, 3.00, 0.75, 0.75, 0.00],  # y0
#      [1, 0.00, 0.75, 0.75, 0.00],  # y1
#      [0.75, 0.75, 0.75, 0.75, 0.75],  # y2
#      [0.00, 0.00, 0.00, 0.75, 0.00]]  # y3

# PlotHeatMap(z, x, y, 'HeatMap Title')