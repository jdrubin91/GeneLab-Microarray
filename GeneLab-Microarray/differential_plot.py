__author__ = 'Jonathan Rubin'

import plotly.plotly as py
import plotly.graph_objs as go
import os, config

def MA_plotly(diffExp_file):
    sig = list()
    sig_text = list()
    non_sig = list()
    non_sig_text = list()
    pval_cut = 0.01
    with open(diffExp_file) as F:
        header = F.readline().strip('\n').split('\t')
        fc_index = [i for i in range(len(header)) if 'FC' in header[i]][0]+1
        exp_index = [i for i in range(len(header)) if 'Exp' in header[i]][0]+1
        p_index = [i for i in range(len(header)) if 'adj' in header[i]][0]+1
        for line in F:
            linelist = line.strip('\n').split('\t')
            FC = linelist[fc_index]
            Exp = linelist[exp_index]
            pval = linelist[p_index]
            if pval < pval_cut:
                sig.append((Exp,FC))
                sig_text.append(linelist[0])
            else:
                non_sig.append((Exp,FC))
                non_sig_text.append(linelist[0])

    trace_comp0 = go.Scatter(
        x=[x for x,y in non_sig],
        y=[y for x,y in non_sig],
        mode='markers',
        marker=dict(size=12,
                    color="navy"
                   ),
        name='Non-Significant Probes',
        text=non_sig_text,
        )

    trace_comp1 = go.Scatter(
        x=[x for x,y in sig],
        y=[y for x,y in sig],
        mode='markers',
        marker=dict(size=12,
                    color="red"
                   ),
        name='Significant Probes',
        text=sig_text,
            )

    data_comp = [trace_comp0, trace_comp1]
    layout_comp = go.Layout(
        title='Microarray MA-Plot',
        hovermode='closest',
        xaxis=dict(
            title='Average Expression',
            ticklen=5,
            zeroline=False,
            gridwidth=2,
        ),
        yaxis=dict(
            title='Log Fold Change',
            ticklen=5,
            gridwidth=2,
            zeroline=True
        ),
    )
    fig_comp = go.Figure(data=data_comp, layout=layout_comp)
    py.iplot(fig_comp, filename='MA-Plot')




if __name__ == "__main__":
    diffExp_file = '/Users/jonathanrubin/Google Drive/NASA/home/batch_out/GLDS-4/microarray/diffExpression.txt'
    MA_plotly(diffExp_file)