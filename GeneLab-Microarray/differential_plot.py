__author__ = 'Jonathan Rubin'

import plotly.plotly as py
import plotly
import plotly.graph_objs as go    
import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
import mpld3
import pygal
from pygal import config
from scipy.stats import gaussian_kde
import os, config

class SliderView(mpld3.plugins.PluginBase):
    """ Add slider and JavaScript / Python interaction. """

    JAVASCRIPT = """
    mpld3.register_plugin("sliderview", SliderViewPlugin);
    SliderViewPlugin.prototype = Object.create(mpld3.Plugin.prototype);
    SliderViewPlugin.prototype.constructor = SliderViewPlugin;
    SliderViewPlugin.prototype.requiredProps = ["idline", "callback_func"];
    SliderViewPlugin.prototype.defaultProps = {}

    function SliderViewPlugin(fig, props){
        mpld3.Plugin.call(this, fig, props);
    };

    SliderViewPlugin.prototype.draw = function(){
      var line = mpld3.get_element(this.props.idline);
      var callback_func = this.props.callback_func;

      var div = d3.select("#" + this.fig.figid);

      // Create slider
      div.append("input").attr("type", "range").attr("min", 0).attr("max", 1).attr("step", 0.01).attr("value", 1)
          .on("change", function() {
              var command = callback_func + "(" + this.value + ")";
              console.log("running "+command);
              var callbacks = { 'iopub' : {'output' : handle_output}};
              var kernel = IPython.notebook.kernel;
              kernel.execute(command, callbacks, {silent:false});
          });

      function handle_output(out){
        //console.log(out);
        var res = null;
        // if output is a print statement
        if (out.msg_type == "stream"){
          res = out.content.data;
        }
        // if output is a python object
        else if(out.msg_type === "pyout"){
          res = out.content.data["text/plain"];
        }
        // if output is a python error
        else if(out.msg_type == "pyerr"){
          res = out.content.ename + ": " + out.content.evalue;
          alert(res);
        }
        // if output is something we haven't thought of
        else{
          res = "[out type not implemented]";  
        }

        // Update line data
        line.data = JSON.parse(res);
        line.elements()
          .attr("d", line.datafunc(line.data))
          .style("stroke", "black");

       }

    };
    """

    def __init__(self, line, callback_func):
        self.dict_ = {"type": "sliderview",
                      "idline": mpld3.utils.get_id(line),
                      "callback_func": callback_func}

def MA_plotly(diffExp_file):
    sig = list()
    sig_text = list()
    non_sig = list()
    non_sig_text = list()
    pval_cut = 0.1
    with open(diffExp_file) as F:
        header = F.readline().strip('\n').split('\t')
        fc_index = [i for i in range(len(header)) if 'FC' in header[i]][0]+1
        exp_index = [i for i in range(len(header)) if 'Exp' in header[i]][0]+1
        p_index = [i for i in range(len(header)) if 'adj' in header[i]][0]+1
        for line in F:
            linelist = line.strip('\n').split('\t')
            FC = linelist[fc_index]
            Exp = linelist[exp_index]
            pval = float(linelist[p_index])
            if pval < pval_cut:
                sig.append((Exp,FC))
                sig_text.append(linelist[0].strip('"'))
            else:
                non_sig.append((Exp,FC))
                non_sig_text.append(linelist[0].strip('"'))


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
    plotly.offline.plot(fig_comp, filename='/Users/jonathanrubin/Google Drive/NASA/home/batch_out/GLDS-4/microarray/MA-Plot.html',
        auto_open=False,)

def updateSlider(pval_cut):
    return pval_cut

def MA_mpld3(diffExp_file):
    foldChange = list()
    averageExpression = list()
    adjustedPvalue = list()
    geneName = list()
    pval_cut = 0.1
    with open(diffExp_file) as F:
        header = F.readline().strip('\n').split('\t')
        fc_index = [i for i in range(len(header)) if 'FC' in header[i]][0]+1
        exp_index = [i for i in range(len(header)) if 'Exp' in header[i]][0]+1
        p_index = [i for i in range(len(header)) if 'adj' in header[i]][0]+1
        for line in F:
            linelist = line.strip('\n').split('\t')
            geneName.append(linelist[0].strip('"'))
            foldChange.append(float(linelist[fc_index]))
            averageExpression.append(float(linelist[exp_index]))
            adjustedPvalue.append(float(linelist[p_index]))
            
    css = """
        table
        {
          border-collapse: collapse;
        }
        th
        {
          color: #ffffff;
          background-color: #000000;
        }
        td
        {
          background-color: #cccccc;
          color: #ffffff;
        }
        table, th, td
        {
          font-family:Arial, Helvetica, sans-serif;
          border: 1px solid white;
          text-align: right;
          opacity: 0.9;
        }
        """

    F = plt.figure(figsize=(18,8))
    ax0 = F.add_subplot(111)
    ax0.grid(color='black', linestyle='dashed')
    x=averageExpression
    y=foldChange
    xy = np.vstack([x,y])
    z = gaussian_kde(xy)(xy)
    # scatter = ax0.scatter(x=averageExpression,
    #     y=foldChange,
    #     color=[(1.0,0.0,0.0,1.0) if pval < pval_cut else (0.0,0.0,0.5,0.3) for pval in adjustedPvalue],
    #     s=[60 if pval < pval_cut else 40 for pval in adjustedPvalue],
    #     edgecolor="")
    scatter = ax0.scatter(x=x,y=y,c=z,s=100,edgecolor="")
    sigx = [x for x,pval in zip(averageExpression,adjustedPvalue) if pval < pval_cut]
    sigy = [y for y,pval in zip(foldChange,adjustedPvalue) if pval < pval_cut]
    sigscatter = ax0.scatter(sigx,sigy,c='r',s=100,edgecolor="")
    ax0.tick_params(axis='y', which='both', left='on', right='off', labelleft='on')
    ax0.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='on')
    ax0.set_title("Microarray MA-Plot", size=25)
    ax0.set_ylabel("Log Fold Change", size=18)
    ax0.set_xlabel("Average Expression", size=18)
    # labels = ['<div style="margin: 10px;background-color: #ff0000;border: 1px solid white;opacity: 0.9"><p style="font-family: arial;color: #ffffff">Gene: {gene}<br />({x}, {y})</p></div>'.format(gene=gene,x="%.2f" % x,y="%.2f" % y) 
    #     if pval < pval_cut else '<div style="margin: 10px;background-color: #0000ff;border: 1px solid white;opacity: 0.9"><p style="font-family: arial;color: #ffffff">Gene: {gene}<br />({x}, {y})</p></div>'.format(gene=gene,x="%.2f" % x,y="%.2f" % y) 
    #     for (gene,x,y,pval) in zip(geneName,averageExpression,foldChange,adjustedPvalue)]

    sig = """<table>
    <tr>
        <th></th>
        <th>{gene}</th>
    </tr>
    <tr>
        <td style="background-color: #000000;">x</td>
        <td style="background-color: #ff0000;">{x}</td>
    </tr>
    <tr>
        <td style="background-color: #000000;">y</td>
        <td style="background-color: #ff0000;">{y}</td>
    </tr>
    <tr>
        <td style="background-color: #000000;">p</td>
        <td style="background-color: #ff0000;">{pval}</td>
    </tr>
    </table>"""

    non_sig = """<table>
    <tr>
        <th></th>
        <th>{gene}</th>
    </tr>
    <tr>
        <td style="background-color: #000000;">x</td>
        <td style="background-color: #0000ff;">{x}</td>
    </tr>
    <tr>
        <td style="background-color: #000000;">y</td>
        <td style="background-color: #0000ff;">{y}</td>
    </tr>
    <tr>
        <td style="background-color: #000000;">p</td>
        <td style="background-color: #0000ff;">{pval}</td>
    </tr>
    </table>"""

    labels = [sig.format(gene=gene,x="%.3f" % x,y="%.3f" % y,pval="%.3f" % pval) 
        if pval < pval_cut else non_sig.format(gene=gene,x="%.3f" % x,y="%.3f" % y,pval="%.3f" % pval) 
        for (gene,x,y,pval) in zip(geneName,averageExpression,foldChange,adjustedPvalue)]

    tooltip = mpld3.plugins.PointHTMLTooltip(scatter, labels=labels, css=css)
    tooltip2 = mpld3.plugins.PointHTMLTooltip(sigscatter, labels=labels, css=css)
    mpld3.plugins.connect(F, tooltip, tooltip2, SliderView(scatter, callback_func="updateSlider"))


    # ax1 = plt.subplot(gs[1])
    # columns = ['Gene','AvgExp','LogFC','AdjPval']
    # cell_text = [[w,"%.3f" % x,"%.3f" % y,"%.3f" % z] for w,x,y,z in zip(geneName,averageExpression,foldChange,adjustedPvalue) if z < pval_cut]
    # the_table = ax1.table(cellText=cell_text,loc='top right',cellLoc='left',colLabels=columns,bbox=[0,0,1,1])
    # the_table.auto_set_font_size(False)
    # the_table.set_fontsize(14)
    # cellDict=the_table.get_celld()
    # for i in range(len(cell_text)+1):
    #     cellDict[(i,0)].set_width(0.35)

    
    # ax1.set_xlim([0,1])
    # ax1.set_ylim([0,1])
    # ax1.grid(color='black', linestyle='solid')
    # i = np.linspace(0,1,len(cell_text)+1)
    # k=0
    # for i in np.linspace(0.95,0,len(cell_text)+1):
    #     if i == 0.95:
    #         ax1.text(0.5, i, ' '.join(columns), size=16, ha='center')
    #     else:
    #         ax1.text(0.5, i, ' '.join(cell_text[k]), size=16, ha='center')
    #         k+=1
            
    # ax1.tick_params(axis='y', which='both', left='off', right='off', labelleft='off')
    # ax1.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='off')

    plt.savefig('/Users/jonathanrubin/Google Drive/NASA/home/batch_out/GLDS-4/microarray/MA-Plot_mpld3.png')
    mpld3.save_html(F,'/Users/jonathanrubin/Google Drive/NASA/home/batch_out/GLDS-4/microarray/MA-Plot_mpld3.html')
    with open('/Users/jonathanrubin/Google Drive/NASA/home/batch_out/GLDS-4/microarray/MA-Plot_mpld3.html','a') as html_file:
        html_file.write('<a style="font-size: 20" href="this is a link">This is a link</a>')

    plt.close(F)


def MA_pygal(diffExp_file):
    foldChange = list()
    averageExpression = list()
    adjustedPvalue = list()
    geneName = list()
    pval_cut = 0.1
    with open(diffExp_file) as F:
        header = F.readline().strip('\n').split('\t')
        fc_index = [i for i in range(len(header)) if 'FC' in header[i]][0]+1
        exp_index = [i for i in range(len(header)) if 'Exp' in header[i]][0]+1
        p_index = [i for i in range(len(header)) if 'adj' in header[i]][0]+1
        for line in F:
            linelist = line.strip('\n').split('\t')
            geneName.append(linelist[0].strip('"'))
            foldChange.append(float(linelist[fc_index]))
            averageExpression.append(float(linelist[exp_index]))
            adjustedPvalue.append(float(linelist[p_index]))

    xy_chart = pygal.XY(stroke=False,title="Microarray MA-Plot",x_title="Log Fold Change",y_title="Average Expression",show_legend=False)
    xy_chart.add('Sig',[{'value':(x,y),'label':g} for g,x,y,p in zip(geneName,averageExpression,foldChange,adjustedPvalue) if p < pval_cut])
    xy_chart.add('Non-Sig',[{'value':(x,y),'label':g} for g,x,y,p in zip(geneName,averageExpression,foldChange,adjustedPvalue) if p >= pval_cut])
    xy_chart.render_to_file('/Users/jonathanrubin/Google Drive/NASA/home/batch_out/GLDS-4/microarray/MA-Plot_pygal.svg')

if __name__ == "__main__":
    diffExp_file = '/Users/jonathanrubin/Google Drive/NASA/home/batch_out/GLDS-4/microarray/diffExpression.txt'
    # MA_plotly(diffExp_file)
    MA_mpld3(diffExp_file)
    # MA_pygal(diffExp_file)

