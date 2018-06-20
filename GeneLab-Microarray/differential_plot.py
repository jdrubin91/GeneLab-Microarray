__author__ = 'Jonathan Rubin'

import plotly.plotly as py
import plotly
import plotly.graph_objs as go    
import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
import math
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

def updateSlider(pval_cut):
    return pval_cut

def MA_mpld3(diffExp_file):
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


    foldChange = list()
    averageExpression = list()
    adjustedPvalue = list()
    geneName = list()
    log10pval = list()
    scattersigx = list()
    scattersigy = list()
    volcanosigx = list()
    volcanosigy = list()
    scatterlabels = list()
    scattersiglabels = list()
    volcanolabels = list()
    volcanosiglabels = list()
    cell_text = list()
    pval_cut = 0.1
    with open(diffExp_file) as F:
        header = F.readline().strip('\n').split('\t')
        fc_index = [i for i in range(len(header)) if 'FC' in header[i]][0]+1
        exp_index = [i for i in range(len(header)) if 'Exp' in header[i]][0]+1
        p_index = [i for i in range(len(header)) if 'adj' in header[i]][0]+1
        for line in F:
            linelist = line.strip('\n').split('\t')
            gene = linelist[0].strip('"')
            fc = float(linelist[fc_index])
            exp = float(linelist[exp_index])
            pval = float(linelist[p_index])
            geneName.append(gene)
            foldChange.append(fc)
            averageExpression.append(exp)
            adjustedPvalue.append(pval)
            if pval < pval_cut:
                cell_text.append([gene,exp,fc,pval])
                scattersigx.append(exp)
                scattersigy.append(fc)
                scattersiglabels.append(sig.format(gene=gene,x="%.3f" % exp,y="%.3f" % fc,pval="%.3f" % pval))
                try:
                    l10p = -math.log(pval,10)
                    log10pval.append(l10p)
                    volcanosigy.append(l10p)
                    volcanosigx.append(fc)
                    volcanosiglabels.append(sig.format(gene=gene,x="%.3f" % fc,y="%.3f" % -math.log(pval,10),pval="%.3f" % pval))
                except ValueError:
                    print "Error: Zero adjusted p-value encountered, cannot display in volcano plot.."
                    
            else:
                log10pval.append(-math.log(pval,10))
                scatterlabels.append(non_sig.format(gene=gene,x="%.3f" % exp,y="%.3f" % fc,pval="%.3f" % pval))
                volcanolabels.append(non_sig.format(gene=gene,x="%.3f" % fc,y="%.3f" % -math.log(pval,10),pval="%.3f" % pval))

            
    

    F = plt.figure(figsize=(18,8))
    gs = gridspec.GridSpec(1, 2, width_ratios=[2, 1])
    ax0 = F.add_subplot(gs[0])
    ax0.grid(color='black', linestyle='dashed')
    xy = np.vstack([averageExpression,foldChange])
    z = gaussian_kde(xy)(xy)
    scatter = ax0.scatter(x=averageExpression,y=foldChange,c=z,s=100,edgecolor="")
    sigscatter = ax0.scatter(scattersigx,scattersigy,c='r',s=100,edgecolor="")
    ax0.tick_params(axis='y', which='both', left='on', right='off', labelleft='on')
    ax0.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='on')
    ax0.set_title("Microarray MA Plot", size=25)
    ax0.set_ylabel("Log Fold Change", size=18)
    ax0.set_xlabel("Average Expression", size=18)

    ax1 = F.add_subplot(gs[1])
    volcano = ax1.scatter(x=foldChange,y=log10pval,s=100,edgecolor="",color='navy')
    sigvolcano = ax1.scatter(volcanosigx,volcanosigy,c='r',s=100,edgecolor="")
    ax1.tick_params(axis='y', which='both', left='on', right='off', labelleft='on')
    ax1.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='on')
    ax1.set_title("Microarray Volcano Plot", size=25)
    ax1.set_ylabel("-Log10 P-value", size=18)
    ax1.set_xlabel("Log Fold Change", size=18)
    ax1.grid(color='black', linestyle='dashed')
    ax1.set_ylim(bottom=0)


    F.subplots_adjust(left=0.05,right=0.95,hspace = 0.05, wspace = 0.15)

    tooltip = mpld3.plugins.PointHTMLTooltip(scatter, labels=scatterlabels[:int(0.5*len(scatterlabels))], css=css)
    tooltip2 = mpld3.plugins.PointHTMLTooltip(sigscatter, labels=scattersiglabels, css=css)
    tooltip3 = mpld3.plugins.PointHTMLTooltip(volcano, labels=volcanolabels[:int(0.5*len(volcanolabels))], css=css)
    tooltip4 = mpld3.plugins.PointHTMLTooltip(sigvolcano, labels=volcanosiglabels, css=css)
    # mpld3.plugins.connect(F, tooltip, tooltip2, tooltip3, tooltip4, SliderView(scatter, callback_func="updateSlider"))
    mpld3.plugins.connect(F, tooltip, tooltip2, tooltip3, tooltip4)


    plt.savefig('/Users/jonathanrubin/Google Drive/NASA/home/batch_out/GLDS-4/microarray/MA-Plot_mpld3.png')
    mpld3.save_html(F,'/Users/jonathanrubin/Google Drive/NASA/home/batch_out/GLDS-4/microarray/MA-Plot_mpld3.html')

    with open('/Users/jonathanrubin/Google Drive/NASA/home/batch_out/GLDS-4/microarray/SignificantGenes.html','w') as sigGenes_file:
        sigGenes_file.write("""<!DOCTYPE html>
            <html>
            <head>
            <title>List of Significant Genes</title>
            <style>
            table {
                font-family: arial, sans-serif;
                border-collapse: collapse;
                width: 100%;
            }

            td, th {
                border: 1px solid #dddddd;
                text-align: left;
                padding: 8px;
            }

            tr:nth-child(even) {
                background-color: #dddddd;
            }
            </style>
            </head>
            <body style="width: 1300px; overflow:scroll">
                <h1>List of Significant Genes</h1>
            <div>
                <div style="float: middle; width: 1300px; overflow:scroll; padding-bottom:25px; padding-top:25px">
                    <table> 
                        <tr>
                            <th>Gene</th>
                            <th>Average Expression</th> 
                            <th>Log Fold Change</th>
                            <th>Adjusted P-value</th>
                        </tr>""")
        for row in cell_text:
            name,exp,fc,pval = row
            sigGenes_file.write("""
                        <tr>
                            <td>"""+name+"""</td>
                            <td>"""+str(exp)+"""</td>
                            <td>"""+str(fc)+"""</td>
                            <td>"""+str(pval)+"""</td>
                        </tr>""")
        sigGenes_file.write("""        </table>
                </div>
            </div>
            
            </body>
            </html>""")


    with open('/Users/jonathanrubin/Google Drive/NASA/home/batch_out/GLDS-4/microarray/MA-Plot_mpld3.html','a') as html_file:
        html_file.write('<a style="font-size: 20" href="./SignificantGenes.html">List of Significant Genes</a>')

    plt.close(F)

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
        auto_open=False)

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

