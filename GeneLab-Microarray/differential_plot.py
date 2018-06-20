__author__ = 'Jonathan Rubin'

import plotly.plotly as py
import plotly
import plotly.graph_objs as go    
import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
import mpld3
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
            

    fig, ax = plt.subplots(figsize=(18,8))
    # scatter = ax.scatter(x=averageExpression,
    #     y=foldChange,
    #     c=['r' if pval < pval_cut else 'b' for pval in adjustedPvalue],
    #     alpha=[1 if pval < pval_cut else 0.3 for pval in adjustedPvalue],
    #     edgecolor="")
    scatter = plt.scatter(x=averageExpression,
        y=foldChange,
        color=[(1.0,0.0,0.0,1.0) if pval < pval_cut else (0.0,0.0,0.5,0.1) for pval in adjustedPvalue],
        s=[60 if pval < pval_cut else 40 for pval in adjustedPvalue],
        edgecolor="")
    ax.grid(color='white', linestyle='solid')
    ax.set_title("Microarray MA-Plot", size=25)
    ax.set_ylabel("Log Fold Change", size=18)
    ax.set_xlabel("Average Expression", size=18)
    # labels = [gene+" ("+str(x)+","+str(y)+")" for (gene,x,y) in zip(geneName,foldChange,averageExpression)]
    labels = ['<div style="margin: 10px;background-color: #ff0000;border: 1px solid white;opacity: 0.9"><p style="font-family: arial;color: #ffffff">Gene: {gene}<br />({x}, {y})</p></div>'.format(gene=gene,x="%.2f" % x,y="%.2f" % y) 
        if pval < pval_cut else '<div style="margin: 10px;background-color: #0000ff;border: 1px solid white;opacity: 0.9"><p style="font-family: arial;color: #ffffff">Gene: {gene}<br />({x}, {y})</p></div>'.format(gene=gene,x="%.2f" % x,y="%.2f" % y) 
        for (gene,x,y,pval) in zip(geneName,averageExpression,foldChange,adjustedPvalue)]

    # tooltip = mpld3.plugins.PointLabelTooltip(scatter, labels=labels)
    tooltip = mpld3.plugins.PointHTMLTooltip(scatter, labels=labels)
    mpld3.plugins.connect(fig, tooltip)
    mpld3.plugins.connect(fig,SliderView(scatter, callback_func="updateSlider"))

    mpld3.save_html(fig,'/Users/jonathanrubin/Google Drive/NASA/home/batch_out/GLDS-4/microarray/MA-Plot_mpld3.html')

    plt.close(fig)


if __name__ == "__main__":
    diffExp_file = '/Users/jonathanrubin/Google Drive/NASA/home/batch_out/GLDS-4/microarray/diffExpression.txt'
    # MA_plotly(diffExp_file)
    MA_mpld3(diffExp_file)

