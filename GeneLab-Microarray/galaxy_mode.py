#!/usr/bin/python

__author__ = 'Jonathan Rubin'

import matplotlib as mpl
mpl.use('Agg')
mpl.rcParams['image.cmap'] = 'jet'
mpl.rcParams.update({'figure.autolayout': True})
import os, subprocess, config, math, mpld3, warnings, json, pylab, requests
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde
import scipy
import scipy.cluster.hierarchy as sch
import scipy.spatial.distance as dist
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


#mpld3 hack to make it run with matplotlib v2.2.2
class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)

from mpld3 import _display
_display.NumpyEncoder = NumpyEncoder

#The main command that is run when galaxy mode is specified by user.
def run(counts_table,metadata,diff_analysis,condition1,condition2,padj_cutoff,outliers,html_main,html_folder):
    #Create folder for secondary html output
    if not os.path.exists(html_folder):
        os.makedirs(html_folder)

    #Copies the NASA GeneLab banner to be used in visualization
    copy_logo_command=["scp",os.path.join(config.srcdir,"logos","NASA_GeneLab_banner.jpg"),html_folder]
    subprocess.call(copy_logo_command)
    copy_logo_command=["scp",os.path.join(config.srcdir,"logos","GeneLab_logo.png"),html_folder]
    subprocess.call(copy_logo_command)

    #Performs differential expression via limma using custom R scripts
    limma_output = limma_differential(counts_table,metadata,diff_analysis,condition1,condition2,outliers,html_folder)
    print '\n'.join(limma_output.split('\n')[:-1])

    #Creates MA and Volcano plots to be used in visualization. Also creates the main html file
    sig_genes = differential_visualize(os.path.join(html_folder,'limma_out.txt'),
                        padj_cutoff,condition1,condition2,html_main,html_folder)

    #Gets significant genes based on differential expression to be plotted in a heatmap later
    x,x_full,x_50,row_header,column_header,sig_header = get_matrix(sig_genes,counts_table,limma_output,condition1,condition2)

    #Create a PCA plot of samples based on user inputted counts table
    pca(x_full,column_header,condition1,condition2,html_folder)

    if len(x) != 0:
        #Creates heatmap of significant genes
        row_method = 'average'
        column_method = 'average'
        row_metric = 'cityblock'
        column_metric = 'euclidean'
        color_gradient = 'red_black_green'
        heatmap(x, sig_header, column_header, row_method,
                column_method, row_metric, column_metric,
                color_gradient, html_folder)
    else:
        #Creates heatmap of significant genes
        print "Warning: No significant genes found at padj < " +str(padj_cutoff)+ ". Displaying top 50 genes in heatmap."
        row_method = 'average'
        column_method = 'average'
        row_metric = 'cityblock'
        column_metric = 'euclidean'
        color_gradient = 'red_black_green'
        heatmap(x_50, row_header, column_header, row_method,
                column_method, row_metric, column_metric,
                color_gradient, html_folder)

def pca(x_full,column_header,condition1,condition2,html_folder):
    pca = PCA(n_components=2)
    x_full = np.transpose(x_full)
    principalComponents = pca.fit_transform(x_full)
    PC = np.transpose(principalComponents)
    v1,v2 = [f*100 for f in pca.explained_variance_ratio_]
    v1 = "%.2f" % v1
    v2 = "%.2f" % v2
    c1_x = [x for x,i in zip(PC[0],column_header) if condition1 in i]
    c1_y = [y for y,i in zip(PC[1],column_header) if condition1 in i]
    c2_x = [x for x,i in zip(PC[0],column_header) if condition2 in i]
    c2_y = [y for y,i in zip(PC[1],column_header) if condition2 in i]
    c1_l = [l.split(': ') for l in column_header if condition1 in l]
    c2_l = [l.split(': ') for l in column_header if condition2 in l]
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
          opacity: 0.95;
        }
        """

    C1_html = """<table>
        <tr>
            <th></th>
            <th>{sample}</th>
        </tr>
        <tr>
            <td style="background-color: #000000;">c</td>
            <td style="background-color: #BA3232;">{condition1}</td>
        </tr>
        <tr>
            <td style="background-color: #000000;">x</td>
            <td style="background-color: #BA3232;">{x}</td>
        </tr>
        <tr>
            <td style="background-color: #000000;">y</td>
            <td style="background-color: #BA3232;">{y}</td>
        </tr>
        </table>"""

    C2_html = """<table>
        <tr>
            <th></th>
            <th>{sample}</th>
        </tr>
        <tr>
            <td style="background-color: #000000;">c</td>
            <td style="background-color: #32BABA;">{condition2}</td>
        </tr>
        <tr>
            <td style="background-color: #000000;">x</td>
            <td style="background-color: #32BABA;">{x}</td>
        </tr>
        <tr>
            <td style="background-color: #000000;">y</td>
            <td style="background-color: #32BABA;">{y}</td>
        </tr>
        </table>"""


    label1 = [C1_html.format(sample=sample,condition1=condition1,x="%.3f" % x,y="%.3f" % y) for (sample,condition1),x,y in zip(c1_l,c1_x,c1_y)]
    label2 = [C2_html.format(sample=sample,condition2=condition2,x="%.3f" % x,y="%.3f" % y) for (sample,condition2),x,y in zip(c2_l,c2_x,c2_y)]

    #Here we create the Volcano plot
    F = plt.figure(figsize=(11,7))
    ax1 = F.add_subplot(111)
    pca1 = ax1.scatter(x=c1_x,y=c1_y,c='#BA3232',s=300,edgecolor="")
    plot1 = ax1.plot([], [], "o", color="#BA3232", label=condition1,markersize=15,markeredgecolor="#BA3232")
    pca2 = ax1.scatter(x=c2_x,y=c2_y,c='#32BABA',s=300,edgecolor="")
    plot2 = ax1.plot([], [], "o", color="#32BABA", label=condition2,markersize=15,markeredgecolor="#32BABA")
    cover = ax1.plot([],[], "o", color='w', markersize=16,markeredgecolor='w')
    ax1.tick_params(axis='y', which='both', left='on', right='off', labelleft='on')
    ax1.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='on')
    ax1.set_title("PCA Plot", size=25)
    ax1.set_ylabel("Principal Component 2 (" + v2 + "%)", size=18)
    ax1.set_xlabel("Principal Component 1 (" + v1 + "%)", size=18)
    ax1.grid(color='gray')
    miny = min(c1_y+c2_y)
    maxy = max(c1_y+c2_y)
    yoff = (maxy-miny)/2.0
    minx = min(c1_x+c2_x)
    maxx = max(c1_x+c2_x)
    xoff = (maxx-minx)/2.0
    ax1.set_ylim([miny-yoff,maxy+yoff])
    ax1.set_xlim([minx-xoff,maxx+xoff])
    l = ax1.legend(loc="best", framealpha=0.5, fancybox=True, numpoints=1,fontsize=18)
    colors = ['#BA3232','#32BABA']
    for i, text in enumerate(l.get_texts()):
        text.set_color(colors[i])
    plt.savefig(os.path.join(html_folder,'PCA.png'))
    plt.savefig(os.path.join(html_folder,'PCA.svg'))

    #This is the bulk of the mpld3 code, basically we're creating 'tooltips' which allows us to have the interactive labels
    tooltip1 = mpld3.plugins.PointHTMLTooltip(pca1, labels=label1, css=css)
    tooltip2 = mpld3.plugins.PointHTMLTooltip(pca2, labels=label2, css=css)

    #This connects our plots together
    mpld3.plugins.connect(F,tooltip1,tooltip2)

    #Here we save both a png version of the plot (non-interactive) and the interactive html version of the plot
    mpld3.save_html(F,os.path.join(html_folder,'pca.html'))

    #Close the figure
    plt.close(F)

    #Here we add a 'back' button to the interactive html plots to go back to the original report.
    with open(os.path.join(html_folder,'Pca.html'),'a') as F:
        F.write('<b><a href="index.html">Back</a></b>')


#Get a ndarray matrix with only significant genes and full genes
def get_matrix(sig_genes,counts_table,limma_output,condition1,condition2):
    limma_output = limma_output.split('\n')
    group1_index = [i+1 for i in range(len(limma_output)) if 'group 1:' in limma_output[i]][0]
    group2_index = [i+1 for i in range(len(limma_output)) if 'group 2:' in limma_output[i]][0]
    column_header = limma_output[group1_index].split() + limma_output[group2_index].split()

    x = list()
    x_full = list()
    x_50 = list()
    row_header = list()
    sig_header = list()
    with open(counts_table,'r') as F:
        header = F.readline().strip('\n').split('\t')
        indexes = [header.index(i) for i in column_header]
        counter = 0
        for line in F:
            line = line.strip('\n').split('\t')
            append_list = map(float,[line[index] for index in indexes])
            x_full.append(append_list)
            if counter < 50:
                x_50.append(append_list)
                row_header.append(line[0])
            if line[0] in sig_genes:
                x.append(append_list)
                sig_header.append(line[0])
            counter += 1

    x = np.squeeze(np.asarray(x))
    x_full= np.squeeze(np.asarray(x_full))
    column_header = [i+': '+condition1 for i in limma_output[group1_index].split()] + [i+': '+condition2 for i in limma_output[group2_index].split()]


    return x,x_full,x_50,row_header,column_header,sig_header



#Runs differential expression on two conditions that must be within the inputted metadata
def limma_differential(counts_table,metadata,diff_analysis,condition1,condition2,outliers,html_folder):
    print diff_analysis
    if diff_analysis == 'limma':
        limma_script = os.path.join(config.R_dir,'limmaDiffExp.R')
        if outliers != 'None':
            outliers = '_'.join(outliers.split(','))
            limma_differential_command = ["Rscript", "--vanilla", limma_script, 
                                            "-d", counts_table, 
                                            "-i", metadata, 
                                            "--group1=" + condition1, 
                                            "--group2=" + condition2, 
                                            "-r", outliers,
                                            "-o", os.path.join(html_folder,'limma_out.txt')]
        else:
            limma_differential_command = ["Rscript", "--vanilla", limma_script, 
                                            "-d", counts_table, 
                                            "-i", metadata, 
                                            "--group1=" + condition1,
                                            "--group2=" + condition2,
                                            "-o", os.path.join(html_folder,'limma_out.txt')]
    elif diff_analysis == 'voom':
        limma_script = os.path.join(config.R_dir,'limmaVoom.R')
        if outliers != 'None':
            outliers = '_'.join(outliers.split(','))
            limma_differential_command = ["Rscript", "--vanilla", limma_script, 
                                            "-m", counts_table, 
                                            "-i", metadata, 
                                            "--group1=" + condition1, 
                                            "--group2=" + condition2, 
                                            "-r", outliers,
                                            "-o", os.path.join(html_folder,'limma_out.txt')]
        else:
            limma_differential_command = ["Rscript", "--vanilla", limma_script, 
                                            "-m", counts_table, 
                                            "-i", metadata, 
                                            "--group1=" + condition1,
                                            "--group2=" + condition2,
                                            "-o", os.path.join(html_folder,'limma_out.txt')]
    else:
        print "Error: Could not detect type of differential analysis to be performed. Exiting..."
        sys.exit(1)


    try:
        output = str(subprocess.check_output(limma_differential_command).decode("utf-8"))
    except:
        subprocess.call(limma_differential_command)

    return output


#Saves an interactive html graph based on user condition inputs.
def differential_visualize(diffExp_file,pval_cut,condition1,condition2,html_main,html_folder):
    #In this section of the code, the style of labels is defined. In this case, we are using a table with non_sig hits being blue and sig hits being red
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


    #In this section, the differential expression file is parsed and all necessary lists are generated (to avoid looping through multiple times)
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
    sig_genes = list()
    pval_cut = float(pval_cut)
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
            fc_short = "%.3f" % fc
            exp_short = "%.3f" % exp
            pval_short = "%.3f" % pval
            geneName.append(gene)
            foldChange.append(fc)
            averageExpression.append(exp)
            adjustedPvalue.append(pval)
            if pval < pval_cut:
                sig_genes.append([gene,exp,fc,pval])
                scattersigx.append(exp)
                scattersigy.append(fc)
                scattersiglabels.append(sig.format(gene=gene,x=exp_short,y=fc_short,pval=pval_short))
                try:
                    l10p = -math.log(pval,10)
                    log10pval.append(l10p)
                    volcanosigy.append(l10p)
                    volcanosigx.append(fc)
                    volcanosiglabels.append(sig.format(gene=gene,x=fc_short,y="%.3f" % -math.log(pval,10),pval=pval_short))
                except ValueError:
                    print "Warning: Gene " + gene + " has an adjusted p-value of zero. Cannot display in volcano plot.."
                    
            else:
                log10pval.append(-math.log(pval,10))
                scatterlabels.append(non_sig.format(gene=gene,x=exp_short,y=fc_short,pval="%.3f" % pval))
                volcanolabels.append(non_sig.format(gene=gene,x=fc_short,y="%.3f" % -math.log(pval,10),pval=pval_short))

    #Here we create the main report html page. We will use this template to create two html documents to support a 'back' button on individual plots
    main_report = """<!DOCTYPE html>
<html>
    <head>
        <title>GeneLab-Visualize</title>
        <style>
            * {
                    font-family:Helvetica Neue, Helvetica, sans-serif;
                }
            img {
                    max-width: 100%;
                    max-height: 100%;
                }
            .banner {
                    position: relative;
                    padding: 0px;
                    text-align: center;
                    color: white;
                    width: 100%;
                        }
            .top-left {
                    position: absolute;
                    top: 0px;
                    left: 25px;
                        }
            .centered {
                    position: absolute;
                    top: 50%;
                    left: 20%;
                    transform: translate(-50%, -50%);
                        }
            .row {
                    display: flex; /* equal height of the children */
                }
            .figure {
                    border: 1px solid #ddd;
                    background-color:#ffffff;
                    border-radius: 2px;
                    padding: 5px;
                    width: 48%;
                    margin: 5px;
                }
            .figure:hover {
                    box-shadow: 0 0 5px 2px rgba(0, 140, 186, 0.5);
                        }
            .footer {
                    position: relative;
                    bottom: 0px;
                    height: 120px;
                    padding: 20px;
                    }
        </style
    </head>
    <body style="width:900px; margin:0 auto;">
        <div class="banner">
            <a target="_blank" href="https://genelab.nasa.gov">
                <img src="NASA_GeneLab_banner.jpg" alt="GeneLab banner">
            </a>
        </div>
        <H2>"""+condition1+' vs. '+condition2+"""</H2>
        <H4>There were <a style="font-size: 20" href="Gene_list.html">"""+str(len(scattersigx))+""" significant genes</a> called with p-adj < """+str(pval_cut)+"""</H4>

        <div class="row">
            <div class="figure" style="float:right">
                <a href="PCA.html">
                  <img src="PCA.png" alt="PCA">
                </a>
                <p style="padding-left:20px; padding-right:20px; text-align:justify"><b>Figure 1: <a target="_blank" href="https://en.wikipedia.org/wiki/Principal_component_analysis">PCA</a></b> - Principal component analysis (PCA) is a statistical procedure that uses an orthogonal transformation to convert a set of observations of possibly correlated variables into a set of values of linearly uncorrelated variables called principal components.</p>
            </div>
            <div class="figure" style="float:left">
                <a href="Heatmap.html">
                  <img src="Heatmap.png" alt="Heatmap">
                </a>
                <p style="padding-left:20px; padding-right:20px; text-align:justify"><b>Figure 2: <a target="_blank" href="https://en.wikipedia.org/wiki/Heat_map">Heatmap</a></b> - A heatmap is is a graphical representation of data where the individual values contained in a matrix are represented as colors. Only significant genes are plotted. Values are normalized to unit scale (mean=0, variance=1).</p>
            </div>
        </div>

        <div class="row">
            <div class="figure" style="float:left">
                <a href="MA_plot.html">
                  <img src="MA_plot.png" alt="MA-plot">
                </a>
                <p style="padding-left:20px; padding-right:20px; text-align:justify"><b>Figure 3: <a target="_blank" href="https://en.wikipedia.org/wiki/MA_plot">MA-Plot</a></b> - An MA plot is an application of a Bland-Altman plot for visual representation of genomic data. The plot visualizes the differences between measurements taken in two samples, by transforming the data onto M (log ratio) and A (mean average) scales, then plotting these values. (Red = significant genes)</p>
            </div>

            <div class="figure" style="float:right">
                <a href="Volcano_plot.html">
                  <img src="Volcano_plot.png" alt="Volcano plot">
                </a>
                <p style="padding-left:20px; padding-right:20px; text-align:justify"><b>Figure 4: <a target="_blank" href="https://en.wikipedia.org/wiki/Volcano_plot_(statistics)">Volcano Plot</a></b> - A volcano plot is a type of scatter-plot that is used to quickly identify changes in large data sets composed of replicate data. It plots significance versus fold-change on the y and x axes, respectively. (Red = significant genes)</p>
            </div>
        </div>


        <div style="width:100%;float:left;padding-top:30px"><hr></div>


        <div class="footer" style="float:right">
            <img style="float:right;width:300px" src="GeneLab_logo.png" alt="GeneLab logo">
            <p><a style="color:#ff8c00; text-decoration:none;" target="_blank" href="https://www.nasa.gov">www.nasa.gov</a></p>
        </div>
        <div class="footer" style="float:left">
            <address style="color:gray;">
                Contact:<br>
                &nbsp;<a href="mailto:jonathan.d.rubin@nasa.gov">Jonathan Rubin</a><br>
                &nbsp;<a href="mailto:daniel.e.mattox@nasa.gov">Daniel Mattox</a><br>
                &nbsp;GeneLab<br>
                &nbsp;Nasa Ames Research Center<br>
                &nbsp;Moffett Blvd<br>
                &nbsp;Mountain View, CA 94035<br>
            </address>
        </div>
    </body>

</html>"""
    

    #Here we write the html report file to the main html file and also to an index.html file which will be linked within the generated plot html files.
    with open(html_main,'w') as outfile:
        outfile.write(main_report)
    with open(os.path.join(html_folder,'index.html'),'w') as outfile:
        outfile.write(main_report)


    #In this section the matplotlib figure is initialized and the MA-plot is created
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        F = plt.figure(figsize=(11,7))
        ax0 = F.add_subplot(111)
        ax0.grid(color='gray')
        x = averageExpression
        y = foldChange
        xy = np.vstack([x,y])
        z = gaussian_kde(xy)(xy)
        idx = np.argsort(z)
        x, y, z = [x[i] for i in idx], [y[i] for i in idx], [z[i] for i in idx]
        scatter = ax0.scatter(x=x,y=y,c=z,s=100,edgecolor="")
        sigscatter = ax0.scatter(scattersigx,scattersigy,c='r',s=100,edgecolor="")
        ax0.tick_params(axis='y', which='both', left='on', right='off', labelleft='on')
        ax0.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='on')
        ax0.set_title("MA Plot", size=25)
        ax0.set_ylabel("Log2 Fold Change ("+condition1+" - "+condition2+")", size=18)
        ax0.set_xlabel("Average Expression", size=18)
        plt.savefig(os.path.join(html_folder,'MA_plot.png'))
        plt.savefig(os.path.join(html_folder,'MA_plot.svg'))


        #This is the bulk of the mpld3 code, basically we're creating 'tooltips' which allows us to have the interactive labels
        tooltip1 = mpld3.plugins.PointHTMLTooltip(sigscatter, labels=[x for x,i in zip(scattersiglabels,range(200))], css=css)

        #This connects our plots together
        mpld3.plugins.connect(F, tooltip1)

        #Here we save both a png version of the plot (non-interactive) and the interactive html version of the plot
        mpld3.save_html(F,os.path.join(html_folder,'MA_plot.html'))

        #Close the figure
        plt.close(F)

        #Here we add a 'back' button to the interactive html plots to go back to the original report.
        with open(os.path.join(html_folder,'MA_plot.html'),'a') as F:
            F.write('<b><a href="index.html">Back</a></b>')




        #Here we create the Volcano plot
        F2 = plt.figure(figsize=(11,7))
        ax1 = F2.add_subplot(111)
        x = foldChange
        y = log10pval
        xy = np.vstack([x,y])
        z = gaussian_kde(xy)(xy)
        idx = np.argsort(z)
        x, y, z = [x[i] for i in idx], [y[i] for i in idx], [z[i] for i in idx]
        volcano = ax1.scatter(x=x,y=y,c=z,s=100,edgecolor="")
        sigvolcano = ax1.scatter(volcanosigx,volcanosigy,c='r',s=100,edgecolor="")
        ax1.tick_params(axis='y', which='both', left='on', right='off', labelleft='on')
        ax1.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='on')
        ax1.set_title("Volcano Plot", size=25)
        ax1.set_ylabel("-Log10 P-value", size=18)
        ax1.set_xlabel("Log2 Fold Change ("+condition1+" - "+condition2+")", size=18)
        ax1.grid(color='gray')
        ax1.set_ylim(bottom=0)
        plt.savefig(os.path.join(html_folder,'Volcano_plot.png'))
        plt.savefig(os.path.join(html_folder,'Volcano_plot.svg'))

        #This is the bulk of the mpld3 code, basically we're creating 'tooltips' which allows us to have the interactive labels
        tooltip2 = mpld3.plugins.PointHTMLTooltip(sigvolcano, labels=[x for x,i in zip(volcanosiglabels,range(200))], css=css)

        #This connects our plots together
        mpld3.plugins.connect(F2,tooltip2)

        #Here we save both a png version of the plot (non-interactive) and the interactive html version of the plot
        mpld3.save_html(F2,os.path.join(html_folder,'Volcano_plot.html'))

        #Close the figure
        plt.close(F2)

        #Here we add a 'back' button to the interactive html plots to go back to the original report.
        with open(os.path.join(html_folder,'Volcano_plot.html'),'a') as F:
            F.write('<b><a href="index.html">Back</a></b>')





    #This section of the code creates an html table of significant genes
    with open(os.path.join(html_folder,'Gene_list.html'),'w') as sigGenes_file:
        sigGenes_file.write("""<!DOCTYPE html>
<html>
    <head>
        <title>List of Significant Genes</title>
        <style>
        * {
            font-family:Arial, Helvetica, sans-serif;
            }
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
    <body style="width:900px; margin:0 auto;">
        <h1>List of Significant Genes</h1>
        <b><a href="index.html">Back</a></b>
    <div>
        <div style="float: middle; width: 100%; padding-top:20px">
            <table> 
                <tr>
                    <th>Gene</th>
                    <th>Average Expression</th> 
                    <th>Log Fold Change</th>
                    <th>Adjusted P-value</th>
                </tr>""")
        for row in sig_genes:
            name,exp,fc,pval = row
            sigGenes_file.write("""
                <tr>
                    <td><a target="_blank" href="https://www.ncbi.nlm.nih.gov/search/?term="""+name+'">'+name+"""</a></td>
                    <td>"""+str(exp)+"""</td>
                    <td>"""+str(fc)+"""</td>
                    <td>"""+str(pval)+"""</td>
                </tr>""")
        sigGenes_file.write("""
            </table>
        </div>
    </div>
    
    </body>
</html>""")


    return [gene for (gene,exp,fc,pval),i in zip(sig_genes,range(50))]


#The below code is taken from https://gist.github.com/peterk87/5505691

def heatmap(x, row_header, column_header, row_method,column_method, row_metric, column_metric,color_gradient, html_folder):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        
        """
        This below code is based in large part on the protype methods:
        http://old.nabble.com/How-to-plot-heatmap-with-matplotlib--td32534593.html
        http://stackoverflow.com/questions/7664826/how-to-get-flat-clustering-corresponding-to-color-clusters-in-the-dendrogram-cre
        x is an m by n ndarray, m observations, n genes
        """
        
        ### Define the color gradient to use based on the provided name
        n = len(x[0]); m = len(x)
        if color_gradient == 'red_white_blue':
            cmap=plt.cm.bwr
        if color_gradient == 'red_black_sky':
            cmap=RedBlackSkyBlue()
        if color_gradient == 'red_black_blue':
            cmap=RedBlackBlue()
        if color_gradient == 'red_black_green':
            cmap=RedBlackGreen()
        if color_gradient == 'yellow_black_blue':
            cmap=YellowBlackBlue()
        if color_gradient == 'seismic':
            cmap=plt.cm.seismic
        if color_gradient == 'green_white_purple':
            cmap=plt.cm.PiYG_r
        if color_gradient == 'coolwarm':
            cmap=plt.cm.coolwarm

        ### Scale the max and min colors so that 0 is white/black
        x = StandardScaler().fit_transform(x)
        vmin = x.min()
        vmax = x.max()
        vmax = max([vmax,abs(vmin)])
        vmin = vmax*-1
        norm = mpl.colors.Normalize(vmin, vmax) ### adjust the max and min to scale these colors

        ### Scale the Matplotlib window size
        default_window_height = 7
        default_window_width = 11
        fig = plt.figure(figsize=(default_window_width,default_window_height)) ### could use m,n to scale here
        color_bar_w = 0.015 ### Sufficient size to show
            
        ## calculate positions for all elements
        # ax1, placement of dendrogram 1, on the left of the heatmap
        [ax1_x, ax1_y, ax1_w, ax1_h] = [0.05,0.22,0.2,0.6]   ### The second value controls the position of the matrix relative to the bottom of the view
        width_between_ax1_axr = -0.004
        height_between_ax1_axc = -0.004 ### distance between the top color bar axis and the matrix
        
        # axr, placement of row side colorbar
        [axr_x, axr_y, axr_w, axr_h] = [0.31,0.1,color_bar_w,0.6] ### second to last controls the width of the side color bar - 0.015 when showing
        axr_x = ax1_x + ax1_w + width_between_ax1_axr
        axr_y = ax1_y; axr_h = ax1_h
        width_between_axr_axm = -0.004

        # axc, placement of column side colorbar
        [axc_x, axc_y, axc_w, axc_h] = [0.4,0.63,0.5,color_bar_w] ### last one controls the height of the top color bar - 0.015 when showing
        axc_x = axr_x + axr_w + width_between_axr_axm
        axc_y = ax1_y + ax1_h + height_between_ax1_axc
        height_between_axc_ax2 = -0.004

        # axm, placement of heatmap for the data matrix
        [axm_x, axm_y, axm_w, axm_h] = [0.4,0.9,2.5,0.5]
        axm_x = axr_x + axr_w + width_between_axr_axm
        axm_y = ax1_y; axm_h = ax1_h
        axm_w = axc_w

        # ax2, placement of dendrogram 2, on the top of the heatmap
        [ax2_x, ax2_y, ax2_w, ax2_h] = [0.3,0.72,0.6,0.15] ### last one controls height of the dendrogram
        ax2_x = axr_x + axr_w + width_between_axr_axm
        ax2_y = ax1_y + ax1_h + height_between_ax1_axc + axc_h + height_between_axc_ax2
        ax2_w = axc_w

        # axcb - placement of the color legend
        [axcb_x, axcb_y, axcb_w, axcb_h] = [0.07,0.88,0.18,0.07]

        # Compute and plot top dendrogram
        if column_method != None:
            d2 = dist.pdist(x.T)
            D2 = dist.squareform(d2)
            ax2 = fig.add_axes([ax2_x, ax2_y, ax2_w, ax2_h], frame_on=False)
            Y2 = sch.linkage(D2, method=column_method, metric=column_metric) ### array-clustering metric - 'average', 'single', 'centroid', 'complete'
            Z2 = sch.dendrogram(Y2)
            ind2 = sch.fcluster(Y2,0.7*max(Y2[:,2]),'distance') ### This is the default behavior of dendrogram
            ax2.set_xticks([]) ### Hides ticks
            ax2.set_yticks([])
        else:
            ind2 = ['NA']*len(column_header) ### Used for exporting the flat cluster data
            
        # Compute and plot left dendrogram.
        if row_method != None:
            d1 = dist.pdist(x)
            D1 = dist.squareform(d1)  # full matrix
            ax1 = fig.add_axes([ax1_x+0.005, ax1_y, ax1_w, ax1_h], frame_on=False) # frame_on may be False
            Y1 = sch.linkage(D1, method=row_method, metric=row_metric) ### gene-clustering metric - 'average', 'single', 'centroid', 'complete'
            Z1 = sch.dendrogram(Y1, orientation='right')
            ind1 = sch.fcluster(Y1,0.7*max(Y1[:,2]),'distance') ### This is the default behavior of dendrogram
            ax1.set_xticks([]) ### Hides ticks
            ax1.set_yticks([])
        else:
            ind1 = ['NA']*len(row_header) ### Used for exporting the flat cluster data
            
        # Plot distance matrix.
        axm = fig.add_axes([axm_x, axm_y, axm_w, axm_h])  # axes for the data matrix
        xt = x
        if column_method != None:
            idx2 = Z2['leaves'] ### apply the clustering for the array-dendrograms to the actual matrix data
            xt = xt[:,idx2]
            # ind2 = ind2[:,idx2] ### reorder the flat cluster to match the order of the leaves the dendrogram
        if row_method != None:
            idx1 = Z1['leaves'] ### apply the clustering for the gene-dendrograms to the actual matrix data
            xt = xt[idx1,:]   # xt is transformed x
            # ind1 = ind1[idx1,:] ### reorder the flat cluster to match the order of the leaves the dendrogram
        ### taken from http://stackoverflow.com/questions/2982929/plotting-results-of-hierarchical-clustering-ontop-of-a-matrix-of-data-in-python/3011894#3011894
        im = axm.matshow(xt, aspect='auto', origin='lower', cmap=cmap, norm=norm) ### norm=norm added to scale coloring of expression with zero = white or black
        axm.set_xticks([]) ### Hides x-ticks
        axm.set_yticks([])

        # Add text
        new_row_header=[]
        new_column_header=[]
        for i in range(x.shape[0]):
            if row_method != None:
                if len(row_header)<200: ### Don't visualize gene associations when more than 100 rows
                    if len(row_header) < 20:
                        fontsize=15
                    else:
                        fontsize=200/len(row_header)
                    axm.text(x.shape[1]-0.5, i-0.1, '  '+row_header[idx1[i]],fontsize=fontsize)
                new_row_header.append(row_header[idx1[i]])
            else:
                if len(row_header)<200: ### Don't visualize gene associations when more than 100 rows
                    if len(row_header) < 20:
                        fontsize=8
                    else:
                        fontsize=200/len(row_header)
                    axm.text(x.shape[1]-0.5, i-0.1, '  '+row_header[i],fontsize=fontsize) ### When not clustering rows
                new_row_header.append(row_header[i])
        for i in range(x.shape[1]):
            if column_method != None:
                axm.text(i, -0.55, ' '+column_header[idx2[i]], rotation=315, verticalalignment="top") # rotation could also be degrees
                new_column_header.append(column_header[idx2[i]])
            else: ### When not clustering columns
                axm.text(i, -0.55, ' '+column_header[i], rotation=315, verticalalignment="top")
                new_column_header.append(column_header[i])

        # Plot colside colors
        # axc --> axes for column side colorbar
        # if column_method != None:
        #     axc = fig.add_axes([axc_x, axc_y, axc_w, axc_h])  # axes for column side colorbar
        #     cmap_c = mpl.colors.ListedColormap(['r', 'g', 'b', 'y', 'w', 'k', 'm'])
        #     dc = np.array(ind2, dtype=int)
        #     dc.shape = (1,len(ind2)) 
        #     im_c = axc.matshow(dc, aspect='auto', origin='lower', cmap=cmap_c)
        #     axc.set_xticks([]) ### Hides ticks
        #     axc.set_yticks([])
        
        # Plot rowside colors
        # axr --> axes for row side colorbar
        # if row_method != None:
        #     axr = fig.add_axes([axr_x+0.005, axr_y, axr_w-0.005, axr_h])  # axes for column side colorbar
        #     dr = np.array(ind1, dtype=int)
        #     dr.shape = (len(ind1),1)
        #     cmap_r = mpl.colors.ListedColormap(['r', 'g', 'b', 'y', 'w', 'k', 'm'])
        #     im_r = axr.matshow(dr, aspect='auto', origin='lower', cmap=cmap_r)
        #     axr.set_xticks([]) ### Hides ticks
        #     axr.set_yticks([])

        # Plot color legend
        axcb = fig.add_axes([axcb_x, axcb_y, axcb_w, axcb_h], frame_on=False)  # axes for colorbar
        cb = mpl.colorbar.ColorbarBase(axcb, cmap=cmap, norm=norm, orientation='horizontal')
        cb.set_ticks([vmin,0,vmax])
        axcb.set_title("Normalized Expression")

        #Save figures
        plt.savefig(os.path.join(html_folder,'Heatmap.png'), bbox_inches='tight', dpi=300)
        plt.savefig(os.path.join(html_folder,'Heatmap.svg'))

        #Create html output file
        html_file = """<!DOCTYPE html>
<html>
    <head>
        <title>Heatmap</title>
        <style>
        * {
            font-family:Arial, Helvetica, sans-serif;
            }
        </style>
    </head>
    <body style="width:100%">
        <img style="width:100%" src="Heatmap.png" alt="Heatmap">
        <a href="index.html"><b>Back</b></a>
    </body>
</html>
"""
        with open(os.path.join(html_folder,'Heatmap.html'),'w') as F:
            F.write(html_file)



################# Create Custom Color Gradients #################
#http://matplotlib.sourceforge.net/examples/pylab_examples/custom_cmap.html

def RedBlackSkyBlue():
    cdict = {'red':   ((0.0, 0.0, 0.0),
                       (0.5, 0.0, 0.1),
                       (1.0, 1.0, 1.0)),
    
             'green': ((0.0, 0.0, 0.9),
                       (0.5, 0.1, 0.0),
                       (1.0, 0.0, 0.0)),
    
             'blue':  ((0.0, 0.0, 1.0),
                       (0.5, 0.1, 0.0),
                       (1.0, 0.0, 0.0))
            }

    my_cmap = mpl.colors.LinearSegmentedColormap('my_colormap',cdict,256)
    return my_cmap

def RedBlackBlue():
    cdict = {'red':   ((0.0, 0.0, 0.0),
                       (0.5, 0.0, 0.1),
                       (1.0, 1.0, 1.0)),

             'green': ((0.0, 0.0, 0.0),
                       (1.0, 0.0, 0.0)),
    
             'blue':  ((0.0, 0.0, 1.0),
                       (0.5, 0.1, 0.0),
                       (1.0, 0.0, 0.0))
            }

    my_cmap = mpl.colors.LinearSegmentedColormap('my_colormap',cdict,256)
    return my_cmap

def RedBlackGreen():
    cdict = {'red':   ((0.0, 0.0, 0.0),
                       (0.5, 0.0, 0.1),
                       (1.0, 1.0, 1.0)),
    
             'blue': ((0.0, 0.0, 0.0),
                       (1.0, 0.0, 0.0)),
    
             'green':  ((0.0, 0.0, 1.0),
                       (0.5, 0.1, 0.0),
                       (1.0, 0.0, 0.0))
            }
    
    my_cmap = mpl.colors.LinearSegmentedColormap('my_colormap',cdict,256)
    return my_cmap

def YellowBlackBlue():
    cdict = {'red':   ((0.0, 0.0, 0.0),
                       (0.5, 0.0, 0.1),
                       (1.0, 1.0, 1.0)),
    
             'green': ((0.0, 0.0, 0.8),
                       (0.5, 0.1, 0.0),
                       (1.0, 1.0, 1.0)),
    
             'blue':  ((0.0, 0.0, 1.0),
                       (0.5, 0.1, 0.0),
                       (1.0, 0.0, 0.0))
            }
    ### yellow is created by adding y = 1 to RedBlackSkyBlue green last tuple
    ### modulate between blue and cyan using the last y var in the first green tuple
    my_cmap = mpl.colors.LinearSegmentedColormap('my_colormap',cdict,256)
    return my_cmap

def getColorRange(x):
    """ Determines the range of colors, centered at zero, for normalizing cmap """
    vmax=x.max()
    vmin=x.min()
    if vmax<0 and vmin<0: direction = 'negative'
    elif vmax>0 and vmin>0: direction = 'positive'
    else: direction = 'both'
    if direction == 'both':
        vmax = max([vmax,abs(vmin)])
        vmin = -1*vmax
        return vmax,vmin
    else:
        return vmax,vmin


if __name__ == "__main__":
    #===========TEST HEATMAP===========
    # row_method = 'average'
    # column_method = 'single'
    # row_metric = 'cityblock' #cosine
    # column_metric = 'euclidean'
    # color_gradient = 'red_black_green'
    # row_header = ['gene_' + str(i) for i in range(5)]
    # column_header = ['sample_' + str(i) for i in range(5)]
    # M = [[j for j in np.arange(-1,1,2.0/5.0)] for i in range(5)]
    # x = np.squeeze(np.asarray(M))
    # heatmap(x, row_header, column_header, row_method,
    #         column_method, row_metric, column_metric,
    #         color_gradient, '/Users/jonathanrubin/Google_Drive/NASA/home/galaxy_test/')





    # #==========TEST LIMMA==========
    # counts_table='/Users/jonathanrubin/Google_Drive/NASA/home/processed_GLDS/GLDS-4/microarray/GLDS-4_microarray_normalized-annotated.txt'
    # metadata='/Users/jonathanrubin/Google_Drive/NASA/home/processed_GLDS/GLDS-4/metadata/s_GLDS-4_microarray_metadata.txt'
    # condition1='Spaceflight'
    # condition2='ground'
    # outliers='GSM458594'
    # html_folder='/Users/jonathanrubin/Google_Drive/NASA/'
    # limma_output = limma_differential(counts_table,metadata,condition1,condition2,outliers,html_folder)
    # limma_output = limma_output.split('\n')
    # group1_index = [i+1 for i in range(len(limma_output)) if 'group 1:' in limma_output[i]][0]
    # group2_index = [i+1 for i in range(len(limma_output)) if 'group 2:' in limma_output[i]][0]
    # column_header = limma_output[group1_index].split() + limma_output[group2_index].split()
    # print column_header



    #=========TEST GO=========
    url = 'http://pantherdb.org/webservices/go/overrep.jsp'
    form_data = {
        'input': 'NM_001013370,NM_001040691',
        'species': "HUMAN",
        'ontology': "biological_process",
        'submit': 'submit',
    }
    response = requests.post(url, data=form_data)
    with open('test.html','w') as outfile:
        outfile.write(response.content)
    print response.content



    #=========TEST ALL=========
    # counts_table='/Users/jonathanrubin/Google_Drive/NASA/home/processed_GLDS/GLDS-4/microarray/GLDS-4_microarray_normalized-annotated.txt'
    # metadata='/Users/jonathanrubin/Google_Drive/NASA/home/processed_GLDS/GLDS-4/metadata/s_GLDS-4_microarray_metadata.txt'
    # condition1='Spaceflight'
    # condition2='ground'
    # outliers='None'
    # html_folder='/Users/jonathanrubin/Google_Drive/NASA/home/galaxy_test/'
    # html_main='/Users/jonathanrubin/Google_Drive/NASA/home/galaxy_test/html_main.html'
    # padj_cutoff='0.1'
    # run(counts_table,metadata,condition1,condition2,padj_cutoff,outliers,html_main,html_folder)



