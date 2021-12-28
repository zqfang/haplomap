
import os
import numpy as np
import pandas as pd
from bokeh.plotting import figure, curdoc
from bokeh.models import ColumnDataSource, TableColumn, DateFormatter, DataTable, HTMLTemplateFormatter, CellFormatter
from bokeh.models import ColorBar, LinearColorMapper, LabelSet, Legend
from bokeh.models import FixedTicker, PrintfTickFormatter
from bokeh.models.glyphs import Patches
from bokeh.models.widgets import Select, TextInput, Dropdown, AutocompleteInput, Div
from bokeh.layouts import column, row

from bokeh.palettes import Category10
from bokeh.transform import factor_cmap, linear_cmap, factor_mark
from bokeh.core.enums import MarkerType
from .helpers import gene_expr_order, codon_flag, load_ghmap, get_color, get_expr2
######### global variable
DATA_DIR = "/data/bases/fangzq/Pubmed/test_cases/TEST" # RUN_000.results.txt
gene_expr_order = []
mesh_terms = {}
fcmap = factor_cmap('impact', palette=Category10[6], factors=['Synonymous','Non-Synonymous','Splicing', 'Stop', 'Non-Coding'])
#fmark = factor_mark('impact', list(MarkerType), ['Synonymous','Non-Synonymous','Splicing', 'Stop', 'Non-Coding'])

####### bokeh #######

### setup widgets
# dataset id
dataset = AutocompleteInput(completions=['MPD_1560','RUN_000'], 
                           title="Enter Dataset:", value="RUN_000", width=200,)
# mesh terms options
meshid = Select(title="MeSH:", value="", options=[], width=200,) # need to update dynamically
# others
symbol = TextInput(value = '', title = "GeneName:", width=200,)
pval = TextInput(value = '', title = "Pvalue:",width=200,)
codon = TextInput(value = '', title = "CodonFlag",width=200,)
# message box
message = Div(text="""Gene Expression: """, width=200, height=100)
# gene expression pattern
exprs = Div(text="""Gene Expression: """, width=200, height=800)
### setup data and plots
## data
source = ColumnDataSource(pd.DataFrame())
source_bar = ColumnDataSource(data=dict(strains=[],traits=[], colors=[]))

source_scatter = ColumnDataSource(data=dict(logPvalue=[], 
                                            MESH=[],
                                            GeneName=[],
                                            impact=[]))

## bar plot

p = figure(plot_width=550, plot_height=400, # x_range=strains, 
           title="Dataset", 
           tools='reset,lasso_select,save', output_backend="svg")
p.toolbar.logo = None
p.vbar(x='strains', top='traits', source=source_bar, line_width=0, fill_color='colors', width=0.7)
p.xgrid.grid_line_color = None
#p.y_range.start = 0
p.xaxis.axis_label = "Strains"
p.yaxis.axis_label = "Values"
p.xaxis.major_label_orientation = np.pi/4
# or alternatively:
#p.xaxis.major_label_orientation = "vertical"


## Scatter plot
#fcmap = factor_cmap('colors', palette=Category10[6], factors=['-1','0','1','2','3', '4'], start= -1, end=4)
#fcmap = factor_cmap('impact', palette=Category10[6], factors=list(codon_flag.values()))
# mapper = linear_cmap(field_name='MESH', palette="Viridis256", low=0, high=1)
# color_bar = ColorBar(color_mapper=mapper['transform'],  
#                      width=10, 
#                      major_label_text_font_size="10pt",
#                      location=(0, 0))
# hover tooltips
TOOLTIPS = [
    ("logPval, LitScore", "($x, $y)"),
    ("GeneName", "@GeneName"),
    ("Impact", "@impact")]

t2 = figure(plot_width=500, plot_height=400, 
            tools="hover,pan,box_zoom,reset,lasso_select,save",
            output_backend="svg", 
            tooltips=TOOLTIPS)

t2.add_layout(Legend(), 'below') # put legend outside
t2.toolbar.logo = None
t2.xaxis.axis_label = "Genetic:  - log10 Pvalue"#r"$$\text{Genetic} - \log_{10} \text{Pvalue}$$"
t2.yaxis.axis_label = "Literature Score" #"MeSH"
t2.scatter('logPvalue', 'MESH', 
            size=8, 
            source=source_scatter, 
            fill_color=fcmap,# mapper,
            #marker=fmark,
            legend_field="impact",
            line_color=None, 
            selection_color="orange", #alpha=0.6, 
            nonselection_alpha=0.1, selection_alpha=0.4)
t2.title.text =  "MeSH Terms"
t2.legend.orientation = "horizontal"
#t2.legend.location = "bottom_right"




## Datatable
columns = ['GeneName', 'CodonFlag','Haplotype', 'Pvalue', 'EffectSize', 'FDR',
           'PopPvalue', 'PopFDR', 'Chr', 'ChrStart', 'ChrEnd', 'LitScore'] 
columns = [ TableColumn(field=c, title=c, formatter=HTMLTemplateFormatter() 
                        if c == 'Haplotype' else CellFormatter()) for c in columns ] # skip index
myTable = DataTable(source=source, columns=columns, width =1000, height = 800, editable = False, autosize_mode="fit_viewport")
# show(myTable)


### setup callbacks
def source_update(attr, old, new):
    try:
        selected_index = source.selected.indices[0]
        symbol.value = str(source.data["GeneName"][selected_index])
        pval.value = str(source.data["Pvalue"][selected_index])
        codon.value = codon_flag.get(str(source.data['CodonFlag'][selected_index]))
        source_bar.data.update(colors=get_color(source.data["Pattern"][selected_index]))
        ep, ep2 = get_expr2(source.data['GeneExprMap'][selected_index], gene_expr_order) 
        exprs.text = "GeneExpression:<br>"+ep2
        message.text = "GeneExpression: <br>"+ep
        
    except IndexError:
        pass
    
def mesh_update(attr, old, new):
    mesh = meshid.value
    t2.title.text = mesh_terms.get(mesh) + " : " + mesh
    mesh = "MeSH_" + mesh_terms.get(mesh)
    mesh_columns = [m for m in source.data.keys() if m.startswith("MeSH") ]
    if mesh not in mesh_columns:
        mesh = mesh_columns[0]
        message.text = f"<p>Sorry, input MeSH not found! <br> Selected: {mesh} </p>"
        return
    source_scatter.data.update(MESH=source.data[mesh])
    myTable.source.data.update(LitScore=source.data[mesh]) # datatable
    
     
def data_update(attr, old, new):
    ds = dataset.value
    DATASET = os.path.join(DATA_DIR, '%s.results.txt' % ds)
    if not os.path.exists(DATASET):
        message.text = f"<p> Error: <br> Dataset {ds} not found <br> Please Input a new dataset name.</p>"
        return
    # update new data
    df, headers = load_ghmap(DATASET)
    df = df[df.CodonFlag>=0]
    global gene_expr_order
    global mesh_terms
    global codon_flag

    gene_expr_order = headers[2]
    # update mesh, bar, scatter
    mesh_columns = [m for m in df.columns if m.startswith("MeSH") ]
    dataset_name, codon_flag, gene_expr_order, strains, traits, mesh_terms = headers[:6]
    meshid.options = list(mesh_terms.keys())
    
    x_range = list(range(0, len(strains)))
    source.data = df
    myTable.source.data = df #.to_dict(orient='list') # datatable
    myTable.source.data.update(LitScore=df.loc[:, mesh_columns[0]].to_list())
    source_bar.data.update(strains=x_range, traits=traits, colors=['#1f77b4']*len(strains))
    source_scatter.data.update(logPvalue=df.loc[:,'logPvalue'].to_list(), 
                               MESH=df.loc[:, mesh_columns[0]].to_list(),
                               GeneName=df.loc[:,'GeneName'].to_list(),
                               impact=df.loc[:,'Impact'].to_list())
    p.xaxis.ticker = FixedTicker(ticks=x_range)
    p.xaxis.major_label_overrides = {k: str(v) for k, v in zip(x_range, strains)}
    p.title.text = "Dataset: "+dataset_name[0]
    t2.title.text = mesh_columns[0].split("_")[-1] + ": "+ list(mesh_terms.keys())[0]
    message.text = f"<p> Loaded dataset {ds}</p>"
    
    
## init data and plot
data_update(None, None, None)

## set change 
dataset.on_change('value', data_update)
source.selected.on_change('indices', source_update)
meshid.on_change('value', mesh_update)
    

## set up layout
sb = column(dataset, meshid, symbol, pval, codon, message, exprs)
figs = row(p, t2)
series_t = column(figs, myTable)
layout = row(series_t, sb)
    
# source.selected.on_change('indices', function_source)
curdoc().add_root(layout)
curdoc().title = "HBCGM"