
import os
import numpy as np
import pandas as pd
from bokeh.plotting import figure, curdoc
from bokeh.models import ColumnDataSource, TableColumn, DateFormatter, DataTable, HTMLTemplateFormatter, CellFormatter
from bokeh.models import ColorBar, LinearColorMapper, LabelSet, Legend
from bokeh.models import FixedTicker, RangeSlider, CDSView, BooleanFilter, GroupFilter, CustomJS
from bokeh.models.glyphs import Patches
from bokeh.models.widgets import Select, TextInput, Dropdown, AutocompleteInput, Div, Button
from bokeh.layouts import column, row

from bokeh.palettes import Category10
from bokeh.transform import factor_cmap, linear_cmap, factor_mark
from bokeh.core.enums import MarkerType
from .helpers import gene_expr_order, codon_flag, mesh_terms
from .helpers import load_ghmap, get_color, get_expr, get_datasets
######### global variable
DATA_DIR = "/data/bases/shared/haplomap/MPD_MeSH_Indel"#"/home/fangzq/github/HBCGM/example/PeltzData/" # RUN_000.results.txt


fcmap = factor_cmap('impact', palette=Category10[6], factors=['0','1','2','3','-1'])#['Synonymous','Non-Synonymous','Splicing', 'Stop', 'Non-Coding'])
#fmark = factor_mark('impact', list(MarkerType), ['Synonymous','Non-Synonymous','Splicing', 'Stop', 'Non-Coding'])

####### bokeh #######
### setup data and plots
## data
source = ColumnDataSource(pd.DataFrame())
source_bar = ColumnDataSource(data=dict(strains=[],traits=[], colors=[]))

### setup widgets
# dataset id
DATASETS = get_datasets(DATA_DIR)
dataset = AutocompleteInput(completions=DATASETS, 
                           title="Enter Dataset:", value=DATASETS[0], width=550,)
# mesh terms options
meshid = Select(title="Select MeSH:", value="", options=[], width=300,) # need to update dynamically
# groupfilter options
impact = Select(title="Select Group:", value="", options=[], width=200,) # need to update dynamically
# others
symbol = TextInput(value = '', title = "Gene Name:", width=300,)
pval = TextInput(value = '', title = "P-value:",width=300,)
litscore = TextInput(value = '', title = "MeSH Score:",width=300,)
codon = TextInput(value = '', title = "Impact (Codon Flag)",width=300,)
# message box
message = Div(text="""<h3> Gene Expression: </h3>""", width=300, height=100)
# gene expression pattern
exprs = Div(text="""<h3> Gene Expression: </h3>""", width=300, height=800)

slider = RangeSlider(title="LitScore Range", start=0.0, end=1.0, value=(0.5, 1.0), step=0.01)
# data view
view = CDSView(source=source, filters=[BooleanFilter(), GroupFilter()])

## Datatable
columns = ['GeneName', 'CodonFlag','Haplotype','EffectSize', 'Pvalue', 'FDR',
           'PopPvalue', 'PopFDR', 'Chr', 'ChrStart', 'ChrEnd', 'LitScore'] 
columns = [ TableColumn(field=c, title=c, formatter=HTMLTemplateFormatter() 
                        if c == 'Haplotype' else CellFormatter()) for c in columns ] # skip index                       
myTable = DataTable(source=source, columns=columns, width =1000, height = 600, index_position=0,
                    editable = False, autosize_mode="fit_viewport", view=view, name="DataTable")

# download
button = Button(label="Download Table", button_type="success")


## bar plot

bar = figure(plot_width=550, plot_height=500, # x_range=strains, 
           title="Dataset", 
           #toolbar_location="right",
           tools='pan,reset,lasso_select,save', output_backend="svg", name="Bar")
bar.toolbar.logo = None
bar.vbar(x='strains', top='traits', source=source_bar, line_width=0, fill_color='colors', width=0.7)
bar.xgrid.grid_line_color = None
#p.y_range.start = 0
bar.xaxis.axis_label = "Strains"
bar.yaxis.axis_label = "Values"
bar.xaxis.major_label_orientation = np.pi/4
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

sca = figure(plot_width=550, plot_height=500, 
            tools="pan,box_zoom,reset,lasso_select,save",
            output_backend="svg",  
            #active_drag="lasso_select",
            #toolbar_location="right",
            tooltips=TOOLTIPS, name='Scatter')

sca.add_layout(Legend(), 'above') # put legend outside
sca.legend.orientation = "horizontal"
sca.toolbar.logo = None
sca.xaxis.axis_label = "Genetic:  - log10 Pvalue"# #r"$$\text{Genetic} - \log_{10} \text{Pvalue}$$" # latex
sca.yaxis.axis_label = "Literature Score" #"MeSH"
sca.scatter('logPvalue', 'LitScore', 
            size=8, 
            source=source, 
            #fill_color=fcmap,# mapper,
            #marker=fmark,
            color='CodonColor',
            legend_field="Impact",
            line_color=None, 
            selection_color="#2F0303", #alpha=0.6, 
            view=view,
            nonselection_alpha=0.1, )#selection_alpha=0.4)
sca.title.text =  "MeSH Terms"
#sca.legend.location = "bottom_right"

                   
### setup callbacks
def source_update(attr, old, new):
    try:
        selected_index = source.selected.indices[0]
        symbol.value = str(source.data["GeneName"][selected_index])
        pval.value = str(source.data["Pvalue"][selected_index])
        litscore.value = str(source.data["LitScore"][selected_index])
        codon.value = codon_flag.get(str(source.data['CodonFlag'][selected_index]))
        source_bar.data.update(colors=get_color(source.data["Pattern"][selected_index]))
        ep, ep2 = get_expr(source.data['GeneExprMap'][selected_index], gene_expr_order) 
        exprs.text = "<h3> Gene Expression: </h3>"+ep2
        message.text = "<h3> Gene Expression: </h3>"+ep
        
    except IndexError:
        pass
    
def mesh_update(attr, old, new):
    mesh = meshid.value
    sca.title.text = mesh_terms.get(mesh) + " : " + mesh
    mesh = "MeSH_" + mesh_terms.get(mesh)
    mesh_columns = [m for m in source.data.keys() if m.startswith("MeSH") ]
    if mesh not in mesh_columns:
        mesh = mesh_columns[0]
        message.text = f"<p>Sorry, input MeSH not found! <br> Selected: {mesh} </p>"
        return
    source.data.update(LitScore=source.data[mesh]) # datatable
    myTable.source.data.update(LitScore=source.data[mesh]) # datatable
    
     
def data_update(attr, old, new):
    global gene_expr_order
    global mesh_terms
    global codon_flag

    ds = dataset.value
    DATASET = os.path.join(DATA_DIR, '%s.results.mesh.txt' % ds)
    if not os.path.exists(DATASET):
        message.text = f"<p> Error: <br> Dataset {ds} not found <br> Please Input a new dataset name.</p>"
        return
    # update new data
    df, headers = load_ghmap(DATASET)
    df = df[df.CodonFlag>=0]
    if df.columns.str.startswith("Pop").sum() == 0: 
        # kick out 'PopPvalue', 'PopFDR',
        myTable.columns = [ columns[i] for i in range(len(columns)) if  i not in [6,7] ]   
    else:
        myTable.columns = columns

    # update mesh, bar, scatter
    mesh_columns = [m for m in df.columns if m.startswith("MeSH") ]
    dataset_name, codon_flag, gene_expr_order, strains, traits, mesh_terms = headers[:6]
    # if (dataset_name[0].lower().find("indel") != -1) or (dataset_name[0].lower().find("_sv") != -1):
    #     codon_flag = {'0':'Low','1':'Moderate','2':'High', '-1':'Modifier'}
    x_range = list(range(0, len(strains)))
    # updates
    message.text = f"<p> Loaded dataset {ds}</p>"
    meshid.options = list(mesh_terms.keys())
    impact.options = list(codon_flag.values())
    # need to reset filters again
    myTable.disabled = True 
    view.filters = [BooleanFilter(), GroupFilter()]
    myTable.disabled = False

    source.data = df
    source.data.update(LitScore=df.loc[:, mesh_columns[0]].to_list())
    myTable.source = source # datatable
    source_bar.data.update(strains=x_range, traits=traits, colors=['#A19B9B']*len(strains))

    bar.xaxis.ticker = FixedTicker(ticks=x_range)
    bar.xaxis.major_label_overrides = {k: str(v) for k, v in zip(x_range, strains)}
    bar.title.text = "Dataset: "+dataset_name[0]
    sca.title.text = mesh_columns[0].split("_")[-1] + ": "+ list(mesh_terms.keys())[0]
    # slider_update(attr, old, new)

def slider_update(attr, old, new):
    low, high = slider.value
    # sca.y_range.start = low
    # sca.y_range.end = high
    # view.filters[1] = GroupFilter()
    myTable.disabled = True # FIXME: froze table first, update filter, then unfroze. This trick updates datatable
    booleans = [True if low <y_val < high else False for y_val in source.data['LitScore']]
    view.filters[0] = BooleanFilter(booleans)
    myTable.disabled = False

def impact_update(attr, old, new):
    imp = impact.value   
    myTable.disabled = True 
    view.filters[1] = GroupFilter(column_name='Impact', group=imp)
    myTable.disabled = False
    
## init data and plot
data_update(None, None, None)

## set change 
dataset.on_change('value', data_update)
# datatable
source.selected.on_change('indices', source_update)
meshid.on_change('value', mesh_update)
slider.on_change('value', slider_update)
impact.on_change('value', impact_update)
button.js_on_click(CustomJS(args=dict(source=myTable.source),
                            code=open(os.path.join(os.path.dirname(__file__), "static/js/download.js")).read()))

## set up layout
inputs = row(dataset, meshid, impact, )
# curdoc.add_root(inputs)
figs = row(bar, sca)
#figs = row(robj1, robj2)
figs2 = column(figs, myTable)
# curdoc.add_root(figs2)
robj4 = column(button, symbol, pval, codon, slider, message, exprs)
o = row(figs2, robj4)
layout = column(inputs, o)


    
# source.selected.on_change('indices', function_source)
curdoc().add_root(layout)
curdoc().title = "HBCGM Dashboard"