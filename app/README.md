# Bokeh app for Data Visualization of HBCGM output

This is a bokeh application for exploring HBCGM results interactively and download results.

## Run 

debug  
```shell
## --dev autoreload files 
bokeh serve --show app --allow-websocket-origin=peltz-app-03:5006 --dev app/*.py --log-level=debug
```

deployment  
```
bokeh serve --show app --allow-websocket-origin=peltz-app-03:5006
```

view at: peltz-app-03:5006

**Note**:
HBCGM results must be prioritzation by MeSH terms (GNN model) first, then run this app


## Dev
see guide [here](https://docs.bokeh.org/en/2.4.1/docs/user_guide/server.html)