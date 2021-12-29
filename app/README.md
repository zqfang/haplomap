# Bokeh app for Data Visualization of HBCGM output

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


## Dev
see guide [here](https://docs.bokeh.org/en/2.4.1/docs/user_guide/server.html)