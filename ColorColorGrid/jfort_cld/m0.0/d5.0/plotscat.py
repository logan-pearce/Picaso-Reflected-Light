import pandas as pd
from bokeh.plotting import figure, show, output_file 
from bokeh.layouts import column,row
import numpy as np
from bokeh.palettes import magma as colfun1
from bokeh.palettes import viridis as colfun2
from bokeh.palettes import Spectral11
from bokeh.models import HoverTool
from bokeh.models import LinearColorMapper, LogTicker, ColorBar,LogColorMapper

cols = colfun1(200)
color_mapper = LinearColorMapper(palette=cols, low=0, high=1)


dat01 = pd.read_csv('m0.0x_rfacv0.5-nc_tint150-f0.01-d5.0.cld', header=None,delim_whitespace=True)
dat1 = pd.read_csv('m0.0x_rfacv0.5-nc_tint150-f0.3-d5.0.cld', header=None,delim_whitespace=True)
dat6 = pd.read_csv('m0.0x_rfacv0.5-nc_tint150-f1-d5.0.cld', header=None,delim_whitespace=True)

scat01 = np.flip(np.reshape(dat01[4],(60,196)),0)#[0:10,:]
scat1 = np.flip(np.reshape(dat1[4],(60,196)),0)#[0:10,:]
scat6 = np.flip(np.reshape(dat6[4],(60,196)),0)#[0:10,:]

#scat01 = np.flip(np.reshape(np.log10(dat01[0]+1e-60),(60,196)),0)
#scat1 = np.flip(np.reshape(np.log10(dat1[0]+1e-60),(60,196)),0)
#scat6 = np.flip(np.reshape(np.log10(dat6[0]+1e-60),(60,196)),0)
xr, yr = scat01.shape

f01a = figure(x_range=[150, yr], y_range=[0,xr],
                         x_axis_label='Wavelength Grid', y_axis_label='Pressure Grid, TOA ->',
                         title="Scattering Fsed = 0.01",
                        plot_width=300, plot_height=300)

f1a = figure(x_range=[150, yr], y_range=[0,xr],
                         x_axis_label='Wavelength Grid', y_axis_label='Pressure Grid, TOA ->',
                         title="Scattering Fsed = 0.3",
                        plot_width=300, plot_height=300)


f6a =figure(x_range=[150, yr], y_range=[0,xr],
                         x_axis_label='Wavelength Grid', y_axis_label='Pressure Grid, TOA ->',
                         title="Scattering Fsed = 1",
                        plot_width=300, plot_height=300)

f01a.image(image=[scat01],  color_mapper=color_mapper, x=0,y=0,dh=xr,dw = yr)
f1a.image(image=[scat1],  color_mapper=color_mapper, x=0,y=0,dh=xr,dw = yr)
f6a.image(image=[scat6], color_mapper=color_mapper, x=0,y=0,dh=xr,dw = yr)

color_bar = ColorBar(color_mapper=color_mapper, #ticker=LogTicker(),
                     label_standoff=12, border_line_color=None, location=(0,0))

f01a.add_layout(color_bar, 'right')
output_file('scattering.html')
show(column(f01a,f1a,f6a))


#PLOT OPD


scat01 = np.flip(np.reshape(dat01[2]+1e-60,(60,196)),0)
scat1 = np.flip(np.reshape(dat1[2]+1e-60,(60,196)),0)
scat6 = np.flip(np.reshape(dat6[2]+1e-60,(60,196)),0)
xr, yr = scat01.shape
cols = colfun2(200)[::-1]
color_mapper = LogColorMapper(palette=cols, low=1e-3, high=10)


f01 = figure(x_range=[150, yr], y_range=[0,xr],
                         x_axis_label='Wavelength Grid', y_axis_label='Pressure Grid, TOA ->',
                         title="Optical Depth Fsed = 0.01",
                        plot_width=300, plot_height=300)

f1 = figure(x_range=[150, yr], y_range=[0,xr],
                         x_axis_label='Wavelength Grid', y_axis_label='Pressure Grid, TOA ->',
                         title="Optical Depth Fsed = 0.3",
                        plot_width=300, plot_height=300)


f6 =figure(x_range=[150, yr], y_range=[0,xr],
                         x_axis_label='Wavelength Grid', y_axis_label='Pressure Grid, TOA ->',
                         title="Optical Depth Fsed = 1",
                        plot_width=300, plot_height=300)

f01.image(image=[scat01],  color_mapper=color_mapper, x=0,y=0,dh=xr,dw = yr)
f1.image(image=[scat1],  color_mapper=color_mapper, x=0,y=0,dh=xr,dw = yr)
f6.image(image=[scat6], color_mapper=color_mapper, x=0,y=0,dh=xr,dw = yr)

color_bar = ColorBar(color_mapper=color_mapper, ticker=LogTicker(),
                     label_standoff=12, border_line_color=None, location=(0,0))

f01.add_layout(color_bar, 'right')
output_file('opd.html')
show(row(column(f01a,f1a,f6a),column(f01,f1,f6)))

