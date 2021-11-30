from skyfield.api import load, EarthSatellite

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.widgets import Slider


ts = load.timescale()
global time_of_interest
time_of_interest = ts.utc(2021,11,15,2,53,44)
TWOPI = 2 * np.pi
HALFPI = np.pi / 2
def mark_satellite(sat,plot_time):
    sat_ra,sat_dec,_ = sat.at(plot_time).radec()
    sat_ra = sat_ra.radians
    sat_dec = sat_dec.radians
    return ax.scatter(sat_ra,sat_dec,s=200,color='r',marker='+')

def slider_update(slider_position):
    global time_of_interest
    #Find the new time
    desired_time = list(time_of_interest.utc)
    desired_time[4] = slider_position
    time_of_interest = ts.utc(*desired_time)
    date_slider.valtext.set_text(time_of_interest.utc_strftime("\n%Y-%m-%d\n%H:%M:%S UTC"))
    update_plot(time_of_interest)#And now replot everything.

def update_plot(plot_time):
    global plotted_objects
    for i in plotted_objects:
        i.remove()
    #These functions are all time-dependent and draw everything on the background of stars.
    plotted_objects = [mark_satellite(satellite,plot_time) for satellite in sats]

#Set up the plot frame.
fig,ax = plt.subplots()
plt.xlabel("Right ascension")
xtick_locations = np.linspace(0,TWOPI,25)

plt.ylabel("Declination")
ytick_locations = np.linspace(-HALFPI,HALFPI,13)
ax.set_axisbelow(True)
ax.invert_xaxis()

#Change the axis labels to be standard units - Hours for RA, Degrees for DEC
ax.set_axisbelow(True)
#Adjust things to add the slider
plt.subplots_adjust(bottom=0.25)
#Add the controls to scroll through time
date_slider = Slider(
    ax=plt.axes([0.2, 0.1, 0.6, 0.03],facecolor='k'),
    label="Time Slider",
    valmin=0,
    valmax=59,
    valinit = time_of_interest.utc[4],#Start at the current minute in our hour
    initcolor='none')#Remove the initial value marker (irrelevant when hour changes)
date_slider.on_changed(slider_update)
date_slider.valtext.set_text(time_of_interest.utc_strftime("\n%Y-%m-%d\n%H:%M:%S UTC"))
global plotted_objects
plotted_objects = []
ax.set_title("")
with open('corrected_tles.txt') as tlefile:
    lines = tlefile.readlines()
    #split file into groups of 3 lines
    triplets= [lines[i:i+3] for i in range(0,len(lines),3)]
    #Leave out the first two objects (the main sat, and the rocket booster)
    triplets = triplets[2:]
sats = [EarthSatellite(t[1],t[2],t[0],ts) for t in triplets]
update_plot(time_of_interest)
plt.show()
