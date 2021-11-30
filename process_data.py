from skyfield.api import load, EarthSatellite
import numpy as np

with open('tles.txt') as tlefile:
    lines = tlefile.readlines()
    #split file into groups of 3 lines
    triplets= [lines[i:i+3] for i in range(0,len(lines),3)]
    #Leave out the first two objects (the main sat, and the rocket booster)
    triplets = triplets[2:]

ts = load.timescale()

objects = [EarthSatellite(t[1],t[2],t[0],ts) for t in triplets]

#Identify the location and time of collision. We'll do this by taking the debris
#and, for each timestep, find the center of mass (assuming uniform masses, ugh).
#Then find the moment when the least-squares deviation is minimized
#(debris cloud is smallest).
#We know it was November 15, now let's find out when.
#Create a time object of every second on that day.
time_range = ts.utc(2021,11,15,np.linspace(0,24,86400))

min_variance = 1e20
time_of_min_var = None
for time in time_range:
    sat_positions = [sat.at(time).position.km for sat in objects]
    #sat_positions is now a snapshot array where every row is xyz for a sat.
    #Get variance along each column
    sat_positions = np.array(sat_positions)
    variance_by_axis = sat_positions.var(axis=0)
    #Find overall variance by pythagorean of each individual one.
    overall_variance = sum([v**2 for v in variance_by_axis])
    if overall_variance < min_variance:
        min_variance = overall_variance
        time_of_min_var = time
print(min_variance)
print(time_of_min_var)
