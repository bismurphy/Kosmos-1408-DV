#The TLE's we got from celestrak do not propagate back to yield
#a single origin point for all debris objects. Mean anomalies are off.
#This script seeks to correct that and bring them together. Load up
#all objects, and adjust mean anomaly until we get them all as close as possible
#to the single center of mass.
from skyfield.api import load,EarthSatellite
import numpy as np

#Define a class of satellite that remembers its mean anomaly and can be adjusted.
class meanSat(EarthSatellite):
    def __init__(self,line1,line2,name=None, ts=None):
        self.line1 = line1
        self.line2 = line2
        self.name = name
        self.ts = ts
        super().__init__(line1,line2,name,ts)
    def set_MA(self,new_mean_anomaly):
        self.chosen_mean_anomaly = new_mean_anomaly
        new_line2 = self.line2[:43] + f'{new_mean_anomaly:8.4f}' + self.line2[51:]
        self.line2 = new_line2
        super().__init__(self.line1,self.line2,self.name,self.ts)

#Step 1: Load initial TLEs.
print("Loading TLE's of debris objects")
with open('tles.txt') as tlefile:
    lines = tlefile.readlines()
    #split file into groups of 3 lines
    triplets= [lines[i:i+3] for i in range(0,len(lines),3)]
    #Leave out the first two objects (the main sat, and the rocket booster)
    triplets = triplets[2:]

ts = load.timescale()

objects = [meanSat(t[1],t[2],t[0],ts) for t in triplets]
print("TLE's loaded.")
#Step 2: Identify the best guess for the breakup location.
time_range = ts.utc(2021,11,15,2,53,np.linspace(40,50,100))
print("Identifying time of collision")
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
print(overall_variance)
impact_time = time_of_min_var
print(f"Collision found. {impact_time.utc_strftime()}")
#Step 3: Get the average position of all the objects at that time
locations_at_impact = np.array([sat.at(impact_time).position.m for sat in objects])
mean_location = locations_at_impact.mean(axis=0)
print(mean_location)
#Step 4: For each object, identify the best mean anomaly to get close to that location.
for sat in objects:
    best_MA = None
    mindist = 1e20
    for ma_val in np.linspace(0,360,1000):
        sat.set_MA(ma_val)
        dist = np.linalg.norm(sat.at(impact_time).position.m - mean_location)
        if dist < mindist:
            best_MA = ma_val
            mindist = dist
#Step 5: Rerun the center-finder.
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
print(overall_variance)
