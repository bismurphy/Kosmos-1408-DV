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
        self.line1 = line1.replace("\n","")
        self.line2 = line2.replace("\n","")
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
with open('tles_noDG.txt') as tlefile:
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
print("Variance:" + str(overall_variance))
impact_time = time_of_min_var
print(f"Collision found. {impact_time.utc_strftime()}")
#Step 3: Get the known location of the target sat at that time
intact_line1 = '1 13552U 82092A   21317.92599714  .00002092  00000-0  71807-4 0  9997'
intact_line2 = '2 13552  82.5637 124.8027 0018519 111.6524 248.6690 15.29385168142636'
intact_sat = EarthSatellite(intact_line1,intact_line2)
impact_site = intact_sat.at(impact_time).position.m
print(impact_site)
#Step 4: For each object, identify the best mean anomaly to get close to that location.
for sat in objects:
    best_MA = None
    mindist = 1e20
    for ma_val in np.linspace(0,360,3600):
        sat.set_MA(ma_val)
        dist = np.linalg.norm(sat.at(impact_time).position.m - impact_site)
        if dist < mindist:
            best_MA = ma_val
            mindist = dist
    sat.set_MA(best_MA)
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
print("Variance:" + str(overall_variance))
print(time_of_min_var - impact_time) #This is near 0, so time of impact doesn't change.
#Step 6: Refine mean anomaly. TLE gives 4 digits of precision so let's go with it.
#Previous step found best result to 0.1 degrees, so refine around that.
for sat in objects:
    current_MA = sat.chosen_mean_anomaly
    best_MA = None
    mindist = 1e20
    for ma_val in np.linspace(current_MA -0.1,current_MA + 0.1,1000):
        sat.set_MA(ma_val)
        dist = np.linalg.norm(sat.at(impact_time).position.m - impact_site)
        if dist < mindist:
            best_MA = ma_val
            mindist = dist
    sat.set_MA(best_MA)
#Step 7: Get one last variance analysis.
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
print("Variance:" + str(overall_variance))
#Step 8: Output TLEs
with open('corrected_tles.txt','w') as f:
    for sat in objects:
        f.write(sat.name+"\n")
        f.write(sat.line1+"\n")
        f.write(sat.line2+"\n")
