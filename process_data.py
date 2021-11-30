from skyfield.api import load, EarthSatellite
import numpy as np
import matplotlib.pyplot as plt

print("Loading TLE's of debris objects")
with open('tles.txt') as tlefile:
    lines = tlefile.readlines()
    #split file into groups of 3 lines
    triplets= [lines[i:i+3] for i in range(0,len(lines),3)]
    #Leave out the first two objects (the main sat, and the rocket booster)
    triplets = triplets[2:]

ts = load.timescale()

objects = [EarthSatellite(t[1],t[2],t[0],ts) for t in triplets]
print("TLE's loaded.")
###Identify the location and time of collision. We'll do this by taking the debris
#and, for each timestep, find the center of mass (assuming uniform masses, ugh).
#Then find the moment when the least-squares deviation is minimized
#(debris cloud is smallest).
#We know it was November 15, now let's find out when.
#Create a time object of every second on that day.
#time_range = ts.utc(2021,11,15,np.linspace(0,24,86400))

#the above takes several minutes to run; knowing it's at 2:53 we can go much faster:
time_range = ts.utc(2021,11,15,2,53,np.linspace(0,60,1000))
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
time_of_impact = time_of_min_var
print(f"Collision found. {time_of_impact.utc_strftime()}")

#Now load up the last TLE of the intact satellite:
intact_line1 = '1 13552U 82092A   21317.92599714  .00002092  00000-0  71807-4 0  9997'
intact_line2 = '2 13552  82.5637 124.8027 0018519 111.6524 248.6690 15.29385168142636'
intact_sat = EarthSatellite(intact_line1,intact_line2)
#Establish the velocity that the sat had at the moment of impact
pre_collision_velocity = intact_sat.at(time_of_impact).velocity.m_per_s
#Now get the velocity of every debris object at the moment of ejection
post_collision_velocities = [deb.at(time_of_impact).velocity.m_per_s for deb in objects]
#Subtract the pre velocity from each of the post velocities (row-wise vector subtraction)
delta_v = post_collision_velocities - pre_collision_velocity

#We now have the delta v values, in the ECI XYZ frame. Convert to orbit frame (prograde, radial, normal)
impact_site = intact_sat.at(time_of_impact).position.m
radial_hat = impact_site / np.linalg.norm(impact_site)

prograde_hat = pre_collision_velocity / np.linalg.norm(pre_collision_velocity)
#Sanity check: For circular orbit this is 0 since radial and prograde are perpendicular.
#print(np.dot(radial_hat,prograde_hat))

rcrossp = np.cross(radial_hat,prograde_hat)
normal_hat = rcrossp / np.linalg.norm(rcrossp)

#Now that we have unit vectors, get the delta-v in that coordinate frame by dotting with each.
orbframe_dv = np.array([[np.dot(dv,radial_hat),np.dot(dv,prograde_hat),np.dot(dv,normal_hat)] for dv in delta_v])
#Check that these magnitudes match the original ones
orig_dv_mag = np.array([np.linalg.norm(x) for x in delta_v])
transform_dv_mag = np.array([np.linalg.norm(x) for x in orbframe_dv])
#print(orig_dv_mag - transform_dv_mag)
#Confirmed, the coord transform preserved magnitudes. Now can get to plotting.

radials,progrades,normals = orbframe_dv.T

print(radials[-5])
print(orbframe_dv[-5])
print(delta_v[-5])
print(post_collision_velocities[-5])
print(objects[-5])
##plt.scatter(progrades,radials)
##plt.title("Prograde DV versus Radial DV")
##plt.xlabel("Prograde velocity change (m/s)")
##plt.ylabel("Radial velocity change (m/s)")
##plt.show()
##plt.scatter(progrades,normals)
##plt.title("Prograde DV versus Normal DV")
##plt.xlabel("Prograde velocity change (m/s)")
##plt.ylabel("Normal velocity change (m/s)")
##plt.show()
##plt.scatter(radials,normals)
##plt.title("Radial DV versus Normal DV")
##plt.xlabel("Radial velocity change (m/s)")
##plt.ylabel("Normal velocity change (m/s)")
##plt.show()

