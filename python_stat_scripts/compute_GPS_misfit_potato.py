#!/bin/python

from io import StringIO
import numpy as np
import numpy.ma as ma
import sys

# This program compares ASPECT's point value velocity output with GPS velocities
# at the same points. It reports the RMS of the angle between the two vectors
# at each point and the RMS of the velocity magnitude.
# It also write to file the angle and velocity magnitude difference at each point.

# read in the file name of the predicted velocities to compare with GPS
data_dir = raw_input("Enter directory name of predicted velocities file: ")
data_file = data_dir+'/point_values.txt'
print 'Comparing ', data_file, ' with GPS.'
print 'Assuming continental plate is C_1'

# set the names of the columns we wish to use from the point_values.txt file
required_names = ['evaluation_point_x', 'evaluation_point_y', 'evaluation_point_z', 'velocity_x', 'velocity_y', 'velocity_z', 'C_1']

# read the point_values.txt file and only select the requested columns
# default dtype is float
aspect_data = np.genfromtxt(data_file, names = True, delimiter = ' ', usecols=(required_names), dtype = float)

# set whether point needs to be considered or not
# this is a 1D array of n_points
aspect_plate_mask = np.genfromtxt('masks.dat', dtype = bool)

aspect_data_masked = ma.masked_array(aspect_data)
aspect_data_masked['velocity_x'].mask = aspect_plate_mask
aspect_data_masked['velocity_y'].mask = aspect_plate_mask
aspect_data_masked['velocity_z'].mask = aspect_plate_mask

# compute the magnitude of the velocity
aspect_vel = np.sqrt(aspect_data_masked['velocity_x']*aspect_data_masked['velocity_x']+aspect_data_masked['velocity_y']*aspect_data_masked['velocity_y']+aspect_data_masked['velocity_z']*aspect_data_masked['velocity_z'])

# the number of points in the file
n_points = float(aspect_vel.shape[0])

# read in the GPS data and set masks
GPS_data = np.genfromtxt('pred_vel_GPS_Nocq_Doubr.dat', names = True, delimiter = ' ', usecols=('velocity0', 'velocity1', 'velocity2'), dtype = float)
GPS_data_masked =  ma.masked_array(GPS_data, mask=aspect_plate_mask)

# compute the magnitude of the velocity of the GPS data
GPS_vel = np.sqrt(GPS_data_masked['velocity0']*GPS_data_masked['velocity0']+GPS_data_masked['velocity1']*GPS_data_masked['velocity1']+GPS_data_masked['velocity2']*GPS_data_masked['velocity2'])

# the number of points in the file
n_points_GPS = float(GPS_vel.shape[0])

# make sure both files have the same number of points
if (not(n_points == n_points_GPS)) :
  sys.exit('ASPECT file and GPS file have a different number of data points: stopping execution.')

# get the number of unmasked points to do statistics
n_unmasked_points = GPS_vel.count()
print '# of unmasked points:', n_unmasked_points

# set up an array that will be filled with point_x, point_y, point_z, angle, diff_magn
diff_output = np.zeros((int(n_points),5))
diff_output[:, 0] = aspect_data_masked['evaluation_point_x']
diff_output[:, 1] = aspect_data_masked['evaluation_point_y']
diff_output[:, 2] = aspect_data_masked['evaluation_point_z']

# compute angle for those data points that lie within the plate (i.e. that are not masked)
angle = np.degrees(np.arccos((aspect_data_masked['velocity_x']*GPS_data_masked['velocity0']+aspect_data_masked['velocity_y']*GPS_data_masked['velocity1']+aspect_data_masked['velocity_z']*GPS_data_masked['velocity2'])/(aspect_vel*GPS_vel)))
angle_sqrd = angle**2
angle_max = np.amax(angle)
angle_min = np.amin(angle)
angle_std = np.std(angle)
angle_mean = np.mean(angle)

# compute the absolute magnitude difference diff_magn for those data points that lie within the plate
diff_vel = np.absolute(aspect_vel - GPS_vel)
diff_vel_sqrd = diff_vel**2
diff_vel_min = np.amin(diff_vel)
diff_vel_max = np.amax(diff_vel)
diff_vel_std = np.std(diff_vel)
diff_vel_mean = np.mean(diff_vel)

# add angle and diff_magn for each point
diff_output[:, 3] = angle
diff_output[:, 4] = diff_vel

# compute RMS angle and print
#RMS_angle = np.sqrt(angle_sqrd.sum()/n_unmasked_points)
#print 'min: ', angle_min, ' max: ', angle_max, ' mean: ', angle_mean, ' standard dev: ', angle_std, ' RMS angle misfit: ', RMS_angle, ' degrees'

# compute RMS magnitude and print
#RMS_diff_vel = np.sqrt(diff_vel_sqrd.sum()/n_unmasked_points)
#print 'min: ', 100.*diff_vel_min, ' max: ', 100.*diff_vel_max, ' mean: ', 100.*diff_vel_mean, ' standard dev: ', 100.*diff_vel_std, ' RMS angle misfit: ', 100.*RMS_diff_vel, ' cm/yr'

# compute RMS angle and print
RMS_angle = np.sqrt(angle_sqrd.sum()/n_unmasked_points)
print angle_sqrd.sum(), n_unmasked_points, RMS_angle
print 'min:', "%.2f" % angle_min, 'max:', "%.2f" % angle_max, 'mean:', "%.2f" % angle_mean, 'std dev:', "%.2f" % angle_std, 'RMS:', "%.2f" % RMS_angle, 'degrees'

# compute RMS magnitude and print
RMS_diff_vel = np.sqrt(diff_vel_sqrd.sum()/n_unmasked_points)
print 'min:', format(100.*diff_vel_min, '0.3f'), 'max:', format(100.*diff_vel_max, '0.3f'), 'mean:', format(100.*diff_vel_mean, '0.3f'), 'std dev:', format(100.*diff_vel_std, '0.3f'), 'RMS', format(100.*RMS_diff_vel, '0.3f'), 'cm/yr'

# print angle and diff_magn for each point
output_file = data_dir+'/GPS_misfit_potato.txt'
np.savetxt(output_file, diff_output, delimiter = ' ')

out = open(output_file, 'a')
out.write ('----------------------\n')
out.write ('Summary:\n')
out.write ('min: ' + str(angle_min) + ' max: ' + str(angle_max) + ' mean: ' + str(angle_mean) + ' std dev: ' + str(angle_std) + ' RMS ' + str(RMS_angle) + ' degrees\n')
out.write ('min: ' +  str(diff_vel_min) + ' max: ' + str(diff_vel_max) + ' mean: ' + str(diff_vel_mean) + ' std dev: ' + str(diff_vel_std) + ' RMS ' + str(RMS_diff_vel) + ' m/yr\n')

