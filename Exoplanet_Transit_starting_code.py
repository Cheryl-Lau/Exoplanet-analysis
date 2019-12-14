# -*- coding: utf-8 -*-

'''
#############################################################
 Yr2 Astro Lab - Exoplanet Analysis: Transit identification
#############################################################

 Good practices of Python:
 1. Do NOT name variables with single letters unless it is an index,
    the name should indicate its meaning or purpose 
 2. Label indices and counters as i,j,k,n,m,l
 3. Try to structure your code with functions - easier to fix 
    if one goes wrong! 

'''

import numpy as np
import matplotlib.pyplot as plt 


filename = 'Kep456.txt'

# For Data binning - set the number of sections 
num_sections = 1000


####################################################################

'''
 The following code is given to start you off 
 Feel free to edit, remove or put them into functions 


 - Data binning -

 This code imports raw data from a file. x_data and y_data are the 
 xy-coordinates of all original data points. (Plotted in blue)

 The code then smooths the transit lightcurve by diving the data 
 into lots of sections. You may set the number of sections above.

 For each section, a median y-value is obtained, and we set the 
 the midpoint of the section as the x-value. If we join these new 
 points in each section, we are forming a 'fit curve'. (Plotted in red)

 section_x_list and section_y_list are the coordinates of this fit
 curve and you may now work with this new set of xy-coordinates. 


 ### IMPORTANT ###
 Zoom in the plot and look at each individual transit. Does the shape 
 of the red curve at transits look approximately like a trapezium? 
 Is the line between tI and tII (or tIII and tIV) straight? 
 Tweak the number of sections until it does. 
'''

def data_quantization(n,section_length,y_data):
    # for each section n, assign a part of raw data into this section 
    section_y_data = []
    for i,point in enumerate(y_data):  
        if i >= section_length*n and i < section_length*(n+1):            
            section_y_data.append(point) 
    section_y_data = np.array(section_y_data)
    return section_y_data

processed_data = []
with open(filename) as f:
    data = f.readlines()   
    for i in range(len(data)):
        entry = data[i].split()
        processed_data.append(entry)
        
x_data = [float(p[0]) for p in processed_data]
y_data = [float(p[1]) for p in processed_data]

section_length = int(len(x_data)/num_sections) 
section_x_list = [] 
section_y_list = []    

for n in range(0,num_sections):
    # loop through the n number of sections 
    # set x-value as the midpoint of data in this section
    section_mid_index = section_length*n + int(section_length/2.)
    section_x_value = x_data[section_mid_index]
    section_x_list.append(section_x_value)
    
    # calling function data_quantization to assign y data into this section
    section_y_data = data_quantization(n,section_length,y_data)
    section_median = np.median(section_y_data)
    section_y_list.append(section_median)

fig1 = plt.figure(figsize=(16,6))
label = filename.replace('.txt','')
ax1 = fig1.add_subplot(111)
ax1.set_title('{} Exoplanet Transit Lightcurve'.format(label))
ax1.set_xlabel('Time (Barycentric Julian Date)')
ax1.set_ylabel('Flux (e−/s)')
ax1.plot(x_data,y_data,'.',color='royalblue',markersize=1) 
ax1.plot(section_x_list,section_y_list,'-',color='red')
plt.show()

# Rename the new set of coordinates as x and y for convenience 
x = section_x_list
y = section_y_list

####################################################################

# You may now use the x- and y-values (red line) to find the transit points





































