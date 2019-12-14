# -*- coding: utf-8 -*-
####################################################################
 ### Yr2 Astro Lab - Exoplanet Analysis: Transit identification ###
####################################################################

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd 


def main(filename,num_sections,grad_filter):
    
    fig1 = plt.figure(figsize=(16,6))
    label = filename.replace('.txt','')
    ax1 = fig1.add_subplot(111)
    ax1.set_title('{} Exoplanet Transit Lightcurve'.format(label))
    ax1.set_xlabel('Time (Barycentric Julian Date)')
    ax1.set_ylabel('Flux (e−/s)')

    # call subroutine to import data from file 
    print('Reading input file',filename)
    data = data_import(filename,ax1)
    
    # call subroutine loop_over_sections
    section_x_values_list,section_median_list = loop_over_sections(num_sections,data,ax1)
    plt.savefig(label+'.png') 
    
    # call subroutine gradients
    start_point_list,end_point_list,start_point_tail_list,end_point_tail_list, \
        = gradients(grad_filter,section_x_values_list,section_median_list)
        
    # call subroutine group_data
    transit_data = group_data(start_point_list,end_point_list,start_point_tail_list,end_point_tail_list)
    
    transit_data_table = pd.DataFrame(transit_data,columns=['tI','tII','tIV','tIII','t0'])       
    print(transit_data_table) # row:[tI,tII,tIV,tIII,t0]; column:transits 

    return 
        
        
        
def data_import(filename,ax1):
    # import data from file as processed_data matrix and plot the graph 
    
    processed_data = []
    with open(filename) as f:
        data = f.readlines()        
        for i in range(len(data)):
            entry = data[i].split()
            processed_data.append(entry)
        
    time = [p[0] for p in processed_data]
    flux = [p[1] for p in processed_data]
    ax1.plot(time,flux,'.',color='royalblue',markersize=1)
    plt.show()
    
    return processed_data 


# Divide data into smaller sections, and for each section obtain the median value
# to eliminate random errors and fluctuations, producing a best-fit curve
      
def loop_over_sections(num_sections,data,ax1):
    
    x_data = [float(d[0]) for d in data]
    y_data = [float(d[1]) for d in data]
    # divide data into n sections 
    # loop over each section n
    section_length = int(len(x_data)/num_sections) 
    section_x_values_list = [] 
    section_median_list = []    
    for n in range(0,num_sections):
        section_mid_index = section_length*n + int(section_length/2.)
        section_x_value = x_data[section_mid_index]
        section_x_values_list.append(section_x_value)
        
        section_y_data = data_quantization(n,section_length,y_data)
        section_median = np.median(section_y_data)
        section_median_list.append(section_median)
        
        if n % 100 == 0:   
            progress = int(n/num_sections*100)
            print('LOADING: {}%'.format(progress))
        
    # plot median fit line 
    ax1.plot(section_x_values_list,section_median_list,'-',markersize=1,color='red')  
    plt.show()        
      
    return (section_x_values_list, section_median_list)


def data_quantization(n,section_length,y_data):
    # called from subroutine loop_over_sections
    # for each section n, assign parts of raw data into this section 
    
    section_y_data = []
    for i,point in enumerate(y_data):  
        if i >= section_length*n and i < section_length*(n+1):            
            current_section_data = point
            section_y_data.append(current_section_data) 
    section_y_data = np.array(section_y_data)
    
    return section_y_data # (len(data)/num_sections)x2 matrix


# Obtain the gradients (diff between two points) along the best-fit curve 
# Regions with sufficiently large increase/decrease in gradient are dips 

def gradients(grad_filter,x_values,y_values):
    # main function for calling subroutines transit_search and transit_check
    
    abs_gradient_distribution = []
    for i in range(0,len(y_values)-1):
        gradient = y_values[i+1] - y_values[i]
        abs_gradient_distribution.append(abs(gradient))           
    # average gradient 
    gradient_sum = 0
    for grad in abs_gradient_distribution:   
        gradient_sum += grad
    mean_abs_gradient = gradient_sum / len(abs_gradient_distribution)
       
    transit_start_point_list = [] #initialize
    transit_end_point_list = [] 
    
    # first run, proceed to find transit
    print('Begin transit_search routine with filtering condition {}'.format(grad_filter))
    transit_start_point_list,transit_end_point_list,transit_start_point_tail_list, \
        transit_end_point_tail_list = transit_search(grad_filter, \
        mean_abs_gradient,x_values,y_values,transit_start_point_list,transit_end_point_list)
    # outputs transits start and end points
    
    # check that results are returned and are valid transits:
        
    while transit_start_point_list and transit_check(transit_start_point_list, \
        transit_end_point_list) is False:
        # if output arrays are NOT empty, AND did NOT pass transit_check test
        # change filtering condition, then repeat process
        grad_filter += 0.5
        print('Error in transit identification - retry with filtering condition {}'.format(grad_filter))
        
        transit_start_point_list,transit_end_point_list,transit_start_point_tail_list, \
            transit_end_point_tail_list = transit_search(grad_filter, \
            mean_abs_gradient,x_values,y_values,transit_start_point_list,transit_end_point_list)
            
    if transit_start_point_list and transit_check(transit_start_point_list, \
    transit_end_point_list) is True:            
        # if output arrays are NOT empty, AND passed transit_check test
        print('----Transits identification successful----')
                                                                        
        final_transit_start_point_list = transit_start_point_list
        final_transit_end_point_list = transit_end_point_list
        final_transit_start_point_tail_list = transit_start_point_tail_list
        final_transit_end_point_tail_list = transit_end_point_tail_list
        final_grad_filter = grad_filter
        
        print('Grad_filter:', final_grad_filter)           
                                            
        return (final_transit_start_point_list,final_transit_end_point_list, \
            final_transit_start_point_tail_list,final_transit_end_point_tail_list)


def transit_search(grad_filter,mean_abs_gradient,x_values,y_values,transit_start_point_list,transit_end_point_list):
    # called from subroutine gradients 
    # Finds the transit start/end (tail) points   ---- tail: transit start/end points at lowest flux
    print('calling transit_search')
    
    gradient_list = []
    for i in range(len(y_values)-1):  
        gradient = y_values[i+1] - y_values[i]
        gradient_list.append(gradient)
   
    transit_start_range_index_list = []
    transit_end_range_index_list = []
    
    for i,grad in enumerate(gradient_list):     
        # if gradient is sufficiently large
        if abs(grad) > mean_abs_gradient * grad_filter:
            
            # get the range of indices of large gradients at start and end of transits
            # from start/end point to tail
            
            if grad >= 0:  # increasing; end of dip
                transit_end_range_index_list.append(i)          
            elif grad < 0:  # decreasing; start of dip
                transit_start_range_index_list.append(i)                   
    # transit_start/end_range_index_list now contains all large increasing/decreasing bits in all transits 
    
    # call subroutine divide_transit_index_range to separate the indices lists into sets for individual transits
    # a set of start_indices and a set of end_indices for each transit
    splitted_transit_start_index_range = divide_transit_index_range(transit_start_range_index_list)
    splitted_transit_end_index_range = divide_transit_index_range(transit_end_range_index_list)
    
    transit_start_point_list = []
    transit_start_point_tail_list = []
    transit_end_point_list = []
    transit_end_point_tail_list = []    
    
    # using the list of indices found above, obtain the transit time values for each set
    
    for transit_start_index in splitted_transit_start_index_range:
        # the indices for start in one transit (one set)
        transit_start_point_index = transit_start_index[0] # first one
        transit_start_point = x_values[transit_start_point_index]
        transit_start_point_tail_index = transit_start_index[-1] #last one
        transit_start_point_tail = x_values[transit_start_point_tail_index+1] # the point on its right
        
        transit_start_point_list.append(float(transit_start_point))
        transit_start_point_tail_list.append(float(transit_start_point_tail))
        
    for transit_end_index in splitted_transit_end_index_range:
        # the indices for end in one transit (one set)
        transit_end_point_index = transit_end_index[-1] # last one
        transit_end_point = x_values[transit_end_point_index+1] # on right
        transit_end_point_tail_index = transit_end_index[0] # first one
        transit_end_point_tail = x_values[transit_end_point_tail_index] 
        
        transit_end_point_list.append(float(transit_end_point))
        transit_end_point_tail_list.append(float(transit_end_point_tail))
    
    return transit_start_point_list,transit_end_point_list,transit_start_point_tail_list, \
        transit_end_point_tail_list

    
def divide_transit_index_range(transit_range_index_list):
    # called from subroutine transit_search, used for both transit_start and _end sets
    # group start/end_range indices which are in consecutive order 
    # a set of index in consecutive order -> belongs to the same transit
    # divides the array for all transits into individual arrays for each transit 

    split = [0] + [i for i in range(1,len(transit_range_index_list)) if transit_range_index_list[i] - \
        transit_range_index_list [i-1] >1] + [None]  
    splitted_transit_range_index_list = [transit_range_index_list[b:e] for (b,e) in [(split[i-1],split[i])  \
        for i in range(1,len(split))]]
    # Ref: https://stackoverflow.com/questions/3149440/python-splitting-list-based-on-missing-numbers-in-a-sequence
    
    return splitted_transit_range_index_list # index of transit start/end range, divided into sets for each dip  


def transit_check(transit_start_point_list,transit_end_point_list):
    # called from subroutine gradients 
    # if the results obtained from subroutine transit_search do not pass the following tests
    # the program repeats with a different filtering condition 
    ### NOTE: this subroutine is NOT suitable for data with less than 3 transits ###
    print('calling transit_check')
    
    # 1. check that the transits are in pairs 
    if len(transit_start_point_list) == len(transit_end_point_list):
        
        # 2. check that the transit start and end points are in alternating order
        # otherwise it is probably not an exoplanet transit (e.g. discontinuities)
        check_box = []
        for i in range(0,len(transit_start_point_list)-1):
            if transit_start_point_list[i] > transit_end_point_list[i] and  \
                transit_end_point_list[i+1] > transit_start_point_list[i]: # end-start-end-start-....
                check_box.append('e-s')
            elif transit_start_point_list[i] < transit_end_point_list[i] and \
                transit_end_point_list[i] < transit_start_point_list[i+1]: # start-end-start-end-....
                check_box.append('s-e')
            else:
                check_box.append('X') 
                
        consistency_check = True
        for check in check_box:
            if check != check_box[0]:
                consistency_check = False 
                
        if consistency_check == True:
            print('valid transits')
            return True 
        elif consistency_check == False:
            print('invalid transits')
            return False 
            
    elif len(transit_start_point_list) != len(transit_end_point_list):
        print('invalid transits')
        return False


def group_data(start_point_list,end_point_list,start_point_tail_list,end_point_tail_list):
    # packs the results from previous subroutines into a nice table 
    # output transit data nx4 matrix where n is number of dip
    grouped_transit_data = [start_point_list,start_point_tail_list,end_point_list,end_point_tail_list]
    transit_data = []
    # for the i-th dip in lightcurve
    for i in range(len(start_point_list)):
        transit_data.append([entry[i] for entry in grouped_transit_data])
      
    transit_midpoint = [(float(dip[2]) + float(dip[0]))/2. for dip in transit_data]    
    transit_data = np.hstack((transit_data, np.atleast_2d(transit_midpoint).T))         

    return transit_data # [tI,tII,tIV,tIII,t0]



filename = 'Kep456.txt'
grad_filter = 0.5
num_sections = 5000

main(filename,num_sections,grad_filter)























