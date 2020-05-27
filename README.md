# exoplanet-analysis

Uni of York physics - Yr2 exoplanet python lab materials 

(Updated 17/10/2019) 

This repo includes codes for data binning and locating the transit points. 
Calculations of physical quantities not included. 


# I got the starting code, what do I do now?

Here's a hint:
Lets say the line between each point on the smoothed coord (already done for you) is a segment
1. Examine the gradient of each segment using a for-loop with a counter (which gives you the index of the segment the loop is currently working on)
2. Set up conditions which filter those with sufficiently large gradients, preferably one for large ups, and one for large drops 
3. Now you've got the index of the segments with large pos gradients 
How can you use this index to locate the two points? (the point at the start and the end of that segment)
Repeat the same for segments with large neg gradients 
4. Using these list of index, find the actual time & flux value of each point 

P.S. That should take no more than ~20 lines, don't overthink! 

# Coding knowkedge required 

Please google the following python tools and learn how to use them: 
1. for-loops with a counter 
2. list 
3. append  (for storing numbers generated from each iteration of a loop into a single list of numbers)
