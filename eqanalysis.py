"""
eqanalysis.py: analyze and plot earthquake data
Authors: Nolan Cassidy
Credits: Jeske Glenn

CIS 210 assignment 7, Fall 2016
"""

# a set of modules that we need to use in the code below
import math
import random
import argparse
from data import *
import turtle
import sys

# constants for the k-means clustering algorithm
# if you change these in your experimentation, you will need to look at 
# all parts of the code that refer to them, as there is some dependence
# on them (such as number of colors used in plotting clusters)
#
# IF YOU DO CHANGE THEM, YOU MUST PUT THEM BACK TO THE ORIGINAL VALUES
# BEFORE SUBMITTING YOUR WORK!!!!!
NO_OF_CLUSTERS = 6
NO_OF_ITERATIONS = 7
COLORS = ["white","gray","wheat","pink","violet","lightskyblue","blue","green","yellow","orange","red"]

def euclid_distance(point1, point2):
    """
    computes the euclidean distance between two points
    Args:
        point1: list of floats, index 0 is longitude, index 1 is latitude
        point2: list of floats, index 0 is longitude, index 1 is latitude
    Returns:
        float, sqrt((x1-x2)**2 + (y1-y2)**2)
    """

    total = 0
    for index in range(2):
        diff = point1[index] - point2[index]
        total += diff * diff

    return math.sqrt(total)

def create_centroids(k, datadict):
    """
    randomly selects 'k' points from 'datadict' as the starting
        centroids for the k-means clustering algorithm
    Args:
        k: int, number of clusters desired
        datadict: list of lists, each contained list represents an EQ event
    Returns:
        list of lists, each contained list is an event to act as the centroid
    """
    centroids = []
    count = 0
    centroid_keys = []

    while count < k:
        rkey = random.randint(1, len(datadict))
        if rkey not in centroid_keys:
            centroids.append(datadict[rkey])
            centroid_keys.append(rkey)
            count += 1

    return centroids

def create_clusters(k, centroids, datadict, iterations):
    """
    k-means clustering algorithm - implementation taken from page 249 of
        ranum and miller text, with some modifications
    Args:
        k: integer, number of clusters
        centroids: list of events, each event is the centroid of its cluster
        datadict: dictionary of all EQ events
        iterations: int, number of clustering iterations to perform
    Returns:
        list of lists: each contained list is the set of indices into 'datadict'
           for events that belong to that cluster
    """
    for iteration in range(iterations):
        #print("****Iteration", iteration, "****")
        clusters = []
        for i in range(k):
            clusters.append([])

        for key in datadict:
            distances = []
            for cl_index in range(k):
                dist = euclid_distance(datadict[key], centroids[cl_index])
                distances.append(dist)
            min_dist = min(distances)
            index = distances.index(min_dist)
            clusters[index].append(key)

        dimensions = 2
        for cl_index in range(k):
            sums = [0]*dimensions
            for key in clusters[cl_index]:
                data_points = datadict[key]
                for ind in range(2):
                    sums[ind] = sums[ind] + data_points[ind]
            for ind in range(len(sums)):
                cl_len = len(clusters[cl_index])
                if cl_len != 0:
                    sums[ind] /= cl_len
            centroids[cl_index] = sums

        #for c in clusters:
            #print("CLUSTER")
            #for key in c:
                #print(datadict[key], end=" ")
            #print()

    return clusters

def read_file(filename):
    """
    read the EQ events from the csv file, 'filename'; any lines starting with
        # are skipped; the longitude, latitude, magnitude, and depth (in miles)
        is extracted from each event record, and stored as a list against its
        record number in a dictionary
    Args:
        filename: string, name of a CSV file containing the EQ data
    Returns:
        dictionary, indexed by integers, each value is a list of floats
            representing an EQ event
    """
    dict = {}
    key = 0

    fd = open(filename, "r")
    for line in fd:
        if line[0] == '#':
            continue		# causes the loop to grab another line
        key += 1
        values = line.rstrip('\n').split(',')
        lat = float(values[7])
        lon = float(values[8])
        mag = float(values[1])
        dep = float(values[10])
        dict[key] = [lon, lat, mag, dep]
    fd.close()
    return dict

# global data for map - if we had ;earmed about classes yet, this would have
# been hidden in a class instance, and the plot_*() functions would be methods
# on that class instance.  for now, these are global variables, and the
# plot functions access them

eq_turtle = None
eq_win = None
# these are the longitudes and latitudes for the Pacific NorthWest map that
# I have provided to you; do not change them!
left_lon = -128.608689
right_lon = -114.084764
top_lat = 51.248522
bot_lat = 38.584004
lon_diff = 0
lat_diff = 0
size_x = 0
size_y = 0
left_x = 0
bot_y = 0

def prepare_turtle():
    """
    Prepares the turtle and the window to plot magnitudes, depths, or clusters
    Args:
        None
    Outputs:
        creates turtle, sets window size, defines remainder of global
        data needed for plot_routines
    """
    global eq_turtle, eq_win
    global left_lon, right_lon, top_lat, bot_lat
    global lon_diff, lat_diff
    global size_x, size_y, left_x, bot_y

    eq_turtle = turtle.Turtle()
    eq_turtle.speed(10)
    eq_win = turtle.Screen()
    eq_win.screensize(655,808)	# number of pixels in the map I have provided
    lon_diff = right_lon - left_lon
    lat_diff = top_lat - bot_lat
    size_x = eq_win.screensize()[0]
    left_x = -size_x/2
    size_y = eq_win.screensize()[1]
    bot_y = -size_y/2
    eq_win.bgpic("PacificNW.gif")	# the map I have provided
    eq_turtle.hideturtle()
    eq_turtle.up()

def xy_calculate(lon, lat):
    """
    compute (x, y) given lon[gitude] and lat[itude]
    Args:
        lon: float, longitude value for point on map
        lat: float, latitude value for point on map
    Returns:
        tuple, corresponding pixel x and y values for use in turtle methods
    """
    global left_lon, right_lon, top_lat, bot_lat
    global lon_diff, lat_diff
    global size_x, size_y, left_x, bot_y

    x = left_x + (lon - left_lon) / lon_diff * size_x
    y = bot_y + (lat - bot_lat) / lat_diff * size_y
    return (x, y)

def plot_clusters(eq_clusters, eq_dict):
    """
    plot the clusters - use turtle.dot() at the appropriate location on the
        map for each event; use a different color for the events in each
        cluster - e.g. for cluster 0, use 'red', for 1, use 'violet' ...
    Args:
        eq_clusters: list of lists, each contained list has the indices for
                     events in that cluster in eq_dict
        eq_dict: list of lists, each contained list represents an EQ event
    Outputs:
        plots all events in a particular cluster as dots on the map
    """
    global eq_turtle
    count = 0
    for i in eq_clusters:
            for n in i:
                if n in eq_dict:
                    lon = eq_dict[n][0]
                    lat = eq_dict[n][1]
                    x,y = xy_calculate(lon,lat)
                    eq_turtle.goto(x,y)
                    eq_turtle.dot(COLORS[count])
            count+=1


def bin_value(value, bounds):
    """
    'bounds' defines a set of bins; this function returns the index of the
        first bin that contains 'value'
    Args:
        value: float, value to place in bin
        bounds: list of floats, bounds[i] is the top value of the bin
                code assumes that bounds is an increasing set of values
    Returns:
        integer, index of smallest value of bounds[] that is >= value
            if value > bounds[-1], returns len(bounds)
    """
    for i in range(len(bounds)):
        if value <= bounds[i]:
            return i
    return len(bounds)

def plot_magnitudes(eq_dict):
    """
    plot the magnitudes - use turtle.dot() at the appropriate location on the
        map for each event; use a different color and size for magnitude
        equivalence classes - e.g. if magnitude of event is <=1, use small dot
        that is 'violet', if between 1 and 2, use slightly larger dot that is
        'blue', ..., if between 9-10, use very large dot that is 'red'
    Args:
        eq_dict: list of lists, each contained list represents an EQ event
    Outputs:
        plots magnitude of all events as dots on the map
    """
    global eq_turtle

    for i in eq_dict:
        mags = int(eq_dict[i][2])
        lon = eq_dict[i][0]
        lat = eq_dict[i][1]
        x, y = xy_calculate(lon, lat)
        eq_turtle.goto(x, y)
        eq_turtle.dot((mags+5)*2,(COLORS[mags]))

def plot_depths(eq_dict):
    """
    plot the depths - use turtle.dot() at the appropriate location on the
        map for each event; use a different color and size for depth
        equivalence classes - e.g. if depth of event is <=1 mile, use a very
        large dot that is 'red', if between 1 and 5, use slightly smaller dot
        that is 'orange', ..., if between 50-100, use a small dot that is
        'violet'
    Args:
        eq_dict: list of lists, each contained list represents an EQ event
    Outputs:
        plots depth of all events as dots on the map
    """
    global eq_turtle

    for i in eq_dict:
        dep = eq_dict[i][3]
        lon = eq_dict[i][0]
        lat = eq_dict[i][1]
        x, y = xy_calculate(lon, lat)
        eq_turtle.goto(x, y)
        if dep <= 1:
            eq_turtle.dot(5, 'violet')
        elif dep <= 5:
            eq_turtle.dot(10, 'light blue')
        elif dep <= 15:
            eq_turtle.dot(15, 'blue')
        elif dep <= 20:
            eq_turtle.dot(20, 'dark blue')
        elif dep <= 25:
            eq_turtle.dot(25, 'yellow')
        elif dep <= 30:
            eq_turtle.dot(30, 'salmon')
        elif dep <= 35:
            eq_turtle.dot(35, 'orange')
        elif dep <= 40:
            eq_turtle.dot(40, 'coral')
        elif dep <= 45:
            eq_turtle.dot(45, 'orange red')
        elif dep <= 50:
            eq_turtle.dot(50, 'red')

def analyze_depths(eq_dict):
    """
    Perform statistical analysis on the depth information in the dictionary
    Args:
        eq_dict: list of lists, each contained list represents an EQ event
    Outputs:
        mean, median, and standard deviation of depth data
        frequency table for the depth data
    """
    median = 0.0
    stdd = 0.0
    deps = []

    for i in eq_dict:                   #creates depth list
        deps.append(eq_dict[i][-1])
    deps.sort()

    n1 = 0  # finds Standard Deviation
    m1 = 0.0
    for i in deps:
        n1 += 1
        delta = i - m1
        m1 += delta / n1

    if len(eq_dict)%2 == 0:             #finds median
        rmid = len(eq_dict)//2
        r = deps[rmid]
        l = deps[rmid-1]
        median = (r+l)/2
    else:
        mdpt = len(eq_dict)//2
        median = deps[mdpt]

    n2 = 0
    m2 = 0.0
    for i in deps:
        n2 += 1
        delta = i*i - m2
        m2 += delta / n2
    N = len(deps)
    fracN = N/(N-1)
    stdd = math.sqrt((fracN*(m2-(m1**2))))

    print("Mean depth = {:.1f} miles".format(m1))
    print("Median depth = {:.1f} miles".format(median))
    print("Median depth = {:.2f} miles".format(stdd))

    count = {}
    for i in sorted(eq_dict):
        if eq_dict[i][-1] in count:
                count[eq_dict[i][-1]]+=1
        else:
            count[eq_dict[i][-1]] = 1
    i_list = sorted(list(count.keys()))
    print("ITEM     FREQUENCY")
    for j in i_list:
        print("{:5}   {:5}".format(j,count[j]))


def analyze_magnitudes(eq_dict):
    """
    Perform statistical analysis on the magnitude information in the dictionary
    Args:
        eq_dict: list of lists, each contained list represents an EQ event
    Outputs:
        mean, median, and standard deviation of magnitude data
        frequency table for the magnitude data
    """
    median = 0.0
    stdd = 0.0
    mags = []

    for i in eq_dict:  # creates depth list
        mags.append(eq_dict[i][2])
    mags.sort()

    n1 = 0  # finds Standard Deviation
    m1 = 0.0
    for i in mags:
        n1 += 1
        delta = i - m1
        m1 += delta / n1
    n2 = 0
    m2 = 0.0

    if len(eq_dict) % 2 == 0:  # finds median
        rmid = len(eq_dict) // 2
        r = mags[rmid]
        l = mags[rmid - 1]
        median = (r + l) / 2
    else:
        mdpt = len(eq_dict) // 2
        median = mags[mdpt]

    for i in mags:
        n2 += 1
        delta = i * i - m2
        m2 += delta / n2
    N = len(mags)
    fracN = N / (N - 1)
    stdd = math.sqrt((fracN * (m2 - (m1 ** 2))))

    print("Mean depth = {:.1f} miles".format(m1))
    print("Median depth = {:.1f} miles".format(median))
    print("Median depth = {:.2f} miles".format(stdd))

    count = {}
    for i in sorted(eq_dict):
        if eq_dict[i][2] in count:
            count[eq_dict[i][2]] += 1
        else:
            count[eq_dict[i][2]] = 1
    i_list = sorted(list(count.keys()))
    print("ITEM     FREQUENCY")
    for j in i_list:
        print("{:5}   {:5}".format(j, count[j]))


def analyze_clusters(eq_clusters, eq_dict):
    """
    Perform statistical analysis on the depth and magnitude information
        for each cluster
    Args:
        eq_clusters: list of lists, each contained list has the indices into
                     eq_dict for events in that cluster
        eq_dict: list of lists, each contained list represents an EQ event
    Outputs:
        for each cluster:
            mean, median, and standard deviation of magnitude data
            mean, median, and standard deviation of depth data
    """
    count = 0
    for i in eq_clusters:           #creats  magnitude and depths list
        deps = []
        mags = []
        for n in i:
            d = eq_dict[n][-1]
            m = eq_dict[n][2]
            deps.append(d)
            mags.append(m)
        deps.sort()
        mags.sort()

        sum = 0
        for m in mags:
            sum += m
        mmean = sum / len(mags)

        if len(mags) % 2 == 0:
            rm = len(mags) // 2
            lm = rm - 1
            rpt = mags[rm]
            lpt = mags[lm]
            mmedian = (rpt + lpt) / 2.0
        else:
            mdpt = len(mags) // 2
            mmedian = mags[mdpt]

        k = 0
        m2 = 0.0
        for j in mags:
            k += 1
            delta = j * j - m2
            m2 += delta / k
        N = len(mags)
        fracN = N / (N - 1)
        x = (fracN * (m2 - (mmean ** 2)))
        mstdd = math.sqrt(x)

        sum = 0
        for j in deps:
            sum += j
        dmean = sum / len(deps)

        if len(deps) % 2 == 0:
            rm = len(deps) // 2
            lm = rm - 1
            rpt = deps[rm]
            lpt = deps[lm]
            dmedian = (rpt + lpt) / 2.0
        else:
            mdpt = len(deps) // 2
            dmedian = deps[mdpt]

        k = 0
        m2 = 0.0
        for j in deps:
            k += 1
            delta = j * j - m2
            m2 += delta / k
        N = len(deps)
        fracN = N / (N - 1)
        x = (fracN * (m2 - (dmean ** 2)))
        dstdd = math.sqrt(x)

        print("Analysis for cluster {}".format(count))
        print("   Analysis of magnitude data")
        print("       Mean value = {:.1f} miles".format(mmean))
        print("       Median value = {:.1f} miles".format(mmedian))
        print("       Standard Deviation = {:.2f} miles".format(mstdd))
        print("   Analysis of depth data")
        print("       Mean value = {:.1f} miles".format(dmean))
        print("       Median value = {:.1f} miles".format(dmedian))
        print("       Standard Deviation = {:.2f} miles".format(dstdd))
        count += 1

def main():
    """
    Interaction if run from the command line.
    Usage:  python3 eqanalysis.py eq_data_file.csv command
    """
    parser = argparse.ArgumentParser(description="Earthquake event file stats")
    parser.add_argument('eq_file', type=str,
                 help='A csv file containing earthquake events, one per line.')
    parser.add_argument('command', type=str,
                 help='One of the following strings: plot analyze')
    parser.add_argument('what', type=str,
                 help='One of the following strings: clusters depths magnitudes')
    args = parser.parse_args()
    eq_file = args.eq_file
    cmd = args.command
    what = args.what
    if cmd != 'plot' and cmd != 'analyze':
        print('Illegal command: {}; must be "plot" or "analyze"'.format(cmd))
        sys.exit(1)
    if what != 'clusters' and what != 'magnitudes' and what != 'depths':
        print('Can only process clusters, magnitudes, or depths')
        sys.exit(1)
    eq_dict = read_file(eq_file)
    prepare_turtle()
    if what == 'clusters':
        eq_centroids = create_centroids(NO_OF_CLUSTERS, eq_dict)
        eq_clusters = create_clusters(NO_OF_CLUSTERS, eq_centroids, eq_dict, NO_OF_ITERATIONS)
    if cmd == 'plot':
        if what == 'clusters':
            plot_clusters(eq_clusters, eq_dict)
        elif what == 'magnitudes':
            plot_magnitudes(eq_dict)
        elif what == 'depths':
            plot_depths(eq_dict)
        print("ALL EVENTS HAVE BEEN PLOTTED")
        eq_win.exitonclick()
    else:
        if what == 'clusters':
            analyze_clusters(eq_clusters, eq_dict)
        elif what == 'magnitudes':
            analyze_magnitudes(eq_dict)
        elif what == 'depths':
            analyze_depths(eq_dict)

if __name__ == "__main__":
    main()
