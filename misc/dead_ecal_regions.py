#Imports
from random import randint
from scipy.stats import poisson
from math import sqrt

# Options (for ECAL Barrel) and parameters 
iphi_min = -180
iphi_max = 180
iphi_width = 0.0175
ieta_min = -85
ieta_max = 85
ieta_width = 0.0175
dead_fraction_of_cells = 0.01
number_of_toys = 1
cone_size = 0.3

number_of_cells = (iphi_max-iphi_min)*(ieta_max-ieta_min)
mean_number_of_dead = int(number_of_cells*dead_fraction_of_cells)
#print number_of_cells

# Utility methods
def cell_id( iphi, ieta ) :
    return (iphi-iphi_min) * (ieta_max-ieta_min) + (ieta-ieta_min)

def cell_coords( cell_id ) :
    return ( cell_id%(iphi_max-iphi_min) + iphi_min, cell_id/(iphi_max-iphi_min) + ieta_min )

def wrap_iphi( iphi ) :
    return (iphi-iphi_min)%(iphi_max-iphi_min)+iphi_min

def delta_iphi( iphi1, iphi2 ) :
    diff = iphi1-iphi2
    half = (iphi_max-iphi_min)/2
    if diff < 0 :
        if abs(diff) < half : return diff
        else : return (half*2)+diff
    else :
        if diff > half : return -1*((half*2)-diff)
        else : return diff
   
def delta_ieta( ieta1, ieta2 ) :
    return ieta1-ieta2

def delta_phi( delta_iphi ) :
    return delta_iphi*iphi_width

def delta_eta( delta_ieta ) :
    return delta_ieta*ieta_width

def delta_r( iphi1, ieta1, iphi2, ieta2 ) :
    d_iphi = delta_iphi(iphi1,iphi2)
    d_ieta = delta_ieta(ieta1,ieta2)
    d_phi = delta_phi(d_iphi)
    d_eta = delta_eta(d_ieta)
    return sqrt( d_phi*d_phi + d_eta*d_eta )

def cells_within_delta_r( id ) :
    data = set()
    coords = cell_coords(id)
    iphi_bins = int(cone_size / iphi_width)+1
    ieta_bins = int(cone_size / ieta_width)+1
    for i in range(-1*iphi_bins+coords[0],iphi_bins+coords[0]) :
        for j in range(-1*ieta_bins+coords[1],ieta_bins+coords[1]) :
            dr = delta_r( wrap_iphi(i), j, coords[0], coords[1] )
            if dr < cone_size : data.add( cell_id(i,j) )
#                print iphi_bins,ieta_bins,i,j,coords,dr,len(data)
    return data

def validation() :
    cntr = 0
    for i in range(iphi_min,iphi_max) :
        for j in range(ieta_min,ieta_max) :
            cntr = cntr + 1
            if cntr > 100 : break
            id = cell_id(i,j)
            coords = cell_coords(id)
            if i != coords[0] or j != coords[1] :
                print i,j,id,coords
#    for i in range(iphi_min*2,iphi_max*2) :
#        print i,wrap_iphi(i)
#    for i in range(-180,180,90) :
#        for j in range(-180,180,90) :
#            print "i,j,delta_iphi,ieta,r:",i,j,delta_iphi(i,0), delta_ieta(j,0), delta_r(i,0,j,0) 
#validation()
        
# Random number of dead cells
rndm = poisson.rvs(mean_number_of_dead,size=number_of_toys)
#print rndm

# Calculate fraction of ECAL Barrel area (using cells as a proxy) that are within deltaR of dead cell 

dead_fractions = []
for i in range(number_of_toys) : 

    dead_cell_ids = set() # guarantee uniqueness
    while len(dead_cell_ids) < mean_number_of_dead:
        iphi = randint(iphi_min,iphi_max)
        ieta = randint(ieta_min,ieta_max)
        id = cell_id(iphi,ieta)
        dead_cell_ids.add(id)
#    print len(dead_cell_ids)

    cells_list = set()
    for j in dead_cell_ids :
        new_cells = cells_within_delta_r(j)
        print len(new_cells)
        #cells_list.union(  )

    dead_fractions.append( float(len(cells_list))/float(number_of_cells) )

print dead_fractions
