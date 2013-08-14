#!/usr/bin/env python2.3
"""
A simple test routine that tests at least some of the TAP_ext package

"""
import unittest

from Numeric import *


class CalcPolygonsTestCase(unittest.TestCase):
    def testCalcPolygons(self):
        pass
    

class check_receptorsTestCase(unittest.TestCase):

    def test_check_receptors(self):
        """

        The Python and C versions of the receptor site hit test should return the same result

        """
        from RandomArray import uniform, seed
        from TAP_ext import check_receptors
        from time import time


        area = 200
        num_LEs = 100
        num_times = 10
        num_sites = 4

        sites = [array([(20,65),(40,35),(70,25),(75,45),(55,50),(45,75),(20,65)],Float)]*num_sites

        # build site bounding boxes
        BBs = []
        for site in sites:
            max_x = site[0,0]
            min_x = site[0,0]
            max_y = site[0,1]
            min_y = site[0,1]	
            max_x = max(max_x, max(site[:,0]))
            min_x = min(min_x, min(site[:,0]))
            max_y = max(max_y, max(site[:,1]))
            min_y = min(min_y, min(site[:,1]))
            BBs.append(array((max_x,min_x,max_y,min_y),Float))


        LEs = uniform(0,area,(num_times,num_LEs,2))

        Hit_Table1 = zeros((num_LEs, num_sites),Int)

        start = time()
        hit_test(LEs,sites,BBs,Hit_Table1,0)
        print "Python version took %.3f seconds"%(time()-start)


        Hit_Table2 = zeros((num_LEs, num_sites),Int)
        start = time()
        check_receptors.hit_test(LEs,sites,BBs,Hit_Table2,0)
        print "c version took %.3f seconds"%(time()-start)

        assert alltrue(equal(Hit_Table1,Hit_Table2)), "Python and C version gave different results"
	

from TAP_ext import NumericExtras as NE

class NumericExtrasTestCase(unittest.TestCase):

    def testFastclip(self):
        print "testing fastclip"
        A = arange(0,10,1,Float)
        B = clip(A, 3, 5)
        NE.fastclip(A, 3, 5)
        assert alltrue(A == B), "fastclip and clip gave different answers"

    def testByteswap(self):
        A = arange(10)
        B = A.copy()
        NE.byteswap(B)
        B = B.byteswapped()

        assert alltrue(A == B), "NE.byteswap and Numeric.array.byteswapped gave different results"
        
    def testChangetypeA(self):
        """
        changetype should fail for non-contiguous arrays
        """
        A = arange(18)
        A.shape = 3,6
        
        B = A[:,3]
        self.assertRaises(ValueError,NE.changetype,B,Float)

    def testChangetypeB(self):
        """
        changetype should fail for arrays the wrong size for the type
        """
        A = arange(25)

        self.assertRaises(ValueError,NE.changetype,A,Float)


    def testChangetypeC(self):
        """
        changetype(m,typecode) should have the same result as:
        m = fromstring(m.tostring(),typecode)
        """       
        A = arange(26)
        B = A.copy()
        NE.changetype(A,Float)
        assert alltrue (A == fromstring(B.tostring(),Float))



## This is the Python version of the check_receptors code, used by the test code above

def hit_test(LEs,sites,BBs,Hit_Table,Start_step):
	"""
	hit_test computes the receptor site hits given a set of LE positions,
	LEs, and the receptor sites, and the bounding boxes of the receptor sites.
	
	LEs is a M X N X 2 NumPy array (of Floats ?)
	    N is the number of LEs (Num_LEs)
		M is the number of timesteps (must be at least 2)
	sites is a list of N X 2 NumPy arrays (of Floats)
	   N is the number of points in a receptor polygon
	BBs is a list of 4 X 1 NumPy arrays (of Floats) of the bounding box of the sites (max_x,min_x,max_y,min_y)
	
	Hit_Table is a NumPy array of Int16 (short) of size (Num_LEs, Num_sites),
	it hold the values of the first timestep that the site was hit by a given LE.
	***Hit_Table is ALTERED by this function!!!***
	
	
	the function returns None
	
	"""
	
	N_LEs = LEs.shape[1]
	N_times = LEs.shape[0]
	N_sites = len(sites)
	
	for T_ind in range(1,N_times): # loop over timesteps
		for LE_ind in range(N_LEs): # loop over LEs
			LE_line = (tuple(LEs[T_ind-1,LE_ind,:]),tuple(LEs[T_ind,LE_ind,:])) # LE-movement segment
			# did the LE move?
			if (LE_line[0] != LE_line[1]):
				# check bounding boxes
				bb_LE = (max(LE_line[0][0],LE_line[1][0]),min(LE_line[0][0],LE_line[1][0]),
				                 max(LE_line[0][1],LE_line[1][1]),min(LE_line[0][1],LE_line[1][1]))
				for site_ind in range(N_sites): # loop over sites
					if BB_check(BBs[site_ind],bb_LE):
					        # do the line cross check
						for segment in map(None,sites[site_ind][:-1],sites[site_ind][1:]):
							if LCross(LE_line,segment):
								if not Hit_Table[LE_ind,site_ind]:
									Hit_Table[LE_ind,site_ind] = Start_step + T_ind
								break
	return None
	
	
	
def BB_check(bb_1, bb_2):
	"""
	bb_1 and bb_2 are two bounding boxes.
	
	Each is a 4 element tuple of :
	(max_x,min_x,max_y,min_y)
	
	BB_check(bb_1, bb_2)
	returns 1 if the two boxes intersect
	returns 0 if the two boxes don't intersect
	"""
	if ( (bb_1[0] > bb_2[1]) and (bb_1[1] < bb_2[0]) and
         (bb_1[2] > bb_2[3]) and (bb_1[3] < bb_2[2]) ):
		return 1
	else:
		return 0	
	
	
def LCross(S1,S2):
	"""
	S1 and S2 are two element tuples of two element tuples of
	x,y coordinates of the two lines:
	
	Routine to check if two line segments intersect
	
	  returns  0 if they don't intersect,   1 if they intersect
	  """
	((px1,py1),(px2,py2)) = S1
	((px3,py3),(px4,py4)) = S2
	
	# First some utility functions:
	def SideOfLineCheck(x1,y1,x2,y2,Px,Py):
		""" Given a line segment x1,y1 to x2,y2
		it checks to see if point Px,Py is to the right
		or to the left of the line segment looking from
		point x1,y1 to point x2,y2.
		If D is positive, then the point Px,Py is to the LEFT of the
		line segment.  If D is negative, P is to the right of segment.
		If D is zero then, P is on the segment
		If D =0 then that means that the point P is on the line
		defined by the two points...they may not be on the segment
		
		The check is done by taking the
		cross product of the vectors x1,y1 to x2,y2
		
		and x1,y1 to Px,Py
		"""
		
		def CrossProduct(x1,y1,x2,y2):
		                # Given vectors x1,y1 and x2,y2
		                # this routine returns the cross product
		                # which is also the determinant
		
			return x1*y2 - y1*x2
			
		dx = x2 - x1
		dy = y2 - y1
		dxp = Px - x1
		dyp = Py - y1
		
		return CrossProduct(dx,dy,dxp,dyp)
		
		# Check to see if point 3 is to the left of segment 1
		
	D1 = SideOfLineCheck(px1,py1,px2,py2,px3,py3)
	
	# Now check if point 4 is to the left of segment 1
	
	D2 = SideOfLineCheck(px1,py1,px2,py2,px4,py4)
	
	# if points 3 and 4 are on the same side of line 1
	# then things don't cross
	
	if(D1*D2 > 0):
		return  0
		# now we need to check the other way...
		
		#Check to see if point 1 is to the left of segment 2
		
	D1 = SideOfLineCheck(px3,py3,px4,py4,px1,py1)
	
	# Now check if point 2 is to the left of segment 2
	
	D2 = SideOfLineCheck(px3,py3,px4,py4,px2,py2)
	
	# if points 1 and 2 are on the same side of line 2 then things don't cross
	
	if(D1*D2 > 0):
		return 0
		#if we get here, the hummers cross
	return 1
	


if __name__ == "__main__":
##    suite()
    unittest.main()
    

    
	
	
	
	
	
	
	
	
	





