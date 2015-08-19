from scipy.stats.distributions import norm
from scipy import spatial

# x = norm(-1,1).rvs(1000)
# x2 = norm(-.3,.2).rvs(1000)
# y = cos(x)+x2


# data = zip(x,y)
# t_ind = tree.query_ball_point([0,0],1)
# d,i = tree.query(data,5)


def plumeKDE(data,mass,out_points,N_neighbors = 20):


	# test for 2D/3D, assume 2D for now 
	dataxy = data

	# call to KDtree
	tree = spatial.KDTree(dataxy)

	tr_dist,tr_indx = tree.query(out_points,N_neighbors)
	# do something special for d=0 when using particle psitions as out_points?




	print tr_dist

	# use nearest-neighbors distances to define lengths of KDE


	# sum up KDEs at input locations (using particle masses to scale KDE)

	return tr_dist,tr_indx

