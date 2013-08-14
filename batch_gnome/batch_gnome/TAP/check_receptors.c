#include <Python.h>
#include "numpy/arrayobject.h"

// These are little macros I use to access array values for specific rank arrays
#define ARRAYVAL0(aType,a) ( *(aType *)(a->data))       
#define ARRAYVAL1(aType,a,i) ( *(aType *)(a->data + (i)*a->strides[0])) 
#define ARRAYVAL2(aType,a,i,j) ( *(aType *)(a->data + (i)*a->strides[0] + (j)*a->strides[1]))   
#define ARRAYVAL3(aType,a,i,j,k) ( *(aType *)(a->data + (i)*a->strides[0] + (j)*a->strides[1] + (k)*a->strides[2]))


/* Function Prototypes */
static double CrossProduct(double,double,double,double);
static double  SideOfLineCheck(double,double,double,double,double,double);
//static short LCross(double,double,double,double,double,double,double,double );
//static short BB_check(double,double,double,double,double,double,double,double);
//static short single_recept_check(double *, int, double[2][2]);
//static int CrossingsTest( double pgon[][2], int, double, double );


/*A diagnostic function*/
// static void print_array(PyArrayObject *array);
// 
// 
// /* The real functions */
// 
// /* A function that prints a 1d-3d NumPy array of doubles */
// static void print_array(PyArrayObject *array){
// 
//   int num_i,num_j,num_k;
//   int i,j,k;
//   
//   if (array->nd == 1){
// 	num_i = array->dimensions[0];
// 	printf("[");
// 	for (i = 0; i < num_i; i++)
// 	  {
// 		printf("%8.4f ",ARRAYVAL1(double,array,i));
// 	  }
// 	printf("]\n");
//   }
//   else if (array->nd == 2){
// 	num_i = array->dimensions[0];
// 	num_j = array->dimensions[1];
// 
// 	printf("[");
// 	for (i = 0; i < num_i; i++)
// 	  {
// 		if (i > 0) printf(" ");
// 		printf("[");
// 		for (j = 0; j < num_j; j++)
// 		  {
// 			printf("%8.4f",ARRAYVAL2(double,array,i,j));
// 			if (j == num_j-1){
// 			  printf("]");
// 			  if (i < num_i-1) printf("\n");
// 			} 
// 			else printf(" ");
// 		  }
// 	  }
// 	printf("]\n");
//   }
//   else if (array->nd == 3){
// 	num_i = array->dimensions[0];
// 	num_j = array->dimensions[1];
// 	num_k = array->dimensions[2];
// 	printf("dimensions of array are: (%i,%i,%i)\n",num_i,num_j,num_k);
// 
// 	printf("[");
// 	for (i = 0; i < num_i; i++){
// 		if (i > 0) printf(" ");
// 		printf("[");
// 		for (j = 0; j < num_j; j++){
// 			if (j > 0) printf("  ");
// 			printf("[");
// 			for (k = 0; k < num_k; k++){
// 				printf("%8.4f",ARRAYVAL3(double,array,i,j,k));
// 				if (k == num_k-1){
// 				  printf("]");
// 				  if (j < num_j-1) printf("\n");
// 				} 
// 				else printf(" ");
// 			  }
// 			if (j == num_j-1){
// 			  printf("]");
// 			  if (i < num_i-1) printf("\n");
// 			}
// 			/* else printf("]");*/
// 		  }
// 	  }
// 	printf("]\n");
//   }
//   else{
// 	printf("An array with more than 3 dimensions was passed into \"print_array\"") ;
//   }
// }
// 
// 
// /*First Define some utility functions*/
// 
// 
// static short single_recept_check(double site[], int num_points, double LE_line[2][2])
// 	 /* site is a num_points X 2 array of doubles, with the lat/long coordinates of the site polygon */
// 	 /* LE_line is a 4 element array of doubles, with the coordinates of the LE at the beginning and end of the time step */
// 
// {
//   int n;
//   double x1, y1, x2, y2;
// 
//   for (n = 0; n < num_points-1; n++)
// 	{
// /* 	  printf("n = %i\n",n); */
// 	  x1 = site[n*2 + 0];
// 	  y1 = site[n*2 + 1];
// 	  x2 = site[(n+1)*2 + 0];
// 	  y2 = site[(n+1)*2 + 1];
// 
// /* 	  printf("x1 = %f, y1 = %f, x2 = %f, y2 = %f\n",x1,y1,x2,y2); */
// 
// 
// 	  if (LCross(x1, y1, x2, y2, LE_line[0][0], LE_line[0][1], LE_line[1][0], LE_line[1][1]))
// 		return 1;
// 	  }
//   return 0;
// }
// 
// 
// 
// static short BB_check(double max_x1,double min_x1,double max_y1,double min_y1,
// 					  double max_x2,double min_x2,double max_y2,double min_y2)
// {
// 
// /* 	bb_1 and bb_2 are two bounding boxes. */
// 
// /* 	Each is a 4 element tuple of : */
// /* 	(max_x,min_x,max_y,min_y) */
// 
// /* 	BB_check(bb_1, bb_2) */
// /* 	returns 1 if the two boxes intersect */
// /* 	returns 0 if the two boxes don't intersect */
// 
//   if ( (max_x1 > min_x2) && (min_x1 < max_x2) &&
// 	   (max_y1 > min_y2) && (min_y1 < max_y2) )
// 	return 1;
//   else
// 	return 0;
// 
// }



static double CrossProduct(double x1, double x2, double y1, double y2)

{
  return (x1*y2 - y1*x2);
}

/* Function to check which side of a line a point is on*/

static double  SideOfLineCheck(double x1,double y1,double x2,double y2,double Px,double Py)
{
  /* Given a line segment x1,y1 to x2,y2
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
  */


  double dx, dy, dxp, dyp;


  dx = x2 - x1;
  dy = y2 - y1;
  dxp = Px - x1;
  dyp = Py - y1;
  
  return CrossProduct(dx,dy,dxp,dyp);
}
// 
// static short LCross(double px1,double py1,double px2,double py2,double px3,double py3,double px4,double py4)
// 
// {
//   /* 	Routine to check if two line segments intersect */
//   
//   /* 	  returns  0 if they don't intersect,   1 if they intersect */
// 
//   double D1,D2;
// 
// 	/* 	 Check to see if point 3 is to the left of segment 1 */
// 	
//   D1 = SideOfLineCheck(px1,py1,px2,py2,px3,py3);
// 	
//   /* 	 Now check if point 4 is to the left of segment 1 */
// 	
//   D2 = SideOfLineCheck(px1,py1,px2,py2,px4,py4);
// 	
//   /* 	 if points 3 and 4 are on the same side of line 1 */
//   /* 	 then things don't cross */
// 
//   if (D1*D2 > 0)
// 	return  0 ;
//   /* 	 now we need to check the other way... */
// 		
//   /* 	Check to see if point 1 is to the left of segment 2 */
// 		
//   D1 = SideOfLineCheck(px3,py3,px4,py4,px1,py1);
// 	
//   /* 	 Now check if point 2 is to the left of segment 2 */
// 	
//   D2 = SideOfLineCheck(px3,py3,px4,py4,px2,py2);
// 	
//   /* 	 if points 1 and 2 are on the same side of line 2 then things don't cross */
// 	
//   if(D1*D2 > 0)
// 	return 0;
//   /* 		if we get here, the hummers cross */
//   return 1;
// }
// 
// 
// /* Now the main loop code: */
// /* called from Python as: hit_test(LE1,LE2,sites,BBs,Hit_Table,Start_step)*/
// 
// static PyObject * check_receptors_hit_test(PyObject *self, PyObject *args)
// {
//   
//   int site_ind, LE_ind, T_ind;
//   int N_sites, N_LEs, N_times;
//   int result;
// 
//   short Start_step = 0; //, Time_step = 0;
// 
//   double LE_line[2][2];
//   double px, py ;
// 
//   double max_x1,min_x1,max_y1,min_y1;
// 
//   PyObject *sites;
//   PyObject *BBs;
// 
//   PyArrayObject *site;
//   PyArrayObject *BB_site;
//   PyArrayObject *LEs;
//   PyArrayObject *Hit_Table;
// 
//   result = 1;
//   
//   /*   printf("I'm starting\n"); */
// 
//   if (!PyArg_ParseTuple(args, "O!O!O!O!h", 
// 						&PyArray_Type, &LEs,
// 						&PyList_Type, &sites,
// 						&PyList_Type, &BBs,
// 						&PyArray_Type, &Hit_Table,
// 						&Start_step))
// 	return NULL;
// 
//   N_sites = PyList_Size(sites);
//   /* Do some input type checking. Note that only the first item in the list is checked.
// 	 It is assumed that the input list is homogenous */
//   /* check sizes*/
//   if (N_sites != PyList_Size(BBs) || N_sites == 0 ) {
// 	PyErr_SetString (PyExc_ValueError,
// 					 "sites and BBs must be list the same length, and at least one item long");
// 	return NULL;
//   }
//  
//   /* check for lists of arrays */
//   if (PyArray_Check(PyList_GetItem(sites, 0))) {
// 	site = (PyArrayObject *) PyList_GetItem(sites, 0);
// 	if (site->descr->type_num != PyArray_DOUBLE) {
// 	  PyErr_SetString (PyExc_ValueError,
// 					   "sites must be a list of NumPy arrays of type Float");
// 	  return NULL;
// 	}
//   }
//   else {
// 	PyErr_SetString (PyExc_ValueError,
// 					 "sites  be a list of NumPy arrays");
// 	return NULL;
//   }
//   if (PyArray_Check(PyList_GetItem(BBs, 0))) {
// 	site = (PyArrayObject *) PyList_GetItem(BBs, 0);
// 	if (site->descr->type_num != PyArray_DOUBLE) {
// 	  PyErr_SetString (PyExc_ValueError,
// 					   "BBs must be a list of NumPy arrays of type Float");
// 	  return NULL;
// 	}
//   }
//   else {
// 	PyErr_SetString (PyExc_ValueError,
// 					 "BBs  be a list of NumPy arrays");
// 	return NULL;
//   }
// 
//   /* check type of LEs */
//   if (LEs->descr->type_num != PyArray_DOUBLE){
// 	PyErr_SetString (PyExc_ValueError,
// 					 "LEs must be a NumPy array of type Float");
// 	return NULL;
//   }
//   else if ((LEs->dimensions[2] != 2)
// 		   || (LEs->dimensions[0] < 2)) {
// 	PyErr_SetString (PyExc_ValueError,
// 					 "LEs must be M X N X 2 NumPy array, with M >=2");
// 	return NULL;
//   }
//   N_times = LEs->dimensions[0];
//   N_LEs = LEs->dimensions[1];
// 	
//   /* check Hit_Table */
// 
// 
//   if (Hit_Table->descr->type_num != PyArray_LONG){
// //	char str[50];
// //	sprintf(str,"%X  ---  %X", PyArray_LONG, Hit_Table->descr->type_num);
// 
// 	PyErr_SetString (PyExc_ValueError, //str);
// 					 "Hit_Table must be a NumPy array of type Int");
// 	return NULL;
//   }
//   else if ((Hit_Table->dimensions[0] != N_LEs) || (Hit_Table->dimensions[1] != N_sites)) {
// 	PyErr_SetString (PyExc_ValueError,
// 					 "Hit_Table must be a Num_LEs X Num_sites NumPy Array");
// 	return NULL;
//   }
//   
//   /*   printf("N_sites = %i\n",N_sites); */
//   /*   printf("N_LEs = %i\n",N_LEs); */
// 
//   /*   printf("LE1 is %i by %i\n",LE1->dimensions[0],LE1->dimensions[1]); */
//   /*   printf("LE2 is %i by %i\n",LE2->dimensions[0],LE2->dimensions[1]); */
// 
//   /*   printf("about to start looping\n"); */
//   /*   return Py_BuildValue("i",1);   */
//   //  Time_step = Start_step;
// 
// 
//   /* Put a Point-in-polygon check for the first time step. */
//   if (Start_step == 0) {
//     //printf("in the first step loop\n");
// 	for (LE_ind = 0; LE_ind < N_LEs; LE_ind++) { /*loop over LEs */
//       px = ARRAYVAL3(double,LEs,0,LE_ind,0);
//       py = ARRAYVAL3(double,LEs,0,LE_ind,1);
//       //printf("LE point is: %g, %g\n",px,py);
//       for (site_ind = 0; site_ind < N_sites; site_ind++){ /* loop over sites*/
//         site = (PyArrayObject *) PyList_GetItem(sites, site_ind);
//         if (CrossingsTest( ( void* )site->data,
//                            site->dimensions[0],px,py)) {
//           //printf("LE: %i is in site: %i\n",LE_ind,site_ind);
//           ARRAYVAL2(long,Hit_Table,LE_ind,site_ind) = 1;
//         }
//       }
//     }
//   }
// 
//   for (T_ind = 1; T_ind < N_times; T_ind++) { /*Loop Over timesteps*/
//     //printf("in the time step loop, step: %i\n",T_ind);
// 	for (LE_ind = 0; LE_ind < N_LEs; LE_ind++) { /*loop over LEs */
// 	  LE_line[0][0] = ARRAYVAL3(double,LEs,T_ind-1,LE_ind,0);
// 	  LE_line[0][1] = ARRAYVAL3(double,LEs,T_ind-1,LE_ind,1);
// 	  LE_line[1][0] = ARRAYVAL3(double,LEs,T_ind,LE_ind,0);
// 	  LE_line[1][1] = ARRAYVAL3(double,LEs,T_ind,LE_ind,1);
// 	  /* has the LE moved? */
// 	  if ((LE_line[0][0] != LE_line[1][0])||(LE_line[0][1] != LE_line[1][1])) {
// 		/* compute bounding box of LE_line */
// 		if (LE_line[0][0]>LE_line[1][0])
// 		  {
// 			max_x1 = LE_line[0][0];
// 			min_x1 = LE_line[1][0];
// 		  }
// 		else
// 		  {
// 			max_x1 = LE_line[1][0];
// 			min_x1 = LE_line[0][0];
// 		  }
// 		if (LE_line[0][1]>LE_line[1][1])
// 		  {
// 			max_y1 = LE_line[0][1];
// 			min_y1 = LE_line[1][1];
// 		  }
// 		else
// 		  {
// 			max_y1 = LE_line[1][1];
// 			min_y1 = LE_line[0][1];
// 		  }
// 		for (site_ind = 0; site_ind < N_sites; site_ind++){ /* loop over sites*/
// 		  site = (PyArrayObject *) PyList_GetItem(sites, site_ind);
// 		  BB_site = (PyArrayObject *) PyList_GetItem(BBs, site_ind);
// 		  /* do the BB check */
// 		  if (BB_check(max_x1,min_x1,max_y1,min_y1,
// 					   ARRAYVAL1(double,BB_site,0),
// 					   ARRAYVAL1(double,BB_site,1),
// 					   ARRAYVAL1(double,BB_site,2),
// 					   ARRAYVAL1(double,BB_site,3))){
// 			if (single_recept_check((double * )site->data, site->dimensions[0], LE_line)){
// 			  if (ARRAYVAL2(long,Hit_Table,LE_ind,site_ind) < 1) {
// 				ARRAYVAL2(long,Hit_Table,LE_ind,site_ind) = T_ind + Start_step;
// 			  }				
// 			}
// 		  }
// 		}
// 	  }
// 	}
//   }  
//   return Py_BuildValue("i",result);
// }
  


/* Now the main loop code for grid receptors: */
/* called from Python as: Grid_hit_test(LEs, grid, Hit_Table, Start_step)*/

static PyObject * check_receptors_Grid_hit_test(PyObject *self, PyObject *args){
  int site_ind, LE_ind, T_ind;
  int N_LEs, N_times;
  int result;
  int temp;
  int i0, j0, i1, j1;
  int i,j;



  /* grid object stuff*/
  //  PyObject *grid;
  double min_long, max_long, min_lat, max_lat, dlat, dlong;
  long  num_lat, num_long; 
  
  long N_sites;

  short Start_step = 0;

  double px,py;
  double px1, py1, px2, py2;
  double s;
  int side;
 
  PyArrayObject *LEs;
  PyArrayObject *Hit_Table;

  result = 1;
  

  if (!PyArg_ParseTuple(args, "O!(ddddll)O!h", 
						&PyArray_Type, &LEs,
						&min_long,
                        &max_long, 
                        &min_lat,
                        &max_lat,
                        &num_lat,
                        &num_long,
						&PyArray_Type, &Hit_Table,
						&Start_step))
	return NULL;

  dlat = (max_lat-min_lat) / num_lat ;
  dlong = (max_long-min_long) / num_long ;

  //  printf("num_lat = %i, num_long = %i, min_lat = %g, max_lat = %g, min_long = %g, max_long  = %g\n",
  //         num_lat, num_long, min_lat, max_lat, min_long, max_long);

  N_sites = num_lat*num_long;

  /* Do some input type checking. Note that only the first item in the list is checked.
	 It is assumed that the input list is homogenous */

  /* check type of LEs */
  if (LEs->descr->type_num != PyArray_DOUBLE){
	PyErr_SetString (PyExc_ValueError,
					 "LEs must be a NumPy array of type Float");
	return NULL;
  }
  else if ((LEs->dimensions[2] != 2)
		   || (LEs->dimensions[0] < 2)) {
	PyErr_SetString (PyExc_ValueError,
					 "LEs must be M X N X 2 NumPy array, with M >=2");
	return NULL;
  }
  N_times = LEs->dimensions[0];
  N_LEs = LEs->dimensions[1];
	
  /* check Hit_Table */


  if (Hit_Table->descr->type_num != PyArray_LONG){
//	char str[50];
//	sprintf(str,"%X  ---  %X", PyArray_LONG, Hit_Table->descr->type_num);

	PyErr_SetString (PyExc_ValueError, //str);
					 "Hit_Table must be a NumPy array of type Int");
	return NULL;
  }
  else if ((Hit_Table->dimensions[0] != N_LEs) || (Hit_Table->dimensions[1] != N_sites)) {
	PyErr_SetString (PyExc_ValueError,
					 "Hit_Table must be a Num_LEs X Num_sites NumPy Array");
	return NULL;
  }
  
  //printf("About to start the real code\n");
  //printf("N_times = %i, N_LEs = %i\n",N_times,N_LEs);
  /*Now the real code!!*/
  for (T_ind = 1; T_ind < N_times; T_ind++) { /*Loop Over timesteps*/
	for (LE_ind = 0; LE_ind < N_LEs; LE_ind++) { /*loop over LEs */
	  px1 = ARRAYVAL3(double,LEs,T_ind-1,LE_ind,0);
	  py1 = ARRAYVAL3(double,LEs,T_ind-1,LE_ind,1);
	  px2 = ARRAYVAL3(double,LEs,T_ind,LE_ind,0);
	  py2 = ARRAYVAL3(double,LEs,T_ind,LE_ind,1);
	  /* has the LE moved, and is it inside the grid (both time steps)*/
	  if ( ((px1 != px2)||(py1 != py2)) &&
           (px1 >= min_long) && (px1 <= max_long) &&
           (py1 >= min_lat)  && (py1 <= max_lat) &&
           (px2 >= min_long) && (px2 <= max_long) &&
           (py2 >= min_lat)  && (py2 <= max_lat) ) {
		/* what box is the first point in? */
        i0 = ((px1 - min_long) / dlong);
        j0 = ((py1 - min_lat ) / dlat) ;
          /* what box is second point in? */
        i1 = ((px2 - min_long) / dlong);
        j1 = ((py2 - min_lat ) / dlat) ;

        if (i0 > i1){
          temp = i0;
          i0 = i1;
          i1 = temp;
        }
        if (j0 > j1){
          temp = j0;
          j0 = j1;
          j1 = temp;
        }
        //printf("i0,i1 = %i,%i,  j0,j1 = %i,%i\n",i0,i1,j0,j1);
        // check the boxes:
        for (i = i0; i < i1+1; i++){
          for (j = j0; j < j1+1; j++){
            //First point
            px = min_long + i*dlong ;
            py = min_lat + j*dlat ;
            site_ind = i*num_lat + j ;
            //printf("i = %i, j = %i, site_ind = %i\n",i,j,site_ind);

            s = SideOfLineCheck(px1,py1,px2,py2,px,py);
            if (s > 0 ){
              side = 1;
            }                      
            else{
              side = 0;
            }
            //if (LE_ind == 4 ) printf("i = %i, j= %i,s = %g\n",i,j,s);
            // This is sort of kludgy. there should be a loop or function call.
            //second point
            px = min_long + (i+1)*dlong ;
            py = min_lat + j*dlat ;
            s = SideOfLineCheck(px1,py1,px2,py2,px,py);
            if ((s > 0) != side) {
              if (ARRAYVAL2(long,Hit_Table,LE_ind,site_ind) < 1) {
                ARRAYVAL2(long,Hit_Table,LE_ind,site_ind) = T_ind + Start_step;
              }
            }
            else {
              // third point
              px = min_long + (i+1)*dlong ;
              py = min_lat  + (j+1)*dlat ;
              s = SideOfLineCheck(px1,py1,px2,py2,px,py);
              if ((s > 0) != side) {
                if (ARRAYVAL2(long,Hit_Table,LE_ind,site_ind) < 1) {
                  ARRAYVAL2(long,Hit_Table,LE_ind,site_ind) = T_ind + Start_step;
                }
              }
              else {
                // fourth point
                px = min_long + i*dlong;
                py = min_lat + (j+1)*dlat;
                s = SideOfLineCheck(px1,py1,px2,py2,px,py);
                if ((s > 0) != side) {
                  if (ARRAYVAL2(long,Hit_Table,LE_ind,site_ind) < 1) {
                    ARRAYVAL2(long,Hit_Table,LE_ind,site_ind) = T_ind + Start_step;
                  }
                }
              }
            }
            //p3 = (min_long + (i+1)*dlong, min_lat + (j+1)*dlat);
            //p4 = (min_long + i*dlong, min_lat + (j+1)*dlat);
            
          }
        }
      }  
    }
  }  
  return Py_BuildValue("i",result);
}
    
  
// static PyObject * check_receptors_single_recept_check(PyObject *self, PyObject *args)
// {
// 
//   double LE_line[2][2];
//   double *site;
//   int num_points;
// 
//   PyArrayObject *site_array;
//   
//   short result;
// 
//   if (!PyArg_ParseTuple(args, "O((dd)(dd))",&site_array,&(LE_line[0][0]),&(LE_line[0][1]),&(LE_line[1][0]),&(LE_line[1][1]) ))
// 	return NULL;
//   if ((site_array->descr->type_num != PyArray_DOUBLE) || (site_array->dimensions[1] != 2))
// 	{
// 	  PyErr_SetString(PyExc_ValueError,
// 					  "Array must be of type Float, and be of shape (N X 2)");
// 	  return NULL;
// 	}
// 
// 
//   site = (double *) site_array->data;
// 
//   num_points = site_array->dimensions[0];
// 
//   result = single_recept_check(site, num_points, LE_line);
// 
//   return Py_BuildValue("i",result);
// 
// }
// 
// 
// static PyObject * check_receptors_LCross(PyObject *self, PyObject *args)
// {
// 
//   double px1,py1,px2,py2,px3,py3,px4,py4;
//   short result;
// 
//   if (!PyArg_ParseTuple(args, "((dd)(dd))((dd)(dd))", &px1, &py1, &px2, &py2, &px3, &py3, &px4, &py4 ))
// 	return NULL;
//   
//   result = LCross(px1,py1,px2,py2,px3,py3,px4,py4);
//   return Py_BuildValue("i",result);
// }
// 
// static PyObject * check_receptors_BB_check(PyObject *self, PyObject *args)
// {
//   double max_x1, min_x1, max_y1, min_y1,max_x2, min_x2, max_y2, min_y2;
// 
//   short result;
// 
//   if (!PyArg_ParseTuple(args, "(dddd)(dddd)", &max_x1, &min_x1, &max_y1, &min_y1, &max_x2, &min_x2, &max_y2, &min_y2 ))
// 	return NULL;
//   
//   result = BB_check(max_x1, min_x1, max_y1, min_y1,max_x2, min_x2, max_y2, min_y2);
//   return Py_BuildValue("i",result);
// }
// 
// static PyObject * check_receptors_point_in_polygon(PyObject *self, PyObject *args)
// {
// 
//   int numverts;
// 
//   double px,py;
// 
//   PyArrayObject *site_array;
//   
//   int result;
// 
//   if (!PyArg_ParseTuple(args, "O(dd)",&site_array,&px,&py)){
// 	return NULL;
//   }
//   if ((site_array->descr->type_num != PyArray_DOUBLE) || (site_array->dimensions[1] != 2))
// 	{
// 	  PyErr_SetString(PyExc_ValueError,
// 					  "Array must be of type Float, and be of shape (N X 2)");
// 	  return NULL;
// 	}
// 
//   numverts = site_array->dimensions[0];
// 
//   result = CrossingsTest((void *) (site_array->data), numverts, px,py);
// 
//   return Py_BuildValue("i",result);
// }
// 
// 
// static PyObject * check_receptors_print_array(PyObject *self, PyObject *args)
// {
// 
//   PyArrayObject *array;
// 
//   if (!PyArg_ParseTuple(args, "O!", &PyArray_Type, &array))
// 	return NULL;
// 
//   print_array(array);
//   
//   return Py_BuildValue("i",1);
// }
// 
// /* The Point in Polygon code:
// 
//    I got this off the web from "Graphics Gems IV"
//    http://www.graphicsgems.org/
// */
// /* ptinpoly.c - point in polygon inside/outside code.
// 
//    by Eric Haines, 3D/Eye Inc, erich@eye.com
// 
//    This code contains the following algorithms:
// 	crossings - count the crossing made by a ray from the test point
// 	crossings-multiply - as above, but avoids a division; often a bit faster
// 	angle summation - sum the angle formed by point and vertex pairs
// 	weiler angle summation - sum the angles using quad movements
// 	half-plane testing - test triangle fan using half-space planes
// 	barycentric coordinates - test triangle fan w/barycentric coords
// 	spackman barycentric - preprocessed barycentric coordinates
// 	trapezoid testing - bin sorting algorithm
// 	grid testing - grid imposed on polygon
// 	exterior test - for convex polygons, check exterior of polygon
// 	inclusion test - for convex polygons, use binary search for edge.
// */
// #define X	0
// #define Y	1
// /* ======= Crossings algorithm ============================================ */
// 
// /* Shoot a test ray along +X axis.  The strategy, from MacMartin, is to
//  * compare vertex Y values to the testing point's Y and quickly discard
//  * edges which are entirely to one side of the test ray.
//  *
//  * Input 2D polygon _pgon_ with _numverts_ number of vertices and test point
//  * _point_, returns 1 if inside, 0 if outside.	WINDING and CONVEX can be
//  * defined for this test.
//  */
// 
// int CrossingsTest(double pgon[][2], int numverts, double tx, double ty )
// {
// 
//   register int	j, yflag0, yflag1, inside_flag, xflag0 ;
//   register double *vtx0, *vtx1 ;
// 
//   vtx0 = pgon[numverts-1] ;
//   /* get test bit for above/below X axis */
//   yflag0 = ( vtx0[Y] >= ty ) ;
//   vtx1 = pgon[0] ;
// 
//   inside_flag = 0 ;
// 
//   for ( j = numverts+1 ; --j ; ) {
// 
// 	yflag1 = ( vtx1[Y] >= ty ) ;
// 	/* check if endpoints straddle (are on opposite sides) of X axis
// 	 * (i.e. the Y's differ); if so, +X ray could intersect this edge.
// 	 */
// 	if ( yflag0 != yflag1 ) {
//       xflag0 = ( vtx0[X] >= tx ) ;
//       /* check if endpoints are on same side of the Y axis (i.e. X's
//        * are the same); if so, it's easy to test if edge hits or misses.
//        */
//       if ( xflag0 == ( vtx1[X] >= tx ) ) {
// 		/* if edge's X values both right of the point, must hit */
// 		if ( xflag0 ) inside_flag = !inside_flag ;
// 
//       } else {
// 		/* compute intersection of pgon segment with +X ray, note
// 		 * if >= point's X; if so, the ray hits it.
// 		 */
// 		if ( (vtx1[X] - (vtx1[Y]-ty)*
//               ( vtx0[X]-vtx1[X])/(vtx0[Y]-vtx1[Y])) >= tx ) {
// 
//           inside_flag = !inside_flag ;
// 		}
//       }
// 	}
// 	/* move to next pair of vertices, retaining info as possible */
// 	yflag0 = yflag1 ;
// 	vtx0 = vtx1 ;
// 	vtx1 += 2 ;
//   }
//   return( inside_flag ) ;
// }





/* The Set up for the Python Module*/	
static PyMethodDef check_receptorsMethods[] = {
//  {"LCross",    check_receptors_LCross,  METH_VARARGS},
//  {"BB_check",    check_receptors_BB_check,  METH_VARARGS},
//  {"hit_test",    check_receptors_hit_test,  METH_VARARGS},
  {"Grid_hit_test",    check_receptors_Grid_hit_test,  METH_VARARGS},
//  {"print_array",    check_receptors_print_array,  METH_VARARGS},
//  {"single_recept_check",    check_receptors_single_recept_check,  METH_VARARGS},
//  {"point_in_polygon",    check_receptors_point_in_polygon,  METH_VARARGS},
  {NULL,      NULL}        /* Sentinel */
};

void initcheck_receptors(void)
{
  (void) Py_InitModule("check_receptors", check_receptorsMethods);
  import_array()
}
