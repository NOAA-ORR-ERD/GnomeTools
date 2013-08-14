// This file was converted from Gennady's original C++ to C

static PyObject *CalcPolygons_GetVolumes(PyObject *self, PyObject *args)
////////////////////////////////////////////////////////////////////////////////
// calculate oil volume for polygon based on LE mapping
//
// input data:
//		GridMatrix  - a map-matrix of polygons with rank [xs,ys] - each ceil is 
//						  a number (int) of the polygon contained this ceil.
//						  It must be a NumPy array
//		BB				- a bonding box (2x2 - x,y geo-coordinates) for map region,
//						  must be a NumPy array
//		LEs			- a NumPy array of LE coordinates - longitude, latitude
//
//		LEsVol		- a NumPy array of oil volumes for correspondent LE
//
//      TrajFlag    - The flag that that LE has in the trajectory (used to indicate beached, off map, etc.
//                    It is a NumPy array of UInt8, of size (N,) 
// output data:
//		PolVol		- a NumPy array with dimension NumPolygons - oil volumes
//
{

	PyArrayObject *GridMatrix;
	PyArrayObject *LEs;
	PyArrayObject *LEsVol;
	PyArrayObject *PolVol;
	PyArrayObject *BB;
	PyArrayObject *TrajFlag;
	//PyListObject  *OutBinLEs;

    int GridSize;   
    int NumLEs;     
    int NumPolygons;

    int xs;
    int ys;

	double xLE,yLE;								// current LE coordinates
	double x0;	// x-start search point
	double y0;	// y-start search point
	double xMax;
	double yMax;
	double dx;	// deltaX - x-step
	double dy;	// deltaY - y-step
	int ix;											// x-index of map array
	int iy;											// y-index of map array
	int iPol;										// index of polygon
	int i;

    PyObject *LENumber ;

    PyObject *EmptyList;
	PyObject *BinLEs;

	// get input/output parameters
   if(!PyArg_ParseTuple(args,"O!O!O!O!O!O!",
								&PyArray_Type, &GridMatrix,
								&PyArray_Type, &BB,
								&PyArray_Type, &LEs,
								&PyArray_Type, &LEsVol,
								&PyArray_Type, &TrajFlag,
								&PyArray_Type, &PolVol))
		return NULL;

	/* do some input checking.  */

	/* check type of GridMatrix */
	if(!PyArray_Check(GridMatrix))
	{
		PyErr_SetString (PyExc_ValueError,
			"GridMatrix must be a NumPy array");
		return NULL;
	}
	else
		if (GridMatrix->descr->type_num != PyArray_SHORT) 
		{
			PyErr_SetString(PyExc_ValueError,
				"Grid Matrix must be a NumPy array of type Int16");
			return NULL;
		}

	/* check type of LEs */
   if(!PyArray_Check(LEs)) 
	{
		PyErr_SetString (PyExc_ValueError,
			"LEs must be a NumPy array");
		return NULL;
	}
	else
		if(LEs->descr->type_num != PyArray_DOUBLE)
		{
			PyErr_SetString (PyExc_ValueError,
				"LEs must be a NumPy array of type Float");
			return NULL;
		}
		else 
			if ((LEs->dimensions[1] != 2) || (LEs->dimensions[0] < 2)) 
			{
				PyErr_SetString (PyExc_ValueError,
					"LEs must be M X 2 NumPy array, with M >=2");
				return NULL;
			}

	/* check type of LEsVol */
	if(!PyArray_Check(LEsVol))
	{
		PyErr_SetString (PyExc_ValueError,
			"LEsVol must be a NumPy array");
		return NULL;
	}
	else
		if (LEsVol->descr->type_num != PyArray_DOUBLE)
		{
			PyErr_SetString(PyExc_ValueError,
				"LEsVol must be a NumPy array of type Float");
			return NULL;
		}

	/* check type of BB array */
	if(!PyArray_Check(BB))
	{
		PyErr_SetString (PyExc_ValueError,
			"BB must be a NumPy array");
		return NULL;
	}
	else
		if(BB->descr->type_num != PyArray_DOUBLE)
		{
			PyErr_SetString(PyExc_ValueError,
				"BB must be a NumPy array of type Float");
			return NULL;
		}
		else
			if((BB->dimensions[0] != 2) || (BB->dimensions[1] != 2)) 
			{
				PyErr_SetString (PyExc_ValueError,
					"BB must be 2 X 2 NumPy array");
				return NULL;
			}


	/* check type of Trajectory Flag */
	if(!PyArray_Check(TrajFlag))
	{
		PyErr_SetString (PyExc_ValueError,
			"Trajectory Flag must be a NumPy array");
		return NULL;
	}
	else
		if (TrajFlag->descr->type_num != PyArray_UBYTE)
		{
			PyErr_SetString(PyExc_ValueError,
				"Trajectory Flag must be a NumPy array of type UBYTE");
			return NULL;
		}

	GridSize	 = PyArray_Size((PyObject*)GridMatrix);
	NumLEs		 = LEs->dimensions[0];
	NumPolygons = PyArray_Size((PyObject*)PolVol);

	// check array sizes
	if(GridSize == 0) 
	{
		PyErr_SetString(PyExc_ValueError,
			"Grid Array must be at least one item long");
		return NULL;
	}

	if(NumLEs != PyArray_Size((PyObject*)LEsVol) ||
		NumLEs != PyArray_Size((PyObject*)TrajFlag)||
		NumLEs == 0)
	{
		PyErr_SetString(PyExc_ValueError,
			"LEsVolumes and Trajectory Flag must be array the same length, and at least one item long");
		return NULL;
	}

	if(NumPolygons == 0)
	{
		PyErr_SetString(PyExc_ValueError,
			"There must be at least one polygon in the polygons array");
		return NULL;
	}

	/* end of input checking */


	xs = GridMatrix->dimensions[0];		// x-dimension of map array
	ys = GridMatrix->dimensions[1];		// y-dimension of map array

	x0 = ARRAYVAL2(double,BB,0,0);	// x-start search point
	y0 = ARRAYVAL2(double,BB,0,1);	// y-start search point
	xMax = ARRAYVAL2(double,BB,1,0);
	yMax = ARRAYVAL2(double,BB,1,1);
	dx = (xMax - x0)/xs;	// deltaX - x-step
	dy = (yMax - y0)/ys;	// deltaY - y-step


	//	FILE* file = fopen("polygons.log","a");
	
//	fprintf(file,"dx=%g dy=%g\n",dx,dy);
//	fprintf(file,"x0=%g y0=%g\n",x0,y0);

	// initialize polygons volume array and dimensions of output array
	for(i=0; i<NumPolygons; i++)
		ARRAYVAL1(double,PolVol,i) = 0.0;


	BinLEs = PyList_New(0);		// output list
	for(i=0; i<NumPolygons; i++){
      EmptyList = PyList_New(0);
      PyList_Append(BinLEs,EmptyList);// initialize output list
      Py_DECREF(EmptyList);
    }


//	fprintf(file,"NumPolygons=%d\n",NumPolygons);
//	fprintf(file,"NumLEs=%d\n",NumLEs);

	// calculate polygons volume
	for(i=0; i<NumLEs; i++)						// for each LE
		if(ARRAYVAL1(double,LEsVol,i) > 0.)	{// if LE's Volume > 0
		
			xLE = ARRAYVAL2(double,LEs,i,0);
			yLE = ARRAYVAL2(double,LEs,i,1);

			if(ARRAYVAL1(int,TrajFlag,i) & 4)// LE is out of map
			{
              // printf("LE (#%i,%.6f,%.6f) is out of map\n",i,xLE,yLE);
              continue;
			}

			ix = (int) (floor((xLE-x0)/dx));
			iy = (int) (floor((yLE-y0)/dy));
			
			if(ix > xs || ix < 0 || iy > ys || iy < 0){
              printf("LE (#%i,%.6f,%.6f) is out of bins region\n",i,xLE,yLE);
             continue;
            }
			iPol = ARRAYVAL2(short,GridMatrix,ix,iy);	// find polygon index corresponds this LE
            
            if(ARRAYVAL1(int,TrajFlag,i) & 2){	// if it is beached
              continue;
            }
            
			if(iPol) { // LE is inside of some polygon
              if((iPol > NumPolygons) || (iPol < 1)){	// polygon number is in limits?
                PyErr_SetString(PyExc_ValueError,
                                "Polygon number is out of range");
                return NULL;
              }
              else
                {
                  ARRAYVAL1(double,PolVol,iPol-1) += ARRAYVAL1(double,LEsVol,i);
                  //PyList_Append(PyList_GetItem(BinLEs,iPol-1),PyInt_FromLong(i));
                  LENumber = PyInt_FromLong(i);
                  PyList_Append(PyList_GetItem(BinLEs,iPol-1),LENumber);
                  Py_DECREF(LENumber);
                }
              //		fprintf(file,"iPol=%d PolVol[i]=%g\n",iPol, ARRAYVAL1(double,PolVol,iPol-1)); 
            }
		}
    

    // fclose(file);

	//for(i=0; i<NumPolygons; i++)
	//	ARRAYVAL1(double,OutBinLEs,i) = PyList_GetItem(BinLEs,i);

	//delete BinLEs

	//return Py_BuildValue("i",0);//BinLEs);
	return Py_BuildValue("N",BinLEs);
}


static PyObject *CalcPolygons_GetArea(PyObject *self, PyObject *args)
////////////////////////////////////////////////////////////////////////////////////
// calculate area of the polygon
//
// input data: Polygon - sequense of (x,y - double) vertex coordinates of polygon.
//					The last vertex is assumed to be as the first - the polygon is closed
//					THE POLYGON MUST NOT BE SELF INTERSECTING!
//					Polygon must be an M X 2 NumPy array with M >= 4
//
// returning data - area of polygon (double)


{
	PyArrayObject *Polygon;	// polygon array
	int N;						// number of vertexes
	int i,j;
	double area = 0.;
    const double pi  = 3.1415926535897931;
	const double mpd = 1852.*60.;	// meters per degree

    double sumy = 0.0;
    double avgy;
    double xi,xj,yi,yj;

   if(!PyArg_ParseTuple(args,"O!",
								&PyArray_Type, &Polygon))
	return NULL;

	/* check type of Polygon */
   if(!PyArray_Check(Polygon)) 
	{
		PyErr_SetString (PyExc_ValueError,
			"Polygon must be a NumPy array");
		return NULL;
	}
	else
		if(Polygon->descr->type_num != PyArray_DOUBLE)
		{
			PyErr_SetString (PyExc_ValueError,
				"Polygon must be a NumPy array of type Float");
			return NULL;
		}
		else 
			if((Polygon->dimensions[1]!=2)||(Polygon->dimensions[0] < 3)) 
			{
				PyErr_SetString (PyExc_ValueError,
					"Polygon must be M X 2 NumPy array, with M >= 3");
				return NULL;
			}

    // calculate area of the polygon in degrees
	N = Polygon->dimensions[0];
    sumy = 0.0;
    for(i=0 ;i < N-1; i++)	// N-1 - because last point of polygon is the same as first one
    {
        j=(i+1)%(N-1);

        xi =  ARRAYVAL2(double,Polygon,i,0);
        xj =  ARRAYVAL2(double,Polygon,j,0);
        yi =  ARRAYVAL2(double,Polygon,i,1);
        yj =  ARRAYVAL2(double,Polygon,j,1);
        
        area += xi*yj;
        area -= yi*xj;
        sumy += yi;
    }

    avgy = sumy / (N-1);
    //printf("Average y is %f\n",avgy);
    
    area =fabs(area)/2.;
    //printf("The area in degrees is: %f\n",area);

    // convert from square degrees to square meters, using a flat earth
    // projection, based on the mid latitude of the polygon. This is
    // accurate enough only for small polygons

    area = area * mpd * mpd * cos(avgy/180.0*pi);

	return Py_BuildValue("d",fabs(area));
}

static PyMethodDef CalcPolygonsMethods[] = 
{
  {"GetVolumes",  CalcPolygons_GetVolumes,  METH_VARARGS},
  {"GetArea",  CalcPolygons_GetArea,  METH_VARARGS},
  {NULL,      NULL}        /* Sentinel */
};


/*
void initCalcPolygons(void)

{
  (void) Py_InitModule("CalcPolygons", CalcPolygonsMethods);
  import_array()
}
*/



















