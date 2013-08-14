
// First include Python and Numeric headers
#include <Python.h>
#include "Numeric/arrayobject.h"

// Now some handy macros I use to access array values for specific rank arrays
#define ARRAYVAL0(aType,a) ( *(aType *)(a->data))	
#define ARRAYVAL1(aType,a,i) ( *(aType *)(a->data + (i)*a->strides[0]))	
#define ARRAYVAL2(aType,a,i,j) ( *(aType *)(a->data + (i)*a->strides[0] + (j)*a->strides[1]))	
#define ARRAYVAL3(aType,a,i,j,k) ( *(aType *)(a->data + (i)*a->strides[0] + (j)*a->strides[1] + (k)*a->strides[2]))	

/* Now include all the modules code */

#include "check_receptors.c"
#include "CalcPolygons.c"


static PyMethodDef noMethods[] = {
  {NULL,      NULL}        /* Sentinel */
};


void initcheck_receptors(void)
{
    (void) Py_InitModule("TAP_ext.check_receptors", check_receptorsMethods);
    /* add module attributes, if any */
    import_array()
}

void initCalcPolygons(void)
{
    (void) Py_InitModule("TAP_ext.CalcPolygons", CalcPolygonsMethods);
    /* add module attributes, if any */
    import_array()
}


void initTAP_ext(void)
{
    PyObject* module;
    PyObject* package = Py_InitModule("TAP_ext", noMethods);
    if(!package) return;

    /* add package attributes, if any */
    /*PyModule_AddStringConstant(package, "foo", "bar");*/

    module = PyImport_AddModule("TAP_ext.check_receptors");
    if(!module) return;
    if(PyModule_AddObject(package, "check_receptors", module))
        return;
    Py_INCREF(module);
    initcheck_receptors();

    module = PyImport_AddModule("TAP_ext.CalcPolygons");
    if(!module) return;
    if(PyModule_AddObject(package, "CalcPolygons", module))
        return;
    Py_INCREF(module);
    initCalcPolygons();

}







