/*=================================================================
 *
 * helper functions for the TOM xmipp wrappers
 *
 *
 * Electron Tomography toolbox of the
 * Max-Planck-Institute for Biochemistry
 * Dept. Molecular Structural Biology
 * 82152 Martinsried, Germany
 * http://www.biochem.mpg.de
 *
 * and
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 *
 * created: 27/08/2007
 * by: Andreas Korinek & Carlos Oscar Sorzano
 *
 *=================================================================*/

#include "mex.h"
#include "matrix.h"
#include "matrix2d.h"
#include "matrix3d.h"

template <class T>
   void getMatrix1D(const mxArray* prhs, Matrix1D<T> &output)
{
    const int  *input_dims = mxGetDimensions(prhs);
    T *indatapt=(T*)mxGetData(prhs);
    output.resize(input_dims[0]);
    FOR_ALL_ELEMENTS_IN_MATRIX1D(output)
       DIRECT_VEC_ELEM(output,i)=*indatapt++;
}


template <class T>
   void getMatrix2D(const mxArray* prhs, Matrix2D<T> &output)
{
    const int  *input_dims = mxGetDimensions(prhs);
    T *indatapt=(T*)mxGetData(prhs);
    output.resize(input_dims[0],input_dims[1]);
    FOR_ALL_ELEMENTS_IN_MATRIX2D(output)
       DIRECT_MAT_ELEM(output,i,j)=*indatapt++;
}

template <class T>
   void getMatrix3D(const mxArray* prhs, Matrix3D<T> &output)
{
    const int  *input_dims = mxGetDimensions(prhs);
    T *indatapt=(T*)mxGetData(prhs);
    output.resize(input_dims[2],input_dims[0],input_dims[1]);
    FOR_ALL_ELEMENTS_IN_MATRIX3D(output)
       DIRECT_VOL_ELEM(output,k,i,j)=*indatapt++;
}


template <class T>
   void setMatrix1D(const Matrix1D<T> &input, mxArray* &plhs)
{
    int output_dims[2];
    output_dims[0] = XSIZE(input);
    output_dims[1] = 1;
    mxClassID typeClass;
    if (typeid(T)==typeid(double))
        typeClass=mxDOUBLE_CLASS;
    else if (typeid(T)==typeid(int))
        typeClass=mxINT32_CLASS;

    if ((plhs = mxCreateNumericArray(1,output_dims,typeClass,mxREAL))==NULL) 
        mexErrMsgTxt("Memory allocation problem in Xmipp converter.\n");
    T *outdata =(T*) mxGetData(plhs);
    FOR_ALL_ELEMENTS_IN_MATRIX1D(input)
        *outdata++=DIRECT_VEC_ELEM(input,i);
}

template <class T>
   void setMatrix2D(const Matrix2D<T> &input, mxArray* &plhs)
{
    int output_dims[2];
    output_dims[0] = XSIZE(input);
    output_dims[1] = YSIZE(input);
    
    mxClassID typeClass;
    if (typeid(T)==typeid(double))
        typeClass=mxDOUBLE_CLASS;
    else if (typeid(T)==typeid(int))
        typeClass=mxINT32_CLASS;
    
    if ((plhs = mxCreateNumericArray(2,output_dims,typeClass,mxREAL))==NULL) 
        mexErrMsgTxt("Memory allocation problem in Xmipp converter.\n");
    T *outdata =(T*) mxGetData(plhs);
    FOR_ALL_ELEMENTS_IN_MATRIX2D(input)
       *outdata++=DIRECT_MAT_ELEM(input,i,j);
}


template <class T>
   void setMatrix3D(const Matrix3D<T> &input, mxArray* &plhs)
{
    int output_dims[3];
    output_dims[0] = XSIZE(input);
    output_dims[1] = YSIZE(input);
    output_dims[2] = ZSIZE(input);
    
    mxClassID typeClass;
    if (typeid(T)==typeid(double))
        typeClass=mxDOUBLE_CLASS;
    else if (typeid(T)==typeid(int))
        typeClass=mxINT32_CLASS;
    
    if ((plhs = mxCreateNumericArray(3,output_dims,typeClass,mxREAL))==NULL) 
        mexErrMsgTxt("Memory allocation problem in Xmipp converter.\n");
    T *outdata =(T*) mxGetData(plhs);
    FOR_ALL_ELEMENTS_IN_MATRIX3D(input)
       *outdata++=DIRECT_VOL_ELEM(input,k,i,j);
}

#define NUMBER_OF_FIELDS (sizeof(field_names)/sizeof(*field_names))
