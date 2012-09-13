/*=================================================================
 *
 * tom_xmipp_psd_enhance is a wrapper to xmipp_psd_enhance
 *
 * The calling syntax is:
 *
 *		psd = tom_xmipp_psd_enhance_wrapper(image,center,take_log,filter_w1,filter_w2,decay_width,mask_w1,mask_w2)
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

/*xmipp includes */
#include "psd_enhance.h"
#include "tom_xmipp_helpers.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[])
{

    /*read image*/
    Matrix2D<double> image;
    getMatrix2D(prhs[0],image);
    
    Prog_Enhance_PSD_Parameters psdParams;
    
    /*read center flag*/
    psdParams.center = (bool) mxGetScalar(prhs[1]);
    
    /*read log flag*/
    psdParams.take_log = (bool) mxGetScalar(prhs[2]);
    
    /*read filter_w1 (Low freq. for band pass filtration,  max 0.5)*/
    psdParams.filter_w1 = (float) mxGetScalar(prhs[3]);
    
    /*read filter_w2 (High freq. for band pass filtration, max 0.5)*/
    psdParams.filter_w2 = (float) mxGetScalar(prhs[4]);
    
    /*read decay_witdh (Decay for the transition bands)*/
    psdParams.decay_width = (float) mxGetScalar(prhs[5]);

    /*mask_w1 (Low freq. for mask, max 0.5)*/
    psdParams.mask_w1 = (float) mxGetScalar(prhs[6]);
    
    /*mask_w2 (High freq. for mask, max 0.5)*/
    psdParams.mask_w2 = (float) mxGetScalar(prhs[7]);
    
    try
    {
        psdParams.apply(image);
    }
    catch (Xmipp_error Xe)
    {
        mexErrMsgTxt(Xe.msg.c_str());
    }
    
    setMatrix2D(image, plhs[0]);
       
}	

