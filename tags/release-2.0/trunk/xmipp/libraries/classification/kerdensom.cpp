/***************************************************************************
 *
 * Authors:     Alberto Pascual Montano (pascual@cnb.uam.es)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.uam.es'
 ***************************************************************************/

//-----------------------------------------------------------------------------
// xmippKerDenSOM.cc
// Implements Smoothly Distributed Kernel Probability Density Estimator Self-Organizing Map
// This is an abstract base class for different variants of the KerDenSOM algorithm
//-----------------------------------------------------------------------------

#include <fstream>

#include "kerdensom.h"

/**
 * Sets the number of training steps
 * Parameter: _nSteps  Number of training steps
 */
void xmippKerDenSOM::nSteps(const unsigned long& _nSteps)
{
    somNSteps = _nSteps;
};

/**
 * Gets the sigma (Kernel Width)
 */
double xmippKerDenSOM::getSigma()
{
    return sigma;
};


//-----------------------------------------------------------------------------

/**
 * Sets the number of deterministic annealing training steps
 * Parameter: _annSteps  Number of steps
 */
void xmippKerDenSOM::setAnnSteps(const unsigned long& _annSteps)
{
    annSteps = _annSteps;
};


//-----------------------------------------------------------------------------

/**
 * Tests the KerDenSOM
 * Parameter: _som        The KerDenSom to test
 * Parameter: _examples   The training set of examples
 */
double xmippKerDenSOM::test(const xmippFuzzyMap& _som, const TS& _examples) const
{

    // Defines verbosity level
    int verbosity = listener->getVerbosity();
    if (verbosity)
    {
        listener->OnReportOperation((string) "\nEstimating quantization error....\n");
        listener->OnInitOperation(_examples.size());
    }


    /* Scan all data entries */
    double qerror = 0.0;
    for (int i = 0; i < _examples.size(); i++)
    {
        SomIn& theBest = _som.fuzzyTest(i); // get the best
        qerror += (double) eDist(theBest, _examples.theItems[i]);
        if (verbosity)
        {
            int tmp = (int)((_examples.size() * 5) / 100);
            if ((tmp == 0) && (i != 0)) tmp = i;
            else tmp = 1;
            if ((i % tmp) == 0)
                listener->OnProgress(i);
        }
    }
    if (verbosity)listener->OnProgress(_examples.size());
    return (qerror / (double) _examples.size());

};


//-----------------------------------------------------------------------------


/**
 * Update Code Vectors
 */
void xmippKerDenSOM::updateV(xmippFuzzyMap* _som, const TS* _examples, const double& _reg)
{
    unsigned j, cc, vv;
    unsigned t2 = 0;  // Iteration index

    // Calculate Temporal scratch values
    for (cc = 0; cc < numNeurons; cc++)
    {
        for (j = 0; j < dim; j++) tmpMap[cc][j] = 0.;
        if (_reg != 0)
            tmpDens[cc] = _reg * _som->getLayout().numNeig(_som, (SomPos) _som->indexToPos(cc));
        else
            tmpDens[cc] = 0.;
        for (vv = 0; vv < numVectors; vv++)
        {
            double tmpU = (double) _som->memb[vv][cc];
            tmpDens[cc] += tmpU;
            for (j = 0; j < dim; j++)
            {
                tmpMap[cc][j] += (double)((double) tmpU * (double)(_examples->theItems[vv][j]));
            }
        }
    }


    // Update Code vectors using a sort of Gauss-Seidel iterative algorithm.
    // Usually 100 iterations are enough.


    double convergence = 1, stopError2, stopError1;
    while ((convergence > 1e-5) && (t2 < 100))
    {
        t2++;
        stopError2 = 0;
        stopError1 = 0;
        for (cc = 0; cc < numNeurons; cc++)
        {
            if (_reg != 0)
                _som->localAve(_som->indexToPos(cc), tmpV);
            for (j = 0; j < dim; j++)
            {
                double tmpU = (tmpMap[cc][j] + tmpV[j] * _reg) / tmpDens[cc];
                stopError1 += fabs((double)(_som->theItems[cc][j]) - tmpU);
                stopError2 += fabs(tmpU);
                _som->theItems[cc][j] = (xmippFeature) tmpU;
            }
        } // for
        convergence = stopError1 / stopError2;
    } // while


}


//-----------------------------------------------------------------------------

// Main iterations
double xmippKerDenSOM::mainIterations(xmippFuzzyMap* _som, const TS* _examples, double& _sigma, const double& _reg)
{
    int verbosity = listener->getVerbosity();
    double stopError = 1e10, alpha, ts2;
    int iter = 0;
    if (somNSteps == 0) return stopError;
    if (verbosity == 1 || verbosity == 3)
        listener->OnInitOperation(somNSteps);
    do
    {
        iter++;
        updateU(_som, _examples, _sigma, alpha);
        ts2 = updateSigmaII(_som, _examples, _reg, alpha);
        stopError = fabs(_sigma / ts2 - 1.);
        _sigma = ts2;
        updateV(_som, _examples, _reg);
        if (verbosity == 1 || verbosity == 3)
            listener->OnProgress(iter);
        if (verbosity >= 2)
        {
            char s[100];
            sprintf(s, "Iteration %d of %d. variation: %g\n", iter, somNSteps, stopError);
            listener->OnReportOperation((string) s);
        }
    }
    while ((stopError > epsilon) && (iter < somNSteps));
    if (verbosity == 1 || verbosity == 3)
        listener->OnProgress(somNSteps);
    return stopError;
}


//-----------------------------------------------------------------------------

// Estimate Sigma Part I
double xmippKerDenSOM::updateSigmaI(xmippFuzzyMap* _som, const TS* _examples)
{
    double r, t = 0;
    int vv, cc, j;

    // Computing Sigma (Part I)
    for (vv = 0; vv < numVectors; vv++)
    {
        for (cc = 0; cc < numNeurons; cc++)
        {
            r = 0.0;
            for (j = 0; j < dim; j++)
                r += ((double)(_examples->theItems[vv][j]) - (double)(_som->theItems[cc][j])) * ((double)(_examples->theItems[vv][j]) - (double)(_som->theItems[cc][j]));
            t += r * (double)(_som->memb[vv][cc]);
        }
    }
    return (double)(t / (double)(numVectors*dim));
}


//-----------------------------------------------------------------------------
/**************** Necessary stuff ***************************************/
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------

/**
 * Special Initialization of Code Vectors
 */
void xmippKerDenSOM::updateV1(xmippFuzzyMap* _som, const TS* _examples)
{
    unsigned j, cc, vv;

    for (cc = 0; cc < numNeurons; cc++)
    {
        for (j = 0; j < dim; j++) tmpMap[cc][j] = 0.;
        tmpDens[cc] = 0.0;
        for (vv = 0; vv < numVectors; vv++)
        {
            double tmpU = (double)(_som->memb[vv][cc]);
            tmpDens[cc] += tmpU;
            for (j = 0; j < dim; j++)
            {
                tmpMap[cc][j] += (double)((double) tmpU * (double)(_examples->theItems[vv][j]));
            }
        }
    }

    for (cc = 0; cc < numNeurons; cc++)
    {
        for (j = 0; j < dim; j++)
        {
            double tmpU = tmpMap[cc][j] / tmpDens[cc];
            _som->theItems[cc][j] = (xmippFeature) tmpU;
        }
    } // for


}

//-----------------------------------------------------------------------------
/**
 * Special Initialization of Membership Matrix (Fuzzy c-means style)
 */
void xmippKerDenSOM::updateU1(xmippFuzzyMap* _som, const TS* _examples)
{

    double auxProd, auxDist, tmp;
    int k, j, i;

    // Update Membership matrix
    for (k = 0; k < numVectors; k++)
    {
        auxProd = 1;
        for (j = 0; j < numNeurons; j++)
            auxProd *= (double) eDist(_som->theItems[j], _examples->theItems[k]);

        if (auxProd == 0.)
        { // Apply k-means criterion (Data-CB) must be > 0
            for (j = 0; j < numNeurons; j ++)
                if (eDist(_som->theItems[j], _examples->theItems[k]) == 0.) _som->memb[k][j] = 1.0;
                else _som->memb[k][j] =  0.0;
        }
        else
        {
            for (i = 0; i < numNeurons; i ++)
            {
                auxDist = 0;
                for (j = 0; j < numNeurons; j ++)
                {
                    tmp = eDist(_som->theItems[i], _examples->theItems[k]) /
                          eDist(_som->theItems[j], _examples->theItems[k]);
                    auxDist += pow(tmp, 2);
                } // for j
                _som->memb[k][i] = (xmippFeature) 1.0 / auxDist;
            } // for i
        } // if auxProd
    } // for k
}

//-----------------------------------------------------------------------------

/**
 * Special Initialization of the Us
 */
void xmippKerDenSOM::initU(xmippFuzzyMap* _som)
{
    unsigned j, cc, vv;

    // Take random samples
    randomize_random_generator();
    for (vv = 0; vv < numVectors; vv++)
    {
        double t = 0.;
        for (cc = 0; cc < numNeurons; cc++)
        {
            _som->memb[vv][cc] = (xmippFeature) rnd_unif();
            t += _som->memb[vv][cc];
        }
        for (cc = 0; cc < numNeurons; cc++)
            _som->memb[vv][cc] /= (xmippFeature) t;
    }
}


//-----------------------------------------------------------------------------


/**
 * Determines the Random Approximation of GVC
 * (Determines the optimal Regularization Factor)
 */
double xmippKerDenSOM::randApproxGVC(const TS* _examples, const xmippFuzzyMap* _som, double _dataSD, double _reg)
{
    unsigned j, vv, cc;
    double num, den, r;
    xmippVector VV;
    VV.resize(dim, 0.0);
    xmippFuzzyMap tmpSOM(*_som);
    TS tmpTS(*_examples);
    num = 0;

    for (vv = 0; vv < numVectors; vv++)
    {
        for (j = 0; j < dim; j++) VV[j] = 0.0;
        for (cc = 0; cc < numNeurons; cc++)
        {
            for (j = 0; j < dim; j++)
                VV[j] += (_examples->theItems[vv][j] - _som->theItems[cc][j]) * _som->memb[vv][cc];
        }
        for (j = 0; j < dim; j++)
            num += VV[j] * VV[j];
    }

    init_random_generator();
    for (vv = 0; vv < numVectors; vv++)
    {
        for (j = 0; j < dim; j++)
            tmpTS.theItems[vv][j] += rnd_gaus() * _dataSD;
    }
    updateV(&tmpSOM, &tmpTS, _reg);
    den = 0.0;

    init_random_generator();
    for (vv = 0; vv < numVectors; vv++)
    {
        for (j = 0; j < dim; j++) VV[j] = 0.0;
        for (cc = 0; cc < numNeurons; cc++)
        {
            for (j = 0; j < dim; j++)
                VV[j] += (tmpTS.theItems[vv][j] - tmpSOM.theItems[cc][j] - _examples->theItems[vv][j] + _som->theItems[cc][j]) * _som->memb[vv][cc];
        }
        r = 0.;
        for (j = 0; j < dim; j++)
            r += VV[j] * rnd_gaus() * _dataSD;
        den += r * r;
    }
    if (den != 0)
        return (double) num / den;
    else
        return 0.;
}

//-----------------------------------------------------------------------------
/**************** Visualization methods ***************************************/
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------

/**
 * Shows the data
 */

void xmippKerDenSOM::showX(const TS* _ts)
{
    cout << "Data (1..nd, 1..nv) "  << endl;
    for (int i = 0; i < _ts->size(); i++)
    {
        for (int j = 0; j < _ts->theItems[0].size(); j++)
        {
            cout << i + 1 << "  " << j + 1 << "  " << _ts->theItems[i][j] << endl;
        }
    }
}

//-----------------------------------------------------------------------------

/**
 * Shows the code vectors
 */

void xmippKerDenSOM::showV(xmippFuzzyMap* _som)
{
    cout << "Code vectors (1..ni, 1..nj, 1..nv) "  << endl;
    for (int i = 0; i < _som->size(); i++)
    {
        int tmpj = _som->indexToPos(i).first;
        int tmpi = _som->indexToPos(i).second;
        for (int j = 0; j < _som->theItems[0].size(); j++)
            cout << tmpi + 1 << "  " << tmpj + 1 << "  " << j + 1 << "  " << _som->theItems[i][j] << endl;
    }
}

//-----------------------------------------------------------------------------
/**
 * Shows the U ("Membership")
 */

void xmippKerDenSOM::showU(xmippFuzzyMap* _som, const TS* _ts)
{
    cout << " Memberships (1..nd,1..ni,1..nj)" << endl;
    for (int i = 0; i <  _ts->size(); i++)
    {
        for (int j = 0; j < _som->size(); j++)
        {
            int tmpj = _som->indexToPos(j).first;
            int tmpi = _som->indexToPos(j).second;
            cout << i + 1 << "  " << tmpi + 1 << "  " << tmpj + 1 << "  " << _som->memb[i][j] << endl;
        }
    }
}

//-----------------------------------------------------------------------------
/**
 * Prints the Us ("Membership")
 */

void xmippKerDenSOM::printV(xmippFuzzyMap* _som, const TS* _ts, FileName& _fname)
{

    FILE* F = fopen(_fname.c_str(), "w");
    if (F == NULL)
    {
        return;
    }

    fprintf(F, "%d %s %d %d gaussian\n", dim, _som->layout().c_str(), _som->width(), _som->height());
    if (_ts->isNormalized())
    {
        if (_ts->getNormalizationInfo().size() != _som->theItems[0].size())
        {
            ostrstream msg;
            msg << "Normalization information does not coincide with codebook structure";
            throw runtime_error(msg.str());
        }
        for (unsigned it = 0; it < _som->size(); it++)
        {
            for (unsigned i = 0; i < _som->theItems[0].size(); i++)
            {
                if (!isnan(_som->theItems[it][i]))
                    _som->theItems[it][i] = _som->theItems[it][i] * _ts->getNormalizationInfo()[i].sd + _ts->getNormalizationInfo()[i].mean;
                fprintf(F, "%g ", _som->theItems[it][i]);
            }
            fprintf(F, "\n");
        }
    }
    else
    {
        for (int i = 0; i < _som->size(); i++)
        {
            for (int j = 0; j < _som->theItems[0].size(); j++)
                fprintf(F, "%g ", _som->theItems[i][j]);
            fprintf(F, "\n");
        }
    }

    fclose(F);

}

//-----------------------------------------------------------------------------
