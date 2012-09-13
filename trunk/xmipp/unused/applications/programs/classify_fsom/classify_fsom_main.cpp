/***************************************************************************
 *
 * Authors:     Alberto Pascual Montano (pascual@cnb.csic.es)
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
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/


/* Part of this code were developed by Lorenzo Zampighi and Nelson Tang
   of the department of Physiology of the David Geffen School of Medicine,
   University of California, Los Angeles
*/

// To avoid problems with long template names
#pragma warning(disable:4786)

#include <fstream>

#include <classification/fuzzy_som.h>

/* Prototypes -============================================================= */

void Usage(char **argv);

/* Main function -============================================================= */

main(int argc, char** argv)
{


    /* Input Parameters ======================================================== */

    FileName       fn_in;     // input file
    FileName       fn_out;    // output file
    FileName       cb_in = "";    // Code vectors input file
    FileName       fn_algo_in = ""; // input algorithm file
    FileName       tmpN;  // Temporary variable
    double         eps = 1e-7; // Stopping criteria
    unsigned       iter = 1000; // Iteration number
    unsigned       verb = 0; // Verbosity level
    bool           norm = 1; // Normalize?
    unsigned       xdim;  // X-dimension (-->)
    unsigned       ydim;  // Y-dimension
    double         m0 = 2.0; // Initial m
    double         m1 = 1.01; // Final m
    double         reg;  // Regularization (smoothness) parameter
    std::string    layout = "RECT"; // topology (layout)
    unsigned       annSteps = 1000; // Deterministic Annealing steps
    bool           saveClusters = false;    // Save clusters in separate files
    bool           use_rand_cvs = false; // NT: flag to truly randomize codevectors or not

    /* Parameters ============================================================== */
    try
    {

        if (checkParameter(argc, argv, "-i"))
            fn_in = getParameter(argc, argv, "-i");
        else
        {
            Usage(argv);
            exit(EXIT_FAILURE);
        }

        if (checkParameter(argc, argv, "-o"))
            fn_out = getParameter(argc, argv, "-o");
        else
        {
            Usage(argv);
            exit(EXIT_FAILURE);
        }

        if (checkParameter(argc, argv, "-cvin"))
            cb_in = getParameter(argc, argv, "-cvin");


        if (checkParameter(argc, argv, "-xdim"))
            xdim = textToInteger(getParameter(argc, argv, "-xdim"));
        else
        {
            Usage(argv);
            exit(EXIT_FAILURE);
        }

        if (checkParameter(argc, argv, "-ydim"))
            ydim = textToInteger(getParameter(argc, argv, "-ydim"));
        else
        {
            Usage(argv);
            exit(EXIT_FAILURE);
        }

        if (checkParameter(argc, argv, "-hexa"))
        {
            if (checkParameter(argc, argv, "-rect"))
            {
                std::cout << "Error: you can not define two topologies" << std::endl;
                exit(EXIT_FAILURE);
            }
            layout = "HEXA";
        }
        else if (checkParameter(argc, argv, "-rect"))
            layout = "RECT";

        m0 =  textToFloat(getParameter(argc, argv, "-m0", "2.0"));
        m1 =  textToFloat(getParameter(argc, argv, "-m1", "1.01"));
        reg =  textToFloat(getParameter(argc, argv, "-reg", "0.5"));

        eps = textToFloat(getParameter(argc, argv, "-eps", "1e-7"));
        iter = textToInteger(getParameter(argc, argv, "-iter", "1000"));
        verb = textToInteger(getParameter(argc, argv, "-verb", "0"));

        if (checkParameter(argc, argv, "-norm"))
            norm = true;
        else norm = false;

        annSteps = textToInteger(getParameter(argc, argv, "-steps", "1000"));

        if (checkParameter(argc, argv, "-saveclusters"))
            saveClusters = true;
        else saveClusters = false;

        if (checkParameter(argc, argv, "-randomcodevectors"))
            use_rand_cvs = true;
        else use_rand_cvs = false;

        if (argc == 1)
        {
            Usage(argv);
        }

    }
    catch (XmippError XE)
    {
        std::cout << XE;
        Usage(argv);
    }


    /* Some validations ===================================================== */


    if (iter < 1)
    {
        std::cerr << argv[0] << ": invalid value for iter (must be > 1): " << iter << std::endl;
        exit(EXIT_FAILURE);
    }

    if (verb < 0 || verb > 2)
    {
        std::cerr << argv[0] << ": invalid value for verbosity (must be between 0 and 2): " << verb << std::endl;
        exit(EXIT_FAILURE);
    }

    if (m0 <= 1)
    {
        std::cerr << argv[0] << ": invalid value for m0 (must be > 1): " << m0 << std::endl;
        exit(EXIT_FAILURE);
    }

    if (m1 <= 1)
    {
        std::cerr << argv[0] << ": invalid value for m1 (must be > 1): " << m1 << std::endl;
        exit(EXIT_FAILURE);
    }

    if ((annSteps != 0) && (m0 <= m1))
    {
        std::cerr << argv[0] << ": invalid value for m0 and m1 (m0 must be > m1) " << std::endl;
        exit(EXIT_FAILURE);
    }

    if ((annSteps < 0) || (annSteps == 1))
    {
        std::cerr << argv[0] << ": invalid value for annSteps (must be > 1): " << annSteps << std::endl;
        exit(EXIT_FAILURE);
    }

    if (reg < 0)
    {
        std::cerr << argv[0] << ": invalid value for smoothness parameter (must be > 0): " << reg << std::endl;
        exit(EXIT_FAILURE);
    }

    if (xdim < 1)
    {
        std::cerr << argv[0] << ": invalid value for xdim (must be > 1): " << xdim << std::endl;
        exit(EXIT_FAILURE);
    }

    if (ydim < 1)
    {
        std::cerr << argv[0] << ": invalid value for ydim (must be > 1): " << ydim << std::endl;
        exit(EXIT_FAILURE);
    }

    if (reg == 0) reg = 1.23456789e-7;


    /* Shows parameters ===================================================== */

    std::cout << std::endl << "Parameters used: " << std::endl;
    std::cout << "Input data file : " << fn_in << std::endl;
    std::cout << "Output code vector file : " << fn_out << ".cod" << std::endl;
    if (cb_in != "")
        std::cout << "Input code vectors file name : " << cb_in << std::endl;
    else if (use_rand_cvs)
        std::cout << "Using randomized code vectors" << std::endl;
    std::cout << "Horizontal dimension (Xdim) = " << xdim << std::endl;
    std::cout << "Vertical dimension (Ydim) = " << ydim << std::endl;
    if (layout == "HEXA")
        std::cout << "Hexagonal topology " << std::endl;
    else
        std::cout << "Rectangular topology " << std::endl;
    std::cout << "Initial Fuzzy Constant (m0) = " << m0 << std::endl;
    std::cout << "Final Fuzzy Constant (m1) = " << m1 << std::endl;
    std::cout << "Regularization Constant (reg) = " << reg << std::endl;
    std::cout << "Deterministic annealing steps = " << annSteps << std::endl;
    std::cout << "Total number of iterations = " << iter << std::endl;
    std::cout << "Stopping criteria (eps) = " << eps << std::endl;
    std::cout << "verbosity level = " << verb << std::endl;
    if (norm)
        std::cout << "Normalize input data" << std::endl;
    else
        std::cout << "Do not normalize input data " << std::endl;


    /* Open training vector ================================================= */

    std::cout << std::endl << "Reading file " << fn_in << "....." << std::endl;

    std::ifstream inStream(fn_in.c_str());
    if (!inStream)
    {
        std::cerr << argv[0] << ": can't open file " << fn_in << std::endl;
        exit(EXIT_FAILURE);
    }

    ClassicTrainingVectors ts(0, true);
    try
    {
        inStream >> ts;
    }
    catch (std::exception& e)
    {
        std::cerr << argv[0] << ": can't read file " << fn_in  << " because " << e.what() << std::endl;
        exit(EXIT_FAILURE);
    }


    /* Real stuff ============================================================== */


    try
    {

        if (norm)
        {
            std::cout << "Normalizing....." << std::endl;
            ts.normalize();        // Normalize input data
        }

        FuzzyMap *myMap;

        if (cb_in != "")
        {
            std::cout << "Reading fuzzy codevectors file " << cb_in << "....." << std::endl;
            std::ifstream codeStream(cb_in.c_str());
            if (!codeStream)
            {
                std::cerr << argv[0] << ": can't open file " << cb_in << std::endl;
                exit(EXIT_FAILURE);
            }
            myMap = new FuzzyMap(codeStream, ts.size(), true);
        }
        else
            myMap = new FuzzyMap(layout, xdim, ydim, ts, use_rand_cvs);


        FuzzySOM *thisSOM;
        if (fn_algo_in == "")
        {
            thisSOM = new FuzzySOM(m0, m1, annSteps, reg, eps, iter);    // Creates FSOM Algorithm
        }
        else
        {
            std::cout << "Reading algorithm file " << fn_algo_in << "....." << std::endl << std::endl;
            std::ifstream algoStream(fn_algo_in.c_str());
            if (!algoStream)
            {
                std::cerr << argv[0] << ": can't open file " << fn_algo_in << std::endl;
                exit(EXIT_FAILURE);
            }
        }

        TextualListener myListener;     // Define the listener class
        myListener.setVerbosity() = verb;     // Set verbosity level
        thisSOM->setListener(&myListener);       // Set Listener

        if (cb_in != "")
        {
            if (ts.isNormalized())
            {
                std::cout << "Normalizing code vectors....." << std::endl;
                myMap->Normalize(ts.getNormalizationInfo());       // normalize code vectors
            }
            thisSOM->train(*myMap, ts);                   // Train algorithm
        }
        else
            thisSOM->train(*myMap, ts);                     // Train algorithm

        // Test algorithm
        double dist = thisSOM->test(*myMap, ts);
        std::cout << std::endl << "Quantization error : " <<  dist << std::endl;

        // Calculates functional value
        double functional, fidelity, penalty;
        functional = thisSOM->functional(ts, *myMap, m1, reg, fidelity, penalty);
        std::cout << "Functional : " <<  functional << " (fidelity = " << fidelity << " penalty = " << penalty << " )" << std::endl << std::endl;


        // Classifying
        std::cout << "Classifying....." << std::endl;
        myMap->classify(&ts);

        // Calibrating
        std::cout << "Calibrating....." << std::endl;
        myMap->calibrate(ts);

        /*******************************************************
            Saving all kind of Information
        *******************************************************/

        std::cout << "Saving algorithm information as " << fn_out << ".inf ....." << std::endl;
        tmpN = fn_out.c_str() + (std::string) ".inf";
        std::ofstream infS(tmpN.c_str());
        infS << "Fuzzy SOM algorithm" << std::endl << std::endl;
        infS << "Input data file : " << fn_in << std::endl;
        if (cb_in != "")
            infS << "Input code vectors file : " << cb_in << std::endl;
        else if (use_rand_cvs)
            infS << "Using randomized code vectors" << std::endl;
        infS << "Code vectors output file : " << fn_out <<  ".cod" << std::endl;
        infS << "Algorithm information output file : " << fn_out <<  ".inf" << std::endl;
        infS << "Number of feature vectors: " << ts.size() << std::endl;
        infS << "Number of variables: " << ts.theItems[0].size() << std::endl;
        infS << "Horizontal dimension (Xdim) = " << xdim << std::endl;
        infS << "Vertical dimension (Ydim) = " << ydim << std::endl;
        if (layout == "HEXA")
            infS << "Hexagonal topology " << std::endl;
        else
            infS << "Rectangular topology " << std::endl;
        if (norm)
            infS << "Input data normalized" << std::endl;
        else
            infS << "Input data not normalized" << std::endl;
        infS << "Initial Fuzzy constant (m0) = " << m0 << std::endl;
        infS << "Final Fuzzy constant (m1) = " << m1 << std::endl;
        infS << "Smoothness factor (reg) = " << reg << std::endl;
        infS << "Deterministic annealing steps = " << annSteps << std::endl;
        infS << "Total number of iterations = " << iter << std::endl;
        infS << "Stopping criteria (eps) = " << eps << std::endl;
        infS << "Quantization error : " <<  dist << std::endl;
        infS << "Functional : " <<  functional << " (fidelity = " << fidelity << " penalty = " << penalty << " )" << std::endl << std::endl;
        infS.flush();

        // assign data to clusters according to fuzzy threshold
        if (saveClusters)
        {
            std::cout << "Saving neurons assigments ....." << std::endl;
            for (unsigned i = 0; i < myMap->size(); i++)
            {
                tmpN = fn_out.c_str() + (std::string) "."  + integerToString(i);
                std::ofstream cStream(tmpN.c_str());
                for (int j = 0; j < myMap->classifAt(i).size(); j++)
                    cStream << myMap->classifAt(i)[j] << std::endl;
                cStream.flush();
            }
        }

        // save .vs file to be compatible with SOM_PAK
        std::cout << "Saving visual file as " << fn_out << ".vs ....." << std::endl;
        tmpN = fn_out.c_str() + (std::string) ".vs";
        std::ofstream vsStream(tmpN.c_str());
        vsStream << ts.theItems[0].size() << " " << myMap->layout() << " " << myMap->width() << " " << myMap->height() << " gaussian" << std::endl;
        for (int i = 0; i < ts.size(); i++)
        {
            int j = myMap->fuzzyWinner(i);
            vsStream << myMap->indexToPos(j).first << " " << myMap->indexToPos(j).second << " " << myMap->memb[i][j] << " " << ts.theTargets[i] << std::endl;
        }
        vsStream.flush();


        // save .his file (Histogram)
        std::cout << "Saving code vectors histogram file as " << fn_out << ".his ....." << std::endl;
        tmpN = fn_out.c_str() + (std::string) ".his";
        std::ofstream hisStream(tmpN.c_str());
        myMap->printHistogram(hisStream);
        hisStream.flush();

        // save .err file (Average Quantization Error)
        std::cout << "Saving code vectors average quantization error file as " << fn_out << ".err ....." << std::endl;
        tmpN = fn_out.c_str() + (std::string) ".err";
        std::ofstream errStream(tmpN.c_str());
        myMap->printQuantError(errStream);
        errStream.flush();

        if (norm)
        {
            std::cout << "Denormalizing code vectors....." << std::endl;
            myMap->unNormalize(ts.getNormalizationInfo()); // de-normalize codevectors
        }

        std::cout << "Saving code vectors as " << fn_out << ".cod ....." << std::endl;
        tmpN = fn_out.c_str() + (std::string) ".cod";
        std::ofstream codS(tmpN.c_str());
        codS << *myMap;
        codS.flush();


        std::cout << std::endl;

        delete myMap;
        delete thisSOM;


    }
    catch (const std::exception& e)
    {
        std::cout << e.what() << std::endl;
    }
    return 0;
}


/* ------------------------------------------------------------------------- */
/* Help Message for this Program                                             */
/* ------------------------------------------------------------------------- */
void Usage(char **argv)
{
    printf(
        "\nUsage: %s [Purpose and Parameters]"
        "\nPurpose: Fuzzy Self-Organizing Feature Map"
        "\n"
        "\nParameter Values: (note space before value)"
        "\n"
        "\n    -i      file_in           Input data file"
        "\n    -o      file_out          Base name for output data files "
        "\n    -cvin   file_in           Codevectors input file"
        "\n    -saveclusters           save clusters in separate files (Default = No)"
        "\n    -xdim   H-dimension       Horizontal size of the map"
        "\n    -ydim   V-dimension       Vertical size of the map"
        "\n    -hexa               Hexagonal topology"
        "\n    -rect               Rectangular topology (default)"
        "\n    -steps  steps           Deterministic annealing steps (default = 1000)"
        "\n    -m0     Initial m         Initial Fuzzy constant (default = 2)"
        "\n    -m1     Final m        Final Fuzzy constant (default = 1.02)"
        "\n    -reg    smoothness        Smoothness factor (default = 0.5)"
        "\n    -eps    Epsilon        Stopping criteria (default = 1e-7)"
        "\n    -iter   iterations        Number of iterations (default = 1000)"
        "\n    -norm                   Normalize training data (default)"
        "\n    -verb   verbosity         Information level while running: "
        "\n             0: No information (default)"
        "\n             1: Progress bar"
        "\n             2: Code vectors change between iterations"
        "\n    -randomcodevectors        Use truly randomized codevectors (default: FALSE)"
        "\n      \n"
        , argv[0]);
}
