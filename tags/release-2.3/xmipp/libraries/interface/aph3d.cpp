/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
 *              Roberto Marabini (roberto@mipg.upenn.edu)
 *
/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
 *
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

#include "aph3d.h"

#include <data/args.h>

#include <fstream>

#define VERBOSE
//#define DEBUG
// APH =====================================================================
void APHFile3D::read_from_prepmklcf(const FileName &fn)
{
    std::ifstream  fh_aph;
    int            line_no = 1;
    int            hmax = 0, kmax = 0, lmax = 0, hmin = 0, kmin = 0, lmin = 0;
    std::string    line;
    int            h, k, l;

    // Empties current APH File
    clear();
    // Open file
    fh_aph.open(fn.c_str(), std::ios::in);
    if (!fh_aph)
        REPORT_ERROR(1601, "aphFile::read: File " + fn + " not found");

    // Read first line and skip it
//   fh_aph.peek();
    getline(fh_aph, line);
    try
    {
        getline(fh_aph, line);
        Reduce_phase_accuracy = textToFloat(line.c_str() + 47);
        getline(fh_aph, line);
        a = textToFloat(line.c_str() + 47);
        getline(fh_aph, line);
        b = textToFloat(line.c_str() + 47);
        getline(fh_aph, line);
        gamma = textToFloat(line.c_str() + 47);
        getline(fh_aph, line);
        c = textToFloat(line.c_str() + 47);
        getline(fh_aph, line);
        resolution = textToFloat(line.c_str() + 47);
        getline(fh_aph, line);
        Amp_Scale = textToFloat(line.c_str() + 47);
    }
    catch (...)
    {
        REPORT_ERROR(1601, "aph3DFile::read: Wrong first line in file " + fn);
    }

//#define DEBUG_read
#ifdef DEBUG_read
    std::cout << "Reduce_phase_accuracy: " << Reduce_phase_accuracy << std::endl;
    std::cout << "a: " << a << std::endl;
    std::cout << "b: " << b << std::endl;
    std::cout << "gamma: " << gamma << std::endl;
    std::cout << "c: " << c << std::endl;
    std::cout << "resolution: " << resolution << std::endl;
    std::cout << "Amp_Scale: " << Amp_Scale << std::endl;
#endif
#undef DEBUG_read

    // look for the begining of the data
    while (!fh_aph.eof())
    {
        try
        {
            getline(fh_aph, line);
            if (std::string::npos !=
                line.find("   H   K   L      A      P     FOM*100         REJECTS"))
                break;
        }
        catch (Xmipp_error)
        {
            std::cout << "3Daph File reading error an error\n";
        }
    }/* while */

    // look for maximun and minimum
    line_no = 1;

    while (!fh_aph.eof())
    {
        try
        {
            getline(fh_aph, line);
            if (std::string::npos != line.find("                                 "))
                continue;
            else if (std::string::npos != line.find("MKLCF FILE COMPLETED"))
                break;
            h = textToInteger(firstToken(line));
            k = textToInteger(nextToken());
            l = textToInteger(nextToken());
            hmax = XMIPP_MAX(hmax, h);
            kmax = XMIPP_MAX(kmax, k);
            lmax = XMIPP_MAX(lmax, l);
            hmin = XMIPP_MIN(hmin, h);
            kmin = XMIPP_MIN(kmin, k);
            lmin = XMIPP_MIN(lmin, l);
        }
        catch (Xmipp_error)
        {
            std::cout << "aph File: Line " << line_no << " is skipped due to an error\n";
        }
        line_no++;
    }/* while */
//Space Group
    Space_Group = 0;
    while (!fh_aph.eof())
    {
        try
        {
            getline(fh_aph, line);
            if (std::string::npos !=
                line.find(" * Space Group ="))
            {
                Space_Group = textToInteger(line.c_str() + 16);
                break;
            }
        }
        catch (Xmipp_error)
        {
            std::cout << "3Daph File reading error an error\n";
        }
    }/* while */

    if (Space_Group == 0)
        REPORT_ERROR(1601, "aphFile::read: File " + fn + " Space Group not found");

#define DEBUG_max
#ifdef DEBUG_max
    std::cout << "hmax: " << hmax << " kmax: " << kmax << " lmax: " << lmax << std::endl;
    std::cout << "hmin: " << hmin << " kmin: " << kmin << " lmin: " << lmin << std::endl;
#endif
#undef DEBUG_max

    // Ask for memory
    spots_abs.initZeros(lmax - lmin + 1, kmax - kmin + 1, hmax - hmin + 1);
    STARTINGX(spots_abs) = hmin;
    STARTINGY(spots_abs) = kmin;
    STARTINGZ(spots_abs) = lmin;
    spots_arg.initZeros(lmax - lmin + 1, kmax - kmin + 1, hmax - hmin + 1);
    STARTINGX(spots_arg) = hmin;
    STARTINGY(spots_arg) = kmin;
    STARTINGZ(spots_arg) = lmin;
    FOM.initZeros(lmax - lmin + 1, kmax - kmin + 1, hmax - hmin + 1);
    STARTINGX(FOM) = hmin;
    STARTINGY(FOM) = kmin;
    STARTINGZ(FOM) = lmin;

    // Read each line (again) and copy values to the matrices
    fh_aph.close();
    fh_aph.open(fn.c_str(), std::ios::in);
    line_no = 1;

    // look for the begining of the data
    while (!fh_aph.eof())
    {
        try
        {
            getline(fh_aph, line);
            if (std::string::npos !=
                line.find("   H   K   L      A      P     FOM*100         REJECTS"))
                break;
        }
        catch (Xmipp_error)
        {
            std::cout << "3Daph File reading error an error\n";
        }
    }/* while */

    while (!fh_aph.eof())
    {
        try
        {
            getline(fh_aph, line);
            if (std::string::npos != line.find("                                 "))
                continue;
            else if (std::string::npos != line.find("MKLCF FILE COMPLETED"))
                break;
            h     = textToInteger(firstToken(line));
            k     = textToInteger(nextToken());
            l     = textToInteger(nextToken());
            spots_abs(l, k, h)  = textToFloat(nextToken());
            spots_arg(l, k, h)  = textToFloat(nextToken());
            FOM(l, k, h)       = textToInteger(nextToken());
            switch (Space_Group)
            {
            case(1):
                            break;
            case(90)://P4212
                            if (h > k || h < 0 || l < 0)
                    {
                        std::cerr << "\nHORROR reflection outside the assymetric unit\n"
                        << "(h,k,l)=" << h << " " << k << " " << l << std::endl;
                        exit(1);
                        break;
                    }
            }//switch end
        }
        catch (Xmipp_error XE)
        {
            std::cout << XE;
            std::cout << "aph File: Line " << line_no << " is skipped due to an error\n";
        }
        line_no++;
    }/* while */
    // Close file
    fh_aph.close();

    fn_aph = fn;

}/*  APHFile2D::read */


/* ------------------------------------------------------------------------- */
void APHFile3D::clear()
{
    a = 0.;
    b = 0.;
    gamma = 0.;
    c = 0.;
    resolution = 0.;
    Amp_Scale = 0.;
    spots_abs.clear();
    spots_arg.clear();
    FOM.clear();
    Space_Group = 0;
    FOM.clear();
    fn_aph = "";
} /*clear*/
