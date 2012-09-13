/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
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

#include <data/progs.h>
#include <data/args.h>
#include <data/docfile.h>
#include <data/geometry.h>

class Shift_Scale_parameters: public Prog_parameters
{
public:
    Matrix1D<double> shift;
    Matrix1D<double> scale;
    bool             wrap;
    bool             store_in_header;
    bool             center_mass;
    DocFile          DF_shifts;
    DocFile          DF_scales;
    int              colX_shift;
    int              colX_scale;
    bool             Docfile;


    void read(int argc, char **argv)
    {
        Prog_parameters::read(argc, argv);
        int i_shift = paremeterPosition(argc, argv, "-shift");
        int i_scale = paremeterPosition(argc, argv, "-scale");
        center_mass = checkParameter(argc, argv, "-center_mass");

        if (i_shift == -1 && i_scale == -1 && !center_mass)
            REPORT_ERROR(1, "Shift_Scale:: Cannot find -shift or -scale");
        else if (ABS(i_shift - i_scale) <= 1 && !center_mass)
            REPORT_ERROR(1, "Shift_cale: Not enough parameters after -shift or -scale");
        Docfile = checkParameter(argc, argv, "-colX_shift") ||
                  checkParameter(argc, argv, "-colX_scale");
        if (Docfile)
        {
            if (i_shift > 0)
            {
                DF_shifts.read(argv[i_shift+1]);
                colX_shift = textToInteger(getParameter(argc, argv, "-colX_shift"));
                // colX_shift -=3;
                // if(colX_shift<0)
                //   REPORT_ERROR(1,"colX_shift must be no less than 3");
            }
            else colX_shift = -1;

            if (i_scale > 0)
            {
                DF_scales.read(argv[i_scale+1]);
                colX_scale = textToInteger(getParameter(argc, argv, "-colX_scale", "5"));
                colX_scale -= 3;
                if (colX_scale < 0)
                    REPORT_ERROR(1, "colX_scale must be no less than 3");
            }
            else colX_scale = -1;
        }
        else
        {
            int my_dim=0;
            if (checkParameter(argc, argv, "-shift"))
            {
                shift = getVectorParameter(argc, argv, "-shift", -1);
                my_dim = XSIZE(shift);
            }
            if (checkParameter(argc, argv, "-scale"))
            {
                scale = getVectorParameter(argc, argv, "-scale", -1);
                my_dim = XSIZE(scale);
            }

            if (!checkParameter(argc, argv, "-shift") && !center_mass)
                shift.initZeros(my_dim);
            if (!checkParameter(argc, argv, "-scale") && !center_mass)
            {
                scale.resize(my_dim);
                scale.initConstant(1);
            }
        }
        wrap = !checkParameter(argc, argv, "-dont_wrap");
        store_in_header = checkParameter(argc, argv, "-store_in_header");
    }

    void show()
    {
        Prog_parameters::show();
        if (wrap) std::cout << "Wrapping image/volume\n";
        else      std::cout << "Not wrapping image/volume\n";
        if (store_in_header) std::cout << "Storing the shift in header\n";
        else                 std::cout << "Shifting image/volume\n";
        if (XSIZE(shift) > 1)
            std::cout << "Shift: " << shift.transpose() << std::endl;
        else if (DF_shifts.name() != "")
        {
            std::cout << "Shift docfile: " << DF_shifts.name() << std::endl;
            std::cout << "colX_shift:  " << colX_shift << std::endl;
        }
        if (XSIZE(scale) > 1)
            std::cout << "Scale: " << scale.transpose() << std::endl;
        else if (DF_scales.name() != "")
        {
            std::cout << "Scale: docfile: "       << DF_scales.name() << std::endl;
            std::cout << "colX_scale:  " << colX_scale << std::endl;
        }
        if (center_mass) std::cout << "Moving center of mass to origin\n";
    }

    void usage()
    {
        Prog_parameters::usage();
        std::cerr << "   -shift \"[<x>,<y>[,<z>]]\" : Shift by (x,y,z) for volumes, (x,y) for images\n"
        << "   -scale \"[<x>,[<y>,<z>]]\" : Scale by (x,y,z)\n"
        << "   -shift <DocFile>         : Shifts are stored in a Docfile\n"
        << "   -scale <DocFile>         : Scales are stored in a Docfile (may be the same\n"
        << "                              Docfile used for shifts\n"
        << "  [-colX_shift <col>]       : Column with  the X shift\n"
        << "                              First column in the DocFile with data is number 0.\n"
        << "  [-colX_scale <col>]       : Column with the scale information\n"
        << "  [-center_mass]            : Move the center of mass to the origin\n"
        << "  [-dont_wrap]              : By default, the image is wrapped\n"
        << "  [-store_in_header]        : Do not shift, but store the shift in the header\n"
        << "                              Not applicable to volumes\n";
    }
};

bool process_img(ImageXmipp &img, const Prog_parameters *prm)
{
    Shift_Scale_parameters *eprm = (Shift_Scale_parameters *) prm;
    Matrix2D<double> A(3, 3);
    A.initIdentity();

    if (eprm->DF_shifts.name() != "")
    {
        eprm->shift.resize(2);
        XX(eprm->shift) = eprm->DF_shifts(eprm->colX_shift);
        YY(eprm->shift) = eprm->DF_shifts(eprm->colX_shift + 1);
        eprm->DF_shifts.next_data_line();
    }
    else if (eprm->Docfile)
    {
        eprm->shift.resize(2);
        eprm->shift.initConstant(0.);
    }
    else if (eprm->center_mass)
    {
        img().centerOfMass(eprm->shift);
        eprm->shift = -eprm->shift;
    }

    A(0, 2) = XX(eprm->shift);
    A(1, 2) = YY(eprm->shift);
    if (eprm->DF_scales.name() != "")
    {
        eprm->scale.resize(2);
        XX(eprm->scale) = eprm->DF_scales(eprm->colX_scale);
        YY(eprm->scale) = eprm->DF_scales(eprm->colX_scale + 1);
        eprm->DF_scales.next_data_line();
    }
    else if (eprm->Docfile || eprm->center_mass)
    {
        eprm->scale.resize(2);
        eprm->scale.initConstant(1.);
    }
    A(0, 0) = XX(eprm->scale);
    A(1, 1) = YY(eprm->scale);
    if (!eprm->store_in_header) img().selfApplyGeometryBSpline(A, 3, IS_NOT_INV, eprm->wrap);
    else                        img.set_originOffsets(XX(eprm->shift), YY(eprm->shift));
    return true;
}

bool process_vol(VolumeXmipp &vol, const Prog_parameters *prm)
{
    Matrix2D<double> A(4, 4);
    A.initIdentity();
    Shift_Scale_parameters *eprm = (Shift_Scale_parameters *) prm;
    if (eprm->DF_shifts.name() != "")
    {
        eprm->shift.resize(3);
        XX(eprm->shift) = eprm->DF_shifts(eprm->colX_shift);
        YY(eprm->shift) = eprm->DF_shifts(eprm->colX_shift + 1);
        ZZ(eprm->shift) = eprm->DF_shifts(eprm->colX_shift + 2);
        eprm->DF_shifts.next_data_line();
    }
    else if (eprm->Docfile)
    {
        eprm->shift.resize(3);
        eprm->shift.initConstant(1.);
    }
    else if (eprm->center_mass)
    {
        vol().centerOfMass(eprm->shift);
        eprm->shift = -eprm->shift;
    }
    A(0, 3) = XX(eprm->shift);
    A(1, 3) = YY(eprm->shift);
    A(2, 3) = ZZ(eprm->shift);

    if (eprm->DF_scales.name() != "")
    {
        eprm->scale.resize(3);
        XX(eprm->scale) = eprm->DF_scales(eprm->colX_scale);
        YY(eprm->scale) = eprm->DF_scales(eprm->colX_scale + 1);
        ZZ(eprm->scale) = eprm->DF_scales(eprm->colX_scale + 2);
        eprm->DF_scales.next_data_line();
    }
    else if (eprm->Docfile || eprm->center_mass)
    {
        eprm->scale.resize(3);
        eprm->scale.initConstant(1.);
    }
    A(0, 0) = XX(eprm->scale);
    A(1, 1) = YY(eprm->scale);
    A(2, 2) = ZZ(eprm->scale);
    vol().selfApplyGeometryBSpline(A, 3, IS_NOT_INV, eprm->wrap);
    return true;
}

int main(int argc, char **argv)
{
    Shift_Scale_parameters prm;
    SF_main(argc, argv, &prm, (void*)&process_img, (void*)&process_vol);
}

/* Menus ------------------------------------------------------------------- */
/*Colimate:
   PROGRAM Shift {
      url="http://www.cnb.uam.es/~bioinfo/NewXmipp/Applications/Src/Shift/Help/shift.html";
      help="Shift volumes and images";
      OPEN MENU menu_shift;
      COMMAND LINES {
 + usual: xmipp_shift
               #include "prog_line.mnu"
                -shift "["$X","$Y[","$Z]"]"
               [-dont_wrap]
      }
      PARAMETER DEFINITIONS {
        #include "prog_vars.mnu"
            $X  {type=float; label="Shift X ";}
            $Y  {type=float; label="Shift Y ";}
        OPT($Z) {type=float; label="Shift Z ";}
        OPT(-dont_wrap) {label="Do not wrap";}
      }
   }

   MENU menu_shift {
      #include "prog_menu.mnu"
      "Shift parameters"
      $X
      $Y
      OPT($Z)
      OPT(-dont_wrap)
   }
*/
