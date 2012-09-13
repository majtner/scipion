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

#include "multidimensional_array.h"

template<>
std::ostream& operator<<(std::ostream& ostrm,
    const MultidimArray< std::complex<double> >& v)
{
    if (v.xdim == 0)
        ostrm << "NULL Array\n";
    else
        ostrm << std::endl;

    for (int k = STARTINGZ(v); k <= FINISHINGZ(v); k++)
    {
        ostrm << "Slice No. " << k << std::endl;
        for (int i = STARTINGY(v); i <= FINISHINGY(v); i++)
        {
            for (int j = STARTINGX(v); j <= FINISHINGX(v); j++)
            {
                ostrm << DIRECT_VOL_ELEM(v, k, i, j) << ' ';
            }
            ostrm << std::endl;
        }
    }

    return ostrm;
}
