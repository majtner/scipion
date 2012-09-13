/***************************************************************************
 *
 * Authors:    Sjors Scheres           scheres@cnb.csic.es (2004)
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

#include <data/args.h>
#include <data/funcs.h>
#include <data/metadata.h>
#include <data/image.h>
#include <data/filters.h>
#include <data/mask.h>

#include <data/symmetries.h>
#include <data/projection.h>
#include "directions.h"
#include "symmetrize.h"

#include <vector>

/**@defgroup BreakSym break_symmetry (Break symmetry of a set of reference volumes)
   @ingroup ReconsLibraryPrograms */
//@{
/** Break_Sym parameters. */
class Prog_Break_Sym_prm
{

public:
    // Various filenames
    FileName fn_sel, fn_vol, fn_sym, fn_root, fn_mask, fn_iter;
    // Selfile with experimental images or reference volumes
    MetaData SF, SFvol;
    // Symmetry list
    SymList SL;
    // Verbosity flag
    int verb;
    // dimension
    int dim;
    // Number of volumes to process
    int Nvols;
    // Iteration numbering
    int Niter, istart;
    // Convergence check
    double eps;
    // Radius for masking of volume
    double mask_radius;
    // File handler for the history file
    std::ofstream fh_hist;
    // Reference volumes
    std::vector< MultidimArray<double> > vols;
    // Mask
    Image<double> mask;

public:

    /// Read additional arguments for 3D-process from command line
    void read(int argc, char **argv);

    /// Usage
    void usage();

    /// Show
    void show();

    /// Project the reference volume in evenly sampled directions
    void process_one_image(Image<double> &img, int &opt_vol, int &opt_sym, double &maxcorr);

    /// Process various images
    void process_selfile(MetaData &SF, std::vector<MetaData> &SFout, double &avecorr);

    /// reconstruction by (weighted ART) or WBP
    void reconstruction(int argc, char **argv, MetaData &SF, int iter, int volno);

};
//@}
