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

#include <XmippData/xmippProgs.hh>
#include <XmippData/xmippArgs.hh>
#include <XmippData/xmippGeometry.hh>

class Statis_parameters: public Prog_parameters {
public:
   ImageXmipp  sumI, sumI2;
   VolumeXmipp sumV, sumV2;
   int         nI, nV;
public:
   Statis_parameters() {nI=nV=0;}
   void final_process();
};

void process_img(ImageXmipp &img, const Prog_parameters *prm) {
   Statis_parameters *eprm=(Statis_parameters *) prm;
   eprm->sumI().resize(img()); eprm->sumI2().resize(img());
   FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(img()) {
      MULTIDIM_ELEM(eprm->sumI(), i) += MULTIDIM_ELEM(img(),i);
      MULTIDIM_ELEM(eprm->sumI2(),i) += MULTIDIM_ELEM(img(),i)*MULTIDIM_ELEM(img(),i);
   }
   eprm->nI++;
}

void process_vol(VolumeXmipp &vol, const Prog_parameters *prm) {
   Statis_parameters *eprm=(Statis_parameters *) prm;
   eprm->sumV().resize(vol()); eprm->sumV2().resize(vol());
   FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(vol()) {
      MULTIDIM_ELEM(eprm->sumV(), i) += MULTIDIM_ELEM(vol(),i);
      MULTIDIM_ELEM(eprm->sumV2(),i) += MULTIDIM_ELEM(vol(),i)*MULTIDIM_ELEM(vol(),i);
   }
   eprm->nV++;
}

void Statis_parameters::final_process() {
   FileName fn_root=fn_in.without_extension();
   if (nI!=0) {
      FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(sumI()) {
         MULTIDIM_ELEM(sumI(), i) /= nI;
         MULTIDIM_ELEM(sumI2(),i) /= nI;
         MULTIDIM_ELEM(sumI2(),i) -= MULTIDIM_ELEM(sumI(),i)*MULTIDIM_ELEM(sumI(),i);
         MULTIDIM_ELEM(sumI2(),i) = sqrt(ABS(MULTIDIM_ELEM(sumI2(),i)));
      }
      sumI.write(fn_root+".med.xmp");
      sumI2.write(fn_root+".sig.xmp");
   }
   if (nV!=0) {
      FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(sumV()) {
         MULTIDIM_ELEM(sumV(), i) /= nV;
         MULTIDIM_ELEM(sumV2(),i) /= nV;
         MULTIDIM_ELEM(sumV2(),i) -= MULTIDIM_ELEM(sumV(),i)*MULTIDIM_ELEM(sumV(),i);
         MULTIDIM_ELEM(sumV2(),i) = sqrt(ABS(MULTIDIM_ELEM(sumV2(),i)));
      }
      sumV.write(fn_root+".med.vol");
      sumV2.write(fn_root+".sig.vol");
   }
}

int main (int argc, char **argv) {
   Statis_parameters prm;
   prm.each_image_produces_an_output=FALSE;
   // Set default action for application of header transformation
   prm.apply_geo=TRUE;
   SF_main(argc, argv, &prm, (void*)&process_img, (void*)&process_vol);
   prm.final_process();
}

/* ------------------------------------------------------------------------- */
/* Menu                                                                      */
/* ------------------------------------------------------------------------- */
/*Colimate:
   PROGRAM Statis {
      url="http://www.cnb.uam.es/~bioinfo/NewXmipp/Applications/Src/Statis/Help/statis.html";
      help="Compute the average and standard deviation of a set of images";
      OPEN MENU menu_selfile;
      COMMAND LINES {
        + usual: xmipp_statis $SELFILE_IN
      }
   }
*/
