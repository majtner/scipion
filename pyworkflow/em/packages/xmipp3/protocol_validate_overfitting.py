# **************************************************************************
# *
# * Authors:     Roberto Marabini (roberto@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************

from pyworkflow.protocol.params import (PointerParam, FloatParam, NumericListParam,
                                        StringParam, BooleanParam, LEVEL_ADVANCED)
from pyworkflow.em.data import Volume
from pyworkflow.em.protocol import ProtReconstruct3D
from pyworkflow.em.packages.xmipp3.convert import writeSetOfParticles
from pyworkflow.utils import getFloatListFromValues
from pyworkflow.utils.path import cleanPattern, cleanPath
import os
import xmipp

class XmippProtValidateOverfitting(ProtReconstruct3D):
    """    
    Check how the FSC changes with the number of projections used for 3D reconstruction. This method
    has been proposed by B. Heymann at ***
    """
    _label = 'validate overfitting'
    
    #--------------------------- DEFINE param functions --------------------------------------------   
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputParticles', PointerParam, pointerClass='SetOfParticles',pointerCondition='hasAlignmentProj',
                      label="Input particles",  
                      help='Select the input images from the project.')     
        form.addParam('symmetryGroup', StringParam, default='c1',
                      label="Symmetry group", 
                      help='See [[Xmipp Symmetry][http://www2.mrc-lmb.cam.ac.uk/Xmipp/index.php/Conventions_%26_File_formats#Symmetry]] page '
                           'for a description of the symmetry format accepted by Xmipp') 
        form.addParam('numberOfParticles', NumericListParam, default="100 200 500 1000 2000 5000", expertLevel=LEVEL_ADVANCED,
                      label="Number of particles") 
        form.addParam('maxRes', FloatParam, default=-1, expertLevel=LEVEL_ADVANCED,
                      label="Maximum resolution (A)",  
                      help='Maximum resolution (in Angstrom) to consider \n'
                           'in Fourier space (default Nyquist).\n'
                           'Param *--maxres* in Xmipp.') 
        form.addParam('pad', FloatParam, default=2, expertLevel=LEVEL_ADVANCED,
                      label="Padding factor")

        form.addParallelSection(threads=4, mpi=1)

    #--------------------------- INSERT steps functions --------------------------------------------

    def _createFilenameTemplates(self):
        """ Centralize how files are called for iterations and references. """
        myDict = {
            'input_xmd': self._getExtraPath('input_particles.xmd')
            }
        self._updateFilenamesDict(myDict)

    def _insertAllSteps(self):
        self._createFilenameTemplates()
        self._insertFunctionStep('convertInputStep')
        numberOfParticles=getFloatListFromValues(self.numberOfParticles.get())
        fractionCounter=0
        for number in numberOfParticles:
            if number<self.inputParticles.get().getSize():
                self._insertFunctionStep('reconstructionStep',number,fractionCounter)
                fractionCounter+=1
        #self._insertFunctionStep('createOutputStep')
        
    
    #--------------------------- STEPS functions --------------------------------------------
    def convertInputStep(self):
        particlesMd = self._getFileName('input_xmd')
        imgSet = self.inputParticles.get()
        writeSetOfParticles(imgSet, particlesMd)

    def reconstructionStep(self, numberOfImages, fractionCounter):
        fnRoot=self._getExtraPath("fraction%02d"%fractionCounter)
        Ts=self.inputParticles.get().getSamplingRate()
        for i in range(0,2):
            fnImgs=fnRoot+"_images_%02d.xmd"%i
            self.runJob("xmipp_metadata_utilities","-i %s -o %s --operate random_subset %d"%\
                        (self._getFileName('input_xmd'),fnImgs,numberOfImages),numberOfMpi=1)
        
            params =  '  -i %s' % fnImgs
            params += '  -o %s' % fnRoot+"_%02d.vol"%i
            params += ' --sym %s' % self.symmetryGroup.get()
            params += ' --max_resolution %0.3f' % self.maxRes.get()
            params += ' --padding %0.3f' % self.pad.get()
            params += ' --thr %d' % self.numberOfThreads.get()
            params += ' --sampling %f' % Ts
            self.runJob('xmipp_reconstruct_fourier', params)

        self.runJob('xmipp_resolution_fsc', "--ref %s -i %s -o %s --sampling_rate %f"%\
                    (fnRoot+"_00.vol",fnRoot+"_01.vol",fnRoot+"_fsc.xmd",Ts), numberOfMpi=1)
        cleanPattern(fnRoot+"_0?.vol")
        cleanPattern(fnRoot+"_images_0?.xmd")
        
        mdFSC = xmipp.MetaData(fnRoot+"_fsc.xmd")
        for id in mdFSC:
            fscValue = mdFSC.getValue(xmipp.MDL_RESOLUTION_FRC,id)
            maxFreq = mdFSC.getValue(xmipp.MDL_RESOLUTION_FREQREAL,id)
            if fscValue<0.5:
                break
        fh = open(fnRoot+"_freq.txt","w")
        fh.write("%f\n"%maxFreq)
        fh.close()
        cleanPath(fnRoot+"_fsc.xmd")
        
    def createOutputStep(self):
        imgSet = self.inputParticles.get()
        volume = Volume()
        volume.setFileName(self._getFileName('output_volume'))
        volume.setSamplingRate(imgSet.getSamplingRate())
        
        self._defineOutputs(outputVolume=volume)
        self._defineSourceRelation(self.inputParticles, volume)
    
    #--------------------------- INFO functions -------------------------------------------- 
    def _summary(self):
        """ Should be overriden in subclasses to 
        return summary message for NORMAL EXECUTION. 
        """
        msg=[]
        msg.append("Fraction of particles: "+self.fractionOfParticles.get())
        return msg
    
    #--------------------------- UTILS functions --------------------------------------------
