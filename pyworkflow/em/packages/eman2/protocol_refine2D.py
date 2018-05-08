# **************************************************************************
# *
# * Authors:     Tomas Majtner (tmajtner@cnb.csic.es)
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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import os
import re
from os.path import exists
from glob import glob
import pyworkflow.em as em
from pyworkflow.em.data import SetOfClasses2D
from pyworkflow.em.packages.eman2.eman2 import getEmanProgram, validateVersion
from pyworkflow.em.packages.eman2.convert import createEmanProcess
from pyworkflow.protocol.params import (PointerParam, FloatParam, IntParam, EnumParam,
                                        StringParam, BooleanParam)
from pyworkflow.utils.path import cleanPattern, makePath, createLink
from convert import rowToAlignment
from pyworkflow.em.packages.xmipp3.convert import xmippToLocation
import pyworkflow.em.metadata as md
                               
class EmanProtRefine2D(em.ProtRefine3D):
    _label = 'refine 2D'

    def _createFilenameTemplates(self):
        """ Centralize the names of the files. """
        
        myDict = {
                  'partSet': 'sets/inputSet.lst',
                  'partFlipSet': 'sets/inputSet__ctf_flip.lst',
                  'data_scipion': self._getExtraPath('data_scipion_it%(iter)02d.sqlite'),
                  'projections': self._getExtraPath('projections_it%(iter)02d_%(half)s.sqlite'),
                  'classes': self._getExtraPath('r2d_01/classes_%(iter)02d.hdf'),
                  'classesEven': self._getExtraPath('refine_%(run)02d/classes_%(iter)02d_even.hdf'),
                  'classesOdd': self._getExtraPath('refine_%(run)02d/classes_%(iter)02d_odd.hdf'),
                  'cls': 'refine_%(run)02d/cls_result_%(iter)02d',
                  'clsEven': self._getExtraPath('refine_%(run)02d/cls_result_%(iter)02d_even.hdf'),
                  'clsOdd': self._getExtraPath('refine_%(run)02d/cls_result_%(iter)02d_odd.hdf'),
                  'angles': self._getExtraPath('projectionAngles_it%(iter)02d.txt'),
                  'mapEven': self._getExtraPath('refine_%(run)02d/threed_%(iter)02d_even.hdf'),
                  'mapOdd': self._getExtraPath('refine_%(run)02d/threed_%(iter)02d_odd.hdf'),
                  'mapFull': self._getExtraPath('refine_%(run)02d/threed_%(iter)02d.hdf'),
                  'mapEvenUnmasked': self._getExtraPath('refine_%(run)02d/threed_even_unmasked.hdf'),
                  'mapOddUnmasked': self._getExtraPath('refine_%(run)02d/threed_odd_unmasked.hdf'),
                  'fscUnmasked': self._getExtraPath('refine_%(run)02d/fsc_unmasked_%(iter)02d.txt'),
                  'fscMasked': self._getExtraPath('refine_%(run)02d/fsc_masked_%(iter)02d.txt'),
                  'fscMaskedTight': self._getExtraPath('refine_%(run)02d/fsc_maskedtight_%(iter)02d.txt'),
                  }
        self._updateFilenamesDict(myDict)
    
    def _createIterTemplates(self, currRun):
        """ Setup the regex on how to find iterations. """
        self._iterTemplate = self._getFileName('mapFull', run=currRun,
                                               iter=1).replace('threed_01',
                                                               'threed_??')
        # Iterations will be identify by threed_XX_ where XX is the iteration
        #  number and is restricted to only 2 digits.
        self._iterRegex = re.compile('threed_(\d{2,2})')
    
    #--------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('doContinue', BooleanParam, default=False,
                      label='Continue from a previous run?',
                      help='If you set to *Yes*, you should select a previous'
                           'run of type *%s* class and most of the input parameters'
                           'will be taken from it.' % self.getClassName())
        form.addParam('inputParticles', PointerParam, label="Input particles",
                      important=True, pointerClass='SetOfParticles',
                      condition='not doContinue', allowsNull=True,
                      help='Select the input particles.\n')  
        form.addParam('input3DReference', PointerParam,
                      pointerClass='Volume', allowsNull=True,
                      label='Initial 3D reference volume:',
                      condition='not doContinue',
                      help='Input 3D reference reconstruction.\n')
        form.addParam('numberOfIterations', IntParam, default=2,
                      label='Number of iterations:',
                      help='Set the number of iterations. Iterative '
                           'reconstruction improves the overall normalization '
                           'of the 2D images as they are inserted into the '
                           'reconstructed volume, and allows for the '
                           'exclusion of the poorer quality'
                           'images.')
        form.addParam('symmetry', StringParam, default='c1',
                      condition='not doContinue',
                      label='Symmetry group',
                      help='Set the symmetry; if no value is given then the '
                           'model is assumed to have no symmetry. \n'
                           'Choices are: i(n), c(n), d(n), tet, icos, or oct.\n'
                           'See http://blake.bcm.edu/emanwiki/EMAN2/Symmetry'
                           'for a detailed descript of symmetry in Eman.')
        form.addParam('resol', FloatParam, default='10.0',
                      label='Resolution of this refinement run (A):',
                      help='Target resolution in A of this refinement run.'
                           'Usually works best in at least two steps'
                           '(low/medium) resolution, then final resolution)'
                           'when starting with a poor starting model.'
                           'Usually 3-4 iterations is sufficient.')
        form.addParam('molMass', FloatParam, default='500.0', 
                      label='Molecular mass of the specimen (kDa):',
                      help='Approximate molecular mass of the particle, in kDa.'
                           'This is used to runnormalize.bymass. Due to'
                           'resolution effects, not always the true mass.')
        form.addParam('doBreaksym', BooleanParam, default=False,
                       label='Do not impose symmetry?',
                       help='If set True, reconstruction will be asymmetric '
                            'with *Symmetry group* parameter specifying a '
                            'known pseudosymmetry, not an imposed symmetry.')
        form.addParam('useE2make3d', BooleanParam, default=False,
                       label='use e2make3d?',
                       help='Use the traditional e2make3d program instead of '
                            'the new e2make3dpar program.')
        form.addParam('classKeep', FloatParam, default='0.9',
                      label='Fraction of particles to use in final average:',
                      help='The fraction of particles to keep in each class,'
                           'based on the similarity score.')
        form.addParam('m3dKeep', FloatParam, default='0.8',
                      label='Fraction of class-averages to use in 3-D map:',
                      help='The fraction of slices to keep in reconstruction.')
        form.addParam('useSetsfref', BooleanParam, default=True,
                       label='Use the setsfref option in class averaging?',
                       help='This matches the filtration of the class-averages '
                            'to the projections for easier comparison. May '
                            'also improve convergence.')
        form.addParam('doAutomask', BooleanParam, default=False,
                       label='Do automask to the class-average?',
                       help='This will apply an automask to the class-average '
                            'during iterative alignment for better accuracy. '
                            'The final class averages are unmasked.')
        form.addParam('doThreshold', BooleanParam, default=False,
                       label='Apply threshold before project the volume?',
                       help='Applies a threshold to the volume just before '
                            'generating projections. A sort of aggressive '
                            'solvent flattening for the reference.')
        form.addParam('m3dPostProcess', EnumParam,
                      choices=['None', 'filter.highpass.autopeak',
                               'filter.highpass.butterworth',
                               'filter.highpass.gauss',
                               'filter.highpass.tanh',
                               'filter.highpassl.tophat',
                               'filter.lowpass.autob',
                               'filter.lowpass.butterworth',
                               'filter.lowpass.gauss',
                               'filter.lowpass.randomphase',
                               'filter.lowpass.tanh',
                               'filter.lowpass.tophat'],
                      label="Mode to Fourier method:", default=0,
                      display=EnumParam.DISPLAY_COMBO)
        form.addParallelSection(threads=4, mpi=0)
    
    #--------------------------- INSERT steps functions ------------------------
    def _insertAllSteps(self):        
        self._createFilenameTemplates()
        self._createIterTemplates(self._getRun())
        from pyworkflow.em.packages.eman2.convert import writeSetOfParticles
        partSet = self._getInputParticles()
        partAlign = partSet.getAlignment()
        storePath = self._getExtraPath("particles")
        makePath(storePath)
        writeSetOfParticles(partSet, storePath, alignType=partAlign)
        args = "--iter=6 --naliref=7 --nbasisfp=5 --input=particles/mic_000000.hdf --ncls=5 --simcmp=ccc --simalign=rotate_translate_flip --simaligncmp=ccc --simralign=refine --simraligncmp=ccc --classcmp=ccc --classalign=rotate_translate_flip --classaligncmp=ccc --classralign=refine --classraligncmp=ccc --classiter=2 --classkeep=1.5 --classnormproc=normalize.edgemean --classaverager=mean --normproj --classkeepsig --parallel=thread:4"
        program = getEmanProgram('e2refine2d.py')
        self.runJob(program, args, cwd=self._getExtraPath())

        set2D = self._createSetOfClasses2D(self.inputParticles.get())
        self._fillClassesFromLevel(set2D)
        result = {'outputClasses': set2D}
        self._defineOutputs(**result)
        self._defineSourceRelation(self.inputParticles, set2D)

        self._insertFunctionStep('createOutputStep')

    def _fillClassesFromLevel(self, clsSet):
        self._loadClassesInfo(self._getExtraPath("r2d_01/classes_06.hdf"))
        iterator = md.SetMdIterator(self.inputParticles.get().getFileName(),
                                    sortByLabel=md.MDL_ITEM_ID,
                                    updateItemCallback=self._updateParticle,
                                    skipDisabled=True)
        clsSet.classifyItems(updateItemCallback=iterator.updateItem,
                             updateClassCallback=self._updateClass)

    def _updateParticle(self, item, row):
        item.setClassId(row.getValue(md.MDL_REF))
        item.setTransform(rowToAlignment(row, '2D'))

    def _updateClass(self, item):
        classId = item.getObjId()
        if classId in self._classesInfo:
            index, fn, _ = self._classesInfo[classId]
            item.setAlignment2D()
            rep = item.getRepresentative()
            rep.setLocation(index, fn)
            rep.setSamplingRate(self.inputParticles.get().getSamplingRate())

    def _loadClassesInfo(self, filename):
        """ Read some information about the produced 2D classes
        from the metadata file.
        """
        self._classesInfo = {}  # store classes info, indexed by class id
        mdClasses = md.MetaData(filename)
        for classNumber, row in enumerate(md.iterRows(mdClasses)):
            index, fn = xmippToLocation(row.getValue(md.MDL_IMAGE))
            self._classesInfo[classNumber + 1] = (index, fn, row.clone())
    
    #--------------------------- STEPS functions -------------------------------
    def createLinkSteps(self):
        continueRun = self.continueRun.get()
        prevPartDir = continueRun._getExtraPath("particles")
        currPartDir = self._getExtraPath("particles")
        runN = self._getRun() - 1
        prevRefDir = continueRun._getExtraPath("refine_%02d" % runN)
        currRefDir = self._getExtraPath("refine_%02d" % runN)
        prevSetsDir = continueRun._getExtraPath("sets")
        currSetsDir = self._getExtraPath("sets")

#         createLink(prevInfoDir, currInfoDir)
        createLink(prevPartDir, currPartDir)
        createLink(prevRefDir, currRefDir)
        createLink(prevSetsDir, currSetsDir)
    
    def convertImagesStep(self):
        from pyworkflow.em.packages.eman2.convert import writeSetOfParticles
        partSet = self._getInputParticles()
        partAlign = partSet.getAlignment()
        storePath = self._getExtraPath("particles")
        makePath(storePath)
        writeSetOfParticles(partSet, storePath, alignType=partAlign)
        if partSet.hasCTF():
            program = getEmanProgram('e2ctf.py')
            acq = partSet.getAcquisition()
            
            args = " --voltage %3d" % acq.getVoltage()
            args += " --cs %f" % acq.getSphericalAberration()
            args += " --ac %f" % (100 * acq.getAmplitudeContrast())
            if not partSet.isPhaseFlipped():
                args += " --phaseflip"
            args += " --computesf --apix %f --allparticles --autofit --curdefocusfix --storeparm -v 8" % (partSet.getSamplingRate())
            self.runJob(program, args, cwd=self._getExtraPath())
        
        program = getEmanProgram('e2buildsets.py')
        args = " --setname=inputSet --allparticles --minhisnr=-1"
        self.runJob(program, args, cwd=self._getExtraPath())
    
    def refineStep(self, args):
        """ Run the EMAN program to refine a volume. """
        if not self.doContinue:
            cleanPattern(self._getExtraPath('refine_01'))
        program = getEmanProgram('e2refine_easy.py')
        self.runJob(program, args, cwd=self._getExtraPath())
    
    def createOutputStep(self):
        partSet = self._getInputParticles()
        sets2D = self._createSetOfClasses2D(partSet)
        from pyworkflow.em.packages.eman2.convert import readSetOfClasses2D
        readSetOfClasses2D(sets2D, self._getPath('../034712_EmanProtRefine2D/extra/r2d_01/classes_06.hdf'))
        outputSet = {'outputClasses': sets2D}
        # self._defineOutputs(outputVolumes=volumes)
        self._defineOutputs(**outputSet)
        self._defineSourceRelation(self._getInputParticlesPointer(), sets2D)
        #self._defineTransformRelation(self._getInputParticlesPointer(), newPartSet)
    
    #--------------------------- INFO functions -------------------------------------------- 
    def _validate(self):
        errors = []
        validateVersion(self, errors)

        particles = self._getInputParticles()
        samplingRate = particles.getSamplingRate()

        if self.resol <  2 * samplingRate:
            errors.append("\nTarget resolution is smaller than nyquist limit.")
        
        if not self.doContinue:
            self._validateDim(particles, self.input3DReference.get(), errors,
                              'Input particles', 'Reference volume')

        return errors
    
    def _summary(self):
        summary = []
        if not hasattr(self, 'outputVolume'):
            summary.append("Output volumes not ready yet.")
        else:
            inputSize = self._getInputParticles().getSize()
            outputSize = self.outputParticles.getSize()
            diff = inputSize - outputSize
            if diff > 0:
                summary.append("Warning!!! There are %d particles "
                               "belonging to empty classes." % diff)
        return summary
    
    #--------------------------- UTILS functions --------------------------------------------
    def _prepareParams(self):
        args1 = "--input=%(imgsFn)s --model=%(volume)s"
        args2 = self._commonParams()
        
        volume = os.path.relpath(self.input3DReference.get().getFileName(), self._getExtraPath()).replace(":mrc","")
        params = {'imgsFn': self._getParticlesStack(),
                  'volume': volume,
                  }
        
        args = args1 % params + args2
        return args
    
    def _prepareContinueParams(self):
        args1 = "--startfrom=refine_%02d" % (self._getRun() - 1)
        args2 = self._commonParams()
        args = args1 + args2
        return args
    
    def _commonParams(self):
        args = " --targetres=%(resol)f --speed=%(speed)d --sym=%(sym)s --iter=%(numberOfIterations)d"
        args += " --mass=%(molMass)f --apix=%(samplingRate)f --classkeep=%(classKeep)f"
        args += " --m3dkeep=%(m3dKeep)f --parallel=thread:%(threads)d --threads=%(threads)d"
        
        samplingRate = self._getInputParticles().getSamplingRate()
        params = {'resol': self.resol.get(),
                  'speed': int(self.getEnumText('speed')),
                  'numberOfIterations': self.numberOfIterations.get(),
                  'sym': self.symmetry.get(),
                  'molMass': self.molMass.get(),
                  'samplingRate': samplingRate,
                  'classKeep': self.classKeep.get(),
                  'm3dKeep': self.m3dKeep.get(),
                  'threads': self.numberOfThreads.get()
                  }
        args = args % params
         
        if self.doBreaksym:
            args += " --breaksym"
        if self.useE2make3d:
            args += " --m3dold"
        if self.useSetsfref:
            args += " --classrefsf"
        if self.doAutomask:
            args += " --classautomask"
        if self.doThreshold:
            args += " --prethreshold"
        if self.m3dPostProcess.get() > 0:
            args += " --m3dpostprocess=%s" % self.getEnumText('m3dPostProcess')
        return args
    
    def _getRun(self):
        if not self.doContinue:
            return 1
        else:
            files = sorted(glob(self.continueRun.get()._getExtraPath("refine*")))
            if files:
                f = files[-1]
                refineNumber = int(f.split("_")[-1]) + 1
            return refineNumber
    
    def _getBaseName(self, key, **args):
        """ Remove the folders and return the file from the filename. """
        return os.path.basename(self._getFileName(key, **args))
    
    def _getParticlesStack(self):
        if not self._getInputParticles().isPhaseFlipped() and self._getInputParticles().hasCTF():
            return self._getFileName("partFlipSet")
        else:
            return self._getFileName("partSet")

    
    def _createItemMatrix(self, item, rowList):
        if rowList[1] == 1:
            item.setTransform(rowToAlignment(rowList[2:], alignType=em.ALIGN_PROJ))
        else:
            setattr(item, "_appendItem", False)
    
    def _getIterNumber(self, index):
        """ Return the list of iteration files, give the iterTemplate. """
        result = None
        files = sorted(glob(self._iterTemplate))
        if files:
            f = files[index]
            s = self._iterRegex.search(f)
            if s:
                result = int(s.group(1)) # group 1 is 3 digits iteration number
                
        return result
    
    def _lastIter(self):
        return self._getIterNumber(-1)

    def _firstIter(self):
        return self._getIterNumber(0) or 1
    
    def _getIterData(self, it):
        data_sqlite = self._getFileName('data_scipion', iter=it)
        if not exists(data_sqlite):
            iterImgSet = em.SetOfParticles(filename=data_sqlite)
            iterImgSet.copyInfo(self._getInputParticles())
            self._fillDataFromIter(iterImgSet, it)
            iterImgSet.write()
            iterImgSet.close()
        
        return data_sqlite
    
    def _getInputParticlesPointer(self):
        if self.doContinue:
            self.inputParticles.set(self.continueRun.get().inputParticles.get())
        return self.inputParticles
    
    def _getInputParticles(self):
        return self._getInputParticlesPointer().get()
    
    def _fillDataFromIter(self, imgSet, iterN):
        numRun = self._getRun()
        self._execEmanProcess(numRun, iterN)
        initPartSet = self._getInputParticles()
        imgSet.setAlignmentProj()
        partIter = iter(initPartSet.iterItems(orderBy=['_micId', 'id'],
                                              direction='ASC'))
        
        imgSet.copyItems(partIter)
    
    def _execEmanProcess(self, numRun, iterN):
        clsFn = self._getFileName("cls", run=numRun, iter=iterN)
        classesFn = self._getFileName("classes", run=numRun, iter=iterN)
        angles = self._getFileName('angles', iter=iterN)
        
        if not exists(angles) and exists(self._getFileName('clsEven', run=numRun, iter=iterN)):
            proc = createEmanProcess(args='read %s %s %s %s'
                                     % (self._getParticlesStack(), clsFn, classesFn,
                                        self._getBaseName('angles', iter=iterN)),
                                        direc=self._getExtraPath())
            proc.wait()
