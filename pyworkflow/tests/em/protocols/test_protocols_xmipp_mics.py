# **************************************************************************
# *
# * Authors:    Laura del Cano (ldelcano@cnb.csic.es)
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

import unittest, sys
from pyworkflow.em import *
from pyworkflow.tests import *
from pyworkflow.em.packages.xmipp3 import *


# Some utility functions to import micrographs that are used
# in several tests.
class TestXmippBase(unittest.TestCase):
    
    @classmethod
    def runImportMicrograph(cls, pattern, samplingRate, voltage, scannedPixelSize, magnification, sphericalAberration):
        """ Run an Import micrograph protocol. """
        # We have two options: passe the SamplingRate or the ScannedPixelSize + microscope magnification
        if not samplingRate is None:
            cls.protImport = ProtImportMicrographs(samplingRateMode=0, pattern=pattern, samplingRate=samplingRate, magnification=magnification, 
                                                   voltage=voltage, sphericalAberration=sphericalAberration)
        else:
            cls.protImport = ProtImportMicrographs(samplingRateMode=1, pattern=pattern, scannedPixelSize=scannedPixelSize, 
                                                   voltage=voltage, magnification=magnification, sphericalAberration=sphericalAberration)
            
        cls.proj.launchProtocol(cls.protImport, wait=True)
        # check that input micrographs have been imported (a better way to do this?)
        if cls.protImport.outputMicrographs is None:
            raise Exception('Import of micrograph: %s, failed. outputMicrographs is None.' % pattern)
        return cls.protImport
        
    @classmethod
    def runImportMicrographBPV1(cls):
        """ Run an Import micrograph protocol. """
        pattern = getInputPath('Micrographs_BPV1', '*.mrc')
        return cls.runImportMicrograph(pattern, samplingRate=1.237, voltage=300, sphericalAberration=2, scannedPixelSize=None, magnification=56000)
    
    @classmethod
    def runFakedPicking(cls, mics, pattern):
        """ Run a faked particle picking. Coordinates already existing. """
        coordsFolder = getInputPath(pattern)
        cls.protPP = XmippProtParticlePicking(importFolder=coordsFolder, runMode=1)                
        cls.protPP.inputMicrographs.set(mics)               
        cls.proj.launchProtocol(cls.protPP, wait=True)
        # check that faked picking has run ok
        if cls.protPP.outputCoordinates is None:
            raise Exception('Faked particle picking: %s, failed. outputCoordinates is None.' % coordsFolder)
        return cls.protPP


class TestImportMicrographs(TestXmippBase):
    """This class check if a set of micrographs is imported properly"""
    @classmethod
    def setUpClass(cls):
        setupProject(cls)
    
    def testImport1(self):
        pattern = getInputPath('Micrographs_BPV1', '*.mrc')
        samplingRate = None
        scannedPixelSize = 7
        magnification = 56000
        voltage = 300
        sphericalAberration = 2
        
        protImport = self.runImportMicrograph(pattern, samplingRate=samplingRate, 
                                              scannedPixelSize=scannedPixelSize, 
                                              magnification=magnification, voltage=voltage, 
                                              sphericalAberration=sphericalAberration)
        if protImport.isFailed():
            raise Exception(protImport.getError())
        
        m = protImport.outputMicrographs.getAcquisition()
        # Check that sampling rate on output micrographs is equal to 
        self.assertEquals(protImport.outputMicrographs.getScannedPixelSize(), scannedPixelSize, "Incorrect ScannedPixelSize on output micrographs.")
        self.assertEquals(m.getMagnification(), magnification, "Incorrect Magnification on output micrographs.")
        self.assertEquals(m.getVoltage(), voltage, "Incorrect Voltage on output micrographs.")
        self.assertEquals(m.getSphericalAberration(), sphericalAberration, "Incorrect SphericalAberration on output micrographs.")

    def testImport2(self):
        pattern = getInputPath('Micrographs_BPV3', '*.mrc')
        samplingRate = 2.56
        scannedPixelSize = 7
        magnification = 56000
        voltage = 400
        sphericalAberration = 2.5
        
        protImport = self.runImportMicrograph(pattern, samplingRate=samplingRate, 
                                              scannedPixelSize=scannedPixelSize, 
                                              magnification=magnification, voltage=voltage, 
                                              sphericalAberration=sphericalAberration)
        m = protImport.outputMicrographs.getAcquisition()
        # Check that sampling rate on output micrographs is equal to 
        self.assertEquals(protImport.outputMicrographs.getSamplingRate(), samplingRate, "Incorrect SamplingRate on output micrographs.")
        self.assertEquals(m.getVoltage(), voltage, "Incorrect Voltage on output micrographs.")
        self.assertEquals(m.getSphericalAberration(), sphericalAberration, "Incorrect Spherical aberration on output micrographs.")


class TestXmippPreprocessMicrographs(TestXmippBase):
    """This class check if the preprocessing micrographs protocol in Xmipp works properly."""
    @classmethod
    def setUpClass(cls):
        setupProject(cls)
        cls.protImport = cls.runImportMicrographBPV1()
    
    def testDownsampling(self):
        # test downsampling a set of micrographs
        downFactor = 2
        protDown = XmippProtPreprocessMicrographs(doDownsample=True, downFactor=downFactor)
        protDown.inputMicrographs.set(self.protImport.outputMicrographs)
        self.proj.launchProtocol(protDown, wait=True)
        
        # check that output micrographs have double sampling rate than input micrographs
        self.assertEquals(protDown.outputMicrographs.getSamplingRate(), self.protImport.outputMicrographs.getSamplingRate()*downFactor, "Micrographs uncorrectly downsampled")
    
    def testCrop(self):
        # test crop on a set of micrographs
        cropPixels = 100
        protCrop = XmippProtPreprocessMicrographs(doCrop=True, cropPixels=cropPixels)
        protCrop.inputMicrographs.set(self.protImport.outputMicrographs)
        self.proj.launchProtocol(protCrop, wait=True)


class TestXmippCTFEstimation(TestXmippBase):
    """This class check if the protocol to determine the CTF in Xmipp works properly."""
    @classmethod
    def setUpClass(cls):
        setupProject(cls)
            
    def doCTF(self, pattern):
        #First, import a set of micrographs
        protImport = self.runImportMicrograph(pattern, samplingRate=3.711, voltage=300, sphericalAberration=2, scannedPixelSize=None, magnification=56000)
        
        # Now estimate CTF on the downsampled micrographs
        print "Performing CTF..."   
        protCTF = XmippProtCTFMicrographs()                
        protCTF.inputMicrographs.set(protImport.outputMicrographs)        
        self.proj.launchProtocol(protCTF, wait=True)
        
        #self.assertTrue(protCTF.outputMicrographs.hasCTF(), "CTF estimation has not been performed.")
        #self.assertEqual(protCTF.outputMicrographs._xmippMd.get(),protCTF._getPath("micrographs.xmd"), "Xmipp md not set on output.")
        self.assertIsNotNone(protCTF.outputCTF, "SetOfCTF has not been produced.") 
        
    def test_Micrographs_BPV1_Down3(self):
        self.doCTF(pattern = getInputPath('Micrographs_BPV1_Down3', '*.mrc'))
        
    def test_Micrographs_BPV3_Down3(self):
        self.doCTF(pattern = getInputPath('Micrographs_BPV3_Down3', '*.mrc')) 


class TestXmippAutomaticPicking(TestXmippBase):
    """This class check if the protocol to pick the micrographs automatically in Xmipp works properly."""
    @classmethod
    def setUpClass(cls):
        setupProject(cls)    
        pattern = getInputPath('Micrographs_BPV3_Down3', '*.mrc')
        protImport = cls.runImportMicrograph(pattern, samplingRate=1.237, voltage=300, sphericalAberration=2, scannedPixelSize=None, magnification=56000)       
        pattern = getInputPath('Micrographs_BPV2_Down3', '*.mrc')
        cls.protImport_other = cls.runImportMicrograph(pattern, samplingRate=1.237, voltage=300, sphericalAberration=2, scannedPixelSize=None, magnification=56000)        
        cls.protPP = cls.runFakedPicking(protImport.outputMicrographs, 'Picking_XmippBPV3_Down3_Super')

    def testAutomaticPicking(self):
        print "Run automatic particle picking"
        protAutomaticPP = XmippParticlePickingAutomatic()
        protAutomaticPP.xmippParticlePicking.set(self.protPP)
        self.proj.launchProtocol(protAutomaticPP, wait=True)
        
        self.assertIsNotNone(protAutomaticPP.outputCoordinates, "There was a problem with the automatic particle picking")
        
    def testAutomaticPickingOther(self):
        print "Run automatic particle picking"
        protAutomaticPP = XmippParticlePickingAutomatic()
        protAutomaticPP.xmippParticlePicking.set(self.protPP)
        protAutomaticPP.inputMicrographs.set(self.protImport_other.outputMicrographs)
        protAutomaticPP.micsToPick.set(1)
        self.proj.launchProtocol(protAutomaticPP, wait=True)
        
        self.assertIsNotNone(protAutomaticPP.outputCoordinates, "There was a problem with the automatic particle picking")


class TestXmippExtractParticles(TestXmippBase):
    """This class check if the protocol to extract particles in Xmipp works properly."""
    SAME_AS_PICKING = 1
    ORIGINAL = 0
    OTHER = 2
    
    @classmethod
    def setUpClass(cls):

        setupProject(cls)    
        pattern = getInputPath('Micrographs_BPV3_Down3', '*.mrc')
        
        protImport = cls.runImportMicrograph(pattern, samplingRate=3.711, voltage=300, sphericalAberration=2, scannedPixelSize=None, magnification=56000)       
        cls.proj.launchProtocol(protImport, wait=True)
            
        cls.protCTF = XmippProtCTFMicrographs()                
        cls.protCTF.inputMicrographs.set(protImport.outputMicrographs)     
        cls.proj.launchProtocol(cls.protCTF, wait=True)        
        
        cls.protPP = cls.runFakedPicking(cls.protCTF.inputMicrographs.get(), 'Picking_XmippBPV3_Down3_Super')

    def testExtractSameAsPicking(self):
        print "Run extract particles with downsampling factor equal to the one at picking"
        protExtract = XmippProtExtractParticles(boxSize=171, downsampleType=self.SAME_AS_PICKING, 
                                                doFlip=False)
        
        protExtract.inputMicrographs.set(self.protImport.outputMicrographs)
        protExtract.inputCoordinates.set(self.protPP.outputCoordinates)
        protExtract.setObjLabel("extract-same as picking")
        self.proj.launchProtocol(protExtract, wait=True)
        self.assertIsNotNone(protExtract.outputParticles, "There was a problem with the extract particles")
        
    def testExtractOriginal(self):
        print "Run extract particles with downsampling factor equal to the original micrographs"
        protExtract = XmippProtExtractParticles(boxSize=256, downsampleType=self.ORIGINAL, doFlip=False)
        protExtract.inputCoordinates.set(self.protPP.outputCoordinates)
        protExtract.inputMicrographs.set(self.protImport.outputMicrographs)
        protExtract.setObjLabel("extract-original")
        self.proj.launchProtocol(protExtract, wait=True)
        
        self.assertIsNotNone(protExtract.outputParticles, "There was a problem with the extract particles")

    def testExtractOther(self):
        print "Run extract particles with downsampling factor equal to other"
        protExtract = XmippProtExtractParticles(boxSize=110, downsampleType=self.OTHER, downFactor=2,doFlip=False)
        protExtract.inputCoordinates.set(self.protPP.outputCoordinates)
        protExtract.inputMicrographs.set(self.protImport.outputMicrographs)
        protExtract.setObjLabel("extract-other")
        self.proj.launchProtocol(protExtract, wait=True)
        
        self.assertIsNotNone(protExtract.outputParticles, "There was a problem with the extract particles")
        
    def testExtractCTF(self):
        print "Run extract particles with CTF"#        
        protExtract = XmippProtExtractParticles(boxSize=250, downsampleType=self.ORIGINAL,doFlip=True)
        protExtract.inputCoordinates.set(self.protPP.outputCoordinates)
        protExtract.inputMicrographs.set(self.protCTF.inputMicrographs.get())
        protExtract.ctfRelations.set(self.protCTF.outputCTF)
        protExtract.setObjLabel("extract-ctf")
        self.proj.launchProtocol(protExtract, wait=True)
        
        self.assertIsNotNone(protExtract.outputParticles, "There was a problem with the extract particles") 
        
        protEmx1 = ProtEmxExport()
        protEmx1.setObjLabel("emx export coordinates")
        protEmx1.inputSet.set(self.protPP.outputCoordinates)
        self.proj.launchProtocol(protEmx1, wait=True)
        
        protEmx2 = ProtEmxExport()
        protEmx2.setObjLabel("emx export particles")
        protEmx2.inputSet.set(protExtract.outputParticles)
        self.proj.launchProtocol(protEmx2, wait=True)


if __name__ == "__main__":
    if len(sys.argv) > 1:
        className = sys.argv[1]
        cls = globals().get(className, None)
        if cls:
            suite = unittest.TestLoader().loadTestsFromTestCase(cls)
            unittest.TextTestRunner(verbosity=2).run(suite)
        else:
            print "Test: '%s' not found." % className
    else:
        unittest.main()