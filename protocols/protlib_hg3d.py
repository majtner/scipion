#!/usr/bin/env xmipp_python
'''
#/***************************************************************************
# * Authors:     Javier Vargas (jvargas@cnb.csic.es)
# *
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
# *  e-mail address 'xmipp@cnb.csic.es'
# ***************************************************************************
'''
# This library contains some common utilities 
# for all particles related protocols: Extract, Import
from protlib_base import *
from xmipp import MetaData, MetaDataInfo
from protlib_utils import runJob, runShowJ

class ProtHG3DBase(XmippProtocol):
    '''This class will serve as base for init volume related protocols'''
    def __init__(self, protocolName, scriptname, project):
        XmippProtocol.__init__(self, protocolName, scriptname, project)        
        self.Import = 'from protlib_hg3d import *; '
        self.Xdim=MetaDataInfo(self.Classes)[0]
        
    def defineSteps(self):
        self.insertStep('createDir',path=self.ExtraDir)
        self.insertStep("linkAcquisitionInfo",InputFile=self.Classes,dirDest=self.WorkingDir)
        fnOutputReducedClass = self.extraPath("reducedClasses")
    
        # Low pass filter and resize        
        self.MaxFreq = float(self.MaxFreq)
        self.Ts = float(self.Ts)
        K = 0.25*(self.MaxFreq/self.Ts)
        if K<1:
            K=1
        self.Xdim2 = self.Xdim/K
        if (self.Xdim2 < 32):
            self.Xdim2 = 32
            K = self.Xdim/self.Xdim2
            
        freq = self.Ts/self.MaxFreq
        self.Ts = K*self.Ts

        self.insertRunJobStep("xmipp_transform_filter","-i %s -o %s.stk --fourier low_pass %f --save_metadata_stack %s.xmd"
                                                %(self.Classes,fnOutputReducedClass,freq,fnOutputReducedClass))
        self.insertRunJobStep("xmipp_image_resize","-i %s.stk --dim %d %d" %(fnOutputReducedClass,self.Xdim2,self.Xdim2))

    def summary(self):
        message=[]
        message.append("Input images: [%s]"%self.Classes)
        if self.InitialVolume!="":
            message.append("Initial volume: [%s]"%self.InitialVolume)
        message.append("Symmetry: "+self.SymmetryGroup)
        return message

