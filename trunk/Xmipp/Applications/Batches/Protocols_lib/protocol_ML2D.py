#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
# Protocol for Xmipp-based 2D alignment and classification,
# using maximum-likelihood principles, according to:
# - Scheres et al. (2005) J.Mol.Biol 348, 139-149
# - Scheres et al. (2005) Bioinformatics, 21(suppl2), ii243-244
#
# Example use:
# ./protocol_ML2D.py
#
# Author: Sjors Scheres, March 2007
#
#------------------------------------------------------------------------------------------------
# {section} Global parameters
#------------------------------------------------------------------------------------------------
# Selfile with the input images (relative path from ProjectDir):
InSelFile="all_images.sel"
# Working subdirectory:
WorkingDir="ML3ref"
# Delete working subdirectory if it already exists?
""" The directory will not be deleted when only visualizing! 
"""
DoDeleteWorkingDir=True
# Batch submission command (use "" to launch without batch submission):
""" This will depend on your queueing system., ask your system administrator...

    Examples: LaunchJobCommand=\"bsub -q 1day\"
      or, if you do not use a queueing system: LaunchJobCommand=\"\"
"""
LaunchJobCommand="" 
# {expert} Root directory name for this project:
ProjectDir="/home2/bioinfo/scheres/work/protocols"
# {expert} Directory name for logfiles:
""" All logfiles will be stored in $ProjectDir/$LogDir
"""
LogDir="Logs"
#------------------------------------------------------------------------------------------------
# {section} MLalign2D parameters
#------------------------------------------------------------------------------------------------
# Perform 2D maximum-likelihood refinement?
DoML2D=False
# Number of references (or classes) to be used:
NumberOfReferences=3
# Also include mirror transformation in the alignment?
"""  Including the mirror transformation is useful if your particles have a handedness
     and may fall either face-up or face-down on the grid
"""
DoMirror=False
# Use the fast version of this algorithm?
""" See Scheres et al., Bioinformatics, 21 (Suppl. 2), ii243-ii244
"""
DoFast=True
# {expert} Additional xmipp_MLalign2D parameters:
ExtraParamsMLalign2D=""
#------------------------------------------------------------------------------------------------
# {section} Parallelization issues
#------------------------------------------------------------------------------------------------
# Use multiple processors in parallel? (see Expert options)
DoParallel=False
# Number of processors to use:
MyNumberOfCPUs=10
# {expert} A list of all available CPUs (the MPI-machinefile):
""" Depending on your system, your standard script to launch MPI-jobs may require this
"""
MyMachineFile="/home2/bioinfo/scheres/machines.dat"
# {expert} Standard script to launch MPI-jobs:
""" This will also depend on your system...
    The simplest script consists of the following two lines:

    #!/usr/bin/env sh
    `which mpirun` -np MyNumberOfCPUs -machinefile MyMachineFile ~/machines.dat \

    Note that the next line with the xmipp_mpi_MLalign2D command will be
    generated automatically, and the variables MyNumberOfCPUs and MyMachineFile
    will be replaced by the corresponding values given here above.
    More scripts for different batch systems can be found at:
    [Wiki link]
"""
ParallelScript="/home2/bioinfo/scheres/submit_mpi_job.sh"
#------------------------------------------------------------------------------------------------
# {section} Analysis of results
#------------------------------------------------------------------------------------------------
# Visualize class averages and show logfile? (Perform this only after job completion!)
DoVisualizeML2D=False
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
# {end-of-header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
#
class ML2D_class:

    #init variables
    def __init__(self,
                 InSelFile,
                 WorkingDir,
                 DoDeleteWorkingDir,
                 ProjectDir,
                 LogDir,
                 DoML2D,
                 NumberOfReferences,
                 DoMirror,
                 DoFast,
                 DoVisualizeML2D,
                 ExtraParamsMLalign2D,
                 DoParallel,
                 MyNumberOfCPUs,
                 MyMachineFile,
                 LaunchJobCommand,
                 ParallelScript):
	     
        import os,sys,shutil
        scriptdir=os.path.expanduser('~')+'/scripts/'
        sys.path.append(scriptdir) # add default search path
        import log

        self.WorkingDir=WorkingDir
        self.ProjectDir=ProjectDir
        self.InSelFile=self.ProjectDir+'/'+str(InSelFile)
        self.NumberOfReferences=NumberOfReferences
        self.DoMirror=DoMirror
        self.DoFast=DoFast
        self.ExtraParamsMLalign2D=ExtraParamsMLalign2D
        self.DoParallel=DoParallel
        self.MyNumberOfCPUs=MyNumberOfCPUs
        self.MyMachineFile=MyMachineFile
        self.LaunchJobCommand=LaunchJobCommand
        self.ParallelScript=ParallelScript

        # Setup logging
        self.log=log.init_log_system(self.ProjectDir,
                                     LogDir,
                                     sys.argv[0],
                                     self.WorkingDir)
                
        # Delete working directory if it exists, make a new one, and go there
        if (DoDeleteWorkingDir and DoML2D): 
            if os.path.exists(self.WorkingDir):
                shutil.rmtree(self.WorkingDir)
        if not os.path.exists(self.WorkingDir):
            os.makedirs(self.WorkingDir)

        # Execute MLalign2D in the working directory
        os.chdir(self.WorkingDir)
        if (DoML2D):
            self.execute_MLalign2D()

        if (DoVisualizeML2D):
            self.visualize_ML2D()

        # Return to parent dir
        os.chdir(os.pardir)


    def execute_MLalign2D(self):
        import os
        import launch_parallel_job
        print '*********************************************************************'
        print '*  Executing MLalign2D program :' 
        params= ' -i '    + str(self.InSelFile) + \
                ' -o '    + str(self.WorkingDir) + \
                ' -nref ' + str(self.NumberOfReferences)
        if (self.DoFast):
            params+= ' -fast '
        if (self.DoMirror):
            params+= ' -mirror '
        # By default, do not write offsets for 2D alignments...
        # This will affect -restart, but that falls outside the scope of the protocol anyway
        params+=' '+self.ExtraParamsMLalign2D
        params+=' -dont_write_offsets'

        launch_parallel_job.launch_job(self.DoParallel,
                                       "xmipp_MLalign2D",
                                       "xmipp_mpi_MLalign2D",
                                       params,
                                       self.ParallelScript,
                                       self.LaunchJobCommand,
                                       self.log,
                                       self.MyNumberOfCPUs,
                                       self.MyMachineFile,
                                       self.WorkingDir,
                                       False)

    def visualize_ML2D(self):
        # Visualize class averages:
        import glob
        import os
        selfiles=glob.glob(str(self.WorkingDir)+'_it?????.sel')
        if len(selfiles)==0:
            print "No selfiles yet. Visualize after job completion..."
        else:
            lastselfile=selfiles[-1]
            command='xmipp_show -sel '+str(lastselfile)+ ' & '
            print '* ',command
            self.log.info(command)
            os.system(command)
            # print logfile to screen:
            logfiles=glob.glob('ml2d_it?????.log')
            lastlogfile=logfiles[-1]
            fh=open(lastlogfile,'r')
            loglines=fh.readlines()
            print "Logfile "+str(lastlogfile)+": "
            for line in loglines:
                print line[:-1]

    def close(self):
        message='Done!'
        print '*',message
        print '*********************************************************************'

#		
# Main
#     
if __name__ == '__main__':

    # create ML2D_class object

    ML2D=ML2D_class(InSelFile,
                    WorkingDir,
                    DoDeleteWorkingDir,
                    ProjectDir,
                    LogDir,
                    DoML2D,
                    NumberOfReferences,
                    DoMirror,
                    DoFast,
                    DoVisualizeML2D,
                    ExtraParamsMLalign2D,
                    DoParallel,
                    MyNumberOfCPUs,
                    MyMachineFile,
                    LaunchJobCommand,
                    ParallelScript)
    # close 
    ML2D.close()

