/***************************************************************************
 *
 * Authors:     Jose Roman Bilbao (jrbcast@ace.ual.es)
 *    Roberto Marabini (roberto@cnb.csic.es)
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
#include "mpi_reconstruct_fourier.h"

/** Empty constructor */
//ProgMPIRecFourier::ProgMPIRecFourier()
//{
//    node = NULL;
//}

/*  constructor ------------------------------------------------------- */
ProgMPIRecFourier::ProgMPIRecFourier(int argc, char *argv[])
{
    this->read(argc, argv);
}

/* constructor providing an MpiNode
 * this is usefull for using this programs from others
 */
ProgMPIRecFourier::ProgMPIRecFourier(MpiNode * node)
{
    this->setNode(node);
}

/* Special way of reading to sync all nodes */
void ProgMPIRecFourier::read(int argc, char** argv)
{
    XmippMpiProgram::read(argc, argv);
    ProgRecFourier::read(argc, argv);
}

/* Usage ------------------------------------------------------------------- */
void ProgMPIRecFourier::defineParams()
{
    ProgRecFourier::defineParams();
    addParamsLine("  [--mpi_job_size <size=10>]    : Number of images sent to a cpu in a single job ");
    addParamsLine("                                : 10 may be a good value");
    addParamsLine("                                : if  -1 the computer will put the maximum");
    addParamsLine("                                : posible value that may not be the best option");
}

/* Read parameters --------------------------------------------------------- */
void ProgMPIRecFourier::readParams()
{
    ProgRecFourier::readParams();
    mpi_job_size=getIntParam("--mpi_job_size");
}

/* Pre Run PreRun for all nodes but not for all works */
void ProgMPIRecFourier::preRun()
{
    if (nProcs < 2)
        REPORT_ERROR(ERR_ARG_INCORRECT,"This program cannot be executed in a single working node");

    if (node->isMaster())
    {
        show();
        SF.read(fn_sel);

        //Send verbose level to node 1
        MPI_Send(&verbose, 1, MPI_INT, 1, TAG_SETVERBOSE, *node->comm);
    }
    else
    {
        produceSideinfo();
        SF.firstObject();
    }

    //leer sel file / dividir por mpi_job_size
    numberOfJobs=ceil((double)SF.size()/mpi_job_size);

    //only one node will write in the console
    if (node->rank == 1 )
    {
        // Get verbose status
        MPI_Recv(&verbose, 1, MPI_INT, 0, TAG_SETVERBOSE, *node->comm, &status);

        //use threads for volume inverse fourier transform, plan is created in setReal()
        //only rank=1 makes inverse Fourier trnasform
        transformerVol.setThreadsNumber(numThreads);

        //#define DEBUG
#ifdef DEBUG

        std::cerr << "SF.ImgNo() mpi_job_size "
        << SF.ImgNo() << " "
        << mpi_job_size
        << std::endl;
        std::cerr << "numberOfJobs: " << numberOfJobs << std::endl <<std::endl;
#endif
#undef DEBUGDOUBLE

    }
}
/* Run --------------------------------------------------------------------- */
void ProgMPIRecFourier::run()
{
    preRun();
    struct timeval start_time, end_time;
    long int total_usecs;
    double total_time;

    // Real workers, rank=0 is the master, does not work
    nProcs = nProcs - 1;

    if( numberOfJobs < nProcs )
    {
        if( node->isMaster() )
        {
            std::cerr << "\nReducing the number of MPI workers from " <<
            nProcs << " to " <<
            numberOfJobs << std::endl;

            std::cerr << std::flush;
        }

        nProcs = numberOfJobs;

        // Unused nodes are removed from the MPI communicator
        node->active = (node->rank <= numberOfJobs);
        node->updateComm();
    }

    if (node->isMaster())
    {
        gettimeofday(&start_time,NULL);


        if ( verbose )
            init_progress_bar(numberOfJobs);

        int FSC=numberOfJobs/2;

        for (int i=0;i<numberOfJobs;i++)
        {

            //#define DEBUG
#ifdef DEBUG
            std::cerr << "master-recv  i=" << i << std::endl;
            std::cerr << "numberOfJobs: " << numberOfJobs << std::endl <<std::endl;
#endif
#undef DEBUG

            MPI_Recv(0, 0, MPI_INT, MPI_ANY_SOURCE, TAG_FREEWORKER,
                     *node->comm, &status);

            if ( status.MPI_TAG != TAG_FREEWORKER )
                REPORT_ERROR(ERR_ARG_INCORRECT,"Unexpected TAG, please contact developers");

            //#define DEBUG
#ifdef DEBUG

            std::cerr << "master-send i=" << i << std::endl;
#endif
#undef DEBUG

            MPI_Send(&i,
                     1,
                     MPI_INT,
                     status.MPI_SOURCE,
                     TAG_WORKFORWORKER,
                     *node->comm);

            if( i == FSC && fn_fsc != "" )
            {
                // sending every worker COLLECT_FOR_FSC
                for ( int worker = 1 ; worker <= nProcs ; worker ++ )
                {
                    MPI_Recv(0,
                             0,
                             MPI_INT,
                             MPI_ANY_SOURCE,
                             TAG_FREEWORKER,
                             *node->comm,
                             &status);

                    MPI_Send( 0,
                              0,
                              MPI_INT,
                              status.MPI_SOURCE,
                              TAG_COLLECT_FOR_FSC,
                              *node->comm);
                }
            }

            if (verbose)
                progress_bar(i);
        }

        // Wait for all processes to finish processing current jobs
        // so time statistics are correct
        for ( int i = 1 ; i <= nProcs ; i ++ )
        {
            MPI_Recv(0,
                     0,
                     MPI_INT,
                     MPI_ANY_SOURCE,
                     TAG_FREEWORKER,
                     *node->comm,
                     &status);
        }

        gettimeofday(&end_time,NULL);

        total_usecs = (end_time.tv_sec-start_time.tv_sec) * 1000000 + (end_time.tv_usec-start_time.tv_usec);
        total_time=(double)total_usecs/(double)1000000;

        if (verbose > 0)
            std::cout << "\n\nProcessing time: " << total_time << " secs." << std::endl;

        // Start collecting results
        for ( int i = 1 ; i <= nProcs ; i ++ )
        {
            MPI_Send(0,
                     0,
                     MPI_INT,
                     i,
                     TAG_TRANSFER,
                     *node->comm );
        }

    }
    else if( node->active )
    {
        // Select only relevant part of selfile for this rank
        // job number
        // job size
        // aux variable
        double * fourierVolume = (double *)VoutFourier.data;
        double * fourierWeights = FourierWeights.data;

        sizeout = MULTIDIM_SIZE(FourierWeights);

        barrier_init( &barrier, numThreads+1);
        pthread_mutex_init( &workLoadMutex, NULL );
        statusArray = NULL;
        th_ids = (pthread_t *)malloc(numThreads * sizeof(pthread_t));
        th_args = (ImageThreadParams *)malloc(numThreads * sizeof(ImageThreadParams));

        for ( int nt = 0 ; nt < numThreads ; nt++ )
        {
            th_args[nt].parent=this;
            th_args[nt].myThreadID = nt;
            th_args[nt].selFile = new MetaData(SF);
            pthread_create((th_ids+nt),NULL,processImageThread,(void*)(th_args+nt));
        }

        while (1)
        {
            int jobNumber;

            //#define DEBUG
#ifdef DEBUG

            std::cerr << "slave-send TAG_FREEWORKER rank=" << node->rank << std::endl;
#endif
     #undef DEBUG
            //I am free
            MPI_Send(0, 0, MPI_INT, 0, TAG_FREEWORKER, *node->comm);
            MPI_Probe(0, MPI_ANY_TAG, *node->comm, &status);

            if (status.MPI_TAG == TAG_COLLECT_FOR_FSC)
            {
                //If I  do not read this tag
                //master will no further process
                //a posibility is a non-blocking send
                MPI_Recv(0, 0, MPI_INT, 0, TAG_COLLECT_FOR_FSC, *node->comm, &status);

                if( node->rank == 1 )
                {
                    // Reserve memory for the receive buffer
                    double * recBuffer = (double *) malloc (sizeof(double)*BUFFSIZE);
                    int receivedSize;
                    double * pointer;
                    pointer = fourierVolume;
                    int currentSource;

                    if ( nProcs > 2 )
                    {
                        // Receive from other workers
                        for ( int i = 2 ; i <= nProcs ; i++)
                        {
                            MPI_Recv(0,0, MPI_INT, MPI_ANY_SOURCE, TAG_FREEWORKER,
                                     *node->comm, &status);

                            currentSource = status.MPI_SOURCE;

                            pointer = fourierVolume;

                            while (1)
                            {
                                MPI_Probe( currentSource, MPI_ANY_TAG, *node->comm, &status );

                                if ( status.MPI_TAG == TAG_FREEWORKER )
                                {
                                    MPI_Recv(0,0, MPI_INT, currentSource, TAG_FREEWORKER, *node->comm, &status );
                                    break;
                                }

                                MPI_Recv( recBuffer,
                                          BUFFSIZE,
                                          MPI_DOUBLE,
                                          currentSource,
                                          MPI_ANY_TAG,
                                          *node->comm,
                                          &status );

                                MPI_Get_count( &status, MPI_DOUBLE, &receivedSize );

                                for ( int i = 0 ; i < receivedSize ; i ++ )
                                {
                                    pointer[i] += recBuffer[i];
                                }

                                pointer += receivedSize;
                            }

                            pointer = fourierWeights;

                            while (1)
                            {
                                MPI_Probe( currentSource, MPI_ANY_TAG, *node->comm, &status );

                                if ( status.MPI_TAG == TAG_FREEWORKER )
                                {
                                    MPI_Recv( 0,0,MPI_INT,currentSource,TAG_FREEWORKER, *node->comm,&status );

                                    break;
                                }

                                MPI_Recv( recBuffer,
                                          BUFFSIZE,
                                          MPI_DOUBLE,
                                          currentSource,
                                          MPI_ANY_TAG,
                                          *node->comm,
                                          &status );

                                MPI_Get_count( &status, MPI_DOUBLE, &receivedSize );

                                for ( int i = 0 ; i < receivedSize ; i ++ )
                                {
                                    pointer[i] += recBuffer[i];
                                }

                                pointer += receivedSize;
                            }
                        }
                    }
                    free( recBuffer );

                    Image<double> auxVolume1;
                    auxVolume1().alias( FourierWeights );
                    auxVolume1.write((std::string)fn_fsc + "_1_Weights.vol");

                    Image< std::complex<double> > auxFourierVolume1;
                    auxFourierVolume1().alias( VoutFourier );
                    auxFourierVolume1.write((std::string) fn_fsc + "_1_Fourier.vol");


                    // Normalize global volume and store data
                    finishComputations(FileName((std::string) fn_fsc + "_split_1.vol"));

                    Vout().initZeros(volPadSizeZ, volPadSizeY, volPadSizeX);

                    transformerVol.setReal(Vout());
                    Vout().clear();
                    transformerVol.getFourierAlias(VoutFourier);
                    FourierWeights.initZeros(VoutFourier);
                    VoutFourier.initZeros();
                }
                else
                {
                    MPI_Send( 0,0,MPI_INT,1,TAG_FREEWORKER, *node->comm );

                    sendDataInChunks( fourierVolume, 1, 2*sizeout, BUFFSIZE, *node->comm );

                    MPI_Send( 0,0,MPI_INT,1,TAG_FREEWORKER, *node->comm );

                    sendDataInChunks( fourierWeights, 1, sizeout, BUFFSIZE, *node->comm);

                    MPI_Send( 0,0,MPI_INT,1,TAG_FREEWORKER,*node->comm);

                    Vout().initZeros(volPadSizeZ, volPadSizeY, volPadSizeX);
                    transformerVol.setReal(Vout());
                    Vout().clear();
                    transformerVol.getFourierAlias(VoutFourier);
                    FourierWeights.initZeros(VoutFourier);
                    VoutFourier.initZeros();
                }
            }
            else if (status.MPI_TAG == TAG_TRANSFER)
            {
                //If I  do not read this tag
                //master will no further process
                MPI_Recv(0, 0, MPI_INT, 0, TAG_TRANSFER, *node->comm, &status);
#ifdef DEBUG

                std::cerr << "Wr" << node->rank << " " << "TAG_STOP" << std::endl;
#endif

                if ( node->rank == 1 )
                {
                    // Reserve memory for the receive buffer
                    double * recBuffer = (double *) malloc (sizeof(double)*BUFFSIZE);
                    int receivedSize;
                    double * pointer;
                    pointer = fourierVolume;
                    int currentSource;

                    gettimeofday(&start_time,NULL);

                    if ( nProcs > 1 )
                    {
                        // Receive from other workers

                        for ( int i = 0 ; i <= (nProcs-2) ; i++)
                        {
                            MPI_Recv(0,0, MPI_INT, MPI_ANY_SOURCE, TAG_FREEWORKER,
                                     *node->comm, &status);

                            currentSource = status.MPI_SOURCE;

                            pointer = fourierVolume;

                            while (1)
                            {
                                MPI_Probe( currentSource, MPI_ANY_TAG, *node->comm, &status );

                                if ( status.MPI_TAG == TAG_FREEWORKER )
                                {
                                    MPI_Recv(0,0, MPI_INT, currentSource, TAG_FREEWORKER, *node->comm, &status );

                                    break;
                                }

                                MPI_Recv( recBuffer,
                                          BUFFSIZE,
                                          MPI_DOUBLE,
                                          currentSource,
                                          MPI_ANY_TAG,
                                          *node->comm,
                                          &status );

                                MPI_Get_count( &status, MPI_DOUBLE, &receivedSize );

                                for ( int i = 0 ; i < receivedSize ; i ++ )
                                {
                                    pointer[i] += recBuffer[i];
                                }

                                pointer += receivedSize;
                            }

                            pointer = fourierWeights;

                            while (1)
                            {
                                MPI_Probe( currentSource, MPI_ANY_TAG, *node->comm, &status );

                                if ( status.MPI_TAG == TAG_FREEWORKER )
                                {
                                    MPI_Recv( 0,0,MPI_INT,currentSource,TAG_FREEWORKER, *node->comm,&status );

                                    break;
                                }

                                MPI_Recv( recBuffer,
                                          BUFFSIZE,
                                          MPI_DOUBLE,
                                          currentSource,
                                          MPI_ANY_TAG,
                                          *node->comm,
                                          &status );

                                MPI_Get_count( &status, MPI_DOUBLE, &receivedSize );

                                for ( int i = 0 ; i < receivedSize ; i ++ )
                                {
                                    pointer[i] += recBuffer[i];
                                }

                                pointer += receivedSize;
                            }
                        }
                    }

                    free( recBuffer );
                    gettimeofday(&end_time,NULL);

                    total_usecs = (end_time.tv_sec-start_time.tv_sec) * 1000000 + (end_time.tv_usec-start_time.tv_usec);
                    total_time=(double)total_usecs/(double)1000000;
                    if (verbose > 0)
                        std::cout << "Transfers time: " << total_time << " secs." << std::endl;

                    if( fn_fsc != "" )
                    {

                        Image<double> auxVolume2;
                        auxVolume2().alias( FourierWeights );
                        auxVolume2.write((std::string)fn_fsc + "_2_Weights.vol");

                        Image< std::complex<double> > auxFourierVolume2;
                        auxFourierVolume2().alias( VoutFourier );
                        auxFourierVolume2.write((std::string) fn_fsc + "_2_Fourier.vol");


                        // Normalize global volume and store data
                        finishComputations(FileName((std::string) fn_fsc + "_split_2.vol"));

                        Vout().initZeros(volPadSizeZ, volPadSizeY, volPadSizeX);
                        transformerVol.setReal(Vout());
                        Vout().clear();
                        transformerVol.getFourierAlias(VoutFourier);
                        FourierWeights.initZeros(VoutFourier);
                        VoutFourier.initZeros();

                        //int x,y,z;

                        //FourierWeights.getDimension(y,x,z);
                        gettimeofday(&start_time,NULL);

                        auxVolume2.sumWithFile((std::string) fn_fsc + "_1_Weights.vol");

                        gettimeofday(&end_time,NULL);
                        total_usecs = (end_time.tv_sec-start_time.tv_sec) * 1000000 + (end_time.tv_usec-start_time.tv_usec);
                        total_time=(double)total_usecs/(double)1000000;
                        if (verbose > 0)
                            std::cout << "SumFile1: " << total_time << " secs." << std::endl;

                        auxVolume2.sumWithFile((std::string) fn_fsc + "_2_Weights.vol");

                        //VoutFourier.getDimension(y,x,z);
                        auxFourierVolume2.sumWithFile((std::string) fn_fsc + "_1_Fourier.vol");
                        auxFourierVolume2.sumWithFile((std::string) fn_fsc + "_2_Fourier.vol");

                        //remove temporary files
                        remove(((std::string) fn_fsc + "_1_Weights.vol").c_str());
                        remove(((std::string) fn_fsc + "_2_Weights.vol").c_str());
                        remove(((std::string) fn_fsc + "_1_Fourier.vol").c_str());
                        remove(((std::string) fn_fsc + "_2_Fourier.vol").c_str());
                        gettimeofday(&end_time,NULL);
                        total_usecs = (end_time.tv_sec-start_time.tv_sec) * 1000000 + (end_time.tv_usec-start_time.tv_usec);
                        total_time=(double)total_usecs/(double)1000000;
                        if (verbose > 0)
                            std::cout << "SumFile: " << total_time << " secs." << std::endl;

                        /*Save SUM
                                                    //this is an image but not an xmipp image
                                                    auxFourierVolume.write((std::string)fn_fsc + "_all_Fourier.vol",
                                                            false,VDOUBLE);
                                                    auxVolume.write((std::string)fn_fsc + "_all_Weights.vol",
                                                            false,VDOUBLE);
                        */
                    }

                    // Normalize global volume and store data
                    gettimeofday(&start_time,NULL);
                    finishComputations(fn_out);
                    gettimeofday(&end_time,NULL);
                    int i=0;

                    total_usecs = (end_time.tv_sec-start_time.tv_sec) * 1000000 + (end_time.tv_usec-start_time.tv_usec);
                    total_time=(double)total_usecs/(double)1000000;
                    if (verbose > 0)
                    {
                        std::cout << "Weighting time: " << total_time << " secs." << std::endl;
                        std::cout << "Execution completed successfully\n" << std::endl;
                    }
                    break;
                }
                else
                {
                    MPI_Send( 0,0,MPI_INT,1,TAG_FREEWORKER, *node->comm );

                    sendDataInChunks( fourierVolume, 1, 2 * sizeout, BUFFSIZE, *node->comm);

                    MPI_Send( 0,0,MPI_INT,1,TAG_FREEWORKER, *node->comm );

                    sendDataInChunks( fourierWeights, 1, sizeout, BUFFSIZE, *node->comm);

                    MPI_Send( 0,0,MPI_INT,1,TAG_FREEWORKER, *node->comm);

                    break;
                }
            }
            else if (status.MPI_TAG == TAG_WORKFORWORKER)
            {
                //get the job number
                MPI_Recv(&jobNumber, 1, MPI_INT, 0, TAG_WORKFORWORKER, *node->comm, &status);
                //LABEL
                //(if jobNumber == -1) break;
                threadOpCode=PROCESS_IMAGE;

                int min_i, max_i;

                min_i = jobNumber*mpi_job_size;
                max_i = min_i + mpi_job_size - 1;

                if ( max_i >= SF.size())
                    max_i  = SF.size()-1;

                processImages( min_i, max_i );
            }
            else
            {
                std::cerr << "3) Received unknown TAG I quit" << std::endl;
                exit(0);
            }
        }
    }

    // Kill threads used on workers
    if ( node->active && !node->isMaster() )
    {
        threadOpCode = EXIT_THREAD;
        barrier_wait( &barrier );

        for ( int nt=0; nt<numThreads; nt++)
        {
            pthread_join(*(th_ids+nt),NULL);
        }
        barrier_destroy( &barrier );
    }
}

int  ProgMPIRecFourier::sendDataInChunks( double * pointer, int dest, int totalSize, int buffSize, MPI_Comm comm )
{
    double * localPointer = pointer;

    int numChunks = ceil((double)totalSize/(double)buffSize);
    int packetSize;
    int err=0;

    for ( int i = 0 ; i < numChunks ; i ++ )
    {
        if ( i == (numChunks-1))
            packetSize = totalSize-i*buffSize;
        else
            packetSize = buffSize;

        if ( (err = MPI_Send( localPointer, packetSize,
                              MPI_DOUBLE, dest, 0, comm ))
             != MPI_SUCCESS )
        {
            break;
        }

        localPointer += packetSize;
    }

    return err;
}





