/***************************************************************************
 *
 * Authors: Sjors Scheres (scheres@cnb.csic.es)
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

#include <data/image.h>
#include <data/fftw.h>
#include <data/args.h>

#include <vector>

class sort_junk_parameters
{

public:
    double cutoff;
    double avg_mean, avg_sig, avg_min, avg_max;
    double sig_mean, sig_sig, sig_min, sig_max;
    double avg_nlowpix, avg_nhighpix, sig_nlowpix, sig_nhighpix;
    double avg_nradlow, avg_nradhigh, sig_nradlow, sig_nradhigh;
    double sig_sigquad, avg_sigquad, sig_meanquad, avg_meanquad;
    double avg_ampl2_1, avg_ampl2_2, avg_ampl2_3, avg_ampl2_4;
    double sig_ampl2_1, sig_ampl2_2, sig_ampl2_3, sig_ampl2_4;
    std::vector<FileName> names;
    Matrix1D<double> zscore;
    std::vector<std::vector<double> > values;
    FileName fn_out;

public:
    sort_junk_parameters()
    {}

    void usage()
    {
        std::cout  << " A sorting program for identifying junk particles \n"
        << " Parameters:\n"
        << " -i <selfile>            : Selfile with (normalized) input images\n"
        << " [-o <root=\"sort_junk\">] : Output rootname\n"
        << " [-train <selfile>]      : Train on selfile with good particles\n"
        << " [-zcut <float=-1>]       : Cut-off for Z-scores (negative for no cut-off) \n"
        ;
    }

    void process_selfile(SelFile &SF, bool do_means, bool do_values)
    {

        FileName fn;
        ImageXmipp img;
        int imgno, dim, nr_imgs;
        double mean, stddev, minval, maxval;
        double nlowpix, nhighpix, nradlow, nradhigh;
        double sum_quadsig, sum2_quadsig, sum_quadmean, sum2_quadmean;
        double ampl2_1, ampl2_2, ampl2_3, ampl2_4;
        std::vector<double> dum;
        Matrix2D<double> Mrad, ampl2;
        Matrix2D<std::complex<double> > IMG;
        Matrix1D<int> radial_count, center(2);
        Matrix1D<double> rmean_ampl2;
        center.initZeros();

        if (do_means && !do_values)
            std::cerr << " Processing training set ..." << std::endl;
        else
            std::cerr << " Sorting particle set ..." << std::endl;

        SF.ImgSize(dim, dim);
        Mrad.resize(dim, dim);
        Mrad.setXmippOrigin();
        FOR_ALL_ELEMENTS_IN_MATRIX2D(Mrad)
        {
            Mrad(i, j) = sqrt((double)(i * i + j * j));
        }
        nr_imgs = SF.ImgNo();
        init_progress_bar(nr_imgs);
        int c = XMIPP_MAX(1, nr_imgs / 60);

        if (do_means)
        {
            avg_mean = sig_mean = 0.;
            avg_sig = sig_sig = 0.;
            avg_min = sig_min = 0.;
            avg_max = sig_max = 0.;
            avg_sigquad = sig_sigquad = 0.;
            avg_meanquad = sig_meanquad = 0.;
            avg_nlowpix = sig_nlowpix = 0.;
            avg_nhighpix = sig_nhighpix = 0.;
            avg_nradlow = sig_nradlow = 0.;
            avg_nradhigh = sig_nradhigh = 0.;
            avg_ampl2_1 = sig_ampl2_1 = 0.;
            avg_ampl2_2 = sig_ampl2_2 = 0.;
            avg_ampl2_3 = sig_ampl2_3 = 0.;
            avg_ampl2_4 = sig_ampl2_4 = 0.;
        }

        imgno = 0;
        while (!SF.eof())
        {
            fn = SF.NextImg();
            if (fn=="") break;
            img.read(fn);
            img().setXmippOrigin();

            // Overall statistics
            img().computeStats(mean, stddev, minval, maxval);

            // Number of low or high-valued pixels
            nhighpix = nlowpix = 0.;
            nradhigh = nradlow = 0.;
            FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(img())
            {
                if (dMij(img(), i, j) > stddev) nhighpix += 1.;
                if (dMij(img(), i, j) < -stddev) nlowpix += 1.;
                if (dMij(img(), i, j) > stddev) nradhigh += dMij(Mrad, i, j);
                if (dMij(img(), i, j) < -stddev) nradlow += dMij(Mrad, i, j);
            }

            // Quadrant statistics
            quadrant_stats(img, sum_quadsig, sum2_quadsig, sum_quadmean, sum2_quadmean);

            // Fourier-space stats
            FourierTransform(img(), IMG);
            FFT_magnitude(IMG, ampl2);
            CenterFFT(ampl2, true);
            ampl2 *= ampl2;
            rmean_ampl2.initZeros();
            radialAverage(ampl2, center, rmean_ampl2, radial_count);
            ampl2_1 = dVi(rmean_ampl2, 1);
            ampl2_2 = dVi(rmean_ampl2, 2);
            ampl2_3 = dVi(rmean_ampl2, 3);
            ampl2_4 = dVi(rmean_ampl2, 4);

            if (do_means)
            {
                avg_mean += mean;
                sig_mean += mean * mean;
                avg_sig += stddev;
                sig_sig += stddev * stddev;
                avg_min += minval;
                sig_min += minval * minval;
                avg_max += maxval;
                sig_max += maxval * maxval;
                avg_nhighpix += nhighpix;
                sig_nhighpix += nhighpix * nhighpix;
                avg_nlowpix += nlowpix;
                sig_nlowpix += nlowpix * nlowpix;
                avg_nradhigh += nradhigh;
                sig_nradhigh += nradhigh * nradhigh;
                avg_nradlow += nradlow;
                sig_nradlow += nradlow * nradlow;
                avg_meanquad += sum2_quadmean;
                sig_meanquad += sum2_quadmean * sum2_quadmean;
                avg_sigquad += sum2_quadsig;
                sig_sigquad += sum2_quadsig * sum2_quadsig;
                avg_ampl2_1 += ampl2_1;
                sig_ampl2_1 += ampl2_1 * ampl2_1;
                avg_ampl2_2 += ampl2_2;
                sig_ampl2_2 += ampl2_2 * ampl2_2;
                avg_ampl2_3 += ampl2_3;
                sig_ampl2_3 += ampl2_3 * ampl2_3;
                avg_ampl2_4 += ampl2_4;
                sig_ampl2_4 += ampl2_4 * ampl2_4;
            }
            if (do_values)
            {
                values.push_back(dum);
                names.push_back(fn);
                values[imgno].push_back(mean);
                values[imgno].push_back(stddev);
                values[imgno].push_back(minval);
                values[imgno].push_back(maxval);
                values[imgno].push_back(nhighpix);
                values[imgno].push_back(nlowpix);
                values[imgno].push_back(nradhigh);
                values[imgno].push_back(nradlow);
                values[imgno].push_back(sum2_quadsig);
                values[imgno].push_back(sum2_quadmean);
                values[imgno].push_back(ampl2_1);
                values[imgno].push_back(ampl2_2);
                values[imgno].push_back(ampl2_3);
                values[imgno].push_back(ampl2_4);
            }

            if (imgno % c == 0) progress_bar(imgno);
            imgno++;

        }
        progress_bar(nr_imgs);

        if (do_means)
        {
            // Finish  average and standard deviation calculations
            avg_mean /= (double)nr_imgs;
            avg_sig /= (double)nr_imgs;
            avg_min /= (double)nr_imgs;
            avg_max /= (double)nr_imgs;
            avg_nhighpix /= (double)nr_imgs;
            avg_nlowpix /= (double)nr_imgs;
            avg_sigquad /= (double)nr_imgs;
            avg_meanquad /= (double)nr_imgs;
            avg_nradhigh /= (double)nr_imgs;
            avg_nradlow /= (double)nr_imgs;
            avg_ampl2_1 /= (double)nr_imgs;
            avg_ampl2_2 /= (double)nr_imgs;
            avg_ampl2_3 /= (double)nr_imgs;
            avg_ampl2_4 /= (double)nr_imgs;
            sig_mean = sqrt(sig_mean / (double)nr_imgs - avg_mean * avg_mean);
            sig_sig = sqrt(sig_sig / (double)nr_imgs - avg_sig * avg_sig);
            sig_min = sqrt(sig_min / (double)nr_imgs - avg_min * avg_min);
            sig_max = sqrt(sig_max / (double)nr_imgs - avg_max * avg_max);
            sig_nhighpix = sqrt(sig_nhighpix / (double)nr_imgs - avg_nhighpix * avg_nhighpix);
            sig_nlowpix = sqrt(sig_nlowpix / (double)nr_imgs - avg_nlowpix * avg_nlowpix);
            sig_sigquad = sqrt(sig_sigquad / (double)nr_imgs - avg_sigquad * avg_sigquad);
            sig_meanquad = sqrt(sig_meanquad / (double)nr_imgs - avg_meanquad * avg_meanquad);
            sig_nradhigh = sqrt(sig_nradhigh / (double)nr_imgs - avg_nradhigh * avg_nradhigh);
            sig_nradlow = sqrt(sig_nradlow / (double)nr_imgs - avg_nradlow * avg_nradlow);
            sig_ampl2_1 = sqrt(sig_ampl2_1 / (double)nr_imgs - avg_ampl2_1 * avg_ampl2_1);
            sig_ampl2_2 = sqrt(sig_ampl2_2 / (double)nr_imgs - avg_ampl2_2 * avg_ampl2_2);
            sig_ampl2_3 = sqrt(sig_ampl2_3 / (double)nr_imgs - avg_ampl2_3 * avg_ampl2_3);
            sig_ampl2_4 = sqrt(sig_ampl2_4 / (double)nr_imgs - avg_ampl2_4 * avg_ampl2_4);
        }

        if (do_values)
        {
            zscore.resize(nr_imgs);
            for (imgno = 0; imgno < nr_imgs; imgno++)
            {
		if (sig_mean > 0.)
		    values[imgno][0] = ABS(avg_mean - values[imgno][0]) / sig_mean;
		else
		    values[imgno][0] = 0.;
		if (sig_sig > 0.)
		    values[imgno][1] = ABS(values[imgno][1] - avg_sig) / sig_sig;
		else
		    values[imgno][1] = 0.;
		if (sig_min > 0.)
		    values[imgno][2] = ABS(values[imgno][2] - avg_min) / sig_min;
		else
		    values[imgno][2] = 0.;
		if (sig_max > 0.)
		    values[imgno][3] = ABS(values[imgno][3] - avg_max) / sig_max;
		else
		    values[imgno][3] = 0.;
		if (sig_nhighpix > 0.)
		    values[imgno][4] = ABS(values[imgno][4] - avg_nhighpix) / sig_nhighpix;
		else
		    values[imgno][4] = 0.;
		if (sig_nlowpix > 0.)
		    values[imgno][5] = ABS(values[imgno][5] - avg_nlowpix) / sig_nlowpix;
		else
		    values[imgno][5] = 0.;
		if (sig_nradhigh > 0.)
		    values[imgno][6] = ABS(values[imgno][6] - avg_nradhigh) / sig_nradhigh;
		else
		    values[imgno][6] = 0.;
		if (sig_nradlow > 0.)
		    values[imgno][7] = ABS(values[imgno][7] - avg_nradlow) / sig_nradlow;
		else
		    values[imgno][7] = 0.;
		if (sig_sigquad > 0.)
		    values[imgno][8] = ABS(values[imgno][8] - avg_sigquad) / sig_sigquad;
		else
		    values[imgno][8] = 0.;
		if (sig_meanquad > 0.)
		    values[imgno][9] = ABS(values[imgno][9] - avg_meanquad) / sig_meanquad;
		else
		    values[imgno][9] = 0.;
		if (sig_ampl2_1 > 0.)
		    values[imgno][10] = ABS(values[imgno][10] - avg_ampl2_1) / sig_ampl2_1;
 		else
		    values[imgno][10] = 0.;
		if (sig_ampl2_2 > 0.)
		    values[imgno][11] = ABS(values[imgno][11] - avg_ampl2_2) / sig_ampl2_2;
		else
		    values[imgno][11] = 0.;
		if (sig_ampl2_3 > 0.)
		    values[imgno][12] = ABS(values[imgno][12] - avg_ampl2_3) / sig_ampl2_3;
		else
		    values[imgno][12] = 0.;
		if (sig_ampl2_4 > 0.)
		    values[imgno][13] = ABS(values[imgno][13] - avg_ampl2_4) / sig_ampl2_4;
		else
		    values[imgno][13] = 0.;

                if (cutoff > 0.)
                {
                    if (values[imgno][0] < cutoff) values[imgno][0] = 0.;
                    if (values[imgno][1] < cutoff) values[imgno][1] = 0.;
                    if (values[imgno][2] < cutoff) values[imgno][2] = 0.;
                    if (values[imgno][3] < cutoff) values[imgno][3] = 0.;
                    if (values[imgno][4] < cutoff) values[imgno][4] = 0.;
                    if (values[imgno][5] < cutoff) values[imgno][5] = 0.;
                    if (values[imgno][6] < cutoff) values[imgno][6] = 0.;
                    if (values[imgno][7] < cutoff) values[imgno][7] = 0.;
                    if (values[imgno][8] < cutoff) values[imgno][8] = 0.;
                    if (values[imgno][9] < cutoff) values[imgno][9] = 0.;
                    if (values[imgno][10] < cutoff) values[imgno][10] = 0.;
                    if (values[imgno][11] < cutoff) values[imgno][11] = 0.;
                    if (values[imgno][12] < cutoff) values[imgno][12] = 0.;
                    if (values[imgno][13] < cutoff) values[imgno][13] = 0.;
                }
                zscore(imgno) = values[imgno][0] + values[imgno][1] + values[imgno][2] +
                                values[imgno][3] + values[imgno][4] + values[imgno][5] +
                                values[imgno][6] + values[imgno][7] + values[imgno][8] +
                                values[imgno][9] + values[imgno][10] + values[imgno][11] +
                                values[imgno][12] + values[imgno][13];
            }
        }
    }

    void quadrant_stats(ImageXmipp &img,
                        double &sum_quadsig, double &sum2_quadsig,
                        double &sum_quadmean, double &sum2_quadmean)
    {

        double mean, stddev, minval, maxval;
        Matrix1D<int> corner1(2), corner2(2);

        sum_quadsig = sum2_quadsig = 0.;
        sum_quadmean = sum2_quadmean = 0.;
        XX(corner1) = STARTINGX(img());
        YY(corner1) = STARTINGY(img());
        XX(corner2) = 0;
        YY(corner2) = 0;
        img().computeStats(mean, stddev, minval, maxval, corner1, corner2);
        sum2_quadmean += mean * mean;
        sum_quadmean += mean;
        sum2_quadsig += stddev * stddev;
        sum_quadsig += stddev;
        XX(corner1) = STARTINGX(img());
        YY(corner1) = 0;
        XX(corner2) = 0;
        YY(corner2) = FINISHINGY(img());
        img().computeStats(mean, stddev, minval, maxval, corner1, corner2);
        sum2_quadmean += mean * mean;
        sum_quadmean += mean;
        sum2_quadsig += stddev * stddev;
        sum_quadsig += stddev;
        XX(corner1) = 0;
        YY(corner1) = STARTINGY(img());
        XX(corner2) = FINISHINGX(img());
        YY(corner2) = 0;
        img().computeStats(mean, stddev, minval, maxval, corner1, corner2);
        sum2_quadmean += mean * mean;
        sum_quadmean += mean;
        sum2_quadsig += stddev * stddev;
        sum_quadsig += stddev;
        XX(corner1) = 0;
        YY(corner1) = 0;
        XX(corner2) = FINISHINGX(img());
        YY(corner2) = FINISHINGX(img());
        img().computeStats(mean, stddev, minval, maxval, corner1, corner2);
        sum2_quadmean += mean * mean;
        sum_quadmean += mean;
        sum2_quadsig += stddev * stddev;
        sum_quadsig += stddev;

        sum_quadsig /= 4.;
        sum2_quadsig = sqrt(sum2_quadsig / 4. - sum_quadsig * sum_quadsig);
        sum_quadmean /= 4.;
        sum2_quadmean = sqrt(sum2_quadmean / 4. - sum_quadmean * sum_quadmean);

    }

};
/**************************************************************************
        Main
/**************************************************************************/

int main(int argc, char **argv)
{

    SelFile SF, SFout, SFtrain;
    int imgno, dim, nr_imgs, isort;

    Matrix1D<int> sorted;
    FileName fn, fn_train;
    std::ofstream fh_zsum, fh_zind;

    sort_junk_parameters prm;

    try
    {
        fn = getParameter(argc, argv, "-i");
        SF.read(fn);
        nr_imgs = SF.ImgNo();
        prm.fn_out = getParameter(argc, argv, "-o", "sort_junk");
        fn_train = getParameter(argc, argv, "-train", "");
        if (fn_train != "") SFtrain.read(fn_train);
        prm.cutoff = textToFloat(getParameter(argc, argv, "-zcut", "-1"));
    }
    catch (Xmipp_error XE)
    {
        std::cout << XE;
        prm.usage();
        exit(1);
    }

    if (fn_train != "")
    {
        prm.process_selfile(SFtrain, true, false);
        prm.process_selfile(SF, false, true);
    }
    else
    {
        prm.process_selfile(SF, true, true);
    }

    fh_zsum.open((prm.fn_out + ".sumZ").c_str(), std::ios::out);
    fh_zind.open((prm.fn_out + ".indZ").c_str(), std::ios::out);
    sorted = prm.zscore.indexSort();
    SFout.clear();
    fh_zind << "image : avg, stddev, min, max, nhighpix, nlowpix, nradhigh, nradlow, quadsig, quadmean, ampl2_1, ampl2_2, ampl2_3, ampl2_4 " << std::endl;
    for (imgno = 0; imgno < nr_imgs; imgno++)
    {
        isort = sorted(imgno) - 1;
        SFout.insert(prm.names[isort]);
        fh_zsum << prm.zscore(isort)/14. << "   " << prm.names[isort] << std::endl;
        fh_zind << prm.names[isort]
        << " : " << prm.values[isort][0]
        << " " << prm.values[isort][1]
        << " " << prm.values[isort][2]
        << " " << prm.values[isort][3]
        << " " << prm.values[isort][4]
        << " " << prm.values[isort][5]
        << " " << prm.values[isort][6]
        << " " << prm.values[isort][7]
        << " " << prm.values[isort][8]
        << " " << prm.values[isort][9]
        << " " << prm.values[isort][10]
        << " " << prm.values[isort][11]
        << " " << prm.values[isort][12]
        << " " << prm.values[isort][13] << std::endl;
    }
    fh_zsum.close();
    fh_zind.close();
    fn = prm.fn_out + ".sel";
    SFout.write(fn);

}
