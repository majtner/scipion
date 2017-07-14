/***************************************************************************
 *
 * Authors:    Tomas Majtner            tmajtner@cnb.csic.es (2017)
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

#include "classify_fast_2d.h"

#include <data/xmipp_funcs.h>
#include <data/mask.h>
#include <data/filters.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>


// Read arguments ==========================================================
void ProgClassifyFast2D::readParams()
{
    fnSel = getParam("-i");
    fnOut = getParam("-o");
    K = getIntParam("-k");
    fnClusters = getParam("-c");
}

// Show ====================================================================
void ProgClassifyFast2D::show()
{
    if (verbose==0)
        return;
    std::cerr
    << "Input selfile:           " << fnSel      << std::endl
    << "Output selfile:          " << fnOut      << std::endl
    << "Number of clusters:      " << K          << std::endl
    << "Filename with clusters:  " << fnClusters << std::endl
    ;
}

// Usage ===================================================================
void ProgClassifyFast2D::defineParams()
{
    addUsageLine("Clusters a set of images");
    addParamsLine("  -i <selfile>                  : Selfile containing images to be clustered");
    addParamsLine("  -o <selfile>                  : Output selfile");
    addParamsLine("  -k <int>                      : Number of clusters");
    addParamsLine("  [-c <image=\"clusters.xmd\">] : Filename with clusters");
}


class Point
{
private:
	int id_point, id_cluster;
	std::vector<double> values;
	int total_values;

public:
	Point(int id_point, std::vector<double>& values)
	{
		this->id_point = id_point;
		total_values = values.size();

		for(int i = 0; i < total_values; i++)
			this->values.push_back(values[i]);

		id_cluster = -1;
	}

	int getID()
	{
		return id_point;
	}

	void setCluster(int id_cluster)
	{
		this->id_cluster = id_cluster;
	}

	int getCluster()
	{
		return id_cluster;
	}

	double getValue(int index)
	{
		return values[index];
	}

	int getTotalValues()
	{
		return total_values;
	}

	void addValue(double value)
	{
		values.push_back(value);
	}
};


class Cluster
{
private:
	int id_cluster;
	std::vector<double> central_values;
	std::vector<Point> points;

public:
	Cluster(int id_cluster, Point point)
	{
		this->id_cluster = id_cluster;

		int total_values = point.getTotalValues();

		for(int i = 0; i < total_values; i++)
			central_values.push_back(point.getValue(i));

		points.push_back(point);
	}

	void addPoint(Point point)
	{
		points.push_back(point);
	}

	bool removePoint(int id_point)
	{
		int total_points = points.size();

		for(int i = 0; i < total_points; i++)
		{
			if(points[i].getID() == id_point)
			{
				points.erase(points.begin() + i);
				return true;
			}
		}
		return false;
	}

	double getCentralValue(int index)
	{
		return central_values[index];
	}

	void setCentralValue(int index, double value)
	{
		central_values[index] = value;
	}

	Point getPoint(int index)
	{
		return points[index];
	}

	int getTotalPoints()
	{
		return points.size();
	}

	int getID()
	{
		return id_cluster;
	}
};


class KMeans
{
private:
	int K; // number of clusters
	int total_values, total_points, max_iterations;
	std::vector<Cluster> clusters;

	// return ID of nearest center (uses euclidean distance)
	int getIDNearestCenter(Point point)
	{
		double sum = 0.0, min_dist;
		int id_cluster_center = 0;

		for(int i = 0; i < total_values; i++)
		{
			sum += pow(clusters[0].getCentralValue(i) -
					   point.getValue(i), 2.0);
		}

		min_dist = sqrt(sum);

		for(int i = 1; i < K; i++)
		{
			double dist;
			sum = 0.0;

			for(int j = 0; j < total_values; j++)
			{
				sum += pow(clusters[i].getCentralValue(j) -
						   point.getValue(j), 2.0);
			}

			dist = sqrt(sum);

			if(dist < min_dist)
			{
				min_dist = dist;
				id_cluster_center = i;
			}
		}

		return id_cluster_center;
	}

public:
	KMeans(int K, int total_points, int total_values, int max_iterations)
	{
		this->K = K;
		this->total_points = total_points;
		this->total_values = total_values;
		this->max_iterations = max_iterations;
	}

	std::vector<Cluster> run(std::vector<Point> & points, FileName fnClusters)
	{
        // choose K distinct points for the centers of the clusters
        std::vector<int> prohibited_indexes;

        for(int i = 0; i < K; i++)
        {
            while(true)
            {
                int index_point = rand() % total_points;

                if(find(prohibited_indexes.begin(), prohibited_indexes.end(),
                        index_point) == prohibited_indexes.end())
                {
                    prohibited_indexes.push_back(index_point);
                    points[index_point].setCluster(i);
                    Cluster cluster(i, points[index_point]);
                    clusters.push_back(cluster);
                    break;
                }
            }
        }

        std::fstream savedClusters(fnClusters.c_str());
        if (savedClusters.good())
        {
            std::string line;
            for(int i = 0; i < K; i++)
            {
                std::getline(savedClusters, line);
                std::stringstream ss(line);
                double point_value;
                for(int j = 0; j < total_values; j++)
                {
                    ss >> point_value;
                    clusters[i].setCentralValue(j, point_value);
                }
            }
        }

		int iter = 1;

		while(true)
		{
			bool done = true;

			// associates each point to the nearest center
			for(int i = 0; i < total_points; i++)
			{
				int id_old_cluster = points[i].getCluster();
				int id_nearest_center = getIDNearestCenter(points[i]);

				if(id_old_cluster != id_nearest_center)
				{
					if(id_old_cluster != -1)
						clusters[id_old_cluster].removePoint(points[i].getID());

					points[i].setCluster(id_nearest_center);
					clusters[id_nearest_center].addPoint(points[i]);
					done = false;
				}
			}

			// recalculating the center of each cluster
			for(int i = 0; i < K; i++)
			{
				for(int j = 0; j < total_values; j++)
				{
					int total_points_cluster = clusters[i].getTotalPoints();
					double sum = 0.0;

					if(total_points_cluster > 0)
					{
						for(int p = 0; p < total_points_cluster; p++)
							sum += clusters[i].getPoint(p).getValue(j);
						clusters[i].setCentralValue(j, sum / total_points_cluster);
					}
				}
			}

			if(done == true || iter >= max_iterations)
				break;

			iter++;
		}

        // TO REMOVE OUTLIERS
        double dist, sum, stddev;
        std::vector<double> cluster_point_dist;
        for(int i = 0; i < K; i++)
		{
		    dist = 0.0;
		    int points_orig_total = clusters[i].getTotalPoints();

            for(int p = 0; p < points_orig_total; p++)
            {
                sum = 0.0;
                for(int j = 0; j < total_values; j++)
                {
                    sum += pow(clusters[i].getCentralValue(j) -
                               clusters[i].getPoint(p).getValue(j), 2.0);
                }
                cluster_point_dist.push_back(sqrt(sum));
                dist += sqrt(sum)/points_orig_total;
			}

			for(int p = 0; p < points_orig_total; p++)
            {
                stddev += pow(cluster_point_dist[p] - dist, 2.0);
			}
			stddev = sqrt(stddev / points_orig_total);

            int pp = 0;
			for(int p = 0; p < points_orig_total; p++)
            {
                // Swich this condition for displaying eliminated particles
                if ((cluster_point_dist[p] > (dist + 1.5*stddev)) || (cluster_point_dist[p] < (dist - 1.5*stddev)))
                    clusters[i].removePoint(clusters[i].getPoint(pp).getID());
                else pp++;
			}
        }


        std::ofstream saveClusters;
        saveClusters.open(fnClusters.c_str());

        for(int i = 0; i < K; i++)
        {
            for(int j = 0; j < total_values; j++)
            {
                saveClusters << clusters[i].getCentralValue(j) << " ";
            }
            saveClusters << std::endl;
        }
        saveClusters.close();

        return clusters;
	}
};


std::vector<double> ProgClassifyFast2D::feature_extraction()
{
    double avg, stddev, min, max;
    std::vector<double> feature_vector;

    // LBP
    unsigned char code;
    double center;
    int lbp_hist[256] = {};
    for (int y=1; y<(YSIZE(Iref())-1); y++)
    {
        for (int x=1; x<(XSIZE(Iref())-1); x++)
        {
            code = 0;
            center = DIRECT_A2D_ELEM(Iref(),y,x);
            code |= (DIRECT_A2D_ELEM(Iref(),y-1,x-1) > center) << 7;
            code |= (DIRECT_A2D_ELEM(Iref(),y-1,x  ) > center) << 6;
            code |= (DIRECT_A2D_ELEM(Iref(),y-1,x+1) > center) << 5;
            code |= (DIRECT_A2D_ELEM(Iref(),y,  x+1) > center) << 4;
            code |= (DIRECT_A2D_ELEM(Iref(),y+1,x+1) > center) << 3;
            code |= (DIRECT_A2D_ELEM(Iref(),y+1,x  ) > center) << 2;
            code |= (DIRECT_A2D_ELEM(Iref(),y+1,x-1) > center) << 1;
            code |= (DIRECT_A2D_ELEM(Iref(),y  ,x-1) > center) << 0;
            int code_min = (int) code;
            for (int i=0; i<7; i++)
            {
                unsigned char c = code & 1;  // extract the low bit
                code >>= 1;  // shift right
                code |= (c << 7);  // put the previous low bit in the high bit
                if ((int) code < code_min)
                    code_min = (int) code;
            }
            lbp_hist[code_min]++;
        }
    }

    for (int i=0; i<256; i++)
        feature_vector.push_back(lbp_hist[i]);

    // ENTROPY
/*    int hist[256] = {};
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Iref())
    {
        int idx = floor(((DIRECT_MULTIDIM_ELEM(Iref(),n) - Iref().computeMin()) * 255.0) / (Iref().computeMax() - Iref().computeMin()));
        hist[idx]++;
    }
    double entropy = 0;
    for (int i=0; i<256; i++)
    {
        entropy += std::max(hist[i], 1) * log2((std::max(hist[i], 1)));
    }
    feature_vector.push_back(entropy);*/


    if (XSIZE(masks[0]) < 1)
    {
        for (int i=0; i<(sizeof(masks)/sizeof(*masks)); i++)
        {
            masks[i].resize(Iref());
            masks[i].setXmippOrigin();
        }

        int wave_size = XSIZE(Iref()) / 2;
        int wave_size_step = XSIZE(Iref()) / 32;

//        for (int i=2; i<(sizeof(masks)/sizeof(*masks)); i++)
//        {
//            BinaryCircularMask(masks[0], wave_size);
//            BinaryCircularMask(masks[1], wave_size - wave_size_step);
//
//            // if we want overlapping regions (rings)
//            // BinaryCircularMask(masks[1], wave_size - 2*wave_size_step);
//
//            FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(masks[0])
//                DIRECT_MULTIDIM_ELEM(masks[i],n) = DIRECT_MULTIDIM_ELEM(masks[0],n) -
//                                                   DIRECT_MULTIDIM_ELEM(masks[1],n);
//
//            wave_size -= wave_size_step;
//        }

        for (int i=2; i<(sizeof(masks)/sizeof(*masks)); i++)
        {
            BinaryCircularMask(masks[0], wave_size);
            BinaryCircularMask(masks[1], wave_size - wave_size_step);

            FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(masks[0])
                DIRECT_MULTIDIM_ELEM(masks[i],n) = 2*DIRECT_MULTIDIM_ELEM(masks[1],n) -
                                                   DIRECT_MULTIDIM_ELEM(masks[0],n);

            BinaryCircularMask(masks[1], wave_size - 2*wave_size_step);
            FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(masks[0])
                DIRECT_MULTIDIM_ELEM(masks[i],n) = DIRECT_MULTIDIM_ELEM(masks[i],n) -
                                                   DIRECT_MULTIDIM_ELEM(masks[1],n);

            // if we want overlapping regions (rings)
            wave_size -= wave_size_step;
            // wave_size -= 2*wave_size_step;
        }
    }

    for (int i=2; i<(sizeof(masks)/sizeof(*masks)); i++)
    {
        // HAAR-LIKE FEATURES
        //computeStats_within_binary_mask(masks[i], Iref(), min, max, avg, stddev);
        //feature_vector.push_back(avg);
        //feature_vector.push_back(stddev);

        // RANDOM VALUES
        // feature_vector.push_back((rand()%1000)/1000);
        // feature_vector.push_back((rand()%1000)/1000);
    }

    return feature_vector;
}

void ProgClassifyFast2D::run()
{
    // Read the input metadata
    SF.read(fnSel);
    FileName fnImg, fnClass, fnTemp;
    int itemId = 0;
    MDRow row;
    MetaData MDsummary, MDclass;
    std::vector<double> fv;
    std::vector<Point> points;
    std::vector<Cluster> clusters;
    srand (time(NULL));

    FOR_ALL_OBJECTS_IN_METADATA(SF)
    {
        itemId++;
    	SF.getValue(MDL_IMAGE, fnImg,__iter.objId);
    	Iref.read(fnImg);
    	Iref().setXmippOrigin();
    	CorrelationAux aux;
    	centerImageTranslationally(Iref(), aux);

        fv = feature_extraction();
        Point p(itemId, fv);
        points.push_back(p);
    }

    std::size_t extraPath = fnOut.find_last_of("/");
    KMeans kmeans(K, itemId, fv.size(), 100);
    clusters = kmeans.run(points, fnOut.substr(0,extraPath+1) + fnClusters);
    fnTemp = fnOut.substr(0,extraPath+1) + "temp.xmd";

    std::ifstream f(fnOut.c_str());
    if (f.good())
        rename(fnOut.c_str(), fnTemp.c_str());

    std::size_t classCount;
    for(int i = 0; i < K; i++)
    {
        int total_points_cluster =  clusters[i].getTotalPoints();
        int old_points_cluster;

        size_t ii = MDsummary.addObject();
        MDsummary.setValue(MDL_REF, i+1, ii);

        std::ifstream f(fnTemp.c_str());
        if (f.good())
        {
            fnClass=formatString("classes@%s", fnTemp.c_str());
            MDclass.read(fnClass);
            MDclass.getRow(row, i+1);
            row.getValue(MDL_CLASS_COUNT, classCount);
            MDsummary.setValue(MDL_CLASS_COUNT, (size_t) total_points_cluster+classCount, ii);
        }
        else
            MDsummary.setValue(MDL_CLASS_COUNT, (size_t) total_points_cluster, ii);

        std::ostringstream clusterValues;
        clusterValues << "[";
        for(int j = 0; j < fv.size()-1; j++)
            clusterValues << clusters[i].getCentralValue(j) << ", ";
        clusterValues << clusters[i].getCentralValue(fv.size()-1) << "]";
        MDsummary.setValue(MDL_FAST2D_CENTROID, clusterValues.str(), ii);

        MDsummary.write(formatString("classes@%s", fnOut.c_str()), MD_APPEND);
        MDclass.clear();

        if (f.good())
        {
            fnClass=formatString("class%06d_images@%s", i+1, fnTemp.c_str());
            MDclass.read(fnClass);
        }

        size_t micrographId;
        for(int j = 0; j < total_points_cluster; j++)
        {
            SF.getRow(row, clusters[i].getPoint(j).getID());
            size_t recId = MDclass.addRow(row);
            MDclass.setValue(MDL_REF, i+1, recId);
            MDclass.write(formatString("class%06d_images@%s", i+1, fnOut.c_str()), MD_APPEND);
        }
    }
}