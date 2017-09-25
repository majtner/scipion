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
#include "classify_extract_features.h"
#include <data/filters.h>
#include <data/morphology.h>
#include <numeric>
#include <set>


// Read arguments ==========================================================
void ProgClassifyFast2D::readParams()
{
    fnSel = getParam("-i");
    fnOut = getParam("-o");
    K = getIntParam("-k");
    fnClusters = getParam("-c");
    maxObjects = getIntParam("-m");
}

// Show ====================================================================
void ProgClassifyFast2D::show()
{
    if (verbose==0)
        return;
    std::cerr
    << "Input selfile:                          " << fnSel      << std::endl
    << "Output selfile:                         " << fnOut      << std::endl
    << "Number of clusters:                     " << K          << std::endl
    << "Filename with clusters:                 " << fnClusters << std::endl
    << "Max. number of objects for clustering:  " << maxObjects << std::endl
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
    addParamsLine("  [-m <int=\"0\">]              : Max. number of objects for clustering");
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
	int total_values, total_points, maxIterations, maxObjects;
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
	KMeans(int K, int total_points, int total_values, int maxIterations,
	       int maxObjects)
	{
		this->K = K;
		this->total_points = total_points;
		this->total_values = total_values;
		this->maxIterations = maxIterations;
		this->maxObjects = maxObjects;
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
				int old_cluster = points[i].getCluster();
				int nearest_center = getIDNearestCenter(points[i]);

				if(old_cluster != nearest_center)
				{
					if(old_cluster != -1)
						clusters[old_cluster].removePoint(points[i].getID());

					points[i].setCluster(nearest_center);
					clusters[nearest_center].addPoint(points[i]);
					done = false;
				}
			}

			// recalculating the center of each cluster
			for(int i = 0; i < K; i++)
			{
				for(int j = 0; j < total_values; j++)
				{
					int total_points = clusters[i].getTotalPoints();
					double sum = 0.0;

					if(total_points > 0)
					{
						for(int p = 0; p < total_points; p++)
							sum += clusters[i].getPoint(p).getValue(j);
						clusters[i].setCentralValue(j, sum / total_points);
					}
				}
			}

			if(done == true || iter >= maxIterations)
				break;

			iter++;
		}

        // This code is removing outliers, whose distance from centroid is
        // 1.5*stddev, its efficiency depends strongly on feature extraction
        /*
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
                    sum += pow(clusters[i].getCentralValue(j) -
                               clusters[i].getPoint(p).getValue(j), 2.0);

                cluster_point_dist.push_back(sqrt(sum));
                dist += sqrt(sum) / points_orig_total;
			}

			for(int p = 0; p < points_orig_total; p++)
                stddev += pow(cluster_point_dist[p] - dist, 2.0);

			stddev = sqrt(stddev / points_orig_total);

            int pp = 0;
			for(int p = 0; p < points_orig_total; p++)
            {
                // Swich this condition for taking only eliminated particles
                if ((cluster_point_dist[p] > (dist + 1.5*stddev)) ||
                    (cluster_point_dist[p] < (dist - 1.5*stddev)))
                    clusters[i].removePoint(clusters[i].getPoint(pp).getID());
                else pp++;
			}
        }
        */

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

bool ProgClassifyFast2D::isParticle(size_t id)
{
    int count = XSIZE(Iref()) + YSIZE(Iref()) / 2;
    int hist[256] = {};
    int low_thresh_val, high_thresh_val, val;
    int low_thresh_cnt = 0, high_thresh_cnt = 0;
    typedef std::pair<int, int> pairs;
    std::set<pairs> low_inten_pts, high_inten_pts;

    double m, M;
    Iref().computeDoubleMinMax(m, M);
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Iref())
    {
        val = floor(((DIRECT_MULTIDIM_ELEM(Iref(), n) - m) * 255.0) / (M - m));
        hist[val]++;
    }

    for (low_thresh_val = 0; low_thresh_cnt < count; low_thresh_val++)
        low_thresh_cnt += hist[low_thresh_val];
    for (high_thresh_val = 255; high_thresh_cnt < count; high_thresh_val--)
        high_thresh_cnt += hist[high_thresh_val];

    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(Iref())
    {
        val = floor(((DIRECT_A2D_ELEM(Iref(), i, j) - m) * 255.0) / (M - m));
        if (val < low_thresh_val)
            low_inten_pts.insert(std::make_pair(i, j));
        else if (val > high_thresh_val)
            high_inten_pts.insert(std::make_pair(i, j));
    }

    double dist_high = 0.0, dist_low = 0.0;
    int comparisons = 0;

    for (std::set<pairs>:: iterator it1 = high_inten_pts.begin();
         it1 != --high_inten_pts.end(); ++it1)
    {
        pairs tmp1 = *it1;
        for (std::set<pairs>::iterator it2 = ++it1;
             it2 != high_inten_pts.end(); ++it2)
        {
            pairs tmp2 = *it2;
            dist_high += sqrt(pow(tmp1.first - tmp2.first, 2.0) +
                              pow(tmp1.second - tmp2.second, 2.0));
            comparisons++;
        }
        --it1;
    }
    dist_high = dist_high / comparisons;

    comparisons = 0;
    for (std::set<pairs>:: iterator it1 = low_inten_pts.begin();
         it1 != --low_inten_pts.end(); ++it1)
    {
        pairs tmp1 = *it1;
        for (std::set<pairs>::iterator it2 = ++it1;
             it2 != low_inten_pts.end(); ++it2)
        {
            pairs tmp2 = *it2;
            dist_low += sqrt(pow(tmp1.first - tmp2.first, 2.0) +
                             pow(tmp1.second - tmp2.second, 2.0));
            comparisons++;
        }
        --it1;
    }
    dist_low = dist_low / comparisons;



//    std::vector<double> vars;
//    double mean_all = 0.0;
//    for (int yy = 1; yy <= 4; yy++)
//    {
//        int y_max = YSIZE(Iref()) / 4 * yy;
//        int y_min = YSIZE(Iref()) / 4 * (yy-1);
//        for (int xx = 1; xx <= 4; xx++)
//        {
//            int x_max = XSIZE(Iref()) / 4 * xx;
//            int x_min = XSIZE(Iref()) / 4 * (xx-1);
//
//            double mean = 0.0;
//            double var = 0.0;
//            int count = 0;
//            for (int y = y_min; y < y_max; y++)
//            {
//                for (int x = x_min; x < x_max; x++)
//                {
//                    mean += DIRECT_A2D_ELEM(Iref(),y,x);
//                    count++;
//                }
//            }
//            mean = mean / count;
//            for (int y = y_min; y < y_max; y++)
//            {
//                for (int x = x_min; x < x_max; x++)
//                {
//                    var += (DIRECT_A2D_ELEM(Iref(),y,x) - mean) *
//                           (DIRECT_A2D_ELEM(Iref(),y,x) - mean);
//                }
//            }
//            vars.push_back(var / count);
//            mean_all += mean;
//        }
//    }
//
//    double var_all = 0.0;
//    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(Iref())
//        var_all += (DIRECT_A2D_ELEM(Iref(),i,j) - mean_all) * (DIRECT_A2D_ELEM(Iref(),i,j) - mean_all);
//    var_all = var_all / (XSIZE(Iref()) * YSIZE(Iref()));
//
//    double average = std::accumulate(vars.begin(), vars.end(), 0.0)/vars.size();
//    double variance = 0.0;
//    for (int i = 0; i < vars.size(); i++)
//	    variance += (vars[i] - average) * (vars[i] - average);

    return true;
}

void ProgClassifyFast2D::run()
{
    SF.read(fnSel);
    FileName fnImg, fnClass, fnTemp;
    int itemId = 0, countItems = 0;
    MDRow row;
    MetaData MDsummary, MDclass;
    std::vector<double> fv;
    std::vector<Point> points;
    std::vector<Cluster> clusters;
    srand (time(NULL));
    ProgExtractFeatures features;

    FOR_ALL_OBJECTS_IN_METADATA(SF)
    {
        countItems++;
    	SF.getValue(MDL_IMAGE, fnImg,__iter.objId);
    	Iref.read(fnImg);
    	Iref().setXmippOrigin();

    	Image<double> granulo;
        granulo().resizeNoCopy(Iref());

        // GRANULO DESCRIPTOR
        double m, M;
        Iref().computeDoubleMinMax(m, M);

        for (int N = 1; N < 5; N++)
        {
            int size = N*2 + 1;
            bool struct_elem[size][size];
            for (int y = 0; y < size; y++)
            {
                for (int x = 0; x < size; x++)
                    struct_elem[x][y] = pow(x - N, 2.0) +
                                        pow(y - N, 2.0) <=
                                        pow(N, 2.0);
            }

            double sum = 0.0;
            for (int y = 0; y < YSIZE(Iref()); y++)
            {
                for (int x = 0; x < XSIZE(Iref()); x++)
                {
                    double struct_min = M;
                    for (int yy = y-N; yy <= y+N; yy++)
                    {
                        if (yy < 0 || yy >= YSIZE(Iref())) continue;
                        for (int xx = x-N; xx <= x+N; xx++)
                        {
                            if (xx < 0 || xx >= XSIZE(Iref())) continue;
                            if (struct_elem[xx+N-x][yy+N-y] &&
                                DIRECT_A2D_ELEM(Iref(), yy, xx) < struct_min)
                                struct_min = DIRECT_A2D_ELEM(Iref(), yy, xx);
                        }
                    }
                    DIRECT_A2D_ELEM(granulo(),y,x) = struct_min;
                }
            }

            for (int y = 0; y < YSIZE(Iref()); y++)
            {
                for (int x = 0; x < XSIZE(Iref()); x++)
                {
                    double struct_max = m;
                    for (int yy = y-N; yy <= y+N; yy++)
                    {
                        if (yy < 0 || yy >= YSIZE(Iref())) continue;
                        for (int xx = x-N; xx <= x+N; xx++)
                        {
                            if (xx < 0 || xx >= XSIZE(Iref())) continue;
                            if (struct_elem[xx+N-x][yy+N-y] &&
                                DIRECT_A2D_ELEM(granulo(), yy, xx) > struct_max)
                                struct_max = DIRECT_A2D_ELEM(granulo(), yy, xx);
                        }
                    }
                    sum += struct_max;
                }
            }
            std::cout << sum << " ";
        }
        std::cout << std::endl;

    	CorrelationAux aux;
    	centerImageTranslationally(Iref(), aux);
    	denoiseTVFilter(Iref(), 50);

        if (isParticle(__iter.objId))
        {
            itemId++;
            fv = features.extractLBP(Iref());
            Point p(countItems, fv);
            points.push_back(p);
        }
    }

    std::size_t extraPath = fnOut.find_last_of("/");
    KMeans kmeans(K, itemId, fv.size(), 100, maxObjects);
    clusters = kmeans.run(points, fnOut.substr(0, extraPath+1) + fnClusters);
    fnTemp = fnOut.substr(0, extraPath+1) + "temp.xmd";

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
            MDsummary.setValue(MDL_CLASS_COUNT,
                (size_t) total_points_cluster + classCount, ii);
        }
        else
            MDsummary.setValue(MDL_CLASS_COUNT,
                (size_t) total_points_cluster, ii);

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

        for(int j = 0; j < total_points_cluster; j++)
        {
            SF.getRow(row, clusters[i].getPoint(j).getID());
            size_t recId = MDclass.addRow(row);
            MDclass.setValue(MDL_REF, i+1, recId);
            MDclass.write(formatString("class%06d_images@%s",
                i+1, fnOut.c_str()), MD_APPEND);
        }
    }
}