/*
 * Software License Agreement (BSD License)
 *
 *  Point Cloud Library (PCL) - www.pointclouds.org
 *  Copyright (c) 2010-2011, Willow Garage, Inc.
 *
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above
 *     copyright notice, this list of conditions and the following
 *     disclaimer in the documentation and/or other materials provided
 *     with the distribution.
 *   * Neither the name of the copyright holder(s) nor the names of its
 *     contributors may be used to endorse or promote products derived
 *     from this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 *  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 *  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 *  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 * $Id$
 *
 */

// Zexin: implement Euclidean clustering based on the Bonn Octree (J. Behley et al., ICRA, 2015.)

#ifndef OCTREE_EXTRACT_CLUSTERS_H
#define OCTREE_EXTRACT_CLUSTERS_H

#include <pcl/pcl_base.h>
#include "octree_unibn.hpp"

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/** \brief Decompose a region of space into clusters based on the Euclidean distance between points
  * \param cloud the point cloud message
  * \param tree the spatial locator (e.g., kd-tree) used for nearest neighbors searching
  * \note the tree has to be created as a spatial locator on \a cloud
  * \param tolerance the spatial cluster tolerance as a measure in L2 Euclidean space
  * \param clusters the resultant clusters containing point indices (as a vector of PointIndices)
  * \param min_pts_per_cluster minimum number of points that a cluster may contain (default: 1)
  * \param max_pts_per_cluster maximum number of points that a cluster may contain (default: max int)
  * \ingroup segmentation
  */
template<typename PointT>
void
extractOctreeEuclideanClusters(const pcl::PointCloud <PointT>& cloud,
                               const unibn::Octree<PointT>& tree,
                               float tolerance,
                               std::vector <pcl::PointIndices>& clusters,
                               unsigned int min_pts_per_cluster = 1,
                               unsigned int max_pts_per_cluster = (std::numeric_limits<int>::max)());

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/** \brief @b EuclideanClusterExtraction represents a segmentation class for cluster extraction in an Euclidean sense.
  * \author Radu Bogdan Rusu
  * \ingroup segmentation
  */
template<typename PointT>
class OctreeEuclideanClusterExtraction : public pcl::PCLBase<PointT> {
    typedef pcl::PCLBase <PointT> BasePCLBase;

public:
    typedef pcl::PointCloud <PointT> PointCloud;
    typedef typename PointCloud::Ptr PointCloudPtr;
    typedef typename PointCloud::ConstPtr PointCloudConstPtr;

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /** \brief Empty constructor. */
    OctreeEuclideanClusterExtraction() :
            cluster_tolerance_(0),
            min_pts_per_cluster_(1),
            max_pts_per_cluster_(std::numeric_limits<int>::max()) {};

    /** \brief Set the spatial cluster tolerance as a measure in the L2 Euclidean space
      * \param[in] tolerance the spatial cluster tolerance as a measure in the L2 Euclidean space
      */
    inline void
    setClusterTolerance(double tolerance) {
        cluster_tolerance_ = tolerance;
    }

    /** \brief Get the spatial cluster tolerance as a measure in the L2 Euclidean space. */
    inline double
    getClusterTolerance() const {
        return (cluster_tolerance_);
    }

    /** \brief Set the minimum number of points that a cluster needs to contain in order to be considered valid.
      * \param[in] min_cluster_size the minimum cluster size
      */
    inline void
    setMinClusterSize(int min_cluster_size) {
        min_pts_per_cluster_ = min_cluster_size;
    }

    /** \brief Get the minimum number of points that a cluster needs to contain in order to be considered valid. */
    inline int
    getMinClusterSize() const {
        return (min_pts_per_cluster_);
    }

    /** \brief Set the maximum number of points that a cluster needs to contain in order to be considered valid.
      * \param[in] max_cluster_size the maximum cluster size
      */
    inline void
    setMaxClusterSize(int max_cluster_size) {
        max_pts_per_cluster_ = max_cluster_size;
    }

    /** \brief Get the maximum number of points that a cluster needs to contain in order to be considered valid. */
    inline int
    getMaxClusterSize() const {
        return (max_pts_per_cluster_);
    }

    /** \brief Cluster extraction in a PointCloud given by <setInputCloud (), setIndices ()>
      * \param[out] clusters the resultant point clusters
      */
    void
    extract(std::vector <pcl::PointIndices>& clusters);

protected:
    // Members derived from the base class
    using BasePCLBase::input_;
    using BasePCLBase::indices_;
    using BasePCLBase::initCompute;
    using BasePCLBase::deinitCompute;

    /** \brief The spatial cluster tolerance as a measure in the L2 Euclidean space. */
    double cluster_tolerance_;

    /** \brief The minimum number of points that a cluster needs to contain in order to be considered valid (default = 1). */
    int min_pts_per_cluster_;

    /** \brief The maximum number of points that a cluster needs to contain in order to be considered valid (default = MAXINT). */
    int max_pts_per_cluster_;

    /** \brief Class getName method. */
    virtual std::string getClassName() const { return ("OctreeEuclideanClusterExtraction"); }

};

/** \brief Sort clusters method (for std::sort).
  * \ingroup segmentation
  */
inline bool
comparePointClusters(const pcl::PointIndices& a, const pcl::PointIndices& b) {
    return (a.indices.size() < b.indices.size());
}


#ifdef PCL_NO_PRECOMPILE
#include <pcl/segmentation/impl/extract_clusters.hpp>
#endif

#endif  //#ifndef OCTREE_EXTRACT_CLUSTERS_H
