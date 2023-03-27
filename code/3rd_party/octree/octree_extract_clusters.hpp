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

#ifndef IMPL_OCTREE_EXTRACT_CLUSTERS_H
#define IMPL_OCTREE_EXTRACT_CLUSTERS_H

#include "octree_extract_clusters.h"

//////////////////////////////////////////////////////////////////////////////////////////////
template<typename PointT>
void
extractOctreeEuclideanClusters(const pcl::PointCloud <PointT>& cloud,
                               const unibn::Octree<PointT>& tree,
                               float tolerance, std::vector <pcl::PointIndices>& clusters,
                               unsigned int min_pts_per_cluster,
                               unsigned int max_pts_per_cluster) {
    // Create a bool vector of processed point indices, and initialize it to false
    std::vector<bool> processed(cloud.points.size(), false);

    std::vector <uint32_t> nn_indices;
    // Process all points in the indices vector
    for (int i = 0; i < static_cast<int> (cloud.points.size()); ++i) {
        if (processed[i])
            continue;

        std::vector<int> seed_queue;
        int sq_idx = 0;
        seed_queue.push_back(i);

        processed[i] = true;

        while (sq_idx < static_cast<int> (seed_queue.size())) {
            // Search for sq_idx
            tree.template radiusNeighbors<unibn::L2Distance<PointT>>(
                    cloud.points[seed_queue[sq_idx]], tolerance, nn_indices);

            if (nn_indices.size() == 1) { // the search point itself
                sq_idx++;
                continue;
            }

            for (unsigned int nn_index: nn_indices) {
                if (nn_index == -1
                    || nn_index == seed_queue[sq_idx] // Exclude the search point
                    || processed[nn_index])           // Has this point been processed before ?
                    continue;

                // Perform a simple Euclidean clustering
                seed_queue.push_back(nn_index);
                processed[nn_index] = true;
            }

            sq_idx++;
        }

        // If this queue is satisfactory, add to the clusters
        if (seed_queue.size() >= min_pts_per_cluster && seed_queue.size() <= max_pts_per_cluster) {
            pcl::PointIndices r;
            r.indices.resize(seed_queue.size());
            for (size_t j = 0; j < seed_queue.size(); ++j)
                r.indices[j] = seed_queue[j];

            // These two lines should not be needed: (can anyone confirm?) -FF
            std::sort(r.indices.begin(), r.indices.end());
            r.indices.erase(std::unique(r.indices.begin(), r.indices.end()), r.indices.end());

            r.header = cloud.header;
            clusters.push_back(r);   // We could avoid a copy by working directly in the vector
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////

template<typename PointT>
void
OctreeEuclideanClusterExtraction<PointT>::extract(std::vector <pcl::PointIndices>& clusters) {
    if (!initCompute() ||
        (input_ != 0 && input_->points.empty()) ||
        (indices_ != 0 && indices_->empty())) {
        clusters.clear();
        return;
    }

    // Initialize the spatial locator
    unibn::Octree<PointT> octree;

    // Send the input dataset to the spatial locator
    octree.initialize(*input_);
    extractOctreeEuclideanClusters(
            *input_, octree, cluster_tolerance_, clusters, min_pts_per_cluster_, max_pts_per_cluster_);

    // Sort the clusters based on their size (largest one first)
    std::sort(clusters.rbegin(), clusters.rend(), comparePointClusters);

    deinitCompute();
}

#endif // IMPL_OCTREE_EXTRACT_CLUSTERS_H
