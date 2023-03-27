/*
 * Software License Agreement (Apache License)
 *
 *  Copyright (C) 2023, Xufei Wang (tjwangxufei@tongji.edu.cn),
 *                      Zexin Yang (zexinyang@tongji.edu.cn),
 *                      Liangliang Nan (liangliang.nan@gmail.com).
 *  All rights reserved.
 *
 *  This file is part of GlobalMatch (https://github.com/zexinyang/GlobalMatch),
 *  which implements the point cloud registration method described in the following paper:
 *  -----------------------------------------------------------------------------------------------------------
 *  GlobalMatch: Registration of forest terrestrial point clouds by global matching of relative stem positions.
 *  Xufei Wang, Zexin Yang, Xiaojun Cheng, Jantien Stoter, Wenbing Xu, Zhenlun Wu, and Liangliang Nan.
 *  ISPRS Journal of Photogrammetry and Remote Sensing. Vol. 197, 71-86, 2023.
 *  -----------------------------------------------------------------------------------------------------------
 *  We kindly ask you to cite the above paper if you use (part of) the code or ideas in your academic work.
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *  http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 */

#ifndef GLOBALMATCH_STEM_MAPPING_H
#define GLOBALMATCH_STEM_MAPPING_H

#include "common.h"
#include <pcl/features/normal_3d.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>

typedef pcl::PointIndices Indices;
typedef pcl::ModelCoefficients Coefficient;
typedef std::vector<Coefficient> CoefficientVector;

class Mapping {
public:
    Mapping() : leaf_size_(0.01),
                search_radius_(0.10),
                verticality_threshold_(0.90),
                min_pts_per_cluster_(200),
                min_dist_between_stems_(2 * search_radius_),
                mesh_ground_(vtkSmartPointer<vtkPolyData>::New()),
                indices_understory_(new Indices),
                cloud_stems_(new Cloud3D) {
    }

    virtual ~Mapping() = default;

    inline void
    setInputCloud(const Cloud3D::ConstPtr& cloud_input) {
        cloud_input_ = cloud_input;
    }

    inline void
    setLeafSize(float leaf_size) {
        leaf_size_ = leaf_size;
    }

    inline void
    setSearchRadius(float search_radius) {
        search_radius_ = search_radius;
    }

    inline void
    setVerticalityThreshold(float verticality_threshold) {
        verticality_threshold_ = verticality_threshold;
    }

    inline void
    setMinClusterSize(int min_cluster_size) {
        min_pts_per_cluster_ = min_cluster_size;
    }

    // TODO: add more setters and getters

    void extract(Cloud3D::Ptr& tree_positions);

private:

    void extractUnderstory();

    void extractStemPoints();

    void extractTreePositions(Cloud3D::Ptr& tree_positions);

    Cloud3D::ConstPtr cloud_input_;
    Indices::Ptr indices_understory_;
    Cloud3D::Ptr cloud_stems_;
    vtkSmartPointer<vtkPolyData> mesh_ground_;

    // Subsampling
    float leaf_size_;
    // Verticality-based filtering
    float search_radius_;
    float verticality_threshold_;
    // Euclidean clustering
    int min_pts_per_cluster_;
    float min_dist_between_stems_;
};


#endif //GLOBALMATCH_STEM_MAPPING_H
