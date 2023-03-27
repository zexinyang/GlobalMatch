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

#ifndef GLOBALMATCH_STEM_MATCHING_H
#define GLOBALMATCH_STEM_MATCHING_H

#include "common.h"
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/common/distances.h>
#include <pcl/PolygonMesh.h>
#include <pcl/registration/transformation_estimation_2D.h>
#include <algorithm>
#include <vector>
#include <set>
#include <fstream>
#include <numeric>

class Matching {
public:
    Matching() : knn_(20),
                 max_stems_for_exhaustive_search_(50),
                 edge_diff_(0.05),
                 stem_positions_src_(new Cloud3D),
                 stem_positions_tgt_(new Cloud3D) {
    }

    inline void
    setPairwiseStemPositions(const Cloud3D::ConstPtr& stem_positions_src,
                             const Cloud3D::ConstPtr& stem_positions_tgt) {
        stem_positions_src_ = stem_positions_src;
        stem_positions_tgt_ = stem_positions_tgt;
    }

    void
    estimateTransformation(Eigen::Matrix4f& transform);

    inline size_t
    getNumberOfMatches() {
        return stem_matches_.size();
    };

    // TODO: add setters and getters

private:
    struct VertexSide {
        int vertex;
        float side;
    };
    typedef std::vector<VertexSide> Triangle;

    void
    constructTriangles(const Cloud3D::ConstPtr& stem_positions,
                       std::vector<Triangle>& triangles) const;

    bool
    satisfyLocalConsistency(const Point3D& feature_src,
                            const Point3D& feature_tgt) const;

    bool
    satisfyGlobalConsistency(const std::vector<int>& pair_initial,
                             const std::vector<int>& pair_candidate);

    void
    localMatching();

    void
    globalMatching();

    /**
     * @brief Accelerate global matching by growing only a limited number of (e.g., 10k)
     * randomly sampled groups of triangle pairs. Use this function if your point cloud contains
     * a large number of trees (e.g., more than 10k trees per point cloud).
     */
    void
    randomGlobalMatching();

    Cloud3D::ConstPtr stem_positions_src_, stem_positions_tgt_;
    std::vector<Triangle> triangles_src_, triangles_tgt_;
    std::vector<std::pair<int, int>> locally_matched_pairs_;
    std::vector<int> globally_matched_pairs_;
    pcl::Correspondences stem_matches_;

    int knn_;
    int max_stems_for_exhaustive_search_;
    float edge_diff_;
};


#endif //GLOBALMATCH_STEM_MATCHING_H
