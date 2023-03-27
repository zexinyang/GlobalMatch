// ======================================================================================
// Copyright 2017 State Key Laboratory of Remote Sensing Science, 
// Institute of Remote Sensing Science and Engineering, Beijing Normal University

// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at

//     http://www.apache.org/licenses/LICENSE-2.0

// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
// ======================================================================================

#define DLL_IMPLEMENT

#include "CSF.h"
#include "XYZReader.h"
#include "Vec3.h"
#include "Rasterization.h"
#include "c2cdist.h"
#include <fstream>


CSF::CSF(int index) {
    params.bSloopSmooth = true;
    params.time_step = 0.65;
    params.dist_min = 0.2;
    params.dist_max = 3.0;
    params.cloth_resolution = 0.5;
    params.rigidness = 2;
    params.interations = 500;

    this->index = index;
}

CSF::CSF() {
    params.bSloopSmooth = true;
    params.time_step = 0.65;
    params.dist_min = 0.2;
    params.dist_max = 3.0;
    params.cloth_resolution = 0.5;
    params.rigidness = 2;
    params.interations = 500;
    this->index = 0;
}

CSF::~CSF() = default;

// Zexin: to be compatible with PCL point clouds
void CSF::setPointCloud(pcl::PointCloud<pcl::PointXYZ>::ConstPtr pcl_cloud) {
    point_cloud.resize(pcl_cloud->size());

    int pointCount = static_cast<int>(pcl_cloud->size());
#if defined(_OPENMP)
#pragma omp parallel for
#endif
    for (int i = 0; i < pointCount; i++) {
        csf::Point las;
        las.x = pcl_cloud->points[i].x;
        las.y = -pcl_cloud->points[i].z;
        las.z = pcl_cloud->points[i].y;
        point_cloud[i] = las;
    }
}

void CSF::setPointCloud(csf::PointCloud& pc) {
    point_cloud.resize(pc.size());
    int pointCount = static_cast<int>(pc.size());
#if defined(_OPENMP)
#pragma omp parallel for
#endif
    for (int i = 0; i < pointCount; i++) {
        csf::Point las;
        las.x = pc[i].x;
        las.y = -pc[i].z;
        las.z = pc[i].y;
        point_cloud[i] = las;
    }
}

void CSF::do_filtering(std::vector<int>& ptid_stem_area,
                       std::vector<int>& ptid_others,
                       bool exportCloth,
                       vtkSmartPointer <vtkPolyData>& mesh) {
    // Terrain
    //std::cout << "[" << this->index << "] Configuring terrain..." << std::endl;
    csf::Point bbMin, bbMax;
    point_cloud.computeBoundingBox(bbMin, bbMax);

    double cloth_y_height = 0.05;

    int clothbuffer_d = 2;
    Vec3 origin_pos(
            bbMin.x - clothbuffer_d * params.cloth_resolution,
            bbMax.y + cloth_y_height,
            bbMin.z - clothbuffer_d * params.cloth_resolution
    );

    int width_num = static_cast<int>(
                            std::floor((bbMax.x - bbMin.x) / params.cloth_resolution)
                    ) + 2 * clothbuffer_d;

    int height_num = static_cast<int>(
                             std::floor((bbMax.z - bbMin.z) / params.cloth_resolution)
                     ) + 2 * clothbuffer_d;

    //std::cout << "[" << this->index << "] Configuring cloth..." << std::endl;
    //std::cout << "[" << this->index << "]  - width: " << width_num << " "
    //          << "height: " << height_num << std::endl;

    Cloth cloth(
            origin_pos,
            width_num,
            height_num,
            params.cloth_resolution,
            params.cloth_resolution,
            0.3,
            9999,
            params.rigidness,
            params.time_step
    );

    //std::cout << "[" << this->index << "] Rasterizing..." << std::endl;
    Rasterization::RasterTerrian(cloth, point_cloud, cloth.getHeightvals());

    double time_step2 = params.time_step * params.time_step;
    double gravity = 0.2;

    //std::cout << "[" << this->index << "] Simulating..." << std::endl;
    cloth.addForce(Vec3(0, -gravity, 0) * time_step2);

    // boost::progress_display pd(params.interations);
    for (int i = 0; i < params.interations; i++) {
        double maxDiff = cloth.timeStep();
        cloth.terrCollision();
        //params.class_threshold / 100
        if ((maxDiff != 0) && (maxDiff < 0.005)) {
            // early stop
            break;
        }
        // pd++;
    }

    if (params.bSloopSmooth) {
        //std::cout << "[" << this->index << "]  - post handle..." << std::endl;
        cloth.movableFilter();
    }

    if (exportCloth)
        cloth.outputMesh(mesh);

    c2cdist c2c(params.dist_min, params.dist_max);
    c2c.calCloud2CloudDist(cloth, point_cloud, ptid_stem_area, ptid_others);
}

void CSF::savePoints(std::vector<int> grp, std::string path) {
    if (path == "") {
        return;
    }

    std::ofstream f1(path.c_str(), std::ios::out);

    if (!f1)
        return;

    for (int i: grp) {
        f1 << std::fixed << std::setprecision(8)
           << point_cloud[i].x << "	"
           << point_cloud[i].z << "	"
           << -point_cloud[i].y << std::endl;
    }

    f1.close();
}
