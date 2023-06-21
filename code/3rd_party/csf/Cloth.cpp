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

#include "Cloth.h"
#include <fstream>


Cloth::Cloth(const Vec3& _origin_pos,
             int _num_particles_width,
             int _num_particles_height,
             double _step_x,
             double _step_y,
             double _smoothThreshold,
             double _heightThreshold,
             int rigidness,
             double time_step)
        : constraint_iterations_(rigidness),
          time_step_(time_step),
          smoothThreshold_(_smoothThreshold),
          heightThreshold_(_heightThreshold),
          origin_pos_(_origin_pos),
          step_x_(_step_x),
          step_y_(_step_y),
          num_particles_width_(_num_particles_width),
          num_particles_height_(_num_particles_height) {
    // I am essentially using this vector as an array with room for
    // num_particles_width_*num_particles_height_ particles
    particles_.resize(num_particles_width_ * num_particles_height_);

    double time_step2 = time_step * time_step;

    // creating particles in a grid of particles from (0,0,0) to
    // (width,-height,0) creating particles in a grid
    for (int i = 0; i < num_particles_width_; i++) {
        for (int j = 0; j < num_particles_height_; j++) {
            Vec3 pos(origin_pos_.f[0] + i * step_x_,
                     origin_pos_.f[1],
                     origin_pos_.f[2] + j * step_y_);

            // insert particle in column i at jth row
            particles_[j * num_particles_width_ + i] = Particle(pos, time_step2);
            particles_[j * num_particles_width_ + i].pos_x = i;
            particles_[j * num_particles_width_ + i].pos_y = j;
        }
    }

    // Connecting immediate neighbor particles with constraints
    // (distance 1 and sqrt(2) in the grid)
    for (int x = 0; x < num_particles_width_; x++) {
        for (int y = 0; y < num_particles_height_; y++) {
            if (x < num_particles_width_ - 1)
                makeConstraint(getParticle(x, y), getParticle(x + 1, y));

            if (y < num_particles_height_ - 1)
                makeConstraint(getParticle(x, y), getParticle(x, y + 1));

            if ((x < num_particles_width_ - 1) && (y < num_particles_height_ - 1))
                makeConstraint(getParticle(x, y), getParticle(x + 1, y + 1));

            if ((x < num_particles_width_ - 1) && (y < num_particles_height_ - 1))
                makeConstraint(getParticle(x + 1, y), getParticle(x, y + 1));
        }
    }

    // Connecting secondary neighbors with constraints (distance 2 and sqrt(4) in the grid)
    for (int x = 0; x < num_particles_width_; x++) {
        for (int y = 0; y < num_particles_height_; y++) {
            if (x < num_particles_width_ - 2)
                makeConstraint(getParticle(x, y), getParticle(x + 2, y));

            if (y < num_particles_height_ - 2)
                makeConstraint(getParticle(x, y), getParticle(x, y + 2));

            if ((x < num_particles_width_ - 2) && (y < num_particles_height_ - 2))
                makeConstraint(getParticle(x, y), getParticle(x + 2, y + 2));

            if ((x < num_particles_width_ - 2) && (y < num_particles_height_ - 2))
                makeConstraint(getParticle(x + 2, y), getParticle(x, y + 2));
        }
    }
}

double Cloth::timeStep() {
    int particleCount = static_cast<int>(particles_.size());
#if defined(_OPENMP)
#pragma omp parallel for
#endif
    for (int i = 0; i < particleCount; i++) {
        particles_[i].timeStep();
    }

    for (int j = 0; j < particleCount; j++) {
        particles_[j].satisfyConstraintSelf(constraint_iterations_);
    }

    double maxDiff = 0;

    for (int i = 0; i < particleCount; i++) {
        if (particles_[i].isMovable()) {
            double diff = fabs(particles_[i].old_pos.f[1] - particles_[i].pos.f[1]);

            if (diff > maxDiff)
                maxDiff = diff;
        }
    }

    return maxDiff;
}

void Cloth::addForce(const Vec3 direction) {
    for (std::size_t i = 0; i < particles_.size(); i++) {
        particles_[i].addForce(direction);
    }
}

void Cloth::terrCollision() {
    int particleCount = static_cast<int>(particles_.size());
#if defined(_OPENMP)
#pragma omp parallel for
#endif
    for (int i = 0; i < particleCount; i++) {
        Vec3 v = particles_[i].getPos();

        if (v.f[1] < heightvals_[i]) {
            particles_[i].offsetPos(Vec3(0, heightvals_[i] - v.f[1], 0));
            particles_[i].makeUnmovable();
        }
    }
}

void Cloth::movableFilter() {
    std::vector <Particle> tmpParticles;

    for (int x = 0; x < num_particles_width_; x++) {
        for (int y = 0; y < num_particles_height_; y++) {
            Particle* ptc = getParticle(x, y);

            if (ptc->isMovable() && !ptc->isVisited) {
                std::queue<int> que;
                std::vector <XY> connected; // store the connected component
                std::vector <std::vector<int>> neibors;
                int sum = 1;
                int index = y * num_particles_width_ + x;

                // visit the init node
                connected.push_back(XY(x, y));
                particles_[index].isVisited = true;

                // enqueue the init node
                que.push(index);

                while (!que.empty()) {
                    Particle* ptc_f = &particles_[que.front()];
                    que.pop();
                    int cur_x = ptc_f->pos_x;
                    int cur_y = ptc_f->pos_y;
                    std::vector<int> neibor;

                    if (cur_x > 0) {
                        Particle* ptc_left = getParticle(cur_x - 1, cur_y);

                        if (ptc_left->isMovable()) {
                            if (!ptc_left->isVisited) {
                                sum++;
                                ptc_left->isVisited = true;
                                connected.push_back(XY(cur_x - 1, cur_y));
                                que.push(num_particles_width_ * cur_y + cur_x - 1);
                                neibor.push_back(sum - 1);
                                ptc_left->c_pos = sum - 1;
                            } else {
                                neibor.push_back(ptc_left->c_pos);
                            }
                        }
                    }

                    if (cur_x < num_particles_width_ - 1) {
                        Particle* ptc_right = getParticle(cur_x + 1, cur_y);

                        if (ptc_right->isMovable()) {
                            if (!ptc_right->isVisited) {
                                sum++;
                                ptc_right->isVisited = true;
                                connected.push_back(XY(cur_x + 1, cur_y));
                                que.push(num_particles_width_ * cur_y + cur_x + 1);
                                neibor.push_back(sum - 1);
                                ptc_right->c_pos = sum - 1;
                            } else {
                                neibor.push_back(ptc_right->c_pos);
                            }
                        }
                    }

                    if (cur_y > 0) {
                        Particle* ptc_bottom = getParticle(cur_x, cur_y - 1);

                        if (ptc_bottom->isMovable()) {
                            if (!ptc_bottom->isVisited) {
                                sum++;
                                ptc_bottom->isVisited = true;
                                connected.push_back(XY(cur_x, cur_y - 1));
                                que.push(num_particles_width_ * (cur_y - 1) + cur_x);
                                neibor.push_back(sum - 1);
                                ptc_bottom->c_pos = sum - 1;
                            } else {
                                neibor.push_back(ptc_bottom->c_pos);
                            }
                        }
                    }

                    if (cur_y < num_particles_height_ - 1) {
                        Particle* ptc_top = getParticle(cur_x, cur_y + 1);

                        if (ptc_top->isMovable()) {
                            if (!ptc_top->isVisited) {
                                sum++;
                                ptc_top->isVisited = true;
                                connected.push_back(XY(cur_x, cur_y + 1));
                                que.push(num_particles_width_ * (cur_y + 1) + cur_x);
                                neibor.push_back(sum - 1);
                                ptc_top->c_pos = sum - 1;
                            } else {
                                neibor.push_back(ptc_top->c_pos);
                            }
                        }
                    }
                    neibors.push_back(neibor);
                }

                if (sum > MAX_PARTICLE_FOR_POSTPROCESSIN) {
                    std::vector<int> edgePoints = findUnmovablePoint(connected);
                    handle_slop_connected(edgePoints, connected, neibors);
                }
            }
        }
    }
}

std::vector<int> Cloth::findUnmovablePoint(std::vector <XY> connected) {
    std::vector<int> edgePoints;

    for (std::size_t i = 0; i < connected.size(); i++) {
        int x = connected[i].x;
        int y = connected[i].y;
        int index = y * num_particles_width_ + x;
        Particle* ptc = getParticle(x, y);

        if (x > 0) {
            Particle* ptc_x = getParticle(x - 1, y);

            if (!ptc_x->isMovable()) {
                int index_ref = y * num_particles_width_ + x - 1;

                if ((fabs(heightvals_[index] - heightvals_[index_ref]) < smoothThreshold_) &&
                    (ptc->getPos().f[1] - heightvals_[index] < heightThreshold_)) {
                    Vec3 offsetVec = Vec3(0, heightvals_[index] - ptc->getPos().f[1], 0);
                    particles_[index].offsetPos(offsetVec);
                    ptc->makeUnmovable();
                    edgePoints.push_back(i);
                    continue;
                }
            }
        }

        if (x < num_particles_width_ - 1) {
            Particle* ptc_x = getParticle(x + 1, y);

            if (!ptc_x->isMovable()) {
                int index_ref = y * num_particles_width_ + x + 1;

                if ((fabs(heightvals_[index] - heightvals_[index_ref]) < smoothThreshold_) &&
                    (ptc->getPos().f[1] - heightvals_[index] < heightThreshold_)) {
                    Vec3 offsetVec = Vec3(0, heightvals_[index] - ptc->getPos().f[1], 0);
                    particles_[index].offsetPos(offsetVec);
                    ptc->makeUnmovable();
                    edgePoints.push_back(i);
                    continue;
                }
            }
        }

        if (y > 0) {
            Particle* ptc_y = getParticle(x, y - 1);

            if (!ptc_y->isMovable()) {
                int index_ref = (y - 1) * num_particles_width_ + x;

                if ((fabs(heightvals_[index] - heightvals_[index_ref]) < smoothThreshold_) &&
                    (ptc->getPos().f[1] - heightvals_[index] < heightThreshold_)) {
                    Vec3 offsetVec = Vec3(0, heightvals_[index] - ptc->getPos().f[1], 0);
                    particles_[index].offsetPos(offsetVec);
                    ptc->makeUnmovable();
                    edgePoints.push_back(i);
                    continue;
                }
            }
        }

        if (y < num_particles_height_ - 1) {
            Particle* ptc_y = getParticle(x, y + 1);

            if (!ptc_y->isMovable()) {
                int index_ref = (y + 1) * num_particles_width_ + x;

                if ((fabs(heightvals_[index] - heightvals_[index_ref]) < smoothThreshold_) &&
                    (ptc->getPos().f[1] - heightvals_[index] < heightThreshold_)) {
                    Vec3 offsetVec = Vec3(0, heightvals_[index] - ptc->getPos().f[1], 0);
                    particles_[index].offsetPos(offsetVec);
                    ptc->makeUnmovable();
                    edgePoints.push_back(i);
                    continue;
                }
            }
        }
    }

    return edgePoints;
}

void Cloth::handle_slop_connected(std::vector<int> edgePoints, std::vector <XY> connected,
                                  std::vector <std::vector<int>> neibors) {
    std::vector<bool> visited;

    for (std::size_t i = 0; i < connected.size(); i++) visited.push_back(false);

    std::queue<int> que;

    for (std::size_t i = 0; i < edgePoints.size(); i++) {
        que.push(edgePoints[i]);
        visited[edgePoints[i]] = true;
    }

    while (!que.empty()) {
        int index = que.front();
        que.pop();

        int index_center = connected[index].y * num_particles_width_ + connected[index].x;

        for (std::size_t i = 0; i < neibors[index].size(); i++) {
            int index_neibor = connected[neibors[index][i]].y * num_particles_width_ + connected[neibors[index][i]].x;

            if ((fabs(heightvals_[index_center] - heightvals_[index_neibor]) < smoothThreshold_) &&
                (fabs(particles_[index_neibor].getPos().f[1] - heightvals_[index_neibor]) < heightThreshold_)) {
                Vec3 offsetVec = Vec3(0, heightvals_[index_neibor] - particles_[index_neibor].getPos().f[1], 0);
                particles_[index_neibor].offsetPos(offsetVec);
                particles_[index_neibor].makeUnmovable();

                if (visited[neibors[index][i]] == false) {
                    que.push(neibors[index][i]);
                    visited[neibors[index][i]] = true;
                }
            }
        }
    }
}

// Zexin
void Cloth::outputMesh(vtkSmartPointer <vtkPolyData>& mesh) {
    vtkSmartPointer <vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    for (auto& particle: particles_)
        points->InsertNextPoint(particle.getPos().f[0], particle.getPos().f[2], -particle.getPos().f[1]);

    /* Create the triangles
     * (see https://github.com/CloudCompare/CloudCompare/blob/
     *      b2eb26a4568767c335d7745578a9726f397cc0eb/plugins/core/Standard/qCSF/src/
     *      Cloth.cpp#L131" )
     */
    vtkSmartPointer <vtkCellArray> triangles = vtkSmartPointer<vtkCellArray>::New();
    for (int x = 0; x < num_particles_width_ - 1; ++x) {
        for (int y = 0; y < num_particles_height_ - 1; ++y) {
            // A ---------- B
            // |            |
            // D ---------- C
            int iA = y * num_particles_width_ + x;
            int iD = iA + num_particles_width_;
            int iB = iA + 1;
            int iC = iD + 1;

            vtkSmartPointer <vtkTriangle>
                    triangle_left = vtkSmartPointer<vtkTriangle>::New(),
                    triangle_right = vtkSmartPointer<vtkTriangle>::New();
            triangle_left->GetPointIds()->SetId(0, iA);
            triangle_left->GetPointIds()->SetId(1, iB);
            triangle_left->GetPointIds()->SetId(2, iD);
            triangle_right->GetPointIds()->SetId(0, iD);
            triangle_right->GetPointIds()->SetId(1, iB);
            triangle_right->GetPointIds()->SetId(2, iC);
            triangles->InsertNextCell(triangle_left);
            triangles->InsertNextCell(triangle_right);
        }
    }

    // Add the geometry and topology to the polydata
    mesh->SetPoints(points);
    mesh->SetPolys(triangles);
}
