#include <omp.h>
#include "common.h"
#include <cmath>
#include <iostream>
#include <vector>
#include <forward_list>
#include <cstdio>

// Apply the force from neighbor to particle
void apply_force(particle_t& particle, particle_t& neighbor) {
    // Calculate Distance
    double dx = neighbor.x - particle.x;
    double dy = neighbor.y - particle.y;
    double r2 = dx * dx + dy * dy;

    // Check if the two particles should interact
    if (r2 > cutoff * cutoff)
        return;
    r2 = fmax(r2, min_r * min_r);
    double r = sqrt(r2);

    // Very simple short-range repulsive force
    double coef = (1 - cutoff / r) / r2 / mass;
    particle.ax += coef * dx;
    particle.ay += coef * dy;
}

// Integrate the ODE
void move(particle_t& p, double size) {
    // Slightly simplified Velocity Verlet integration
    // Conserves energy better than explicit Euler method
    p.vx += p.ax * dt;
    p.vy += p.ay * dt;
    p.x += p.vx * dt;
    p.y += p.vy * dt;

    // Bounce from walls
    while (p.x < 0 || p.x > size) {
        p.x = p.x < 0 ? -p.x : 2 * size - p.x;
        p.vx = -p.vx;
    }

    while (p.y < 0 || p.y > size) {
        p.y = p.y < 0 ? -p.y : 2 * size - p.y;
        p.vy = -p.vy;
    }
}

static int nbinsx;  // number of bins in one dimension (total number of bins = nbinsx^2)
static std::vector<std::vector<std::forward_list<int>>> bins;
static double dxbin; // length&width of each bin
//static int total,nthreadsx,nthreadsy,bins_per_thread_x,bins_per_thread_y;

void init_simulation(particle_t* parts, int num_parts, double size) {
    nbinsx = ((int)((double) size / (cutoff)) + 1);  // bin size greater or equal to cutoff length
    
    /// Allocate memory ///
    for (int ib = 0; ib < nbinsx; ++ib) {
        std::vector<std::forward_list<int>> row;
        bins.push_back(row);
        for (int jb = 0; jb < nbinsx; ++jb) {
            std::forward_list<int> list = {};
            bins[ib].push_back(list);
        }
    }
        
    dxbin = size / (double) nbinsx;

    /// INITIALIZE - compute which bin each particle is in  ///
    for (int i = 0; i < num_parts; ++i) {
        int ib = (int)(parts[i].x / dxbin);
        int jb = (int)(parts[i].y / dxbin);
        bins[ib][jb].emplace_front(i);
    }
}


typedef std::forward_list<int>::iterator IT;


void simulate_one_step(particle_t* parts, int num_parts, double size) {
    int total = omp_get_num_threads();
    int nthreadsx = (int) sqrt((double) total);
    int nthreadsy = total/nthreadsx;
    int bins_per_thread_x = nbinsx / nthreadsx +1;
    int bins_per_thread_y = nbinsx / nthreadsy +1;
    
    int rank = omp_get_thread_num();
    int ib_min = (rank/nthreadsx)*bins_per_thread_x;
    int jb_min = (rank%nthreadsy)*bins_per_thread_y;
       
    // Compute Forces
    for (int ib = ib_min; ib < fmin(nbinsx,ib_min+bins_per_thread_x); ++ib) {
        for (int jb = jb_min; jb < fmin(nbinsx,jb_min+bins_per_thread_y); ++jb) {
            for (int i : bins[ib][jb]) {
                parts[i].ax = parts[i].ay = 0;
                for (int local_i = fmax(0,ib-1); local_i <= fmin(nbinsx-1,ib+1); ++local_i){
                    for (int local_j = fmax(0,jb-1); local_j <= fmin(nbinsx-1,jb+1); ++local_j){
                        for (int j : bins[local_i][local_j]) {
                            apply_force(parts[i], parts[j]);
                        }
                    }
                }
            }
        }
    }
    #pragma omp barrier
    
    // Move Particles
    #pragma omp single
    {
    for (int i = 0; i < num_parts; ++i) {
        int old_ib = (int)(parts[i].x / dxbin);
        int old_jb = (int)(parts[i].y / dxbin);
        move(parts[i], size);
        
        int ib = (int)(parts[i].x / dxbin);
        int jb = (int)(parts[i].y / dxbin);
        if (ib != old_ib || jb != old_jb) {
            bins[old_ib][old_jb].remove(i);
            bins[ib][jb].emplace_front(i);
        }
    }
    }
    //#pragma omp single
    //{    
    //for (int ib = 0; ib < nbinsx; ++ib) {
    //    for (int jb = 0; jb < nbinsx; ++jb) {
    //        IT before = bins[ib][jb].before_begin();
    //        for (IT it = bins[ib][jb].begin(); it != bins[ib][jb].end(); ){
    //            int i = *it;
    //            move(parts[i], size);
    //            int new_ib = (int)(parts[i].x / dxbin);
    //            int new_jb = (int)(parts[i].y / dxbin);
    //            if (ib != new_ib || jb != new_jb) {
    //                it = bins[ib][jb].erase_after( before );
    //                bins[new_ib][new_jb].emplace_front(i);
    //            } 
    //            else{
    //                before = it;
    //                ++it;
    //            }
    //        }
    //        
    //    }
    //}
    //}
    //#pragma omp barrier
}