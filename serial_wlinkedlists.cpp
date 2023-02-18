#include "common.h"
#include <cmath>
#include <iostream>
#include <vector>
#include <forward_list>
#include <iostream>

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

// Linked List class
class Node {
   public:
   int partnum;
   Node* next;
};

static int nbinsx;  // number of bins in one dimension (total number of bins = nbinsx^2)
static std::vector<std::vector<Node*>> bins;        // 2D array of linkedlists, that represents particles in each bin
static std::vector<std::vector<Node*>> newelements; // 2D array of linkedlists, that holds the particles that have changed bins at the end of the time step
static double dxbin; // length&width of each bin


void init_simulation(particle_t* parts, int num_parts, double size) {
    nbinsx = ((int)((double) size / (cutoff)) + 1);  // bin size greater or equal to cutoff length
    
    /// Allocate memory ///
    for (int ib = 0; ib < nbinsx; ++ib) {
        std::vector<Node*> row,row2;
        bins.push_back(row);
        newelements.push_back(row2);
        for (int jb = 0; jb < nbinsx; ++jb) {
            // each bin's LinkedList has a head Node that is just a placeholder; it points to the first particle in the bin
            Node* head = new Node();head->partnum = -1;head->next = NULL;
            bins[ib].push_back(head);
            // The newelements array does not need such placeholders
            newelements[ib].push_back((Node*) NULL);
        }
    }
        
    dxbin = size / (double) nbinsx;

    /// INITIALIZE - compute which bin each particle is in  ///
    for (int i = 0; i < num_parts; ++i) {
        int ib = (int)(parts[i].x / dxbin);
        int jb = (int)(parts[i].y / dxbin);
        // Push i at as the second element of the bins[ib][jb] linked list (remember the head is just a placeholder)
        Node* nd = new Node();
        nd->partnum = i;
        nd->next = bins[ib][jb]->next;
        bins[ib][jb]->next = nd;
    }
}


void simulate_one_step(particle_t* parts, int num_parts, double size) {
    // Compute Forces
    for (int ib = 0; ib < nbinsx; ++ib) {
        for (int jb = 0; jb < nbinsx; ++jb) {
            // iterate over the linked list
            Node* it = bins[ib][jb]->next;
            while(it != NULL){
                int i = it->partnum;
                parts[i].ax = parts[i].ay = 0;
                // iterate over neighboring bins
                for (int local_i = fmax(0,ib-1); local_i <= fmin(nbinsx-1,ib+1); ++local_i){
                    for (int local_j = fmax(0,jb-1); local_j <= fmin(nbinsx-1,jb+1); ++local_j){
                        // iterate over the particles in each neighboring bin
                        Node* it_n = bins[local_i][local_j]->next;
                        while(it_n != NULL){
                            apply_force(parts[i], parts[it_n->partnum]);
                            it_n = it_n->next;
                        }
                    }
                }
                it = it->next;
            }
        }
    }
    
    // Move particles
    for (int ib = 0; ib < nbinsx; ++ib) {
        for (int jb = 0; jb < nbinsx; ++jb) {
            Node* it_before = bins[ib][jb];
            while(it_before->next != NULL){
                int i = it_before->next->partnum;
                move(parts[i], size);
                int new_ib = (int)(parts[i].x / dxbin);
                int new_jb = (int)(parts[i].y / dxbin);   
                // if the particle has changed bins
                if (ib != new_ib || jb != new_jb) {
                    // delete particle from current bin
                    it_before->next = it_before->next->next;
                    // emplace particle in newelements bin
                    Node* nd = new Node();
                    nd->partnum = i;
                    nd->next = newelements[new_ib][new_jb];
                    newelements[new_ib][new_jb] = nd;
                } 
                else{
                    it_before = it_before->next;
                }
            }
            
        }
    }
    
    // Take all the particles that have moved (and have been placed temporarily in newelements) and put them in their new bins
    for (int ib = 0; ib < nbinsx; ++ib) {
        for (int jb = 0; jb < nbinsx; ++jb) {
            Node* newelem = newelements[ib][jb];
            while(newelem != NULL){
                Node* nd = new Node();
                nd->partnum = newelem->partnum;
                nd->next = bins[ib][jb]->next;                
                bins[ib][jb]->next = nd;
                newelem = newelem->next;
            }   
            newelements[ib][jb] = (Node*) NULL;
        }
    }

}
