//
//  Swarm.h
//  GreyWolfOptimizer
//
//  Created by Pablo Reynoso, Zineng Xu, Carless Coll on 10/05/17.
//  Copyright © 2017 Pablo Reynoso. All rights reserved.
//


#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "Wolf.h"

#define n_pack 100

using namespace std;


class Swarm{
    
private:
    
    int n_gen;
    int n_wolves;
    int wolf_dims;
    int dimension;
    int l_bound;
    int u_bound;
    wolf * wolves_pack;
    Wolf W;
    
    
public:
    
    Swarm(void);
    ~Swarm(void);
    void create_pack_wolves(void);
    void print_pack_wolves(void);
    void evaluate_pack_wolves(void);
    void GrayWolfOptimizer(float *, float *, float *);
    void Constructor(int, int, float);
    void store_results_as_file(float, float, float, float);
    
};

Swarm::Swarm(){}

Swarm::~Swarm(){
    for(int i=0; i<this->n_wolves; i++) free(&this->wolves_pack[i]);
}

void Swarm::Constructor(int n_generations, int n_wolves, float lambda){
    W.CNP.setLambda(lambda);
    this->n_gen = n_generations;
    this->n_wolves = n_wolves;
    this->wolf_dims = W.get_space_size();
    this->dimension = W.get_dimension_size();
    this->wolves_pack = (wolf*) malloc(n_wolves*sizeof(wolf));
}

void Swarm::create_pack_wolves(){
    
    char wolf_type [10] = "omega";
    
    for(int i=0; i<this->n_wolves; i++){
        wolf *w = (wolf*) malloc(sizeof(wolf));
        W.create_wolf(w,wolf_type);
        this->wolves_pack[i] = *w;
    }
    
}

void Swarm::print_pack_wolves(){
    for(int i=0; i<this->n_wolves; i++) W.print_wolf(&this->wolves_pack[i]);
}

void Swarm::evaluate_pack_wolves(){
    for(int i=0; i<this->n_wolves; i++) W.evaluate_wolf(&this->wolves_pack[i]);
}

void Swarm::store_results_as_file(float link_density, float distance, float entropy, float lamb){
    
    char full_path [200],file_name [15],file_id [15], file_name2[10];
    
    struct timeval tv;
    gettimeofday(&tv,NULL);
    
    strcpy(full_path,"/Users/usuario/Desktop/CN_Optimization/CN_Optimization/GWO_results/");
    strcpy(file_name,"GWO_CN_measures");
    sprintf(file_name2,"_%f_", lamb);
    sprintf(file_id,"%ld",tv.tv_sec);
    
    strcat(full_path,file_name);
    strcat(full_path,file_id);
    strcat(full_path, file_name2);
    strcat(full_path,".net");
    
    
    FILE *fptr;
    fptr = fopen(full_path, "w+");
    
    fprintf(fptr,"\n%f;%f;%f", link_density, distance, entropy);
    fclose(fptr);
    
}

void Swarm::GrayWolfOptimizer(float *ld_avg, float *d_avg, float *en_avg){
    
    //alpha, beta, delta wolves
    wolf *alpha_wolf = (wolf*) malloc(sizeof(wolf));
    wolf *beta_wolf  = (wolf*) malloc(sizeof(wolf));
    wolf *delta_wolf = (wolf*) malloc(sizeof(wolf));
    
    //mutated wolves
    wolf *mutated_wolves = (wolf*) malloc(3 * sizeof(wolf));
    
    
    //Initialize alpha, beta, delta wolves
    char wolf_type1 [10] = "alpha";
    W.create_wolf(alpha_wolf, wolf_type1);
    
    char wolf_type2 [10] = "beta";
    W.create_wolf(beta_wolf, wolf_type2);
    
    char wolf_type3 [10] = "delta";
    W.create_wolf(delta_wolf, wolf_type3);
    
    
    //Initialize wolves population considering proper problem solution space boundaries.
    this->create_pack_wolves();
    
    //Initialize convergence curve
    float convergence_curve [this->n_gen];
    for(int i=0; i<this->n_gen; i++) convergence_curve[i] = 0.0;
    
    
    //Start Algorithm time
    clock_t begin = clock();
    
    char ans[5];
    for(int l=0; l<this->n_gen; l++){
        
        for(int i=0; i<this->n_wolves; i++){
            
            //Calculate objective function for each search agent
            W.evaluate_wolf(&this->wolves_pack[i]);
            
            //Update Alpha, Beta, and Delta
            if (this->wolves_pack[i].aptitude < alpha_wolf->aptitude){
                
                //Wolves shiftting positions
                W.copy_wolf(beta_wolf,delta_wolf);
                W.copy_wolf(alpha_wolf,beta_wolf);
                W.copy_wolf(&this->wolves_pack[i],alpha_wolf);
                
            }else if(this->wolves_pack[i].aptitude < beta_wolf->aptitude){
                
                //Wolves shiftting positions
                W.copy_wolf(beta_wolf,delta_wolf);
                W.copy_wolf(&this->wolves_pack[i],beta_wolf);
                
            }else if(this->wolves_pack[i].aptitude < delta_wolf->aptitude){
                
                W.copy_wolf(&this->wolves_pack[i],delta_wolf);
                
            }
            
            
        }
        
        // mutate alpha, beta and delta wolves
        W.copy_wolf(alpha_wolf, &mutated_wolves[0]);
        W.copy_wolf(beta_wolf, &mutated_wolves[1]);
        W.copy_wolf(delta_wolf, &mutated_wolves[2]);
        
        W.mutate_wolf(&mutated_wolves[0]);
        W.mutate_wolf(&mutated_wolves[1]);
        W.mutate_wolf(&mutated_wolves[2]);
        
        // select alpha, beta and delta among the mutated and original wolves
        for(int i=0; i<3; i++){
            
            //Calculate objective function for each search agent
            W.evaluate_wolf(&mutated_wolves[i]);
            
            //Update Alpha, Beta, and Delta
            if (mutated_wolves[i].aptitude < alpha_wolf->aptitude){
                
                //Wolves shiftting positions
                W.copy_wolf(beta_wolf,delta_wolf);
                W.copy_wolf(alpha_wolf,beta_wolf);
                W.copy_wolf(&mutated_wolves[i],alpha_wolf);
                
            }else if(mutated_wolves[i].aptitude < beta_wolf->aptitude){
                
                //Wolves shiftting positions
                W.copy_wolf(beta_wolf,delta_wolf);
                W.copy_wolf(&mutated_wolves[i],beta_wolf);
                
            }else if(mutated_wolves[i].aptitude < delta_wolf->aptitude){
                
                W.copy_wolf(&mutated_wolves[i],delta_wolf);
                
            }
        }
        
        double a = 2.0-(2.0*l/this->n_gen); // a decreases linearly fron 2 to 0
        
        double r1;
        double r2;
        double A1[this->wolf_dims];
        double C1[this->wolf_dims];
        double A2[this->wolf_dims];
        double C2[this->wolf_dims];
        double A3[this->wolf_dims];
        double C3[this->wolf_dims];
        
        for(int i=0;i<this->wolf_dims; i++){
            
            r1 = (double)rand()/(double)RAND_MAX;
            r2 = (double)rand()/(double)RAND_MAX;
            A1[i] = 2.0*a*r1 - a;
            C1[i] = 2.0*r2;
            
            r1 = (double)rand()/(double)RAND_MAX;
            r2 = (double)rand()/(double)RAND_MAX;
            A2[i] = 2.0*a*r1 - a;
            C2[i] = 2.0*r2;
            
            r1 = (double)rand()/(double)RAND_MAX;
            r2 = (double)rand()/(double)RAND_MAX;
            A3[i] = 2.0*a*r1 - a;
            C3[i] = 2.0*r2;
            
        }
        
        
        //Update the Position of search agents including omegas
        for(int i=0; i<this->n_wolves; i++){
            
            //Obtaining alpha_wolf positional contribution
            double D_alpha[this->wolf_dims];
            double X1[this->wolf_dims];
            W.directioning_wolf(alpha_wolf,&this->wolves_pack[i],C1,D_alpha); // Eq(3.5)-p1
            W.wolf_direction_contribution(alpha_wolf,D_alpha,A1,X1); // Eq (3.6)-part 1
            
            //Obtaining beta_wolf positional contribution
            double D_beta[this->wolf_dims];
            double X2[this->wolf_dims];
            W.directioning_wolf(beta_wolf,&this->wolves_pack[i],C2,D_beta); // Eq(3.5)-p2
            W.wolf_direction_contribution(beta_wolf,D_alpha,A2,X2); // Eq (3.6)-part 2
            
            //Obtaining delta_wolf positional contribution
            double D_delta[this->wolf_dims];
            double X3[this->wolf_dims];
            W.directioning_wolf(delta_wolf,&this->wolves_pack[i],C3,D_delta); // Eq(3.5)-p3
            W.wolf_direction_contribution(delta_wolf,D_alpha,A3,X3); // Eq (3.6)-part 3
            
            
            W.wolf_optimal_moving(&this->wolves_pack[i],X1,X2,X3);  // Equation (3.7)
            
            
        }
        
        convergence_curve[l] = alpha_wolf->aptitude;
    }
    
    
    clock_t end = clock();
    double time_spent_s = (double)(end - begin) / CLOCKS_PER_SEC;
     
    int wolf_dim = W.get_dimension_size();
    int network_tmp[wolf_dim*wolf_dim];
    
    W.CNP.copy_network(&alpha_wolf->spatial_position[0], network_tmp, true);
    W.CNP.resize_network(network_tmp,"nosl_sl");
    W.CNP.print_network_graph(network_tmp,false);
    W.CNP.store_net_as_file(network_tmp);
      
    *ld_avg += alpha_wolf->link_density;
    *d_avg += alpha_wolf->average_distance; 
    *en_avg += W.CNP.get_entropy(&alpha_wolf->spatial_position[0]);

}
