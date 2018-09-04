//
//  Wolf.h
//  GreyWolfOptimizer
//
//  Created by Pablo Reynoso, Zineng Xu, Carless Coll on 10/05/17.
//  Copyright © 2017 Pablo Reynoso. All rights reserved.
//

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

using namespace std;

#define max_dim 10000 //considering (wolf) network max_nodes = 100.
#define hg_length 10 //group name length

#include "ComplexNetworksProblem.h"

struct wolf{
    
    int spatial_position[max_dim];
    float aptitude;
    char herarchy_group [hg_length];
    float link_density;
    float average_distance;
    
};


class Wolf{
    
private:
    
    int space;
    int dimension;
    bool nosl_halfmatrix;
    
public:
    
    ComplexNetworksProblem CNP;
    Wolf();
    ~Wolf();
    int get_space_size(void);
    int get_dimension_size(void);
    void create_wolf(wolf *, char *);
    void copy_wolf(wolf *, wolf *);
    void evaluate_wolf(wolf *);
    void directioning_wolf(wolf *, wolf *, double *, double *);
    void wolf_direction_contribution(wolf *, double *, double *, double *);
    void wolf_optimal_moving(wolf *, double *, double *, double *);
    void constrain_wolf(double *);
    void print_wolf(wolf *);
    void mutate_wolf(wolf *);
    
    
};

Wolf::Wolf(){
    
    CNP = ComplexNetworksProblem();
    this->dimension = CNP.getNetworkNumNodes();
    this->space = this->dimension * this->dimension;
    this->nosl_halfmatrix = CNP.getIfNoslHalfMatrix();
    if(this->nosl_halfmatrix == true){
        this->space -= this->dimension;
        this->space = this->space / 2;
    }
    
}

Wolf::~Wolf(){}

int Wolf::get_space_size(){
    return this->space;
}

int Wolf::get_dimension_size(){
    return this->dimension;
}

void Wolf::create_wolf(wolf* w1, char *t){
    
    for(int i=0; i<hg_length; i++, t++) *(w1->herarchy_group+i) = *t;
    
    int aux_wolf[this->dimension*this->dimension];
    CNP.generate_network(aux_wolf,this->nosl_halfmatrix);
    for(int i=0; i<this->space; i++) *(w1->spatial_position+i) = aux_wolf[i];
    
    w1->aptitude = notoptimized;
    w1->average_distance = notoptimized;
    w1->link_density = notoptimized;
}

void Wolf::copy_wolf(wolf *w1, wolf *w2){
    
    w2->aptitude = w1->aptitude;
    w2->average_distance = w1->average_distance;
    w2->link_density = w1->link_density;
    for(int i=0; i<this->space; i++) *(w2->spatial_position+i) = *(w1->spatial_position+i);
    
}

void Wolf::evaluate_wolf(wolf *w){
    
    int num_weights = CNP.getWeightsNum();
    
    float weights[num_weights];
    CNP.getObjectivesWeights(weights);
    
    float avg_distance = CNP.get_average_distance(w->spatial_position,this->nosl_halfmatrix);
    float lnk_density = CNP.get_link_density(w->spatial_position,this->nosl_halfmatrix);
    
    w->aptitude = weights[0]*avg_distance + weights[1]*lnk_density;
    
    w->average_distance = avg_distance;
    w->link_density = lnk_density;
    
}


void Wolf::directioning_wolf(wolf *w1, wolf *w2, double *C, double *D){
    for(int i=0;i<this->space; i++) *(D+i) = fabs(*(C+i)*w1->spatial_position[i] - w2->spatial_position[i]);
}

void Wolf::wolf_direction_contribution(wolf *w1, double *D, double *A, double *X){
    for(int i=0; i<this->space; i++) *(X+i) = w1->spatial_position[i] - *(D+i)*(*(A+i));
}

void Wolf::wolf_optimal_moving(wolf *w, double *X1, double *X2, double *X3){
    
    double aux_wolf [this->space];
    for(int i=0; i<this->space; i++) aux_wolf[i] = ((*(X1+i) + *(X2+i) + *(X3+i)) / 3.0);
    
    this->constrain_wolf(aux_wolf);
    for(int i=0;i<this->space; i++) w->spatial_position[i] = (int) *(aux_wolf+i);
    
}

void Wolf::constrain_wolf(double *wolf_optimal_pos){
    
    double pos_space_mean = 0.0;
    for(int i=0;i<this->space; i++) pos_space_mean += *(wolf_optimal_pos+i);
    pos_space_mean /= this->space;
    
    for(int i=0;i<this->space; i++){
        if(*(wolf_optimal_pos+i) > pos_space_mean) *(wolf_optimal_pos+i) = 1;
        else *(wolf_optimal_pos+i) = 0;
    }
    
}


void Wolf::print_wolf(wolf* w){
    
    printf("\n:::::::::::::::::::::::::::");
    printf("\nWolf: Canis Lupus");
    printf("\nHierarchy Group: %s",w->herarchy_group);
    printf("\nPosition:");
    printf("\n[ ");
    for(int i=0; i<this->space; i++) printf("%d ",*(w->spatial_position+i));
    printf("]");
    printf("\nAptitude: %f",w->aptitude);
    printf("\nLink_Density: %f",w->link_density);
    printf("\nAverage_Distance: %f",w->average_distance);
    printf("\n::::::::::::::::::::::::::");
    
}


void Wolf::mutate_wolf(wolf * w){
    
    for (int i = 0; i<this->space; i++){
        if (rand()%100 < 4){
            int connection = *(w->spatial_position+i);
            if(connection == 0) *(w->spatial_position+i)=1;
            if(connection == 1) *(w->spatial_position+i)=0;
        }
    }
    
}
