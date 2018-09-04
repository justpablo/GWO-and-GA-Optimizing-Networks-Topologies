//
//  Cromosome.h
//  Genetic Algorithm
//
//  Created by Pablo Reynoso, Zineng Xu, Carless Coll on 15/05/17.
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

struct cromosome{
    
    int genes[max_dim];
    float aptitude;
    float link_density;
    float average_distance;
    
};


class Cromosome{
    
private:
    
    int genes_size;
    int dimension;
    bool nosl_halfmatrix;
    
    
public:
    
    ComplexNetworksProblem CNP;
    Cromosome ();
    ~Cromosome();
    int get_genes_size(void);
    int get_dimension_size(void);
    void create_cromosome(cromosome *);
    void copy_cromosome(cromosome *, cromosome *);
    void evaluate_cromosome(cromosome *);
    void print_cromosome(cromosome *);
    void genes_exchange(cromosome *, cromosome *, int, int);
    
    
};

Cromosome::Cromosome(){
    
    CNP = ComplexNetworksProblem();
    this->dimension = CNP.getNetworkNumNodes();
    this->genes_size = this->dimension * this->dimension;
    this->nosl_halfmatrix = CNP.getIfNoslHalfMatrix();
    if(this->nosl_halfmatrix == true){
        this->genes_size -= this->dimension;
        this->genes_size = this->genes_size / 2;
    }
    
}

Cromosome::~Cromosome(){}

int Cromosome::get_genes_size(){
    return this->genes_size;
}

int Cromosome::get_dimension_size(){
    return this->dimension;
}

void Cromosome::create_cromosome(cromosome *c1){
    
    
    int aux_cromosome[this->dimension*this->dimension];
    CNP.generate_network(aux_cromosome,this->nosl_halfmatrix);
    for(int i=0; i<this->genes_size; i++) *(c1->genes+i) = aux_cromosome[i];
    
    c1->aptitude = notoptimized;
    c1->average_distance = notoptimized;
    c1->link_density = notoptimized;
    
}

void Cromosome::copy_cromosome(cromosome *c1, cromosome *c2){
    
    c2->aptitude = c1->aptitude;
    c2->average_distance = c1->average_distance;
    c2->link_density = c1->link_density;
    for(int i=0; i<this->genes_size; i++) *(c2->genes+i) = *(c1->genes+i);
    
}

void Cromosome::evaluate_cromosome(cromosome *c){
    
    int num_weights = CNP.getWeightsNum();
    
    float weights[num_weights];
    CNP.getObjectivesWeights(weights);
    
    float avg_distance = CNP.get_average_distance(c->genes,this->nosl_halfmatrix);
    float lnk_density = CNP.get_link_density(c->genes,this->nosl_halfmatrix);
    
    c->aptitude = weights[0]*avg_distance + weights[1]*lnk_density;
    
    c->average_distance = avg_distance;
    c->link_density = lnk_density;
    
}


void Cromosome::print_cromosome(cromosome *c){
    
    printf("\n:::::::::::::::::::::::::::");
    printf("\nCromosome");
    printf("\nGenes:");
    printf("\n[ ");
    for(int i=0; i<this->genes_size; i++) printf("%d ",*(c->genes+i));
    printf("]");
    printf("\nAptitude: %f",c->aptitude);
    printf("\nLink_Density: %f",c->link_density);
    printf("\nAverage_Distance: %f",c->average_distance);
    printf("\n::::::::::::::::::::::::::");
    
}

//Crossover-Swapping: swap genes from different cromosomes
void Cromosome::genes_exchange(cromosome *c1, cromosome *c2, int start, int end){
    
    for(int i=start; i<=end; i++){
        int aux_gene = c1->genes[i];
        c1->genes[i] = c2->genes[i];
        c2->genes[i] = aux_gene;
    }
    
}







