//
//  Genome.h
//  GreyWolfOptimizer
//
//  Created by Pablo Reynoso, Zineng Xu, Carless Coll on 15/05/17.
//  Copyright © 2017 Pablo Reynoso. All rights reserved.
//

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include "Cromosome.h"

#define n_croms 100

class Genome{
    
private:
    
    int cromosomes;
    int generations;
    int genes;
    int dimension;
    int tournaments;
    int cross_splits;
    double cross_prob;
    double mutac_prob;
    cromosome *parents;
    cromosome *children;
    cromosome *best_cromosome;
    Cromosome C;
    
public:
    
    Genome();
    ~Genome();
    void Constructor(int, int, int, int, double, double, float);
    void create_cromosomes(void);
    void print_cromosomes(void);
    void evaluate_cromosomes(void);
    void selection(void);
    void crossover(void);
    void mutation(int);
    void elitism(void);
    void generation_gap(void);
    void GeneticAlgorithm(float *, float *, float *);
    void store_results_as_file(float, float, float, float);
    
    
};

Genome::Genome(){}

Genome::~Genome(){}


void Genome::Constructor(int cromosomes, int generations, int tournaments, int cross_splits, double crossover, double mutation,float lambda){
    
    C = Cromosome();
    C.CNP.setLambda(lambda);
    this->cromosomes = cromosomes;
    this->generations = generations;
    this->cross_prob = crossover;
    this->mutac_prob = mutation;
    this->tournaments = tournaments;
    this->cross_splits = cross_splits;
    this->genes = C.get_genes_size();
    this->dimension = C.get_dimension_size();
    this->parents = (cromosome*) malloc(this->cromosomes*sizeof(cromosome));
    this->children = (cromosome*) malloc(this->cromosomes*sizeof(cromosome));
    this->best_cromosome = (cromosome*) malloc(sizeof(cromosome));
    
}


void Genome::create_cromosomes(void){
    for(int i=0; i<this->cromosomes; i++) C.create_cromosome(&this->parents[i]);
}

void Genome::print_cromosomes(){
    for(int i=0; i<this->cromosomes; i++) C.print_cromosome(&this->parents[i]);
}

void Genome::evaluate_cromosomes(void){
    for(int i=0; i<this->cromosomes; i++) C.evaluate_cromosome(&this->parents[i]);
}


//Selection Mechanism: k-Tournament
void Genome::selection(void){
    
    if(this->tournaments > this->cromosomes) return;
    
    int k_tournament[this->tournaments];
    for(int i=0; i<(this->cromosomes/2)-1; i++){
        
        for(int j=0; j<this->tournaments; j++){
            k_tournament[j] = rand()%this->cromosomes;
        }
        
        int winner = k_tournament[0];
        float winner_apt = this->parents[winner].aptitude;
        
        for(int j=1; j<this->tournaments; j++){
            
            int current_tour = k_tournament[j];
            if(this->parents[current_tour].aptitude < winner_apt){
                winner = current_tour;
                winner_apt = this->parents[current_tour].aptitude;
            }
        }
        C.copy_cromosome(&this->parents[winner],&this->children[i]);
        
    }
    C.copy_cromosome(this->best_cromosome, &this->children[this->cromosomes/2]);
}


//Crossover Mechanism: n-splits, cross_prob
void Genome::crossover(void){
    
    if(this->cross_splits > this->genes) return;
    
    cromosome *father = (cromosome *)malloc(sizeof(cromosome));
    cromosome *mother = (cromosome *)malloc(sizeof(cromosome));
    
    int slice = this->genes / this->cross_splits;
    for(int i=0; i<this->cromosomes/2; i+=2){
        
        C.copy_cromosome(&this->children[i],father);
        C.copy_cromosome(&this->children[i+1],mother);
        
        double prob = (double)rand()/(double)RAND_MAX;
        if(prob <= this->cross_prob){
            
            for(int i=0; i<this->cross_splits; i++){
                
                if(i%2 == 1) continue;
                
                int start = i*slice;
                int end = (i+1)*slice - 1-((i+1)/3);
                C.genes_exchange(father,mother,start,end);
                
            }
            
        }
        
        C.copy_cromosome(father,&this->children[(this->cromosomes/2)+i]);
        C.copy_cromosome(mother,&this->children[(this->cromosomes/2)+i+1]);
        
    }
    
}


//Mutation Mechanism: random-crom-gene, muta_prob
void Genome::mutation(int i){
    
    int crom_to_mute = i; //rand()%this->cromosomes;
    int gene_to_mute = rand()%this->genes;
    
    int exist_gene = this->children[crom_to_mute].genes[gene_to_mute];
    
    if(exist_gene == 0) this->children[crom_to_mute].genes[gene_to_mute] = 1;
    if(exist_gene == 1) this->children[crom_to_mute].genes[gene_to_mute] = 0;
    
}

//Elitism Mechanism: select the best from current parents
void Genome::elitism(void){
    
    int elite = 0;
    int elite_aptitude = this->parents[elite].aptitude;
    for(int i=0; i<this->cromosomes; i++){
        if(this->parents[i].aptitude <= elite_aptitude){
            elite = i;
            elite_aptitude = this->parents[elite].aptitude;
        }
    }
    
    C.copy_cromosome(&this->parents[elite], this->best_cromosome);
    
}

void Genome::generation_gap(void){
    for(int i=0; i<this->cromosomes; i++) C.copy_cromosome(&this->children[i], &this->parents[i]);
}


void Genome::store_results_as_file(float link_density, float distance, float entropy, float lamb){
    
    char full_path [200],file_name [15],file_id [15], file_name2[10];
    
    struct timeval tv;
    gettimeofday(&tv,NULL);
    
    strcpy(full_path,"/Users/usuario/Desktop/CN_Optimization/CN_Optimization/GA_results/");
    strcpy(file_name,"GA_CN_measures");
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

void Genome::GeneticAlgorithm(float *ld_avg, float *d_avg, float *en_avg){
    
    //Initialize convergence curve
    float convergence_curve [this->generations];
    for(int i=0; i<this->generations; i++) convergence_curve[i] = 0.0;
    
    
    //Start Algorithm time
    clock_t begin = clock();
    
    
    this->create_cromosomes();
    for(int i=0; i<this->generations; i++){
        
        this->evaluate_cromosomes();
        this->elitism();
        this->selection();
        this->crossover();
        
        
        double m_prob = (double)rand()/(double)RAND_MAX;
        for(int i=0; i<this->cromosomes; i++){
            if(m_prob <= this->mutac_prob){
                this->mutation(i);
            }
        }
        
        this->generation_gap();
        
    }
    clock_t end = clock();
    double time_spent_s = (double)(end - begin) / CLOCKS_PER_SEC;
    
    
    int network_tmp[this->dimension*this->dimension];
    C.CNP.copy_network(&best_cromosome->genes[0], network_tmp, true);
    C.CNP.resize_network(network_tmp,"nosl_sl");
    C.CNP.print_network_graph(network_tmp,false);
    C.CNP.store_net_as_file(network_tmp);
    
    
    *ld_avg += best_cromosome->link_density;
    *d_avg += best_cromosome->average_distance;
    *en_avg += C.CNP.get_entropy(&best_cromosome->genes[0]);
    
    
}



