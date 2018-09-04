//
//  ComplexNetworksProblem.h
//  GreyWolfOptimizer
//
//  Created by Pablo Reynoso, Zineng Xu, Carles Coll on 05/05/17.
//  Copyright © 2017 Pablo Reynoso. All rights reserved.
//

#include <time.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

#define max_nodes 100
#define num_objectives 2

#define unreachable 99999
#define notoptimized 9999
#define undefined -1

class ComplexNetworksProblem{
    
private:
    
    int network_nodes;
    int potential_network_arcs;
    int num_weights;
    float lambda;
    bool nosl_halfmatrix;
    float weights[num_objectives];
    
public:
    
    ComplexNetworksProblem();
    ~ComplexNetworksProblem();
    int getNetworkNumNodes(void);
    int getPotentialNetworkNumArcs(void);
    int getWeightsNum(void);
    bool getIfNoslHalfMatrix(void);
    void getObjectivesWeights(float *);
    void generate_network(int *, bool);
    void resize_network(int *, char*);
    void copy_network(int *, int *, bool);
    void print_network(int *, bool);
    void print_network_graph(int *, bool);
    float get_link_density(int *, bool);
    float get_average_distance(int *, bool);
    int factorial(int);
    float binomial_coefficient(int, int);
    int dijsktra(int*, int, int);
    void dijsktra_print_path(int *, char *, int);
    void store_net_as_file(int *);
    float get_entropy(int *);
    void setLambda(float);
    
};

ComplexNetworksProblem::ComplexNetworksProblem(){
    
    this->nosl_halfmatrix = true;
    this->num_weights = num_objectives;
    this->network_nodes = 50;
    this->potential_network_arcs = this->network_nodes*this->network_nodes;
    this->lambda = 0.0;
    this->weights[0] = lambda;
    this->weights[1] = 1.0-lambda;
    
}

ComplexNetworksProblem::~ComplexNetworksProblem(){}

void ComplexNetworksProblem::setLambda(float lamb){
    this->lambda = lamb;
    this->weights[0] = lamb;
    this->weights[1] = 1.0 - lamb;
}

int ComplexNetworksProblem::getNetworkNumNodes(){
    return this->network_nodes;
}

int ComplexNetworksProblem::getPotentialNetworkNumArcs(){
    return this->potential_network_arcs;
}

bool ComplexNetworksProblem::getIfNoslHalfMatrix(){
    return this->nosl_halfmatrix;
}

int ComplexNetworksProblem::getWeightsNum(){
    return this->num_weights;
}

void ComplexNetworksProblem::getObjectivesWeights(float *W){
    for(int i=0; i<this->num_weights; i++) *(W+i) = this->weights[i];
}

void ComplexNetworksProblem::generate_network(int *net, bool exist_nosl_halfmat){
    
    //NOTE: ConsideratE generating a Barbasi-Albert network type as initialization
    int num_nodes = this->network_nodes;
    int num_arcs = this->potential_network_arcs;
    
    if(exist_nosl_halfmat == true){
        num_arcs -= num_nodes;
        num_arcs = num_arcs / 2;
        
        for(int i=0; i<num_arcs; i++){
            double r = (double)rand()/(double)RAND_MAX;
            if (r < 0.4) *(net+i) = 1;
            else *(net+i) = 0;
        }
        
    }
    
    
    
}

void ComplexNetworksProblem::resize_network(int *net, char *conversion){
    
    
    int num_nodes = this->network_nodes;
    int num_arcs = this->potential_network_arcs;
    
    int aux_net[num_arcs];
    
    
    if(strcmp(conversion,"nosl_sl") == 0){
        
        for(int i=0; i<(num_arcs-num_nodes)/2; i++) aux_net[i] = *(net+i);
        
        for(int i=0; i<num_arcs; i++) *(net+i) = 0;
        
        int k = 0;
        for(int i=0; i<num_nodes; i++){
            for(int j=i+1; j<num_nodes; j++,k++){
                *(net+(i*num_nodes)+j) = aux_net[k];
                *(net+(j*num_nodes)+i) = aux_net[k];
            }
        }
        
    }
    
    if(strcmp(conversion,"sl_nosl") == 0){
        
        for(int i=0; i<num_arcs; i++) aux_net[i] = *(net+i);
        
        int k = 0;
        for(int i=0; i<num_nodes; i++){
            for(int j=i+1; j<num_nodes; j++){
                *(net+(k++)) = aux_net[(i*num_nodes)+j];
                
            }
        }
        
        
    }
    
    
}


void ComplexNetworksProblem::copy_network(int *net1, int *net2, bool exist_nosl_halfmat){
    
    int num_nodes = this->network_nodes;
    int num_arcs = this->potential_network_arcs;
    
    if(exist_nosl_halfmat == true){
        num_arcs -= num_nodes;
        num_arcs = num_arcs/2;
    }
    
    for(int i=0; i<num_arcs; i++) *(net2+i) = *(net1+i);
    
}


//NOTE: The network adjacency matrix have to be rescaled to a nxn matrix (with selfloops as 0) to be used.
void ComplexNetworksProblem::print_network(int *net, bool exist_nosl_halfmat){
    
    int num_nodes = this->network_nodes;
    int num_arcs = this->potential_network_arcs;
    
    if(exist_nosl_halfmat == false){
        
        printf("\n:::::::::::::::::::::::::::::::::::::::::::::::::::::");
        for(int i=0; i<num_arcs; i++){
            if(i%num_nodes == 0) printf("\n");
            if(*(net+i) == unreachable) printf("X ");
            else printf("%d ",*(net+i));
        }
        printf("\n:::::::::::::::::::::::::::::::::::::::::::::::::::::");
        
    }
    
    
    
    
}

//NOTE: The network adjacency matrix have to be rescaled to a nxn matrix (with selfloops as 0) to be used.
void ComplexNetworksProblem::print_network_graph(int *net, bool exist_nosl_halfmat){
    
    int num_nodes = this->network_nodes;
    int num_arcs = this->potential_network_arcs;
    
    if(exist_nosl_halfmat == false){
        
        printf("\n:::::::::::::::::::::::::::::::::::::::::::::::::::::");
        for(int i=0; i<num_nodes; i++){
            printf("\n( %d ) ",i+1);
            for(int j=i*num_nodes; j<(i*num_nodes)+num_nodes; j++){
                if((j%(num_nodes+1) != 0) && *(net+j) == 1) printf("-[ %d ]",(j%num_nodes)+1);
            }
        }
        printf("\n:::::::::::::::::::::::::::::::::::::::::::::::::::::");
        
        
    }
    
}


float ComplexNetworksProblem::binomial_coefficient(int n, int k=2){
    
    float c_n_k = -1.0;
    if(k <= n) c_n_k = n*(n-1)/k;
    return c_n_k;
}

void ComplexNetworksProblem::dijsktra_print_path(int *prev, char *path, int start){
    
    *(path) = '\0';
    
    int j = 0;
    int i = start;
    while(i != undefined){
        *(path+(j++)) = i+65;
        i = *(prev+i);
    }
    
    printf("\n");
    for(int k=0; k<j; k++) printf("%c", *(path+k));
    
}


//NOTE: The network adjacency matrix have to be rescaled to a nxn matrix (with selfloops as 0) to be used.
int ComplexNetworksProblem::dijsktra(int* network,int source,int target){
    
    int num_nodes = this->network_nodes;
    int num_arcs = this->potential_network_arcs;
    
    int dist[num_nodes];
    int prev[num_nodes];
    int selected[num_nodes];
    char path[num_nodes];
    
    //Initialization of distances to node and previous node for all nodes in network.
    for(int i=0; i<num_nodes; i++){
        
        dist[i] = unreachable;
        prev[i] = undefined;
        selected[i] = 0;
    }
    
    //Initializing the source node selected status and distance.
    selected[source] = 1;
    dist[source] = 0;
    
    bool exist_path = true;
    int start = source;
    while(selected[target] == 0){
        
        int min_dist_i = unreachable;
        int node_i = -2;
        for(int i=0; i<num_nodes; i++){
            
            if(selected[i] == 0){
                
                int d = dist[start] + *(network+(start*num_nodes)+i);
                if(d < dist[i]){
                    
                    dist[i] = d;
                    prev[i] = start;
                }
                if(min_dist_i > dist[i]){
                    
                    min_dist_i = dist[i];
                    node_i = i;
                }
                
                
            }
            
            
        }
        
        if(min_dist_i == unreachable && node_i == -2){
            exist_path = false;
            break;
        }
        
        start = node_i;
        selected[start] = 1;
    }
    
    int path_dist = unreachable;
    if(exist_path == true) path_dist = dist[target];
    
    return path_dist;
    
}



//Normalized Link (Network) Density Fitness Function
float ComplexNetworksProblem::get_link_density(int *net, bool exist_nosl_halfmat){
    
    int num_nodes = this->network_nodes;
    int num_arcs = this->potential_network_arcs;
    
    if(exist_nosl_halfmat == true){
        num_arcs -= num_nodes;
        num_arcs = num_arcs / 2;
    }
    
    float bin_coeff = this->binomial_coefficient(num_nodes, 2);
    
    int actual_arcs = 0;
    for(int i=0; i<num_arcs; i++) actual_arcs += *(net+i);
    
    return (1.0 / bin_coeff) * actual_arcs;
}


//Normalized Average Distance Fitness Function
float ComplexNetworksProblem::get_average_distance(int *net, bool exist_nosl_halfmat){
    
    int num_nodes = this->network_nodes;
    int num_arcs = this->potential_network_arcs;
    int network [num_arcs];
    
    this->copy_network(net,network,exist_nosl_halfmat);
    
    if(exist_nosl_halfmat == true){
        this->resize_network(network,"nosl_sl");
    }
    
    for(int i=0; i<num_arcs; i++) if(*(network+i) == 0) *(network+i) = unreachable;
    
    float bin_coeff = this->binomial_coefficient(num_nodes, 2);
    float d_linear = (num_nodes + 1) / 3;
    
    
    int min_dist_i_j = 0;
    int d = 0;
    for(int i=0; i<num_nodes; i++){
        for(int j=i+1; j<num_nodes; j++){
            min_dist_i_j = this->dijsktra(network, i, j);
            d += min_dist_i_j;
        }
    }
    
    return (1.0 * d / bin_coeff) / d_linear;
}


//NOTE: The network adjacency matrix have to be rescaled to a nxn matrix (with selfloops as 0) to be used.
void ComplexNetworksProblem::store_net_as_file(int *network){
    
    char full_path [200],file_name [15],file_id [15];
    
    struct timeval tv;
    gettimeofday(&tv,NULL);
    
    strcpy(full_path,"/Users/usuario/Desktop/CN_Optimization/CN_Optimization/networks/");
    strcpy(file_name,"GA_CN_");
    sprintf(file_id,"%ld",tv.tv_sec);
    
    strcat(full_path,file_name);
    strcat(full_path,file_id);
    strcat(full_path,".net");
    
    
    FILE *fptr;
    fptr = fopen(full_path, "w+");
    
    fprintf(fptr,"*Vertices %d",this->network_nodes);
    for(int i=0; i<this->network_nodes; i++) fprintf(fptr,"\n%d %d",i+1,i);
    
    fprintf(fptr,"\n\n*Edges");
    
    
    for(int i=0; i<this->network_nodes; i++){
        for(int j=i+1; j<this->network_nodes; j++){
            if(*(network+(i*this->network_nodes+j)) == 1){
                fprintf(fptr,"\n%d %d %d",i+1,j+1,1);
            }
        }
        
    }
    
    fclose(fptr);
    
    
}

float ComplexNetworksProblem::get_entropy(int *net){
    
    int network [this->potential_network_arcs];
    int nodes_degrees[this->potential_network_arcs];
    
    this->copy_network(net,network,true);
    this->resize_network(network,"nosl_sl");
    
    for (int i=0;i<this->potential_network_arcs;i++) nodes_degrees[i]=0;
    
    int max_degree = 0;
    for(int i=0; i<this->network_nodes; i++){
        for(int j=0; j<this->network_nodes; j++){
            nodes_degrees[i] += *(network+i * this->network_nodes+j);
        }
        if(nodes_degrees[i] > max_degree) max_degree=nodes_degrees[i];
        //printf("\n%i node has %i edges", i+1, nodes_degrees[i]);
        
    }
    //printf("\nMax degree %i", max_degree);
    
    float prob_degree[max_degree+1];
    for (int i=0; i<=max_degree; i++) prob_degree[i]=0.0;
    
    for (int i=0; i<this->network_nodes; i++) prob_degree[nodes_degrees[i]] += 1.0/this->network_nodes;
    
    float sum = 0.0;
    float entropy = 0.0;
    for (int i=0; i <=max_degree; i++) {
        //printf("\n%i degree with prob %f", i, prob_degree[i]);
        sum += prob_degree[i];
        if (prob_degree[i] != 0.0) entropy -= prob_degree[i] * log(prob_degree[i]);
        
    }
    return entropy;
}
