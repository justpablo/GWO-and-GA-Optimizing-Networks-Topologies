//
//  main.cpp
//  CN_Optimization
//
//  Created by Pablo Reynoso, Zineng Xu, Carless Coll on 10/06/17.
//  Copyright © 2017 Pablo Reynoso. All rights reserved.
//


#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//#include "Genome.h"
#include "Swarm.h"

using namespace std;

int main(){
    
    float lambda = 0.0;
    float link_densities[51];
    float average_distances[51];
    float entropies[51];
    
    cout<<"Grey Wolf Optimizer begin running...."<<endl;
    
    Swarm *pack_of_packs = new Swarm[51];
    
    srand(17);
    
    for (int i=0; i<51; i++){
        float ld_avg = 0.0;
        float d_avg = 0.0;
        float en_avg = 0.0;
        
        //Swarm S = Swarm(200,15, lambda);
        for (int j=0; j<1; j++){
            printf("\ni = %i   j = %i", i, j);
            pack_of_packs[i].Constructor(200, 15, lambda);
            pack_of_packs[i].GrayWolfOptimizer(&ld_avg, &d_avg, &en_avg);
            
        }
        link_densities[i] = ld_avg/50;
        average_distances[i] = d_avg/50;
        entropies[i] = en_avg/50;
        printf("\nlamb: %f  LD: %f  AD: %f  H: %f", lambda, link_densities[i], average_distances[i], entropies[i]);
        pack_of_packs[i].store_results_as_file(link_densities[i], average_distances[i], entropies[i], lambda);
        lambda += 0.02;
    }
    
    cout<<"Grey Wolf Optimizer finished running."<<endl;
    
    
    /*
    cout<<"Genetic Algorithm begin running...."<<endl;
    
    Genome *pack_of_packs = new Genome[51];
    
    srand(17);
    
    for (int i=0; i < 51; i++){
        float ld_avg = 0.0;
        float d_avg = 0.0;
        float en_avg = 0.0;
        
        //Swarm S = Swarm(200,15, lambda);
        for (int j=0; j <1; j++){
            printf("\ni = %i   j = %i", i, j);
            pack_of_packs[i].Constructor(20,100,5,3,0.3,0.05,lambda);
            pack_of_packs[i].GeneticAlgorithm(&ld_avg, &d_avg, &en_avg);
            
        }
        link_densities[i] = ld_avg/50;
        average_distances[i] = d_avg/50;
        entropies[i] = en_avg/50;
        printf("\nlamb: %f  LD: %f  AD: %f  H: %f", lambda, link_densities[i], average_distances[i], entropies[i]);
        pack_of_packs[i].store_results_as_file(link_densities[i], average_distances[i], entropies[i], lambda);
        lambda += 0.02;
    }
    
    

    cout<<"Genetic Algorithm finished running."<<endl;
    */
    return 0;
    
}
    
