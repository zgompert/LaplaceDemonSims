#include <iostream>
#include <sstream>
#include <fstream>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_sf_gamma.h>
#include <float.h>
#include <math.h>

#include "pred_sims.H"

using namespace std;

// print software usage
void usage(char * name){
  fprintf(stdout,"\n%s version %s\n\n", name, VERSION);
  fprintf(stdout, "Usage: abcfs -g genefile -e envfile -f nefile -t traitfile [options]\n");
  fprintf(stdout, "-g     Infile with allele frequency data\n");
  fprintf(stdout, "-e     Infile with environmental covariate data\n");
  fprintf(stdout, "-f     Infile with varNe estimates\n");
  fprintf(stdout, "-t     Infile with trait genetic arch. estimates\n");
  fprintf(stdout, "-o     Outfile for simulation summary stats. [out_fsabc.txt]\n");
  fprintf(stdout, "-n     Number of simulations [1000]\n");
  fprintf(stdout, "-s     SS to print: 0 = bv, 1 = snp, 2 = both [0]\n");
  // fprintf(stdout, "-w     Lower bnd. on U prior for w [1]\n");
  //fprintf(stdout, "-x     Upper bnd. on U prior for w [8]\n");
  fprintf(stdout, "-p     Prior prob. of non-zero selection by component [1]\n");
  fprintf(stdout, "-a     Mean on N prior for sel. function intercept [0.18]\n");
  fprintf(stdout, "-c     SD on N prior for sel. function intercept [DBL_MIN]\n");
  fprintf(stdout, "-b     Mean on N prior for sel. function slope [-0.81]\n");
  fprintf(stdout, "-d     SD on N prior for sel. function slope [DBL_MIN]\n");
  
  exit(1);
}

// ------ Functions for input and output ---------------

// read input from the infiles
void getdata(string geneFile, string envFile, dataset * data){
  int i, j;
  string line, element;
  ifstream infile;
  istringstream stream;

  // read initial allele freq. data
  infile.open(geneFile.c_str());
  if (!infile){
    cerr << "Cannot open file " << geneFile << endl;
    exit(1);
  }

  // read line with data dimensions
  getline(infile, line);
  stream.str(line);
  stream.clear();
  stream >> element; // number of loci
  data->nLoci = atoi(element.c_str()); 
  stream >> element; // number of populations
  data->nPops = atoi(element.c_str());  

  // dynamic memory allocation for allele frequency data
  data->p0 = gsl_matrix_calloc(data->nLoci, data->nPops);
  
  // read and store allele frequencies
  for(i=0; i<data->nLoci; i++){
    getline(infile, line); // data for one locus
    stream.str(line);
    stream.clear();
    for(j=0; j<data->nPops; j++){
      stream >> element;
      gsl_matrix_set(data->p0,i, j, atof(element.c_str()));
		      
    }
  }
  infile.close();
 
  // read environmental covariate, single column of values
  infile.open(envFile.c_str());
  if (!infile){
    cerr << "Cannot open file " << envFile << endl;
    exit(1);
  }

  getline(infile, line);
  stream.str(line);
  stream.clear();
  stream >> element; // number of gens
  data->nGens = atoi(element.c_str()); 
  stream >> element; // number of populations

  // dynamic memory allocation of env data
  data->envData = gsl_matrix_calloc(data->nGens-1, data->nPops);
  
  // loop through file, one row per gen (-1) one column per pop
  for(i=0; i<(data->nGens-1); i++){
    stream.clear();
    getline(infile, line);
    stream.str(line);
    stream.clear();
    for(j=0; j<data->nPops; j++){
      stream >> element;
      gsl_matrix_set(data->envData, i, j, atof(element.c_str()));
    }
  }
  infile.close();

}

// read existing estimates of Ne
void getne(string neFile, dataset * data){
  int i, j;
  string line, element;
  ifstream infile;
  istringstream stream;

  // read ne data, first line has dimensions, combes by samples, followed by samples of Ne
  infile.open(neFile.c_str());
  if (!infile){
    cerr << "Cannot open file " << neFile << endl;
    exit(1);
  }

  // read line with data dimensions
  getline(infile, line);
  stream.str(line);
  stream.clear();
  stream >> element; // number of posterior samples
  data->nNeSams = atoi(element.c_str());
  stream >> element; // number of pops 
 

  // dynamic memory allocation for existing Ne estimates
  data->ne = gsl_matrix_calloc(data->nNeSams, data->nPops);

  // read and store Ne estimates
  for(i=0; i<data->nNeSams; i++){
    getline(infile, line); 
    stream.str(line);
    stream.clear();
    for(j=0; j<data->nPops; j++){
      stream >> element;
      gsl_matrix_set(data->ne, i, j, atof(element.c_str()));
    }
  }
  infile.close();

}

 // read existing estimates of Ne
void getqtl(string traitFile, dataset * data){
  int i, j;
  string line, element;
  ifstream infile;
  istringstream stream;

  // read ne data, first line has dimensions, combes by samples, followed by samples of Ne
  infile.open(traitFile.c_str());
  if (!infile){
    cerr << "Cannot open file " << traitFile << endl;
    exit(1);
  }

  // read line with data dimensions
  getline(infile, line);
  stream.str(line);
  stream.clear();
  stream >> element; // number of loci
  stream >> element; // 2, pip, beta
 

  // dynamic memory allocation for QTL effects
  data->qtl = gsl_matrix_calloc(data->nLoci, 2);

  // read and store QTL estimates
  for(i=0; i<data->nLoci; i++){
    getline(infile, line); // data for one locus
    stream.str(line);
    stream.clear();
    for(j=0; j<2; j++){
      stream >> element;
      gsl_matrix_set(data->qtl, i, j, atof(element.c_str()));
    }
  }
  infile.close();

}


// wf simulation wrapper function
void runsim(dataset * data, param * params, int x, int ss){

  int i, j;
  
  // sample params from priors
  // QTL effects
  samqtl(data, params);

  // selection
  // uses spike and slab prior
  // average effect
  if(gsl_ran_flat(r,0,1) < params->prsel)
    params->fsA = params->fsAlb + gsl_ran_gaussian(r, params->fsAub);
  else
    params->fsA = 0;
  // environmental effect
  if(gsl_ran_flat(r,0,1) < params->prsel)
    params->fsB =  params->fsBlb + gsl_ran_gaussian(r, params->fsBub);
  else
    params->fsB = 0;


  cout << params->fsA << " " << params->fsB;
  // conduct the simulations
  for(j=0; j<data->nPops; j++){
    simpop(data, params, j);		   
  }

  // compute derived parameters
  compdpar(data, params);

  // calculate and store the summary statistics
  for(i=0; i<data->nLoci; i++){
    if(ss > 0)
      calcss(data, params, i);
  }
    if((ss == 0) | (ss == 2))
      calcssbv(data, params);
  
}

// sample non-zero QTL effects based on their PIPs
void samqtl(dataset * data, param * params){
  int i;
  double pr;

  params->Nqtl = 0;
  gsl_vector_set_zero(params->beta);
  for(i=0; i<data->nLoci; i++){
    pr = gsl_rng_uniform(r);
    if(pr < gsl_matrix_get(data->qtl, i, 0)){ // has non-zero effect
      gsl_vector_set(params->beta, i, gsl_matrix_get(data->qtl, i, 1));
      params->Nqtl++;
    }
    // uncomment to set to model-averaged effect
    //gsl_vector_set(params->beta, i, gsl_matrix_get(data->qtl, i, 1) * gsl_matrix_get(data->qtl, i, 0));
  }
  params->sivec = gsl_vector_calloc(params->Nqtl * data->nPops * (data->nGens-1));
}

// simulation evolution for one population given a selection differential
void simpop(dataset * data, param * params, int j){
  int i, k, x;
  double S, si;
  double p, dp, pprime;
  unsigned int Ne;
  int N;

  // sample posterior value of Ne
  x = gsl_ran_flat(r, 0, data->nNeSams);
  Ne = floor(gsl_matrix_get(data->ne, x, j));
  
  // set initial allele freqs. to p0
  for(i=0; i<data->nLoci; i++){
    p = gsl_matrix_get(data->p0, i, j);
    gsl_matrix_set(params->pp, i, j * data->nGens, p); 
  }
  
  // sims
  for(k=0; k<(data->nGens-1); k++){
    N = 0;
    // calculate the selection differential
    S = params->fsA + params->fsB * gsl_matrix_get(data->envData, k, j);
    gsl_vector_set(params->Svec, j * (data->nGens-1) + k, S);
    
    // calculate/simulate dp by selection and drift for each SNP
    for(i=0; i<data->nLoci; i++){
      if( gsl_vector_get(params->beta, i) != 0){ // SNP has effect on trait
	
	si = gsl_vector_get(params->beta, i) * (S/data->sigma2); // from L&W 5.21
	// note that this is an approximation; this is a good approximation when z ~ N
	// it is a first order approximation, won't work if selection is on variance rather
	// than the mean of the trait
	gsl_vector_set(params->sivec, N + k * params->Nqtl +
		       j * (data->nGens - 1) * params->Nqtl,
		       abs(si));
	dp =  si * gsl_matrix_get(params->pp, i, j * data->nGens + k); // from L&W 5.8b
	N++;
      }
      else{
	dp = 0;
      }

      pprime = dp + gsl_matrix_get(params->pp, i, j * data->nGens + k);
      if(pprime > 1)
	pprime = 1;
      if(pprime < 0)
	pprime = 0;
      
      // sample new pi
      p = gsl_ran_binomial(r, pprime, 2 * Ne)/(2.0 * Ne);
      gsl_matrix_set(params->pp, i, j * data->nGens + k + 1, p); // store allele freq.
      // obs dp
      dp = p - gsl_matrix_get(params->pp, i, j * data->nGens + k);
      gsl_matrix_set(params->dpp, i, j * (data->nGens - 1) + k, dp); // store allele freq. change
    }
  }
}

// Compute derived paramers describing mean and variance of S and si
void compdpar(dataset * data, param * params){

  // mean and sd for S
  params->Smn = gsl_stats_mean(params->Svec->data, params->Svec->stride, params->Svec->size);
  params->Ssd = gsl_stats_sd_m(params->Svec->data, params->Svec->stride, params->Svec->size, params->Smn);

  // mean and sd for si
  params->simn = gsl_stats_mean(params->sivec->data, params->sivec->stride, params->sivec->size);
  params->sisd = gsl_stats_sd_m(params->sivec->data, params->sivec->stride, params->sivec->size, params->simn);

}

// calculate summary statistics, sum dp and cov dp env.
void calcss(dataset * data, param * params, int i){

  int j, k;
  double sumdp = 0;
  double val, envcov = 0;
  double cov = 0;
  
  gsl_vector * sumdpv;
  gsl_vector * envv;
  
  sumdpv = gsl_vector_calloc(data->nPops * (data->nGens - 1));
  envv = gsl_vector_calloc(data->nPops * (data->nGens - 1));
  
  for(j=0; j<data->nPops; j++){
    for(k=0; k<(data->nGens-1); k++){
      sumdp += gsl_matrix_get(params->dpp, i, j * (data->nGens - 1) + k);
      val = gsl_matrix_get(data->envData, k, j);
      gsl_vector_set(envv, j * (data->nGens - 1) + k, val);
      envcov += val;
    }
  }
  gsl_vector_set(params->ssSumdp, i, sumdp);

   
  //compute means
  sumdp = sumdp / (double) (data->nPops * (data->nGens - 1));
  envcov = envcov / (double) (data->nPops * (data->nGens - 1));

  // compute covariance
  for(j=0; j<data->nPops; j++){
    for(k=0; k<(data->nGens-1); k++){
      cov += (gsl_matrix_get(params->dpp, i, j * (data->nGens - 1) + k) - sumdp) *
	(gsl_vector_get(envv, j * (data->nGens - 1) + k) - envcov);
    }
  }
      
  cov = cov / (double) (data->nPops * (data->nGens - 1));
    
  gsl_vector_set(params->ssEnvcov, i, cov);
  
  gsl_vector_free(sumdpv);
  gsl_vector_free(envv);
}

// calculate breeding value level summary statistics
void calcssbv(dataset * data, param * params){

  int j, k, i;
  double sumdpA, mudp = 0;
  double val, envcov = 0;
  double bv = 0;
  double cov = 0;
 
  
  gsl_vector * mudpv;
  gsl_vector * envv;
  
  mudpv = gsl_vector_calloc(data->nPops * (data->nGens - 1));
  envv = gsl_vector_calloc(data->nPops * (data->nGens - 1));
  
  for(j=0; j<data->nPops; j++){
    for(k=0; k<(data->nGens-1); k++){
      sumdpA = 0;
      bv = 0;
      for(i=0; i<data->nLoci; i++){
	sumdpA += gsl_vector_get(params->beta, i) * gsl_matrix_get(params->dpp, i, j * (data->nGens-1) + k);
	bv += gsl_vector_get(params->beta, i) * gsl_matrix_get(params->pp, i, j * (data->nGens-1) + k);
      }
      mudp += 2 * sumdpA;
      gsl_matrix_set(params->ssMuBv, k, j, 2 * bv);
      gsl_vector_set(mudpv, j * (data->nGens - 1) + k, 2 * sumdpA);
      val = gsl_matrix_get(data->envData, k, j);
      gsl_vector_set(envv, j * (data->nGens - 1) + k, val);
      envcov += val;
    }
  }

  params->ssMudpBv = mudp;
  
  //compute means
  mudp = mudp / (double) (data->nPops * (data->nGens - 1));
  envcov = envcov / (double) (data->nPops * (data->nGens - 1));

  // compute covariance
  for(j=0; j<data->nPops; j++){
    for(k=0; k<(data->nGens-1); k++){
      cov += (gsl_vector_get(mudpv, j * (data->nGens - 1) + k) - mudp) *
	(gsl_vector_get(envv, j * (data->nGens - 1) + k) - envcov);
    }
  }
      
  cov = cov / (double) (data->nPops * (data->nGens - 1));
  
  params->ssEnvcovBv = cov; 
  
  gsl_vector_free(mudpv);
  gsl_vector_free(envv);

}
    
    
// write parameter values and summary statistics to the outfile
void writesim(dataset * data, param * params, FILE * OUT, int ss){
  int i, j;
  
  // write param values
  // fprintf(OUT, "%.3f %.3f", params->fsA, params->fsB);
  // write derived params, note si is abs(si)
  // fprintf(OUT, " %.3f %.3f %.4f %.4f", params->Smn, params->Ssd,
  // 	  params->simn, params->sisd);
  
  if((ss == 0) | (ss == 2)){
    for(j=0; j<data->nGens; j++){
      fprintf(OUT, " %.4f", gsl_matrix_get(params->ssMuBv, j, 0)); // assumes only one population
    }
  }
  if(ss > 0){
    // write ss, sum of dp

    for(i=0; i<data->nLoci; i++){
	for(j=0; j<data->nGens; j++){
	  fprintf(OUT, "%.3f ", gsl_matrix_get(params->pp, i, j));
      }
      fprintf(OUT,"\n");
    }

  }
}
