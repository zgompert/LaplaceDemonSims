// file: main.C for psim

// simulates evolution of a polygenic trait by fluctuating selection

// Requires the following input files:
// geneFile: point estimates of allele freqs. nloci rows by npop columns, just the initial freqs (as freqs)
// envFile: environmental covariates, one row per generation, one column per pop, note that the gen t (the final gen.) should be left off
// neFile: post. samples for the variance effective population size, one row per sample, one column per pop
// traitFile: trait gen arch. estimates, one row per SNP, pip followed by beta | lambda = 1
// NOTE: files begin with a row giving the dimensions of the file (row then column)

// Time-stamp: <Friday, 17 August 2012, 14:46 CDT -- zgompert>

#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>
#include <time.h>
#include <getopt.h>
#include <omp.h>
#include "pred_sims.H"

using namespace std;

gsl_rng * r;  /* global state variable for random number generator */

/* ----------------- */
/* Beginning of main */
/* ----------------- */

int main(int argc, char *argv[]) {

  time_t start = time(NULL);
  time_t end;
  int rng_seed = 0;
  int ch = 0;
  int x;
  int nsims = 1000;
  int ss = 0; // 0 = bv ss, 1 = snp ss, 2 = bv and snp ss
  
  string geneFile = "undefined";
  string envFile = "undefined";
  string neFile = "undefined";
  string traitFile = "undefined";
  string outFile = "out_fsabc.txt";
  
  dataset data;
  param params;

  // set defaults
  data.sigma2 = 0.02;
  params.fsAlb = 0.18;// here lb = mean, ub = sd
  params.fsAub = DBL_MIN; // or use 0.2
  params.fsBlb = -0.81;
  params.fsBub = DBL_MIN;
  params.prsel = 1;

  
  // get command line arguments
  if (argc < 2) {
    usage(argv[0]);
  }
  
  while ((ch = getopt(argc, argv, "g:e:f:t:o:n:s:a:c:b:d:p:")) != -1){
    switch(ch){
    case 'g':
      geneFile = optarg;
      break;
    case 'e':
      envFile = optarg;
      break;
    case 'f':
      neFile = optarg;
      break;
    case 't':
      traitFile = optarg;
      break;
    case 'o':
      outFile = optarg;
      break;
    case 'n':
      nsims = atoi(optarg);
      break;
    case 's':
      ss = atoi(optarg);
      break;
    case 'a':
      params.fsAlb = atof(optarg);
      break;
    case 'c':
      params.fsAub = atof(optarg);
      break;
    case 'b':
      params.fsBlb = atof(optarg);
      break;
    case 'd':
      params.fsBub = atof(optarg);
      break;
    case 'p':
      params.prsel = atof(optarg);
      break;
    case '?':
    default:
      usage(argv[0]);
    }
  }
  
  // set up gsl random number generation 
  gsl_rng_env_setup();
  r = gsl_rng_alloc (gsl_rng_default);
  srand(time(NULL));
  rng_seed = rand();
  gsl_rng_set(r, rng_seed); /* seed gsl_rng with output of rand, which
                               was seeded with result of time(NULL) */

  // read infiles, record genotype likelihoods and data dimensions
  cout << "Reading input from files: " << geneFile << " and " << 
    envFile << endl;
  getdata(geneFile, envFile, &data);

  // read Ne from a file
  cout << "Reading posterior samples of Ne from " << neFile << endl;
  getne(neFile,&data);

  // read Ne from a file
  cout << "Reading GWA data " << traitFile << endl;
  getqtl(traitFile,&data);

  // memory allocation for params
  params.beta = gsl_vector_calloc(data.nLoci);
  params.pp = gsl_matrix_calloc(data.nLoci, data.nPops * data.nGens);
  params.dpp =  gsl_matrix_calloc(data.nLoci, data.nPops * (data.nGens-1));
  params.ssSumdp = gsl_vector_calloc(data.nLoci);
  params.ssEnvcov = gsl_vector_calloc(data.nLoci);
  params.Svec = gsl_vector_calloc(data.nPops * (data.nGens-1));
  params.ssMuBv = gsl_matrix_calloc(data.nGens, data.nPops);
  
  // open outfile
  FILE * OUT;
  OUT = fopen(outFile.c_str(), "w");
  
  // run wf sims
  for(x=0; x<nsims; x++){
    cout << "sim " << x << endl;
    runsim(&data, &params, x, ss);
    // write sims
    writesim(&data, &params, OUT, ss);
    gsl_vector_free(params.sivec); // need to free each time b/c allocated by samqtl
    // and depends on the # of qtl sampled
  }

  // close outfile
  fclose(OUT);

  // prints run time
  end = time(NULL);
  cout << "Runtime: " << (end-start)/3600 << " hr " << (end-start)%3600/60 << " min ";
  cout << (end-start)%60 << " sec" << endl;
  return 0;
}
