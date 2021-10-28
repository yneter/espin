#include <random>
#include <vector>
#include <iostream>
#include <iomanip>

typedef std::numeric_limits< double > dbl;
inline double sq(double x) { return x*x; }

struct genotype
{
  std::vector<double> gene;
  double score;
  double fitness;
  double rfitness;
  double cfitness;

  genotype(int nvars) : gene(nvars) {      
  }
};


//    Modified simple GA 
//    Original version by Dennis Cormier/Sita Raghavan/John Burkardt.
//    modified for C++
//  Reference:
//    Zbigniew Michalewicz,
//    Genetic Algorithms + Data Structures = Evolution Programs,
//    Third Edition,
//    Springer, 1996,
//    ISBN: 3-540-60676-9,
//    LC: QA76.618.M53.
//

template <class fitness_finder> class simple_GA { 
   fitness_finder &ffinder;
   std::vector<fitness_finder> ffinder_list;
   int nvars;
   int popsize;
   std::vector<double> upper;
   std::vector<double> lower;
    
   int int_uniform_ab ( int a, int b ) { 
      static thread_local std::mt19937 generator( rand() );
      std::uniform_int_distribution<int> distribution(a, b);
      return distribution(generator);
   }

   double real_uniform_ab ( double a, double b ) { 
      static thread_local std::mt19937 generator( rand() );
      std::uniform_real_distribution<double> distribution(a, b);
      return distribution(generator);
   }

   void Xover ( int one, int two ) {
      //  Select the crossover point.
      int point = int_uniform_ab ( 0, nvars - 1 );
      //  Swap genes in positions 0 through POINT-1.
      for (int i = 0; i < point; i++ ) {
	 double t = population[one].gene[i];
	 population[one].gene[i] = population[two].gene[i];
	 population[two].gene[i] = t;
      }
   }

   void copy_gene(int from, int to) { 
      for (int i = 0; i < nvars; i++ ) {
        population[to].gene[i] = population[from].gene[i];
      }
      population[to].score = population[from].score;
      population[to].fitness = population[from].fitness;      
   }

   std::vector<bool> is_fixed;
   std::vector<bool> is_integer;
   std::vector<bool> fixed_value;

public : 


   double pxover;
   double pmutation;
  //   struct std::vector<genotype> population[popsize+1];
  //   struct std::vector<genotype> newpopulation[popsize+1]; 
   struct std::vector<genotype> population;
   struct std::vector<genotype> newpopulation; 
   double temp;

   simple_GA(fitness_finder &f, int pop_size) : ffinder(f) {
      pxover = 0.8;
      pmutation = 0.1;
      popsize = pop_size;
      ffinder_list.reserve(popsize);
      for (int i = 0; i < popsize; i++) { 
	ffinder_list.push_back( f );
      }

      nvars = f.nvars();
      population.reserve(popsize+1);
      newpopulation.reserve(popsize+1);
      for (int n = 0; n <= popsize; n++) { 
	population.push_back( genotype( nvars ) );
	newpopulation.push_back( genotype( nvars ) );
      }

      is_fixed.resize(nvars);
      is_integer.resize(nvars);
      upper.resize(nvars);
      lower.resize(nvars);
      for (int var = 0; var < nvars; var++) { 
	 is_fixed[var] = false;
	 is_integer[var] = false;
      }
   }

   void mutate (void) { 
      const double a = 0.0;
      const double b = 1.0;
      double lbound;
      double ubound;
      double x;

      for (int i = 0; i < popsize; i++ ) {
	 for (int j = 0; j < nvars; j++ ) {
	    if ( is_fixed[j] ) { 
	       population[i].gene[j] = fixed_value[j];
	    } else { 
	       x = real_uniform_ab (a, b);
	       if (x < pmutation ) {	      
		  lbound = lower[j];
		  ubound = upper[j];
		  if (!is_integer[j]) { 
		     population[i].gene[j] = real_uniform_ab (lbound, ubound);
		  } else { 
population[i].gene[j] = int_uniform_ab ( round(lbound) , round(ubound) );
		  }
	       } 
	    }
	 }
      }
   }
  
   void mutate (double amplitude) { 
      const double a = 0.0;
      const double b = 1.0;
      double lbound;
      double ubound;
      double x;

      for (int i = 0; i < popsize; i++ ) {
	 for (int j = 0; j < nvars; j++ ) {
	    if ( is_fixed[j] ) { 
	       population[i].gene[j] = fixed_value[j];
	    } else { 
	       x = real_uniform_ab (a, b);
	       if ( x < pmutation ) {	      
		  lbound = std::max(lower[j], population[i].gene[j] - amplitude);
		  ubound = std::min(upper[j], population[i].gene[j] + amplitude);
		  population[i].gene[j] = real_uniform_ab (lbound, ubound);
	       }
	    }
	 }
      }
   }

   void crossover (void) {
      const double a = 0.0;
      const double b = 1.0;
      int mem;
      int one;
      int first = 0;
      
      for ( mem = 0; mem < popsize; ++mem ) {
	double x = real_uniform_ab ( a, b );
	
	if ( x < pxover ) {
	  ++first;
	  
	  if ( first % 2 == 0 ) {
	    Xover ( one, mem );
	  } else {
	    one = mem;
	  }
	}
      }
      return;
   }

   void set_as_integer(int var) { 
      is_integer[var] = true;
   }

   void set_as_float(int var) { 
      is_integer[var] = false;
   }

   void fix(int var, double value) { 
     if (var >= 0 && var < nvars) {
        is_fixed[var] = true;
	fixed_value[var] = value;
     }
   }

   void fix_to_gene(double gene[]) { 
      for (int i = 0; i < nvars; i++ ) {
         is_fixed[i] = true;
         fixed_value[i] = gene[i];
      }
   }

   void fix_to_gene(int var, double gene[]) { 
      if (var >= 0 && var < nvars) {
         is_fixed[var] = true;
	 fixed_value[var] = gene[var];
      }
   }

   void unfix(int var) { 
     if (var >= 0 && var < nvars) {
        is_fixed[var] = false;
     }
   }

// 
//  If the best individual from the new population is better than 
//  the best individual from the previous population, then 
//  copy the best from the new population; else replace the 
//  worst individual from the current population with the 
//  best one from the previous generation                     
//  
   void elitist(void) {
     int i;
     double best, worst;
     int best_mem, worst_mem;
//
//
// elitist based on scores and not fitness since scores are temperature independent. 
//   
     best = worst = population[0].score;
     best_mem = worst_mem = 0;

     for (i = 0; i < popsize - 1; ++i) {
        if ( population[i+1].score < population[i].score ) {
	   if ( best <= population[i].score ) {
	      best = population[i].score;
	      best_mem = i;
	   }
	   
	   if ( population[i+1].score <= worst ) {
	      worst = population[i+1].score;
	      worst_mem = i + 1;
	   }
	} else {
	  if ( population[i].score <= worst ) {
	     worst = population[i].score;
	     worst_mem = i;
	  }
	  if ( best <= population[i+1].score ) {
	     best = population[i+1].score;
	     best_mem = i + 1;
	  }
	}
     }

     if ( population[popsize].score <= best ) {
        copy_gene(best_mem, popsize); 
     } else {
        copy_gene(popsize, worst_mem);
     } 
   }


   void evaluate(void) {
      // when we are dealing with random fitnesses we recompute the fitness of the best individual at popsize 
      // not now 
      #pragma omp parallel for
      for (int member = 0; member < popsize; member++ ) { 
	 population[member].score = ffinder_list[member].score(population[member].gene);	 
      }
      double avscore = 0.0;
      for (int member = 0; member < popsize; member++ ) { 
	 avscore += population[member].score/(double) (popsize+1);
      }
      for (int member = 0; member < popsize; member++ ) { 
	// towards minimum
	//	 population[member].fitness = exp ( -(population[member].score - avscore)/temp );
	// towards maximum
	population[member].fitness = exp ( (population[member].score - avscore)/temp );
	//	 population[member].fitness = population[member].score;
      }
   }

   void initialize (void) {
      for (int i = 0; i < nvars; i++ ) {
	 lower[i] = ffinder.lower(i);
	 upper[i] = ffinder.upper(i);
      }
      for (int j = 0; j <= popsize; j++ ) {
	 population[j].fitness = 0;	
	 population[j].score = 0;	
	 population[j].rfitness = 0;
	 population[j].cfitness = 0;
         for (int i = 0; i < nvars; i++ ) {
	    double lbound, ubound;
	    lbound = lower[i];
	    ubound = upper[i];

	    if (is_fixed[i]) { 
	       population[j].gene[i] = fixed_value[i]; 	       
	    } else { 
	       if (!is_integer[i]) { 
		  population[j].gene[i] = real_uniform_ab (lbound, ubound);
	       } else { 
		  population[j].gene[i] = int_uniform_ab ( round(lbound) , round(ubound) );
	       }
	    }
	 }
      }
   }  

   void initial_values(int num, const double *gene) { 
      std::cout << "# setting initial_value for gene " << num << std::endl;
      for (int i = 0; i < nvars; i++ ) {
         population[num].gene[i] = gene[i];
      }
      std::cout << "# with initial score " << ffinder.score(population[num].gene) << std::endl;
   }

   void initial_values(int num, const std::vector<double> &gene) { 
      std::cout << "# setting initial_value for gene " << num << std::endl;
      for (int i = 0; i < nvars; i++ ) {
         population[num].gene[i] = gene[i];
      }
      std::cout << "# with initial score " << ffinder.score(population[num].gene) << std::endl;
   }

   void keep_the_best (void) { 
      int cur_best;
      int mem;

      cur_best = 0;
      population[popsize].fitness = 0;
      
      for ( mem = 0; mem < popsize; mem++ ) {
        if ( population[popsize].fitness < population[mem].fitness ) {
	  cur_best = mem;
	  population[popsize].fitness = population[mem].fitness;
	}
      }
      // 
      //  Once the best member in the population is found, copy the genes.
      //
      copy_gene(cur_best, popsize);
      return;
   }


   void report ( int generation ) {
      double avg;
      double best_val;
      double square_sum;
      double stddev;
      double sum;
      double sum_square;
      double av_score; 

      if ( generation == 0 ) {
	 std::cout << "\n";
	 std::cout << "Value     Generation    Best         Best       Average    Average    Standard \n";
	 std::cout << "Value     number        value        Score      fitness    score      deviation \n";
	 std::cout << "\n";
      }

      sum = 0.0;
      sum_square = 0.0;
      av_score = 0.0;

      for (int i = 0; i < popsize; i++ ) {
	 sum += population[i].fitness;
	 sum_square += population[i].fitness * population[i].fitness;
	 av_score += population[i].score;
      }

      avg = sum / ( double ) popsize;
      av_score /= (double) popsize;
      square_sum = avg * avg * popsize;
      stddev = sqrt ( ( sum_square - square_sum ) / (double) popsize );
      best_val = population[popsize].fitness;
      double best_score = population[popsize].score;

      std::cout << "  " << std::setw(8) << "equal " 
                << "  " << std::setw(8) << generation 
		<< "  " << std::setw(14) << best_val 
		<< "  " << std::setw(14) << best_score
		<< "  " << std::setw(14) << avg 
		<< "  " << std::setw(14) << av_score	
		<< "  " << std::setw(14) << stddev << "\n";

      std::cout << std::flush;
   }

 
   void selector (void) {
      const double a = 0.0;
      const double b = 1.0;
      int i;
      int j;
      int mem;
      double p;
      double sum;
      //
      //  Find the total fitness of the population.
      //
      sum = 0.0;
      for ( mem = 0; mem < popsize; mem++ ) {
	 sum = sum + population[mem].fitness;
      }
      //
      //  Calculate the relative fitness of each member.
      //
      for ( mem = 0; mem < popsize; mem++ ) {
	 population[mem].rfitness = population[mem].fitness / sum;
      }
      // 
      //  Calculate the cumulative fitness.
      //
      population[0].cfitness = population[0].rfitness;
      for ( mem = 1; mem < popsize; mem++ ) {
	 population[mem].cfitness = population[mem-1].cfitness + population[mem].rfitness;
      }
      // 
      //  Select survivors using cumulative fitness. 
      //
      for ( i = 0; i < popsize; i++ ) { 
	 p = real_uniform_ab (a, b);
	 if ( p < population[0].cfitness ) {
	    newpopulation[i] = population[0];      
	 } else {
	    // could use a dichotomic search - 
	    for ( j = 0; j < popsize; j++ ) { 
	       if ( population[j].cfitness <= p && p < population[j+1].cfitness ) {
		  newpopulation[i] = population[j+1]; 
	       }
	    }
	 }
      }
      // 
      //  Overwrite the old population with the new one.
      //
      for ( i = 0; i < popsize; i++ ) {
	 population[i] = newpopulation[i]; 
      }

      return;     
   }

   void step(void) { 
       selector();
       crossover();
       mutate();
       evaluate();
       elitist();
   }

   void print_info(void) { 
      std::cout << "# popsize " << popsize << std::endl;
      std::cout << "# nvars " << nvars << std::endl;
      std::cout << "# pxover " << pxover << std::endl;
      std::cout << "# pmutation " << pmutation << std::endl;
      std::cout << "# temp " << temp << std::endl;
   }

   void print_best(int generation) { 
      std::cout << "# best gene = " << population[popsize].fitness << "\n";
      for (int i = 0; i < nvars; i++ ) {
         std::cout << generation << "   " << i << "    " << population[popsize].gene[i] << "  %" << std::endl;
      }
      std::cout << "# with fitness = " << population[popsize].fitness << "\n";
      std::cout << "# with score = " << population[popsize].score << "\n";
   }

   void copy_best(std::vector<double> &x) { 
      x.insert(x.begin(), std::begin(population[popsize].gene), std::end(population[popsize].gene));
   }

};

/***
struct OptExample { 
   int nvars(void) { return 4 ; }
   double score(const std::vector<double> &v) { 
      double s = 0;
      for (int i = 0; i < nvars(); i++) { 
	s += (50.0 - (v[i] - (double)i)*(v[i] - (double)i));
      }
      return s;
   }

   double upper(int n) { return 10; } 
   double lower(int n) { return -10; } 
};



main()
{
   OptExample f;
   simple_GA<OptExample> ga(f, 10);
   ga.pxover = 0.8;
   ga.pmutation = 0.1;
   ga.initialize();
   ga.temp = 10.0;
   ga.initialize();

   std::vector<double> bestgene;
   double anneal_eps = 1e-4;
   double temp_min = 2e-3;
   std::cout << "# anneal_eps "  << anneal_eps << std::endl;
   std::cout << "# temp_min "  << temp_min << std::endl;

   int maxgens = 1000;
    for (int generation = 0; generation < maxgens; generation++ ) {
       ga.temp *= (1.0 - anneal_eps);
       if (ga.temp < temp_min) ga.temp = temp_min;
       ga.step();
       if (!(generation % 10)) {
	 ga.report(generation);
       }
       if (!(generation % 100)) {
	  ga.print_best(generation);	  
	  ga.copy_best(bestgene);
       }
    }


}
****/
