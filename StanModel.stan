
data {
	int <lower = 0, upper = 1>  DoVariants;       // Boolean for Variants model
	int <lower = 0, upper = 1>  DoVaxVariants;    // Boolean for specific variant's VEs
	int <lower = 0, upper = 1>  DoAge;            // Boolean for age model
	int <lower = 0, upper = 1>  DoVaxAge;         // Boolean for specific age's VEs
	
	int<lower = 1>	NumDatapoints; 	// Number of Rt values (all timepoints x regions) 
	int<lower = 1>	NumDoses; 		  // Number of parameters: Vax doses
	int<lower = 1>  NumVar;         // Number of parameters: SARS-CoV-2 variants
	int<lower = 1>  NumVaxVar;      // Number of parameters: variants considered for the VEs
	int<lower = 1>  NumGroup;       // Number of parameters: age groups
	int<lower = 1>  NumVaxGroup;    // Number of parameters: agr groups considered for the VEs
	int<lower = 1>	NumLTLAs;		    // Number of regions / LTLAs
	int<lower = 1>	NumTimepoints;	// Number of timepoints (weeks)
	int<lower = 1>  NumTrendPar;    // Number of NationalTrend pars (Knots or Free)
	
	int LTLAs[NumDatapoints];  // vector giving LTLA number for each Rt-region combination.
	    
	vector[NumDatapoints]   Timepoints; // Sequence of timepoints
	
	vector <lower = 0> [NumDatapoints] 			                RtVals;   // y: Rt values across all time points and regions (expressed as giant vector) 
	// matrix <lower = 0, upper = 1> [NumDatapoints, NumDoses]	VaxProp;	// x predictor: Binary design matrix. Each row is a region and date combination.
	// instead of the above, try this.
	real<lower = 0, upper = 1> VaxProp [NumDatapoints, NumDoses, NumGroup];
	                                                                  // Each column has the proportion of vax at that time/LTLA with 1, 2, 3 doses
	real<lower = 0, upper = 1> VarProp [NumDatapoints, NumVar];  // x predictor: Each col has the proportion of variants 1, 2, etc. up to NumVar
	//matrix <lower = 0, upper = 1> [NumDatapoints, 1]        AgeProp;  // x predictor: Each age group proportion in each LTLA
	// Two changes here - first, make AgeProp a vector, not a matrix.
	// Second, it only needs a few numbers in it, namely the number of age groups, not the number of data points (i.e. LTLAs x Timepoints)
	vector <lower = 0, upper = 1> [NumGroup] AgeProp;  // x predictor: Each age group proportion in each LTLA
}


parameters {

	real<lower = 0> sigma; 
	
	vector <lower = 0>[NumLTLAs] 	 RegionalScale; 	  	    // Prev alpha, indexed by: i) LTLA; 
	vector <lower = 0>[NumTrendPar] NationalTrend_condensed; // indexed by: i) Knots/Timepoints, 
	real<lower = 0, upper = 1> VaxEffect[NumDoses, NumVaxVar, NumVaxGroup];

	vector <lower = 1>[NumVar-1] VarAdvantage; 
}

transformed parameters{
  
	real FinalRtperVariantTimeRegion 	= 0;
	
	// allocate
  
	vector [NumDatapoints] NationalTrend	= rep_vector(0, NumDatapoints); // To calculate lambda from line or free
	vector [NumDatapoints] RegionalTrends 	= rep_vector(0, NumDatapoints);
	
	vector [NumDatapoints] LogPredictions 		= rep_vector(0, NumDatapoints);

		// Initialise NationalTrend if we are not doing knots: free parameters allowed
	{  
			// initialize NationalTrend to have same values for every LTLA (i.e. repeat NationalTrend_condensed for every LTLA)
			int ind_par = 0; // initialize index
			for (i in 1:NumLTLAs){
		    NationalTrend[(ind_par + 1):(ind_par + NumTimepoints)] = NationalTrend_condensed[1:NumTimepoints];
				ind_par = ind_par + NumTimepoints; // update index
			}
	}
	
	/// I've taken out the various booleans (IncludeIntercept etc.) to save time, but can add back in later.
	for (i in 1:NumDatapoints) {
	  RegionalTrends[i] += NationalTrend[i] * RegionalScale[LTLAs[i]]; // lambda * RegionalScale^T in manuscript. Note use of intercept makes this line akin to IntDim (B) = 2 with column vector of 1s for one column of lambda
	}

	// Danny changes - consolodating loops below
	{ // Bracket to compile
	int VaxVariantIndex;  // set to 1 by default
	int VaxAgeIndex    ;  // set to 1 by default
	
	real ProductOfDoses_PerVariantPerAgeGroup; 
	real WeightedSumEfficacyOverAgeGroups_PerVariant;   
	
	for (TimeRegion in 1:NumDatapoints) 
	{
	  LogPredictions[TimeRegion] = 0; // initialize final Rt to zero for each time point and region.
	  
	  for (Variant in 1:NumVar) // sum over variants
	  {
	      // choose appropriate index for variant vaccine 
	      if (NumVaxVar == NumVar) VaxVariantIndex = Variant; else VaxVariantIndex = 1; // don't need the else as is specified by default above, but you get idea.
        
        // (re-)initialize to 0
        WeightedSumEfficacyOverAgeGroups_PerVariant = 0;   
        
        // Select Variant Advantage parameter 

	      FinalRtperVariantTimeRegion = VarProp[TimeRegion, Variant] * RegionalTrends[TimeRegion]; 
	      if (Variant != 1) FinalRtperVariantTimeRegion = FinalRtperVariantTimeRegion * VarAdvantage[Variant - 1];
	      
	      for (Group in 1:NumGroup)
	      {
  	      // choose appropriate index for age vaccine 
	        if (NumVaxGroup == NumGroup) VaxAgeIndex = Group; else VaxAgeIndex = 1; // don't need the else as is specified by default above, but you get idea.
	        
	        // calculate product (i.e. effect) of all doses (for this variant and age group)
	        ProductOfDoses_PerVariantPerAgeGroup = 1; // (re-)initialize.
  	      for (Dose in 1:NumDoses)
  	      {
  	        ProductOfDoses_PerVariantPerAgeGroup                   *= 
  	            (1 - (VaxProp[TimeRegion,Dose, Group]              *     // note has Group, not index
  	            VaxEffect[Dose, VaxVariantIndex, VaxAgeIndex]));         // note this has indices, not Group/Variant.
  	      }
  	      // add to age group sum
  	      WeightedSumEfficacyOverAgeGroups_PerVariant += AgeProp[Group] * ProductOfDoses_PerVariantPerAgeGroup;
	      }
	      
	      FinalRtperVariantTimeRegion *= WeightedSumEfficacyOverAgeGroups_PerVariant; 

	      LogPredictions[TimeRegion] += FinalRtperVariantTimeRegion; 
	  }
	}
	
	} // Bracket to compile
}

model {
  
  RegionalScale ~ normal(1, 3);
  NationalTrend_condensed ~ normal(1, 3);
	
	sigma 		  ~ std_normal();
	
	for (i in 1:NumDoses)
	  for (j in 1:NumVaxVar)
	    for(k in 1:NumVaxGroup)
		    VaxEffect[i, j, k] ~ beta(1, 1);
	
	VarAdvantage ~ normal(1, 4);
		
	RtVals 			~ normal(LogPredictions, sigma);
}

generated quantities {
	vector[NumDatapoints] log_lik;
	
	for (i in 1:NumDatapoints) {
		log_lik[i] = normal_lpdf(RtVals[i] | LogPredictions[i], sigma); // NOTES: Log of normal function: real normal_lpdf(reals y | reals mu, reals sigma)
	}
}
