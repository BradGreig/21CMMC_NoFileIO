#include "../Parameter_files/INIT_PARAMS.H"
#include "../Parameter_files/ANAL_PARAMS.H"
#include "Variables.h"
#include "bubble_helper_progs.c"
#include "filter.c"
#include "gsl/gsl_sf_erf.h"

/**** Arrays declared and used *****/

void init_21cmMC_arrays() {
    
    xH = (float *) fftwf_malloc(sizeof(float)*HII_TOT_NUM_PIXELS);
    
    deltax_unfiltered_original = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);
    deltax_filtered = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);
    deltax = (float *) malloc(sizeof(float)*HII_TOT_FFT_NUM_PIXELS);
    Fcoll = (float *) malloc(sizeof(float)*HII_TOT_FFT_NUM_PIXELS);
    delta_T = (float *) malloc(sizeof(float)*HII_TOT_NUM_PIXELS);
    
    k_factor = 1.25;
	k_first_bin_ceil = DELTA_K;
	k_max = DELTA_K*HII_DIM;
	// initialize arrays
	// ghetto counting (lookup how to do logs of arbitrary bases in c...)
	NUM_BINS = 0;
	k_floor = 0;
	k_ceil = k_first_bin_ceil;
	while (k_ceil < k_max){
		NUM_BINS++;
		k_floor=k_ceil;
		k_ceil*=k_factor;
	}
    
    p_box = (double *)malloc(sizeof(double)*NUM_BINS);
    k_ave =  (double *)malloc(sizeof(double)*NUM_BINS);
    in_bin_ct = (unsigned long long *)malloc(sizeof(unsigned long long)*NUM_BINS);
    
    deldel_T = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);

}

void destroy_21cmMC_arrays() {
    
    fftwf_free(xH);
    
    fftwf_free(deltax_unfiltered_original);
    fftwf_free(deltax_filtered);
    free(deltax);
    free(Fcoll);
    free(delta_T);
    
    free(p_box);
    free(k_ave);
    free(in_bin_ct);
    
    fftwf_free(deldel_T);
}

struct ReturnParams{
    float x;
    float y;
};

struct ReturnParamsUsefulFloat{
    float x;
};

struct ReturnParamsUsefulInt{
    int x;
    int y;
};

struct ReturnParamsUsefulInt IntFunction() {
    
    struct ReturnParamsUsefulInt result;
    result.x = HII_DIM;
    result.y = VELOCITY_COMPONENT;
    
    return result;
}

struct ReturnParamsUsefulFloat FloatFunction() {
    
    struct ReturnParamsUsefulFloat result;
    result.x = BOX_LEN;
    
    return result;
}

struct ReturnParams drive_21CMMC(fftwf_complex *deltax_unfiltered, float *v, double *MockObservation, double *TelescopeSensitivity, float ForeGroundCut, float ShotNoiseCut, float ModellingUncertainty, float REDSHIFT, float ION_EFF_FACTOR, float MFP, float TVIR_MIN, int PERFORM_PS){

    omp_set_num_threads(NUMCORES);
    
    fftwf_plan plan;
    
    struct ReturnParams DataToReturn;
    
    // Various parameters to be used for the MCMC code
    int short_completely_ionised;
    
    // Other parameters used in the code
    float growth_factor;

    float R, R_begin, pixel_mass, cell_length_factor;
	float ave_N_min_cell, M_MIN;
	int x,y,z, N_min_cell, LAST_FILTER_STEP;
	unsigned long long ion_ct;
	float f_coll_crit, pixel_volume, density_over_mean, erfc_num, erfc_denom, erfc_denom_cell, res_xH;
	float xHI_from_xrays, std_xrays;
	
    double global_xH, global_step_xH, ave_xHI_xrays, ave_den, ST_over_PS, mean_f_coll_st, mean_f_coll_ps, f_coll, ave_fcoll;
	const gsl_rng_type * T;
	gsl_rng * r;

    unsigned long long ct, nonlin_ct;
    
    int i,j,k;
    
    float pixel_x_HI, pixel_deltax, H;
	int n_x, n_y, n_z;
	double dvdx, ave, max_v_deriv, temp;
	float max, maxi, maxj, maxk, maxdvdx, min, mini, minj, mink, mindvdx;
	float k_x, k_y, k_z, k_mag;
	float const_factor, T_rad;
    
    TVIR_MIN = pow(10.,TVIR_MIN);
    
    /**** Perform 'perturb_field.c' *****/
			
	/***************   BEGIN INITIALIZATION   **************************/
			 
	// perform a very rudimentary check to see if we are underresolved and not using the linear approx
	if ((BOX_LEN > DIM) && !EVOLVE_DENSITY_LINEARLY){
		printf("perturb_field.c: WARNING: Resolution is likely too low for accurate evolved density fields\n It Is recommended that you either increase the resolution (DIM/Box_LEN) or set the EVOLVE_DENSITY_LINEARLY flag to 1\n");
	}
    
	// initialize power spectrum
	init_ps();
	growth_factor = dicke(REDSHIFT);
    
    init_21cmMC_arrays();

    memcpy(deltax_unfiltered_original, deltax_unfiltered, sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);

	i=0;
		
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);

	pixel_volume = pow(BOX_LEN/(float)HII_DIM, 3);
	pixel_mass = RtoM(L_FACTOR*BOX_LEN/(float)HII_DIM); 
	f_coll_crit = 1/ION_EFF_FACTOR;
	cell_length_factor = L_FACTOR;
	
	//set the minimum source mass
	if (TVIR_MIN > 0){ // use the virial temperature for Mmin
		if (TVIR_MIN < 9.99999e3) // neutral IGM
			M_MIN = TtoM(REDSHIFT, TVIR_MIN, 1.22);
		else // ionized IGM
			M_MIN = TtoM(REDSHIFT, TVIR_MIN, 0.6);
	}
	else if (TVIR_MIN < 0){ // use the mass
		M_MIN = ION_M_MIN;
	}
	// check for WDM
	if (P_CUTOFF && ( M_MIN < M_J_WDM())){
		printf( "The default Jeans mass of %e Msun is smaller than the scale supressed by the effective pressure of WDM.\n", M_MIN);
		M_MIN = M_J_WDM();
		printf( "Setting a new effective Jeans mass from WDM pressure supression of %e Msun\n", M_MIN);
	}
	   
    for (ct=0; ct<HII_TOT_NUM_PIXELS; ct++){
        xH[ct] = 1;
    }

	// lets check if we are going to bother with computing the inhmogeneous field at all...
	mean_f_coll_st = FgtrM_st(REDSHIFT, M_MIN);
	mean_f_coll_ps = FgtrM(REDSHIFT, M_MIN);
	if (mean_f_coll_st/f_coll_crit < HII_ROUND_ERR){ // way too small to ionize anything...
		printf( "The ST mean collapse fraction is %e, which is much smaller than the effective critical collapse fraction of %e\n I will just declare everything to be neutral\n", mean_f_coll_st, f_coll_crit);
        // find the neutral fraction
        global_xH = 1;
        
        for (ct=0; ct<HII_TOT_NUM_PIXELS; ct++){
            xH[ct] = 1;
        }
		// print out the xH box
	}
		
	plan = fftwf_plan_dft_r2c_3d(HII_DIM, HII_DIM, HII_DIM, (float *)deltax_unfiltered, (fftwf_complex *)deltax_unfiltered, FFTW_ESTIMATE);

	fftwf_execute(plan);
	fftwf_destroy_plan(plan);
    fftwf_cleanup();
	// remember to add the factor of VOLUME/TOT_NUM_PIXELS when converting from
	//  real space to k-space
	// Note: we will leave off factor of VOLUME, in anticipation of the inverse FFT below
        
    for (ct=0; ct<HII_KSPACE_NUM_PIXELS; ct++){
        deltax_unfiltered[ct] /= (HII_TOT_NUM_PIXELS+0.0);
    }
    
    short_completely_ionised = 0;
	// loop through the filter radii (in Mpc)
	erfc_denom_cell=1; //dummy value
	R=fmin(MFP, L_FACTOR*BOX_LEN);
    R_begin = R;
	LAST_FILTER_STEP = 0;
    
	while (!LAST_FILTER_STEP){//(R > (cell_length_factor*BOX_LEN/(HII_DIM+0.0))){
        
        if ((R/DELTA_R_HII_FACTOR) <= (cell_length_factor*BOX_LEN/(float)HII_DIM)){
			LAST_FILTER_STEP = 1;
		}
		memcpy(deltax_filtered, deltax_unfiltered, sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);
        
		HII_filter(deltax_filtered, HII_FILTER, R);
        
        plan = fftwf_plan_dft_c2r_3d(HII_DIM, HII_DIM, HII_DIM, (fftwf_complex *)deltax_filtered, (float *)deltax_filtered, FFTW_ESTIMATE);
		fftwf_execute(plan);
        fftwf_destroy_plan(plan);
        fftwf_cleanup();
        
		// Check if this is the last filtering scale.  If so, we don't need deltax_unfiltered anymore.
		// We will re-read it to get the real-space field, which we will use to se the residual
		// neutral fraction
		ST_over_PS = 0;
		f_coll = 0;
		if (LAST_FILTER_STEP){
			
            memcpy(deltax_unfiltered, deltax_unfiltered_original, sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);
            
            erfc_denom_cell = sqrt( 2*(pow(sigma_z0(M_MIN), 2) - pow(sigma_z0(RtoM(cell_length_factor*BOX_LEN/(float)HII_DIM)), 2) ) );
            if (erfc_denom_cell < 0){  // our filtering scale has become too small
                break;
            }
                
            // renormalize the collapse fraction so that the mean matches ST,
            // since we are using the evolved (non-linear) density field
            for (x=0; x<HII_DIM; x++){
                for (y=0; y<HII_DIM; y++){
                    for (z=0; z<HII_DIM; z++){
                        density_over_mean = 1.0 + *((float *)deltax_unfiltered + HII_R_FFT_INDEX(x,y,z));
                        erfc_num = (Deltac - (density_over_mean-1)) /  growth_factor;
                        f_coll += splined_erfc(erfc_num/erfc_denom_cell);
                        Fcoll[HII_R_FFT_INDEX(x,y,z)] = splined_erfc(erfc_num/erfc_denom_cell);
                    }
                }
            }
            f_coll /= (double) HII_TOT_NUM_PIXELS;
            ST_over_PS = mean_f_coll_st/f_coll;
            
		} // end if last filter step conditional statement
        
		// not the last filter step, and we operating on the density field
		else if (!USE_HALO_FIELD){
            
            erfc_denom = sqrt( 2*(pow(sigma_z0(M_MIN), 2) - pow(sigma_z0(RtoM(R)), 2) ) );
            if (erfc_denom < 0){  // our filtering scale has become too small
                break;
            }
            
            // renormalize the collapse fraction so that the mean matches ST,
            // since we are using the evolved (non-linear) density field
            for (x=0; x<HII_DIM; x++){
                for (y=0; y<HII_DIM; y++){
                    for (z=0; z<HII_DIM; z++){
                        density_over_mean = 1.0 + *((float *)deltax_filtered + HII_R_FFT_INDEX(x,y,z));
                        erfc_num = (Deltac - (density_over_mean-1)) /  growth_factor;
                        f_coll += splined_erfc(erfc_num/erfc_denom);
                        Fcoll[HII_R_FFT_INDEX(x,y,z)] = splined_erfc(erfc_num/erfc_denom);
                    }
                }
            }
            f_coll /= (double) HII_TOT_NUM_PIXELS;
            ST_over_PS = mean_f_coll_st/f_coll;
		}

		/************  MAIN LOOP THROUGH THE BOX **************/

        // now lets scroll through the filtered box
        ave_xHI_xrays = ave_den = ave_fcoll = std_xrays = 0;
        ion_ct=0;
        for (x=0; x<HII_DIM; x++){
            for (y=0; y<HII_DIM; y++){
                for (z=0; z<HII_DIM; z++){
                    if (LAST_FILTER_STEP)
                        density_over_mean = 1.0 + *((float *)deltax_unfiltered + HII_R_FFT_INDEX(x,y,z));
                    else
                        density_over_mean = 1.0 + *((float *)deltax_filtered + HII_R_FFT_INDEX(x,y,z));
					
                    // check for aliasing which can occur for small R and small cell sizes,
                    // since we are using the analytic form of the window function for speed and simplicity
                    if (density_over_mean <= 0){
                    //	    fprintf(LOG, "WARNING: aliasing during filtering step produced density n/<n> of %06.2f at cell (%i, %i, %i)\n Setting to 0\n", density_over_mean, x,y,z);
                        density_over_mean = FRACT_FLOAT_ERR;
                    }
                    
                    f_coll = ST_over_PS * Fcoll[HII_R_FFT_INDEX(x,y,z)];
					
                    // adjust the denominator of the collapse fraction for the residual electron fraction in the neutral medium
                    xHI_from_xrays = 1;
					
                    // check if ionized!
                    if (f_coll > f_coll_crit){ // ionized
						
                        if (FIND_BUBBLE_ALGORITHM == 2) // center method
                            xH[HII_R_INDEX(x, y, z)] = 0;
						
                        else if (FIND_BUBBLE_ALGORITHM == 1) // sphere method
                            update_in_sphere(xH, HII_DIM, R/BOX_LEN, x/(HII_DIM+0.0), y/(HII_DIM+0.0), z/(HII_DIM+0.0));
						
                        else{
                            printf( "Incorrect choice of find bubble algorithm: %i\nAborting...", FIND_BUBBLE_ALGORITHM);
                            fflush(NULL);
                            z=HII_DIM;y=HII_DIM,x=HII_DIM;R=0;
                        }
                    }
					
                    // check if this is the last filtering step.
                    // if so, assign partial ionizations to those cells which aren't fully ionized
                    else if (LAST_FILTER_STEP && (xH[HII_R_INDEX(x, y, z)] > TINY)){
                        
                        f_coll = ST_over_PS * Fcoll[HII_R_FFT_INDEX(x,y,z)];
                        if (f_coll>1) f_coll=1;
                        ave_N_min_cell = f_coll * pixel_mass*density_over_mean / M_MIN; // ave # of M_MIN halos in cell
                        if (ave_N_min_cell < N_POISSON){
                            // the collapsed fraction is too small, lets add poisson scatter in the halo number
                            N_min_cell = (int) gsl_ran_poisson(r, ave_N_min_cell);
                            f_coll = N_min_cell * M_MIN / (pixel_mass*density_over_mean);
                        }
                        
                        if (f_coll>1) f_coll=1;
                        res_xH = xHI_from_xrays - f_coll * ION_EFF_FACTOR;
                        // and make sure fraction doesn't blow up for underdense pixels
                        if (res_xH < 0)
                            res_xH = 0;
                        else if (res_xH > 1)
                            res_xH = 1;
						
                        xH[HII_R_INDEX(x, y, z)] = res_xH;
                    } // end partial ionizations at last filtering step
                } // k
            } // j
        } // i
    
        global_step_xH = 0;
        for (ct=0; ct<HII_TOT_NUM_PIXELS; ct++){
            global_step_xH += xH[ct];
        }
        global_step_xH /= (float)HII_TOT_NUM_PIXELS;
        
        if(global_step_xH==0.0) {
            short_completely_ionised = 1;
            break;
        }
        
        R /= DELTA_R_HII_FACTOR;
	}
	
	// find the neutral fraction
	global_xH = 0;
    
    for (ct=0; ct<HII_TOT_NUM_PIXELS; ct++){
        global_xH += xH[ct];
    }
	global_xH /= (float)HII_TOT_NUM_PIXELS;
    
	// deallocate
	gsl_rng_free (r);

	/**** End of perform 'find_HII_bubbles.c' *****/

	/**** Perform 'delta_T.c' ******/

	/************  BEGIN INITIALIZATION ****************************/

    max = -1e3;
	min = 1e3;
	ave = 0;
	nonlin_ct=0;
    
    double *k_box_ave = (double *)malloc(sizeof(double)*NUM_BINS);
    double *p_box_ave = (double *)malloc(sizeof(double)*NUM_BINS);
    double *p_box_ave_poisson = (double *)malloc(sizeof(double)*NUM_BINS);
    
    double chi_squared;
    
    chi_squared = 0.0;
    
    if(PERFORM_PS==1) {
        
        T_rad = T_cmb*(1+REDSHIFT);
        H = hubble(REDSHIFT);
        const_factor = 27 * (OMb*hlittle*hlittle/0.023) *
        sqrt( (0.15/OMm/hlittle/hlittle) * (1+REDSHIFT)/10.0 );

        memcpy(deltax, deltax_unfiltered_original, sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);

        /************  END INITIALIZATION ****************************/

        // ok, lets fill the delta_T box; which will be the same size as the bubble box
        temp=0;
        
        for (i=0; i<HII_DIM; i++){
            for (j=0; j<HII_DIM; j++){
                for (k=0; k<HII_DIM; k++){
				
                    pixel_deltax = deltax[HII_R_FFT_INDEX(i,j,k)];
                    pixel_x_HI = xH[HII_R_INDEX(i,j,k)];
				
                    if (pixel_x_HI > TINY){
                        temp = pixel_deltax;
                    }
				
                    delta_T[HII_R_INDEX(i,j,k)] = const_factor*pixel_x_HI*(1+pixel_deltax);
								
                    if (max < delta_T[HII_R_INDEX(i,j,k)]){ max = delta_T[HII_R_INDEX(i,j,k)];}
                    if (min > delta_T[HII_R_INDEX(i,j,k)]){ min = delta_T[HII_R_INDEX(i,j,k)];}
                    ave += delta_T[HII_R_INDEX(i,j,k)];
                }
            }
        }
        ave /= (float)HII_TOT_NUM_PIXELS;
			
        // now write out the delta_T box
        if (T_USE_VELOCITIES){
            max = -1;
            min = 1e3;
            ave = 0;
            
            // let's take the derivative in k-space
            plan = fftwf_plan_dft_r2c_3d(HII_DIM, HII_DIM, HII_DIM, (float *)v, (fftwf_complex *)v, FFTW_ESTIMATE);
            fftwf_execute(plan);
            fftwf_destroy_plan(plan);

            for (n_x=0; n_x<HII_DIM; n_x++){
                if (n_x>HII_MIDDLE)
                    k_x =(n_x-HII_DIM) * DELTA_K;  // wrap around for FFT convention
                else
                    k_x = n_x * DELTA_K;
			
                for (n_y=0; n_y<HII_DIM; n_y++){
                    if (n_y>HII_MIDDLE)
                        k_y =(n_y-HII_DIM) * DELTA_K;
                    else
                        k_y = n_y * DELTA_K;
            
                    for (n_z=0; n_z<=HII_MIDDLE; n_z++){
                        k_z = n_z * DELTA_K;
					
                        // take partial deriavative along the line of sight
                        switch(VELOCITY_COMPONENT){
                            case 1:
                                *((fftwf_complex *) v + HII_C_INDEX(n_x,n_y,n_z)) *= k_x*I/(float)HII_TOT_NUM_PIXELS;
                                break;
                            case 3:
                                *((fftwf_complex *) v + HII_C_INDEX(n_x,n_y,n_z)) *= k_z*I/(float)HII_TOT_NUM_PIXELS;
                                break;
                            default:
                                *((fftwf_complex *) v + HII_C_INDEX(n_x,n_y,n_z)) *= k_y*I/(float)HII_TOT_NUM_PIXELS;
                        }
                    }
                }
            }
            
            plan = fftwf_plan_dft_c2r_3d(HII_DIM, HII_DIM, HII_DIM, (fftwf_complex *)v, (float *)v, FFTW_ESTIMATE);
            fftwf_execute(plan);
            fftwf_destroy_plan(plan);
		
            // now add the velocity correction to the delta_T maps
            max_v_deriv = fabs(MAX_DVDR*H);
            for (i=0; i<HII_DIM; i++){
                for (j=0; j<HII_DIM; j++){
                    for (k=0; k<HII_DIM; k++){
					
                        dvdx = v[HII_R_FFT_INDEX(i,j,k)];
					
                        // set maximum allowed gradient for this linear approximation
                        if (fabs(dvdx) > max_v_deriv){
                            if (dvdx < 0) dvdx = -max_v_deriv;
                            else dvdx = max_v_deriv;
//                               nonlin_ct++;
                        }
					
                        delta_T[HII_R_INDEX(i,j,k)] /= (dvdx/H + 1.0);
					
                        if (max < delta_T[HII_R_INDEX(i,j,k)]){
                            maxi = i;
                            maxj = j;
                            maxk = k;
                            maxdvdx = dvdx;
                            max = delta_T[HII_R_INDEX(i,j,k)];
                        }
                        if (min > delta_T[HII_R_INDEX(i,j,k)]){
                            mini = i;
                            minj = j;
                            mink = k;
                            mindvdx = dvdx;
                            min = delta_T[HII_R_INDEX(i,j,k)];
                        }
					
                        ave += delta_T[HII_R_INDEX(i,j,k)];
                    }
                }
            }
            ave /= (HII_TOT_NUM_PIXELS+0.0);
        }

        /******  PRINT OUT THE POWERSPECTRUM  *********/

        for (ct=0; ct<NUM_BINS; ct++){
            p_box[ct] = k_ave[ct] = 0;
            in_bin_ct[ct] = 0;
        }
		
        // fill-up the real-space of the deldel box
        for (i=0; i<HII_DIM; i++){
            for (j=0; j<HII_DIM; j++){
                for (k=0; k<HII_DIM; k++){
                    *((float *)deldel_T + HII_R_FFT_INDEX(i,j,k)) = (delta_T[HII_R_INDEX(i,j,k)]/ave - 1)*VOLUME/(HII_TOT_NUM_PIXELS+0.0);
                    if (DIMENSIONAL_T_POWER_SPEC){
                        *((float *)deldel_T + HII_R_FFT_INDEX(i,j,k)) *= ave;
                    }
                    // Note: we include the V/N factor for the scaling after the fft
                }
            }
        }
        
        // transform to k-space
        plan = fftwf_plan_dft_r2c_3d(HII_DIM, HII_DIM, HII_DIM, (float *)deldel_T, (fftwf_complex *)deldel_T, FFTW_ESTIMATE);
        fftwf_execute(plan);
        fftwf_destroy_plan(plan);
	
        // now construct the power spectrum file

        for (n_x=0; n_x<HII_DIM; n_x++){
            if (n_x>HII_MIDDLE)
                k_x =(n_x-HII_DIM) * DELTA_K;  // wrap around for FFT convention
            else
                k_x = n_x * DELTA_K;
		
            for (n_y=0; n_y<HII_DIM; n_y++){
                if (n_y>HII_MIDDLE)
                    k_y =(n_y-HII_DIM) * DELTA_K;
                else
                    k_y = n_y * DELTA_K;
			
                for (n_z=0; n_z<=HII_MIDDLE; n_z++){
                    k_z = n_z * DELTA_K;
				
                    k_mag = sqrt(k_x*k_x + k_y*k_y + k_z*k_z);
				
                    // now go through the k bins and update
                    ct = 0;
                    k_floor = 0;
                    k_ceil = k_first_bin_ceil;
                    while (k_ceil < k_max){
                        // check if we fal in this bin
                        if ((k_mag>=k_floor) && (k_mag < k_ceil)){
                            in_bin_ct[ct]++;
                            p_box[ct] += pow(k_mag,3)*pow(cabs(deldel_T[HII_C_INDEX(n_x, n_y, n_z)]), 2)/(2.0*PI*PI*VOLUME);
                            // note the 1/VOLUME factor, which turns this into a power density in k-space
						
                            k_ave[ct] += k_mag;
                            break;
                        }
					
                        ct++;
                        k_floor=k_ceil;
                        k_ceil*=k_factor;
                    }
                }
            }
        } // end looping through k box
        
        // From the input power spectrum data (observation, error and k-cuts) compute the chi-squared for this 21 cm PS
        for (ct=1; ct<NUM_BINS; ct++){
            if (in_bin_ct[ct]>0) {
                k_box_ave[ct] = k_ave[ct]/(in_bin_ct[ct]+0.0);
                p_box_ave[ct] = p_box[ct]/(in_bin_ct[ct]+0.0);
                p_box_ave_poisson[ct] = p_box[ct]/(in_bin_ct[ct]+0.0)/sqrt(in_bin_ct[ct]+0.0);
            
                if ((k_box_ave[ct] > ForeGroundCut) && (k_box_ave[ct] < ShotNoiseCut)) {
                    chi_squared = chi_squared + ( p_box_ave[ct] - MockObservation[ct-1] )*( p_box_ave[ct] - MockObservation[ct-1] )/pow(( sqrt(TelescopeSensitivity[ct-1]*TelescopeSensitivity[ct-1] + ModellingUncertainty*ModellingUncertainty*p_box_ave[ct]*p_box_ave[ct]) ),2.);
                }
            
            }
        }
        
        
    }
    
	/****** END POWER SPECTRUM STUFF   ************/
	
	/**** End of perform 'delta_T.c' *****/

    DataToReturn.y = (float)global_xH;
    DataToReturn.x = (float)chi_squared;
   
    destroy_21cmMC_arrays();
    
    return DataToReturn;
}
