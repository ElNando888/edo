// Copyright (C) 2014 Fernando Portela <nando8888@gmail.com>
//  
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software 
// Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.


#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

extern "C" {
#include "RNAstruct.h"
#include "fold_vars.h"
#include "fold.h"
#include "params.h"
#include "part_func.h"
#include "utils.h"
#include "convert_epars.h"
#include "read_epars.h"
#include "MEA.h"
}

#include "edo_fold.h"


#define BC_LEN        7

#define TAIL5P_SS     ".."
#define HAIRPIN_SS    "(((((((....)))))))."
#define TAIL3P_SS     "...................."

#define TAIL5P_NTS    "GG"
#define HAIRPIN_NTS   ".......UUCG.......A"
#define TAIL3P_NTS    "AAAGAAACAACAACAACAAC"



static lab_config daslab_cnf = {
    BC_LEN,
    TAIL5P_SS,
    HAIRPIN_SS,
    TAIL3P_SS,
    TAIL5P_NTS,
    HAIRPIN_NTS,
    TAIL3P_NTS,
    39, 0, 11, 2
};

lab_config* cnf = &daslab_cnf;



bool is_legal( char* seq )
{
#if 0
    if( strstr( seq, "AAAAA" ) ) return false;
#endif
    if( strstr( seq, "CCCC" ) ) return false;
    if( strstr( seq, "GGGG" ) ) return false;

    return true;
}


double score_seq( int s, char* seq )
{
    double score = 0.0;
    if( !is_legal( seq ) ) return score;

    int length = strlen( seq );
    int* ix = get_iindx( length );

    double betaScale = 1.0;
    double kT = ( betaScale*( ( temperature+K0 )*GASCONST ) )/1000.; /* in Kcal */
    model_detailsT md;
    set_model_details( &md );

    char* secstr = strdup( seq );
    secstr[0] = 0;
    fold_constrained = 0;
    paramT* params = get_scaled_parameters( temperature, md );
    double min_en = fold_par( seq, secstr, params, fold_constrained, 0 );
    if( strncmp( secstr + length - cnf->ofs1, 
                 cnf->hairpin_ss, strlen( cnf->hairpin_ss ) ) != 0
        || strcmp( secstr + length - strlen( cnf->tail3p_ss ),
                   cnf->tail3p_ss ) != 0 ) {
        free( params );
        free( secstr );
        free( ix );
        return score;
    }
    #pragma omp atomic update
    num_scored++;

    double pf_scale = exp( -( 1.07*min_en )/kT/length );
    pf_paramT* pf_params = get_boltzmann_factors( temperature, betaScale, md, pf_scale );
    // Either patch fold_vars.h by inserting following at line 166
    //
    // #ifdef _OPENMP
    // #pragma omp threadprivate(iindx)
    // #endif
    //
    // or uncomment this pragma below
    //
    // #pragma omp critical(pf_fold)
    double e = pf_fold_par( seq, NULL, pf_params, 1, fold_constrained, 0 );
    FLT_OR_DBL* ppm = export_bppm();

#define pr_ij(i,j) (i == j? 0.0 : (i < j ? ppm[ix[i]-j] : ppm[ix[j]-i]))

    score = 10.0;
    int i, o;
    for( i = 1; i <= length; i++ ) {
        for( o = 1; o <= 19; o++ ) {
            int j = length - 39 + o;
            double v = pr_ij( i, j );
            score -= v * ( 1.0 - v );
        }
    }

    free( pf_params );
    free( params );
    free( secstr );
    free( ix );
    return score;
}
