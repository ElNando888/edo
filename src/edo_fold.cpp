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
#include "libconfig.h"

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

static const char* daslab_forbidden[] = { "CCCC", "GGGG", NULL };

static lab_config daslab_cnf = {
    BC_LEN,
    TAIL5P_SS,
    HAIRPIN_SS,
    TAIL3P_SS,
    TAIL5P_NTS,
    HAIRPIN_NTS,
    TAIL3P_NTS,
    daslab_forbidden,
    39, 0, 11, 2,
    10., 1.
};

lab_config* cnf = &daslab_cnf;



bool load_config( char* filename )
{
    bool ok = true;
    config_t cfg;
    config_init( &cfg );
    if( config_read_file( &cfg, filename ) == CONFIG_TRUE ) {
        lab_config* lab = (lab_config*) calloc( 1, sizeof( lab_config ) );
        const char* buf;
        if( ok && config_lookup_string( &cfg, "settings.secstr.tail5p",
                                        &buf ) == CONFIG_TRUE ) {
            lab->tail5p_ss = strdup( buf );
        } else {
            ok = false;
        }
        if( ok && config_lookup_string( &cfg, "settings.secstr.hairpin",
                                        &buf ) == CONFIG_TRUE ) {
            lab->hairpin_ss = strdup( buf );
        } else {
            ok = false;
        }
        if( ok && config_lookup_string( &cfg, "settings.secstr.tail3p",
                                        &buf ) == CONFIG_TRUE ) {
            lab->tail3p_ss = strdup( buf );
        } else {
            ok = false;
        }
        if( ok && config_lookup_string( &cfg, "settings.sequences.tail5p",
                                        &buf ) == CONFIG_TRUE ) {
            lab->tail5p_nts = strdup( buf );
        } else {
            ok = false;
        }
        if( ok && config_lookup_string( &cfg, "settings.sequences.hairpin",
                                        &buf ) == CONFIG_TRUE ) {
            lab->hairpin_nts = strdup( buf );
        } else {
            ok = false;
        }
        if( ok && config_lookup_string( &cfg, "settings.sequences.tail3p",
                                        &buf ) == CONFIG_TRUE ) {
            lab->tail3p_nts = strdup( buf );
        } else {
            ok = false;
        }
        if( ok ) {
            // consistency checks
            ok &= strlen( lab->tail5p_ss ) == strlen( lab->tail5p_nts );
            ok &= strlen( lab->hairpin_ss ) == strlen( lab->hairpin_nts );
            ok &= strlen( lab->tail3p_ss ) == strlen( lab->tail3p_nts );
        }
        if( ok ) {
            char key[128];
            int i;
            for( i = 0; /* forever */; i++ ) {
                sprintf( key, "settings.forbidden.[%d]", i );
                if( config_lookup_string( &cfg, key, &buf ) != CONFIG_TRUE) break;
                lab->forbidden = (const char**) realloc( lab->forbidden,
                                                         (i+1) * sizeof( char* ) );
                lab->forbidden[i] = strdup( buf );
            }
            if( lab->forbidden ) {
                lab->forbidden = (const char**) realloc( lab->forbidden,
                                                         (i+1) * sizeof( char* ) );
                lab->forbidden[i] = NULL;
            }
        }
        if( ok ) {
            // populate the other fields in the lab settings
            lab->ofs1 = strlen( lab->hairpin_ss ) + strlen( lab->tail3p_ss );
            lab->ofs2 = strcspn( lab->hairpin_ss, "(" );
            lab->bclen = strspn( lab->hairpin_ss + lab->ofs2, "(" );
            int j = lab->ofs2 + lab->bclen;
            lab->ofs3 = j + strcspn( lab->hairpin_ss + j, ")" );
            lab->ofs4 = strlen( lab->tail5p_ss );
            lab->s_max = round( strlen( lab->hairpin_ss ) + .9 ) * 0.5;
            lab->s_scale = 10.0 / lab->s_max;

            if( verbose ) {
                fprintf( stderr, "%s / / %s / / %s\n", lab->tail5p_nts,
                                                       lab->hairpin_nts,
                                                       lab->tail3p_nts );
                fprintf( stderr, "%s / / %s / / %s\n", lab->tail5p_ss,
                                                       lab->hairpin_ss,
                                                       lab->tail3p_ss );
                fprintf( stderr, "[%d,%d,%d,%d]-[%4.1f,%6.4f]\n",
                                 lab->ofs1, lab->ofs2, lab->ofs3, lab->ofs4,
                                 lab->s_max, lab->s_scale );
                fprintf( stderr, "Configuration successfully loaded from %s\n",
                                 filename );
            }

            // activate
            cnf = lab;

        } else {
            free( lab );
        }
        
    } else {
        ok = false;
        fprintf( stderr, "%s:%d: %s\n", filename, config_error_line( &cfg ),
                                                  config_error_text( &cfg ) );
    }
    
    config_destroy( &cfg );
    return ok;
}


void load_fold_params( char* filename ) {
    read_parameter_file( filename );
}


bool is_legal( char* seq )
{
#if 0
// TODO: these rules are pipeline-specific...
    if( strstr( seq, "AAAAA" ) ) return false;
    if( strstr( seq, "CCCC" ) ) return false;
    if( strstr( seq, "GGGG" ) ) return false;
#endif

    if( cnf->forbidden ) {
        const char** p = cnf->forbidden;
        while( *p ) {
            if( strstr( seq, *p ) ) return false;
            p++;
        }
    }

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

    score = cnf->s_max;
    int i, o;
    for( i = 1; i <= length; i++ ) {
        for( o = 1; o <= strlen( cnf->hairpin_ss ); o++ ) {
            int j = length - cnf->ofs1 + o;
            double v = pr_ij( i, j );
            score -= v * ( 1.0 - v );
        }
    }
    score *= cnf->s_scale;

    free( pf_params );
    free( params );
    free( secstr );
    free( ix );
    return score;
}
