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


#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <limits.h>
#include <unistd.h>
#include <omp.h>

#include <map>
#include <set>
#include <vector>

#include "edo_fold.h"
#include "edo_cmdl.h"

#define BASES "AUGC"
#define SEPARATORS "_-+='*"


int verbose = 0;

float epsilon = 0.01;

int num_seq;
int num_bc;
int num_loc;
int t = 0;

float**        prices[2];    // [][num_seq][num_bc]
int**          bidders[2];   // [][num_seq][num_bc]
long*          alpha[2];     // [][num_seq]
long long*     locations[2]; // [][num_seq]

typedef std::map<long long, float> val_map;
typedef val_map::iterator          val_map_it;
typedef std::vector<val_map>       val_map_array;

char**         sequences = NULL;    //   [num_seq]
int*           hip = NULL;          //   [num_seq]
               // as in Hairpin Insertion Point
val_map_array  values;

char**         seqs = NULL;
int            nseqs = 0;

std::set<long> tabu;

int            num_scored = 0;




void init( void )
{
    int s;
    for( s = 0, num_bc = 1; s < cnf->bclen; s++ ) num_bc *= 4;
    for( s = 0, num_loc = 1; s < cnf->bclen; s++ ) num_loc *= 6;

    int t;
    for( t = 0; t < 2; t++ ) {
        prices[t] = ( float** ) calloc( num_seq, sizeof( float* ) );
        for( s = 0; s < num_seq; s++ )
            prices[t][s] = ( float* ) calloc( num_bc, sizeof( float ) );
        bidders[t] = ( int** ) calloc( num_seq, sizeof( int* ) );
        for( s = 0; s < num_seq; s++ )
            bidders[t][s] = ( int* ) calloc( num_bc, sizeof( int ) );
        alpha[t] = ( long* ) calloc( num_seq, sizeof( long ) );
        locations[t] = ( long long* ) calloc( num_seq, sizeof( long long ) );
    }

    sequences = (char**) calloc( num_seq, sizeof(char*) );
    values.resize( num_seq );

    srand48( time( NULL ) );
}

#if 0
bool neighbor( int i, int j )
{
    int k;
    unsigned int a = locations[t][i];
    unsigned int b = locations[t][j];
    int d = 0;
    for( k = 0; k < 7; k++ ) {
        if( ( a & 0xF ) != ( b & 0xF ) ) {
            d++;
            if( d > 1 ) return false;
        }
        a >>= 4;
        b >>= 4;
    }
    return true;
}
#endif


int int_urn( int from, int to )
{
    return ( ( (int) (drand48()*(to-from+1)) ) + from );
}
 
// Internal bitmap representation
//
// MSB  LSB  nucleobase    MSB    LSB   5'  3'
//    00        C             0011        CG
//    01        U             1100        GC
//    10        A             1001        AU
//    11        G             0110        UA
//                            0111        UG
//                            1101        GU
//
//
// MSB  (64 bits => up to 16 base pairs) LSB
//
//                ---++---++---++---++--  5'
//              ... -dddd-cccc-bbbb-aaaa
//                  ---++---++---++---++  3'
//

char bases[] = "CUAG";
long long pmap[] = { 0x03, 0x0C, 0x09, 0x06, 0x07, 0x0D };


long location_to_barcode( long long loc )
{
    int k;
    long c = 0;
    for( k = 0; k < cnf->bclen; k++ ) {
        c <<= 2;
        c |= ( loc & 0x03 );
        loc >>= 4;
    }
    return c;
}

void location_to_bcseq( long long loc, char* s1, char* s2 )
{
    int k;
    for( k = 0; k < cnf->bclen; k++ ) {
        s1[k  ] = bases[( loc & 0x0C )>>2];
        s2[6-k] = bases[loc & 0x03];
        loc >>= 4;
    }
}

long bcseq_to_barcode( char* seg )
{
    int k;
    long c = 0;
    for( k = 0; k < cnf->bclen; k++ ) {
        c <<= 2;
        c |= seg[k]=='C'? 0 : ( seg[k]=='U'? 1 : ( seg[k]=='A'? 2 : 3 ) );
    }
    return c;
}


bool is_tabu( char* seg )
{
    long c = bcseq_to_barcode( seg );
    return tabu.find( c ) != tabu.end();
}

void compute_value( int i, long long a )
{
    int len = 1 + strlen( sequences[i] ) + cnf->hp_tl;
    char* seq = (char*) malloc( len );
    sprintf( seq, "%s%.*s%s%s%s", cnf->tail5p_nts,
                                  hip[i], sequences[i],
                                  cnf->hairpin_nts, 
                                  sequences[i] + hip[i],
                                  cnf->tail3p_nts );
    location_to_bcseq( a, seq + cnf->tl5p + hip[i] + cnf->hp5p,
                          seq + cnf->tl5p + hip[i] + cnf->hp3p );
    values[i][a] = is_tabu( seq + cnf->tl5p + hip[i] + cnf->hp3p ) ?
                       0.0 : score_seq( i, seq, cnf->tl5p + hip[i] );
    free( seq );
}

void update_local_values( int i )
{
    int k;
    #pragma omp parallel for
    for( k = 0; k < cnf->bclen; k++ ) {
        int j;
        long long mask = 0x0F << ( k*4 );
        for( j = 0; j < 6; j++ ) {
            long long a = ( locations[t][i] & ~mask ) | ( pmap[j] << ( k*4 ) );
            if( values[i].find( a ) != values[i].end() ) continue;
            compute_value( i, a );
        }
    }
}


void do_round( int r, int* num_swap, float* delta )
{
    int i, cnt = 0;
    int this_num_swap;
    if( !num_swap ) num_swap = &this_num_swap;
    (*num_swap) = 0;
    float this_delta;
    if( !delta ) delta = &this_delta;
    (*delta) = 0.0;

    #pragma omp parallel for schedule(dynamic)
    for( i = 0; i < num_seq; i++ ) {
        int j;
        long k;
        int next_t = (t+1)%2;

        if( verbose )
        #pragma omp critical(console)
        {
            fprintf( stderr, "\rRound %3d, %d/%d", r, ++cnt, num_seq );
            fflush( stderr );
        }

        // step 1
        float* old_p = prices[t][i];
        float* new_p = prices[next_t][i];
        int* old_b = bidders[t][i];
        int* new_b = bidders[next_t][i];

        for( k = 0; k < num_bc; k++ ) {
            new_p[k] = old_p[k];
            new_b[k] = old_b[k];
        }

        for( j = 0; j < num_seq; j++ ) {
#if 0 // as long as the data is available in memory, why not use it all?
            if( (j == i) || !neighbor( i, j ) ) continue;
#endif
            for( k = 0; k < num_bc; k++ ) {
                if( ( prices[t][j][k] > new_p[k] ) 
                    || ( prices[t][j][k] == new_p[k] 
                         && bidders[t][j][k] > new_b[k] ) ) {
                    new_p[k] = prices[t][j][k];
                    new_b[k] = bidders[t][j][k];
                }
            }
        }


        k = alpha[t][i];
        if( values[i][locations[t][i]] < epsilon ) {
            new_b[k] = num_seq; // artificial nudge
        }

        // step 2
        if( old_p[k] <= new_p[k] && new_b[k] != i ) {
            // step 3
            #pragma omp atomic update
            (*num_swap)++;
            // 3.a update local values
            update_local_values( i );

            // 3.b choose a new assignment
            float best_g = -100.0;
            long long best_x = -1;
            long best_c = -1;
            val_map_it m;
            for( m = values[i].begin(); m != values[i].end(); m++ ) {
                long long x = m->first;
                long c = location_to_barcode( x );
                if( values[i][x] - new_p[c] > best_g ) {
                    best_g = values[i][x] - new_p[c];
                    best_x = x;
                    best_c = c;
                }
            }
            float next_g = -100.0;
            for( m = values[i].begin(); m != values[i].end(); m++ ) {
                long long x = m->first;
                long c = location_to_barcode( x );
                if( c == best_c ) continue;
                if( values[i][x] - new_p[c] > next_g ) {
                    next_g = values[i][x] - new_p[c];
                }
            }

            alpha[next_t][i] = best_c;
            locations[next_t][i] = best_x;
            new_b[best_c] = i;
            new_p[best_c] += ( best_g - next_g ) + epsilon;
        } else {
            // step 4
            alpha[next_t][i] = k;
            locations[next_t][i] = locations[t][i];
        }

        #pragma omp atomic update
        (*delta) += values[i][locations[next_t][i]] 
                    - prices[next_t][i][alpha[next_t][i]];
    } // for( i

    (*delta) /= num_seq;
}



int load_input( void )
{
    char* line = NULL;
    size_t len = 0;
    ssize_t r;
    do {
        r = getline( &line, &len, stdin );
        if( r < 0 ) break;
        int i, j, k;
        if( line[0] == '-' ) { // pre-selected barcode, put in exclusion list
            i = strcspn( line, BASES );
            j = strspn( line + i, BASES );
            if( j != cnf->bclen ) continue;
            line[i + j] = '\0';
            tabu.insert( bcseq_to_barcode( line+i ) );
            continue;
        }
        // targets
        i = strspn( line, BASES );
        j = strspn( line + i, SEPARATORS );
        if( j == 0 && i == 0 ) continue;
        k = 0;
        if( j ) {
            k = strspn( line + i + j, BASES );
        }
        line[i + j + k] = '\0';
        if( !is_legal( line ) ) {
            fprintf( stderr, "Warning, illegal sequence %s ignored.\n", line );
        } else {
            char* clean = strdup( line );
            if( j ) strcpy( clean + i, line + i + j );
            clean[i + k] = '\0';
            // add to the list
            nseqs++;
            seqs = (char**) realloc( seqs, nseqs * sizeof(char*) );
            seqs[nseqs - 1] = clean;
            hip = (int*) realloc( hip, nseqs * sizeof(int) );
            hip[nseqs - 1] = i;
        }
    } while( r >= 0 );

    free( line );
    return nseqs;
}


void eval_input( void )
{
    char* line = NULL;
    size_t len = 0;
    ssize_t r;
    do {
        r = getline( &line, &len, stdin );
        if( r < 0 ) break;
        int i = strspn( line, BASES );
        int j = strspn( line + i, SEPARATORS );
        if( j == 0 ) continue;
        int k = strspn( line + i + j, BASES );
        if( k != strlen( cnf->hairpin_nts ) ) continue;
        char* clean = strdup( line );
        strcpy( clean + i, line + i + j );
        int l = strspn( line + i + j + k, SEPARATORS );
        int m = 0;
        if( l ) {
            m = strspn( line + i + j + k + l, BASES );
            strcpy( clean + i + k, line + i + j + k + l );
        }
        clean[i + k + m] = '\0';
        double score = score_seq( 0, clean, i );
        fprintf( stdout, "%9.6f\t%s\n", score, clean );
        free( clean );
    } while( r >= 0 );
    free( line );
}


void auction_barcodes( void )
{
    init();

    int j;
    long k;
    for( j = 0; j < num_seq; j++ ) {
        sequences[j] = seqs[j];
    }
    t = 0;

    // choose an "origin" (see populating the arrays, below)
    long long origin = 0;
    for( j = 0; j < cnf->bclen; j++ ) {
        // try to make it a 'good' one: CG/GC closing, and 60% CG/GC pairs
        int v = int_urn( 0, j == 0 || j == cnf->bclen-1 ? 1 : 9 ) % 4;
        origin <<= 4;
        origin |= v==0 ? 0x3 : v==1 ? 0xC : v==2 ? 0x6 : 0x9;
    }
    long org_bc = location_to_barcode( origin );
    
    // populate arrays
    for( j = 0; j < num_seq; j++ ) {
        for( k = 0; k < num_bc; k++ ) {
            prices[0][j][k] = 1.0 + drand48();
            bidders[0][j][k] = j;
        }
        // all agents start from the same location, and are artificially
        // outbid, which forces them to explore the neighborhood and start
        // spreading in all directions
        locations[0][j] = origin;
        alpha[0][j] = org_bc;
        bidders[0][j][alpha[0][j]] = num_seq;
    }

    // do the loop thing...
    int r = 0;
    int swaps;
    do {
        float delta;
        do_round( ++r, &swaps, &delta );
        if( verbose ) 
            fprintf( stderr, "\rRound %3d: %8d PFs computed, "
                             "%5d reassignments\n", r, num_scored, swaps );
        t++;
        t %= 2;
    // ... until no one complains anymore
    } while( swaps > 0 );

    // output
    fprintf( stderr, "\n" );
    for( j = 0; j < num_seq; j++ ) {
        long long a = locations[0][j];
        int len = 1 + strlen( sequences[j] ) + cnf->hp_tl;
        char* seq = (char*) malloc( len );
        sprintf( seq, "%s%.*s_%s_%s%s", cnf->tail5p_nts, 
                                        hip[j], sequences[j],
                                        cnf->hairpin_nts, 
                                        sequences[j]+hip[j],
                                        cnf->tail3p_nts );
        location_to_bcseq( a, seq + cnf->tl5p + hip[j] + 1 + cnf->hp5p,
                              seq + cnf->tl5p + hip[j] + 1 + cnf->hp3p );
        if( verbose ) fprintf( stderr, "%s\t%6.3f\n", seq, values[j][a] );
        fprintf( stdout, "%s\n", seq );
        free( seq );
    }

}


int main( int argc, char** argv )
{
    edo_args_info args_info;
    
    if( edo_cmdline_parser( argc, argv, &args_info ) != 0 )
        exit( 1 );
    
    if( args_info.verbose_given ) verbose = 1;

    if( args_info.lab_settings_given ) {
        if( !load_config( args_info.lab_settings_arg ) ) {
            exit( 1 );
        }
    }

    if( args_info.paramFile_given ) {
        load_fold_params( args_info.paramFile_arg );
    }

    if( args_info.eval_given ) {
        eval_input();
        exit( 0 );
    }

    edo_cmdline_parser_free( &args_info );

    load_input();
    if( nseqs == 0 ) {
        fprintf( stderr, "Not much work to do..." );
        return 0;
    }

    num_seq = nseqs;
    if( verbose )
        fprintf( stderr, "Loaded %d sequences, %ld barcode(s)"
                         " are reserved.\n\n", num_seq, tabu.size() );

    auction_barcodes();
    return 0;
}

