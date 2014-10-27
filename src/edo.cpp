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

#define dbg 0


float epsilon = 0.01;

int num_seq;
int num_bc;
int num_loc;
int t = 0;

float**        prices[2];    // [][num_seq][num_bc]
int**          bidders[2];   // [][num_seq][num_bc]
int*           alpha[2];     // [][num_seq]
unsigned int*  locations[2]; // [][num_seq]

typedef std::map<int, float> val_map;
typedef val_map::iterator    val_map_it;
typedef std::vector<val_map> val_map_array;

char**         sequences;    //   [num_seq]
val_map_array  values;

char**         seqs = NULL;
int            nseqs = 0;

std::set<int>  tabu;

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
        alpha[t] = ( int* ) calloc( num_seq, sizeof( int ) );
        locations[t] = ( unsigned int* ) calloc( num_seq, sizeof( unsigned int ) );
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

// C U A G
char bases[] = "CUAG";
int pmap[] = { 0x03, 0x0C, 0x09, 0x06, 0x07, 0x0D };


int location_to_barcode( int loc )
{
    int k;
    int c = 0;
    for( k = 0; k < cnf->bclen; k++ ) {
        c <<= 2;
        c |= ( loc & 0x03 );
        loc >>= 4;
    }
    return c;
}

void location_to_bcseq( int loc, char* s1, char* s2 )
{
    int k;
    for( k = 0; k < cnf->bclen; k++ ) {
        s1[k  ] = bases[( loc & 0x0C )>>2];
        s2[6-k] = bases[loc & 0x03];
        loc >>= 4;
    }
}

int bcseq_to_barcode( char* seq )
{
    int k;
    int c = 0;
    seq += strlen( seq ) - cnf->ofs1 + cnf->ofs3;
    for( k = 0; k < cnf->bclen; k++ ) {
        c <<= 2;
        c |= seq[k]=='C'? 0 : ( seq[k]=='U'? 1 : ( seq[k]=='A'? 2 : 3 ) );
    }
    return c;
}


bool is_tabu( char* seq )
{
    int c = bcseq_to_barcode( seq );
    return tabu.find( c ) != tabu.end();
}

void compute_value( int i, int a )
{
    int ofs = strlen( sequences[i] );
    char* seq = (char*) malloc( 1 + ofs + cnf->ofs1 + cnf->ofs4 );
    sprintf( seq, "%s%s%s%s", cnf->tail5p_nts, sequences[i],
                              cnf->hairpin_nts, cnf->tail3p_nts );
    location_to_bcseq( a, seq + ofs + cnf->ofs2 + cnf->ofs4,
                          seq + ofs + cnf->ofs3 + cnf->ofs4 );
    values[i][a] = is_tabu( seq ) ? 0.0 : score_seq( i, seq );
    free( seq );
}

void update_local_values( int i )
{
    int k;
    #pragma omp parallel for
    for( k = 0; k < cnf->bclen; k++ ) {
        int j;
        unsigned int mask = 0x0F << ( k*4 );
        for( j = 0; j < 6; j++ ) {
            int a = ( locations[t][i] & ~mask ) | ( pmap[j] << ( k*4 ) );
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

    #pragma omp parallel for
    for( i = 0; i < num_seq; i++ ) {
        int j, k;
        int next_t = (t+1)%2;

        #pragma omp critical(console)
        if( !dbg ) {
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
            // if( (j == i) || !neighbor( i, j ) ) continue;
            for( k = 0; k < num_bc; k++ ) {
                if( ( prices[t][j][k] > new_p[k] ) || ( prices[t][j][k] == new_p[k] && bidders[t][j][k] > new_b[k] ) ) {
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
            int best_x = -1;
            int best_c = -1;
            val_map_it m;
            for( m = values[i].begin(); m != values[i].end(); m++ ) {
                int x = m->first;
                int c = location_to_barcode( x );
                if( values[i][x] - new_p[c] > best_g ) {
                    best_g = values[i][x] - new_p[c];
                    best_x = x;
                    best_c = c;
                }
            }
            float next_g = -100.0;
            for( m = values[i].begin(); m != values[i].end(); m++ ) {
                int x = m->first;
                int c = location_to_barcode( x );
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
        (*delta) += values[i][locations[next_t][i]] - prices[next_t][i][alpha[next_t][i]];
    } // for( i

    (*delta) /= num_seq;
}



int load_file( char* filename )
{
    FILE *f = fopen( filename, "rt" );
    if( f ) {
        char* line = NULL;
        size_t len = 0;
        ssize_t r;
        do {
            r = getline( &line, &len, f );
            if( r >= 0 ) {
                int l = strspn( line, "AUGC" );
                if( l >= cnf->bclen ) {
                    line[l] = 0;
                    if( l >= cnf->ofs1 + cnf->ofs4
                        && strncmp( line, cnf->tail5p_nts, strlen( cnf->tail5p_nts ) )==0
                        && strcmp( line+l-strlen( cnf->tail3p_nts ), cnf->tail3p_nts )==0 ) {
                        // has lab tails, so it's a player-reserved barcode
                        tabu.insert( bcseq_to_barcode( line ) );
                    } else {
                        if( !is_legal( line ) ) {
                            fprintf( stderr, "Warning, illegal sequence %s excluded.\n", line );
                        } else {
                            // add to the list
                            nseqs++;
                            seqs = (char**) realloc( seqs, nseqs * sizeof(char*) );
                            seqs[nseqs - 1] = strdup( line );
                        }
                    }
                }
            }
            free( line );
            line = NULL;
            len = 0;
        } while( r >= 0 );
        fclose( f );
    }

    return nseqs;
}


void eval_file( char* filename )
{
    FILE* fin = fopen( filename, "r" );
    if( fin == NULL ) return;
    char *outname = (char*) malloc( 1 + strlen( filename ) + 6 );
    sprintf( outname, "%s.score", filename );
    FILE* fout = fopen( outname, "w" );

    char* line = NULL;
    size_t len = 0;
    ssize_t r;
    do {
        r = getline( &line, &len, fin );
        if( r >= 0 ) {
            int l = strspn( line, "AUGC" );
            if( l >= 41 ) {
                line[l] = 0;
                double score = score_seq( 0, line );
                fprintf( fout, "%9.6f\t%s\n", score, line );
            }
        }
        free( line );
        line = NULL;
        len = 0;
    } while( r >= 0 );

    fclose( fout );
    free( outname );
    fclose( fin );
}


void auction_barcodes( char* basename )
{
    init();

    int j, k;
    for( j = 0; j < num_seq; j++ ) {
        sequences[j] = seqs[j];
    }
    t = 0;
    
    // populate arrays
    for( j = 0; j < num_seq; j++ ) {
        for( k = 0; k < num_bc; k++ ) {
            prices[0][j][k] = 1.0 + drand48();
            bidders[0][j][k] = j;
        }
        // locations[0][j] = index_to_location( (int)(drand48() * num_loc) );
        locations[0][j] = 0x3339CCC;
        alpha[0][j] = location_to_barcode( locations[0][j] );
        bidders[0][j][alpha[0][j]] = num_seq;
    }

    // do the loop thing
    int r = 0;
    int swaps;
    do {
        float delta;
        do_round( ++r, &swaps, &delta );
        fprintf( stderr, "\rRound %3d: %8d PFs computed, %5d reassignments\n", r, num_scored, swaps );
        t++;
        t %= 2;
    } while( swaps > 0 );

    // output
    char *outname = (char*) malloc( 1 + strlen( basename ) + 4 );
    sprintf( outname, "%s.out", basename );
    FILE* f = fopen( outname, "w" );

    fprintf( stderr, "\n" );
    for( j = 0; j < num_seq; j++ ) {
        int a = locations[0][j];
        int ofs = strlen( sequences[j] );
        char* seq = (char*) malloc( 1 + ofs + cnf->ofs1 + cnf->ofs4 );
        sprintf( seq, "%s%s%s%s", cnf->tail5p_nts, sequences[j],
                                  cnf->hairpin_nts, cnf->tail3p_nts );
        location_to_bcseq( a, seq + ofs + cnf->ofs2 + cnf->ofs4,
                              seq + ofs + cnf->ofs3 + cnf->ofs4 );
        printf( "%s\t%6.3f\n", seq, values[j][a] );
        if( f ) fprintf( f, "%s\n", seq );
        free( seq );
    }

    if( f ) fclose( f );
    free( outname );
}


int main( int argc, char** argv )
{
    if( argc < 2 ) return 1;

    if( strcmp( argv[1], "-e" )==0 ) {
        if( argc > 2 ) {
            eval_file( argv[2] );
        } else {
            fprintf( stderr, "erm... eval what?\n" );
        }
        return 0;
    }

    load_file( argv[1] );
    if( nseqs == 0 ) {
        fprintf( stderr, "Not much work to do..." );
        return 0;
    }

    num_seq = nseqs;
    fprintf( stderr, "Loaded %d sequences, %ld barcode(s) are reserved.\n\n", num_seq, tabu.size() );

    auction_barcodes( argv[1] );    
    
    return 0;
}

