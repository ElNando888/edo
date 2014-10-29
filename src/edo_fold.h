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



typedef struct _lab_config {
    int	  bclen;
    const char*  tail5p_ss;
    const char*  hairpin_ss;
    const char*  tail3p_ss;
    const char*  tail5p_nts;
    const char*  hairpin_nts;
    const char*  tail3p_nts;
    // forbidden subsequences
    const char** forbidden;
    // various precomputed offsets
    int   ofs1, ofs2, ofs3, ofs4;
    // scores
    double s_max, s_scale;
} lab_config;


extern int verbose;
extern lab_config* cnf;
extern int num_scored;

bool load_config( char* filename );
void load_fold_params( char* filename );
bool is_legal( char* seq );
double score_seq( int s, char* seq, int ip );
