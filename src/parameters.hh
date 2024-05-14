/*
 * CLARK, CLAssifier based on Reduced K-mers.
 */

/*
   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

   Copyright @The Regents of the University of California. All rights reserved.

*/
/*
 * @author: Rachid Ounit, Ph.D.
 * @project: CLARK, Metagenomic and Genomic Sequences Classification project.
 * @note: C++ IMPLEMENTATION supported on latest Linux and Mac OS.
 *
 */


#ifndef PARAMETERS_HH
#define PARAMETERS_HH

#define VERSION "1.2.6.1"
#define YEARS "2013-2019"

#define SB              4       
#define DBCTRESH	4
#define LHTSIZE 	57777779
#define HTSIZE  	1610612741
#define NBN		4
#define SFACTORMAX 	30
#define MAXTS		65535
#define MAXTSSM		16383
#define WEIGHT		31
#define LENGTH          31

typedef uint64_t      T64;
typedef uint32_t      T32;
typedef uint16_t      T16;

#endif
