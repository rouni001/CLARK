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

   Copyright 2013-2019, Rachid Ounit <clark.ucr.help at gmail.com>
 */

/*
 * @author: Rachid Ounit, Ph.D Candidate.
 * @project: CLARK, Metagenomic and Genomic Sequences Classification project.
 * @note: C++ IMPLEMENTATION supported on latest Linux and Mac OS.
 *
 */

#include <cstring>
#include <vector>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include "./kmersConversion.hh"
#include "./contiguousToSpaced_hh.hh"
#include "./parameters.hh"

using namespace std;

int main(int argc, const char** argv)
{
	if (argc != 5)
	{
		cerr << argv[0] << " <DatabaseFilename: TSK file> <weight> <len> <setting:T295,...>"<< endl;
		return 1;
	}
	size_t w  = atoi(argv[2]);
	size_t len = atoi(argv[3]);
	string setting(argv[4]);
	size_t t_b = log(HTSIZE)/log(4.0);
        size_t max16 = t_b + 8;
        size_t max32 = t_b + 16;

	string file(argv[1]);
	string sfile = file.substr(0,file.size()-3);
	if (max16 >= w)
	{
		contiguousTospaced<T32,T16> converter(w, setting, len);
		converter.populate(sfile.c_str());
		converter.write(sfile.c_str());
		exit(0);
	}
	if (max32 >= w)
        {
                contiguousTospaced<T32,T32> converter(w, setting, len);
                converter.populate(sfile.c_str());
                converter.write(sfile.c_str());
                exit(0);
        }
	return 1;
}

