/*
 * 	CLARK, CLAssifier based on Reduced K-mers.
 */

/*
 * 	This program is free software: you can redistribute it and/or modify
 *   	it under the terms of the GNU General Public License as published by
 *      the Free Software Foundation, either version 3 of the License, or
 *      (at your option) any later version.
 *   
 *      This program is distributed in the hope that it will be useful,
 *      but WITHOUT ANY WARRANTY; without even the implied warranty of
 *      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *      GNU General Public License for more details.
 *           
 *      You should have received a copy of the GNU General Public License
 *      along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *              
 *      Copyright @The Regents of the University of California. All rights reserved.
 *
 */
/*
 * 		@author: Rachid Ounit, Ph.D.
 *      @project: CLARK, Metagenomic and Genomic Sequences Classification project.
 *      @note: C++ IMPLEMENTATION supported on latest Linux and Mac OS.
 *     
 */


#include <cstdlib>
#include <fstream>
#include <vector>
#include <string>
#include <iostream>

#include "file.hh"
using namespace std;

int main(int argc, char** argv)
{
	if (argc != 2)
	{
		cerr << "Usage: " << argv[0] << " <File: assembly_summary.txt> " << endl;
		exit(1);
	}

	FILE * fd =  fopen(argv[1],"r");
	if (fd == NULL)
	{	cerr << "Failed to open " << argv[1] << endl; 
		exit(1);
	}
	string line = "";
	vector<char> sep;
	sep.push_back('/');
	vector<string> ele;
	while (getLineFromFile(fd, line))
	{
		ele.clear();
		getElementsFromLine(line,sep,ele);
		if (ele.size() < 1)
		{	continue; }
		cout << "wget " << line << "/" << ele[ele.size()-1] <<"_genomic.fna.gz" << endl;
	}
	fclose(fd);
	return 0;
}
