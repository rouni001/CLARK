/*
 *	CLARK, CLAssifier based on Reduced K-mers.
 */

/*
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *    Copyright 2013-2019, Rachid Ounit <clark.ucr.help at gmail.com>
 *                                   
 */

/*
 *    @author: Rachid Ounit, Ph.D Candidate.
 *    @project: CLARK, Metagenomic and Genomic Sequences Classification project.
 *    @note: C++ IMPLEMENTATION supported on latest Linux and Mac OS.
 *   
 */

#include <cstdlib>
#include <iostream>
#include <vector>
#include <cstring>
#include <iomanip>
#include <stdint.h>
#include <time.h>
#include "./file.hh"
using namespace std;

int main(int argc, const char** argv)
{
	if (argc != 3)
	{
		std::cerr << argv[0] << " <multi-fasta file> <foldername>"<< std::endl;
		exit(1);
	}
	FILE * fd = fopen(argv[1], "r");
	if (fd == NULL)
	{
		std::cerr << "Failed to open file: " << argv[1] << std::endl;
                exit(1);
	}
	
	string line;
	vector<char> sep;
	vector<string> ele;
	sep.push_back(' ');
	sep.push_back('\t');
	sep.push_back('>');

	FILE * sfd = NULL;
	char tab[10000];
	int t = 1;
	srand(time(NULL));
	int r = rand();
	while (getLineFromFile(fd, line))
	{
		if (line[0] =='>')
		{
			if (sfd != NULL)
			{	fclose(sfd);	}
			ele.clear();
			getElementsFromLine(line, sep, ele);
			sprintf(tab,"%s/%s.%i.%i.fa", argv[2],ele[0].c_str(), r, t);
			sfd = fopen(tab,"w");
			t++;
		}
		fprintf(sfd, "%s\n", line.c_str());
	}
	fclose(fd);
	return 0;
}

