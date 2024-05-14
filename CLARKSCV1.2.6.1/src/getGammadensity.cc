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

#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <cstring>
#include <vector>
#include "./file.hh"
using namespace std;

#define T 50

void hbar(const size_t& _scale, const double& _ratio)
{
	size_t cpt = (size_t) (((double) _scale)*_ratio);
	for(size_t t = 0; t < cpt ; t++)
		cout << "*";
	cout << endl;
}

#define MINVAL 0.0

int main(int argc, char** argv)
{
	string line;
	vector<char> sep;
        vector<string> ele;
	sep.push_back(',');
	sep.push_back('\r');
	vector<size_t> counts(T+1,0);
	double total = 0;
	
	for(size_t t = 1; t < argc; t++)
	{
		FILE * fd = fopen(argv[t],"r");
		if (fd == NULL)
		{
			cerr << "Failed to open " << argv[t] << ". This file will be ignored for the calculations." << endl;
			continue;
		}
		std::cerr << "\rProcessing file: " << argv[t]<< "\t\t\t" ;
		getLineFromFile(fd,line);
		ele.clear();
		getElementsFromLine(line, sep, ele);
		if (ele.size() < 6 || ele[ele.size()-6] != "Gamma")
		{
			cerr << argv[t] << " does not contain Gamma scores. This file will be ignored for the calculations." << endl;
                        fclose(fd);
			continue;
		}
		while (getLineFromFile(fd, line))
		{
			ele.clear();
                	getElementsFromLine(line, sep, ele);
			if (ele.size() < 6)
			{	std::cerr << "Error with assignment: "<< line << std::endl;	
				continue;
			}
			double f = atof(ele[ele.size()-6].c_str());
			if (f > 0)
			{
				size_t idx = (size_t) (((f - MINVAL)*((double) T)));
				if (idx < 0)
				{	std::cerr << "Error with: "<< line << std::endl; exit(1);	}
				if (idx >= T)
				{	idx = T;	}
				counts[idx]++;
				total++;
			}
		}
		fclose(fd);
	}
	if (total == 0)
	{
		cerr << "No data found." << endl;
		exit(1);
	}
        cerr << total << " assignments with Gamma score found."  << endl;                                                                                  
	cout << "Interval     \t    Count  \t Density  \tHistogram"<< endl;
	for(double t = 0; t < T+1; t++)
	{
		if (t < T)
		{
			double r = ((double) counts[t])/((double) total);
			//cout << "["<< t/(T)+MINVAL <<"," << (t+1)/(T)+MINVAL<<"]\t"<< counts[t] << "  \t" ;
			//cout << setprecision(2) << fixed << r*100.0 << "\t";
			printf("[%4.2f,%4.2f[\t%10.f\t%10.2f\t", t/(T)+MINVAL,  (t+1)/(T)+MINVAL, (double) counts[t],r*100.0);
			hbar(100,r);
			continue;
		}
		double r = ((double) counts[t])/((double) total);
		//cout << setprecision(2) << fixed;
                //cout << "    [>=1]   \t"<< counts[t] << " \t" ;
		//cout << setprecision(2) << fixed << r*100.0 << "\t";
		printf("  [>=1]       \t%10.f\t%10.2f\t", (double) counts[t],r*100.0);
                hbar(100,r);
	}
	return 0;
}
