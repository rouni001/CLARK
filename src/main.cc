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

#include<iostream>
#include<cstdlib>
#include<stdio.h>
#include<cmath>
#include<fstream>

#include "./CLARK_hh.hh"
#include "./parameters.hh"
#define MAXK 32

using namespace std;

void printUsage()	
{
	cout << endl;
	cout << "CLARK -- ``CLAssifier based on Reduced K-mers'', UCR CS&E " << endl ;
	cout << endl;
	cout << "./CLARK -k <kmerSize> -t <minFreqTarget> -o <minFreqObject> -T <fileTargets> -D <directoryDB/> -O <fileObjects> -R <fileResults> -m <mode> -n <numberofthreads> ...\n\n" << endl;
	cout << "Definitions of parameters (cf. README file, for details):" << endl;
	cout << endl;
	cout << "-k <kmerSize>,       \t k-mer length:\tinteger, >= 2 and <= 32 (for CLARK only). The default value is 31." << endl;
	cout << "-t <minFreqTarget>,  \t minimum of k-mer frequency in targets (for CLARK only):\tinteger, >=0." << endl;
	cout << "-o <minFreqtObject>, \t minimum of k-mer frequency in objects  (for CLARK only):\tinteger, >=0." << endl;
	cout << "-T <fileTargets>,    \t filename of the targets definition:\t text." << endl;
	cout << "-D <directoryDB/>,   \t directory name for the database (to load/save database files):\t text." << endl;
	cout << "-O <fileObjects>,    \t filename of objects (or list of objects):\t text." << endl;
	cout << "-P <file1> <file2>,  \t filenames of paired-end reads:\t texts." << endl;
	cout << "-R <fileResults>,    \t filename to store results (or corresponding list of results file):\t text.\n";
	cout << "-m <mode>,           \t mode of execution: 0 (full), 1 (default), 2 (express) and 3 (spectrum).\n";
	cout << "                     \t For CLARK-S, the full and default mode are the same.\n";
	cout << "-n <numberofthreads>,\t number of threads:\tinteger >= 1.\n";
	cout << "--long,              \t to indicate that the objects files contains very long/large sequences (e.g.,\n";
	cout << "                     \t long contigs from genome assembly, long sequencing reads from Nanopore or Pacbio, etc.)\n";
	cout << "                     \t and allocate more memory accordingly. This option is only for running the full mode or \n";
	cout << "                     \t running CLARK-S, in the case of sequences up to 50 Mbp.\n";
	cout << "--tsk,               \t to request a detailed creation of the database (target specific k-mers files). This option is no more supported." << endl;
	cout << "--ldm,               \t to request the loading of the database by memory mapped-file (in multithreaded mode, multiple parallel threads are requested)." << endl;
	cout << "--kso,               \t to request a preliminary k-spectrum analysis of each object (for mode 3 only)." << endl;
	cout << "--extended,          \t to request an extended output of the full mode (for CLARK only)." << endl;
	cout << "-g <iteration>,      \t gap or number of non-overlapping k-mers to pass for the database creation (for CLARK-l only). The default value is 4." << endl;
	cout << "-s <factor>,         \t sampling factor value in the default mode (for CLARK/CLARK-S only). " << endl;
	cout << endl;
	cout << "--help,              \t to print help/options." << endl;
	cout << "--version,           \t to print the version info." << endl;
	cout << "Version: " << VERSION << ", Contact/Help: Rachid Ounit, rouni001 at cs.ucr.edu" << endl;
	cout << endl;
	return;
}

int main(int argc, char** argv)
{
	if (argc == 2)
	{
		string val(argv[1]);
		if (val == "--help" || val == "--HELP" )
		{       printUsage();   return 0; }
		if (val == "--version" || val == "--VERSION")
		{
			cout << "Version: " << VERSION << ", Contact/Help: Rachid Ounit, rouni001 at cs.ucr.edu" << endl;	
			return 0; 
		}
	}
	if (argc < 6)
	{		
		cerr << "To run " << argv[0] << ", at least four  parameters are necessary:\n" ;
		cerr << "filename of the targets definition, directory of database, filename for objects, filename for results."<< endl;
		printUsage();
		return -1;
	}
	size_t	k 		= LENGTH, w = 0, mode = 1, cpu = 1, iterKmers = 0;
	ITYPE minT 		= 0, minO = 0, sfactor = 0;
	bool cLightDB 		= false, spacedK = false, ldm = false, tsk = false, kso= false, ext = false, isReduced = false, longR =  false;
	int i_targets	 	= -1, i_objects = -1, i_objects2 = -1, i_folder=-1, i_results =-1;
	std::vector<std::string> DSS;

	if (WEIGHT != LENGTH) 
	{
		spacedK = true;
		mode = 0;
		string s1("T295"), s2("T58570"), s3("T38570");
		DSS.push_back(s1);
		DSS.push_back(s2);
		DSS.push_back(s3);
	}

	for(size_t i = 1; i < argc ; i++)
	{
		string val(argv[i]);
		if (val ==  "-k")
		{	if (++i >= argc) {cerr << "Please specify the k-mer length!"<< endl; exit(1);	}
			k =  atoi(argv[i]); 
			if (k <= 1 || k > MAXK) { cerr <<"The k-mer length should be in [2,"<< MAXK << "]." << endl; exit(1);}
			continue;}
		if (val ==  "-t")
		{
			if (++i >= argc) {cerr << "Please specify the minimum frequency (targets)!"<< endl; exit(1);    }
			minT =  atoi(argv[i]);
			if (minT >= 65536) { cerr <<"The min k-mer frequency should be in [0,65535]." << endl; exit(1);}
			continue;}
		if (val ==  "-o")
		{
			if (++i >= argc) {cerr << "Please specify the minimum frequency (objects)!"<< endl; exit(1);    }
			minO =  atoi(argv[i]);
			if ( minO >= 65536) { cerr <<"The min k-mer frequency should be in [0,65535]."  << endl; exit(1);}
			continue;}
		if (val ==  "-m")
		{
			if (++i >= argc) {cerr << "Please specify the mode!"<< endl; exit(1);    }
			mode =  atoi(argv[i]);
			if (mode > 3) {	cerr <<"The mode of execution should be 0 (full), 1 (default), 2 (express) or 3 (colum-based)." << endl; exit(1);}
			if (mode != 3) {	kso = false;	}
			continue;
		}
		if (val ==   "-n")
		{
			if (++i >= argc) {cerr << "Please specify the number of threads!"<< endl; exit(1);    }
			cpu =  atoi(argv[i]);
			if (cpu< 1) { cerr <<"The number of threads should be higher than 0." << endl; exit(1);}
			continue;}
		if (val == "--long")
		{
			longR = true; continue;}
		if (val ==   "--ldm")
		{
			ldm = true; continue;}
		if (val ==  "--tsk")
		{
			cerr << "The option 'tsk' is no more supported.\n";
			continue;
		}
		if (val ==   "--kso")
		{	ext = true;
			kso = true; continue;}
		if (val ==  "--extended")
		{
			ext = true; continue;}
		if (val ==  "-T")
		{
			if (++i >= argc) {cerr << "Please specify the targets!"<< endl; exit(1);    }
			i_targets =  i;
			if (!validFile(argv[i])) { cerr <<"Failed to find/read the file of the targets definition: " << argv[i] << endl; exit(1);}
			continue;}
		if (val ==  "-O")
		{
			if (++i >= argc) {cerr << "Please specify the objects!"<< endl; exit(1);    }
			i_objects =  i;
			if (!validFile(argv[i])) { cerr <<"Failed to find/read the filename of objects: " << argv[i] << endl; exit(1);}
			continue;}
		if (val ==  "-P")
		{
			if (i+2 >= argc) {cerr << "Please specify the paired-end reads!"<< endl; exit(1);    }
			i++;
			i_objects = i;
			i_objects2 = i+1;
			if (!validFile(argv[i++])) { cerr <<"Failed to find/read " << argv[i-1] << endl; exit(1);}
			if (!validFile(argv[i]))   { cerr <<"Failed to find/read " << argv[i] << endl; exit(1);}
			continue;}
		if (val ==   "-D")
		{
			if (++i >= argc) {cerr << "Please specify the database directory!"<< endl; exit(1);    }
			i_folder =  i;
			if (!validFile(argv[i])) { cerr <<"Failed to find/read the directory:  " << argv[i] << endl; exit(1);}
			continue;}
		if (val ==   "-R")
		{
			if (++i >= argc) {cerr << "Please specify where to store results!"<< endl; exit(1);    }
			i_results =  i;
			continue;}
		if (val == "-g")
		{
			if (++i >= argc) {cerr << "Please specify a gap value!"<< endl; exit(1);    }
			iterKmers =  atoi(argv[i]);
			if (iterKmers < 4) {cerr << "Beware! With a gap value lower than 4, the RAM usage can be higher than 4 GB."<< endl; }
			if (iterKmers > 50) {cerr << "The gap value should be lower than 50."<< endl; exit(1);	}
			continue;}
		if (val == "-s")
		{
			if (++i >= argc) {cerr << "Please specify a sampling factor value!"<< endl; exit(1);    }
			sfactor =  atoi(argv[i]);
			if (sfactor < 1 || sfactor > SFACTORMAX) 
			{	cerr << "The sampling factor value should be in the interval [1,"<< SFACTORMAX <<"]."<< endl; 
				exit(1);    }
			continue;}
		cerr << "Failed to recognize option: " << val << endl;
		exit(1);
	}

	if (HTSIZE == LHTSIZE)
	{
		cLightDB = true; 
		if (iterKmers == 0)
		{ 	iterKmers = 4 ;}
		k = 27;
		w = 27; 
		sfactor = 1;
	}
	else
	{
		iterKmers = 0;
	}
	if (spacedK)
        {
                k = LENGTH;
                w = WEIGHT;
        }
	else
	{
		w = k;
	}
	if (sfactor == 0)
	{
		sfactor = (mode == 1)?2:1;
	}
	if (k == 0)
	{
		cerr << "Please specify a k-mer length: -k <integer>" << endl;
		exit(1);
	}
	if ( i_targets < 0 || i_folder < 0 || i_objects < 0 ||  i_results < 0)
	{
		cerr << "Failed to run " << argv[0] << ": at least four  parameters are necessary" ;
		cerr << ": file of targets, directory of database, file of objects, file for results."<< endl;

		printUsage();
		exit(1);
	}
	if (kso && mode != 3)
	{
		cerr << "Please, the option '--kso' is only for the spectrum mode."<< endl;
		exit(1);
	}
	const char * objects 	= i_objects > 0 ? argv[i_objects] : NULL;
	const char * objects2 	= i_objects2 > 0 ? argv[i_objects2] : NULL;
	const bool paired 	= i_objects2 > 0;

	string folder(argv[i_folder]);
	if (folder[folder.size()-1] != '/')
	{	folder.push_back('/');	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	size_t t_b = log(HTSIZE)/log(4.0);
	size_t max16 = t_b + 8;
	size_t max32 = t_b + 16;

	if (w <= max16)
	{
		// Use 2Bytes to store each discriminative k-mer
		CLARK<T16> classifier(k, argv[i_targets], folder.c_str(), w, DSS, minT, tsk, cLightDB, spacedK, iterKmers, cpu, sfactor, longR, ldm);
		if (paired)
		{	classifier.run(objects, objects2, argv[i_results], mode, minO, kso, ext, true);	}
		else
		{	classifier.run(objects, argv[i_results], mode, minO, kso, ext, true); 	}
		exit(0);
	}
	if (w <= max32)
	{
		// Use 4Bytes to store each discriminative k-mer
		CLARK<T32> classifier(k, argv[i_targets], folder.c_str(), w, DSS, minT, tsk, cLightDB, spacedK, iterKmers, cpu, sfactor, longR, ldm);
		if (paired)
                {       classifier.run(objects, objects2, argv[i_results], mode, minO, kso, ext, true); }
                else
                {       classifier.run(objects, argv[i_results], mode, minO, kso, ext, true);   }  
		exit(0);
	}
	if (w <= MAXK)
	{
		// Use 8Bytes to store each discriminative k-mer
		CLARK<T64> classifier(k, argv[i_targets], folder.c_str(), w, DSS, minT, tsk, cLightDB, spacedK, iterKmers, cpu, sfactor, longR, ldm);
		if (paired)
                {       classifier.run(objects, objects2, argv[i_results], mode, minO, kso, ext, true); }
                else
                {       classifier.run(objects, argv[i_results], mode, minO, kso, ext, true);   }  
		exit(0);
	}
	std::cout <<"This version of CLARK does not support k-mer length strictly higher than " << MAXK << std::endl;
	exit(-1);

}

