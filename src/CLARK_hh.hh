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

#ifndef CLARK_HH
#define CLARK_HH

#include <cstring>
#include <vector>
#include <iostream>
#include <cstdlib>
#include <sstream>
#include <time.h>
#include "./dataType.hh"
#include "./HashTableStorage_hh.hh"
#include "./spacedKmer.hh"
#include "./FileHandlerQ.hh"
#include "./FileHandlerA.hh"
#include "./FileHandler.hh"
#include "./HashTop.hh"

#define MAXRSIZE	10000
#define MAXNBREADS	10000
#define MAXREADLEN	500000
#define MAXSCAFLEN	1000000
#define MAXLONGLEN	50000000
#define MINCFSP		0.75
#define MINGMSP		0.06
#define	MXNMLEN		10000

static std::string clarkSizeToString(const size_t value)
{
	std::ostringstream out;
	out << value;
	return out.str();
}

static std::string clarkJoinPath(const std::string& folder, const std::string& name)
{
	if (folder.empty() || folder[folder.size()-1] == '/')
	{
		return folder + name;
	}
	return folder + "/" + name;
}

// CuD fine profiling: RDTSC (x86) reads the CPU's time-stamp counter
// directly (no syscall), giving a much lower-overhead per-call timer
// than clock_gettime/gettimeofday. Falls back to clock_gettime on
// non-x86 targets (Apple Silicon, ARM), where the accumulated value is
// nanoseconds rather than cycles.
#if defined(__x86_64__) || defined(__i386__)
#include <x86intrin.h>
#define CUD_HAVE_RDTSC 1
#else
#define CUD_HAVE_RDTSC 0
#endif

static inline uint64_t cudRdtsc()
{
#if CUD_HAVE_RDTSC
	return __rdtsc();
#else
	struct timespec ts;
	clock_gettime(CLOCK_MONOTONIC, &ts);
	return (uint64_t) ts.tv_sec * 1000000000ULL + (uint64_t) ts.tv_nsec;
#endif
}

// Calibrates cycles-per-second for cudRdtsc() by bracketing a short sleep
// with CLOCK_MONOTONIC (ground truth) and cudRdtsc() reads. On the
// non-x86 fallback, cudRdtsc() already returns nanoseconds, so this
// returns 1e9 (i.e., a no-op conversion factor).
static inline double cudCalibrateTscHz()
{
#if CUD_HAVE_RDTSC
	struct timespec ts0, ts1;
	clock_gettime(CLOCK_MONOTONIC, &ts0);
	uint64_t c0 = cudRdtsc();
	struct timespec sleepFor;
	sleepFor.tv_sec = 0;
	sleepFor.tv_nsec = 50000000L; // 50 ms
	nanosleep(&sleepFor, NULL);
	uint64_t c1 = cudRdtsc();
	clock_gettime(CLOCK_MONOTONIC, &ts1);
	double elapsedSec = (ts1.tv_sec - ts0.tv_sec) + (ts1.tv_nsec - ts0.tv_nsec) / 1e9;
	if (elapsedSec <= 0.0)
	{
		return 1e9;
	}
	return (double) (c1 - c0) / elapsedSec;
#else
	return 1e9;
#endif
}

template <typename HKMERr>
class CLARK
{
	private:
		std::vector< std::pair< std::string, std::string > >	m_targetsID;
		std::vector< std::string >				m_labels;
		std::vector< std::string >				m_labels_c;
		std::vector< std::string >				m_targetsName;

		// Using spectrum files
		std::vector< ObjectData >		m_objectsData;
		// Using spaced seeds
		std::vector< spacedKmer >		m_DSS;
		const size_t				m_weight;
		// Using fasta/fastq files
		ITYPE					m_nbObjects;
		std::vector< std::string>		m_objectsName;
		std::vector< ITYPE>                     m_objectsNorm;
		bool					m_isFastaFile;

		// choice 0: Full, 1:Default, 2: Fast, 3: Spectrum
		size_t					m_mode; 
		std::vector< std::vector<ILBL> >        m_targetsBest;
		std::vector< std::vector<size_t> > 	m_seqSNames;
		std::vector< std::vector<size_t> >      m_seqENames;
		std::vector< std::vector<size_t> >    	m_readsSPos; // Starting 
		std::vector< std::vector<size_t> >    	m_readsEPos; // Ending
		std::vector< std::vector<size_t> >    	m_readsLength;
		std::vector< size_t >			m_posReads;

		// Options for loading db				
		bool					m_isLightLoading;
		bool					m_isSpacedLoading;
	
		const char*				m_folder;
		bool					m_isSpecificTargetsPresent;

		// Dictionary k-mers <---> TargetsID
		EHashtable<HKMERr, bigElement> *	m_centralHt;

		// Tables for storing results
		std::vector< std::vector<ITYPE> >	m_ResultsCentral;
		std::vector< ILBL >			m_resultsFast;
		std::vector< std::string >		m_scoresLines;

		size_t					m_nbCPU;
		const size_t				m_kmerSize;
		const uint8_t				m_k;
		const ITYPE				m_minCountTarget;
		uint64_t				m_iterKmers;
		ITYPE					m_minCountObject;
		bool					m_useWeight;
		bool 					m_spectrumAnalysis;
		bool					m_isPaired;
		bool					m_isExtended;
		bool					m_isLongSequence;

		// Tables for storing temp results in default mode
		std::vector< std::vector <ITYPE> >	m_ITables;
		std::vector< std::vector <ILBL > >	m_Indexes;
		std::vector< std::vector <ITYPE> >	m_resultTargets;

		// Tables storing common and repetitive values for reading sequences
		int					m_Letter[256];
		int					m_table[256];
		int					m_rTable[256];
		uint64_t		    		m_pTable[4];
		std::vector< std::vector<uint64_t> >	m_powerTable;
		int	 				m_separators[256];

		// CuD profiling: phase timings, opt-in via CLARK_CUD_PROFILE=1
		double					m_buildTimeSec;
		double					m_loadTimeSec;
		struct timeval				m_matchEndTime;

		// CuD fine profiling: per-thread RDTSC cycle/call accumulators for
		// k-mer matching, hits update, and classification, opt-in via
		// CLARK_CUD_PROFILE_FINE=1 (adds per-k-mer timer overhead, unlike
		// the phase timings above).
		bool					m_fineProfile;
		double					m_tscHz;
		std::vector<uint64_t>			m_matchCycles;
		std::vector<uint64_t>			m_hitsCycles;
		std::vector<uint64_t>			m_classifyCycles;
		std::vector<uint64_t>			m_matchCalls;
		std::vector<uint64_t>			m_hitsCalls;
		std::vector<uint64_t>			m_classifyCalls;

	public:
		CLARK(const size_t& 			_kmerLength,
				const char* 		_filesName,
				const char* 		_folderName,
				const size_t&		_weight,
				const std::vector<std::string>& _DSS,
				const ITYPE& 		_minCountT,
				const bool&		_creatingkmfiles,
				const bool&		_isLightLoading,
				const bool&		_SpacedLoading,
				const uint64_t&		_iterKmers,
				const size_t& 		_nbCPU,
				const ITYPE&		_samplingFactor,
				const bool&		_isLongSequence,
				const bool& 		_mmapLoading 	= false
		     ); 

		~CLARK();

		void runSimple(const char* 		_fileTofilesname, 
				const char* 		_fileResult,
				const size_t& 		_mode, 
				const ITYPE& 		_minCountO, 
				const bool& 		_spectrumAnalysis = false,
				const bool&             _useWeight      = true
			      );

		void run(const char*			_filesToObjects,
				const char* 		_fileToResults,
				const size_t&   	_mode,
				const ITYPE& 		_minCountO,
				const bool& 		_spectrumAnalysis = false,
				const bool&		_isExtended     = false,	
				const bool&             _useWeight      = true
			);

		void run(const char*			_pairedfile1,
				const char*		_pairedfile2,
                                const char*             _fileToResults,
                                const size_t&           _mode,
                                const ITYPE&            _minCountO,
                                const bool&             _spectrumAnalysis = false,
                                const bool&             _isExtended     = false,
                                const bool&             _useWeight      = true
                        );

		void clear();

	private:
		void loadComputeObjectsSpectrumData();

		void createTargetFilesNames(std::vector< std::string >& 	_filesHT, 
				std::vector< std::string >& 			_filesHTC
				) const;

		void loadSpecificTargetSets(const std::vector<std::string>& 	_filesHT, 
				const std::vector<std::string>& 		_filesHTC, 
				const size_t& 					_sizeMotherHT,
				const ITYPE&            			_samplingFactor, 
				const bool&  					_mmapLoading
				);

		size_t makeSpecificTargetSets(const std::vector<std::string>& 	_filesHT, 
				const std::vector<std::string>& 		_filesHTC
				) const;

		bool getObjectsDataSpectrum(FILE * 				fileToScore
				);

		void getObjectsDataComputeFull(const bool&			isfasta,
				const char*					filename,
				const char* 					_fileResult
				);

		void getObjectsDataCompute(const uint8_t *			_map,
				const size_t& 					nb,
				FILE*             				_fout
				);

		void getObjectsDataComputeFast(const uint8_t *			_map,
				const size_t& 					nb,
				FILE*                   		        _fout
				);

		void getObjectsDataComputeFastLight(const uint8_t *  		_map,
				const size_t&                  			 nb,
				FILE*            			   	_fout
				);

		void getObjectsDataComputeFastSpaced(const bool&		isfasta,
				const char*                                     filename,
				const char*  
				);

		bool getTargetsData(const char* 				_filesName, 
				std::vector< std::string >&		 	_filesHT, 
				std::vector< std::string >&		 	_filesHTC,
				const bool&					_creatingkmfiles,
				const ITYPE&					_samplingfactor = 1
				);

		void print(const bool& 						_creatingkmfiles = false,
				const ITYPE&                                    _samplingfactor = 1
			  ) const;

		void printExtendedResultsHeader(const char* 			_fileResult,
				const bool& 					_sf = false
				) const;

		void printExtendedResults(const char*				 _fileResult
				) const;

		void printExtendedSResults(const char*                            _fileResult
				) const;

		void printExtendedSFResults(const char*                            _fileResult
				) const;

		void printSpeedStats(const struct timeval& 			_requestEnd,
				const struct timeval& 				_requestStart,
				const char* 					_fileResult
				) const;

		void printCuDProfile(const struct timeval& 			_requestStart,
				const struct timeval& 				_requestEnd
				) const;

		void printCuDProfileFine() const;

		// CuD profiling: split out of the hot loop as separate (non-inlined)
		// symbols so `perf report`/`perf annotate` can attribute samples to
		// "hits update" vs "classification" vs the k-mer matching call
		// (EHashtable::queryElement) without adding per-k-mer timer overhead.
		void updateHits(ITYPE* 					_resultTargets,
				ITYPE* 						_iTable,
				ILBL* 						_idx,
				size_t& 					_iSize,
				const ILBL& 					_h,
				const ITYPE& 					_token
				) const __attribute__((noinline));

		void classifyBest(const ITYPE& 				_score,
				const ILBL& 					_h,
				ILBL& 						_opt_h,
				ITYPE& 						_s_best
				) const __attribute__((noinline));

		void getdbName(char * 						_dbname,
				const int& 					_htID  = 0    	
			      ) const;

		std::string getdbNameString(const int& 			_htID = 0
				) const;
};

#endif
/*
 * @author: Rachid Ounit, Ph.D Candidate.
 * @date: 04/01/2014
 *
 */
#ifdef _OPENMP
#include <omp.h>
#endif

#include <string.h>
#include <fstream>
#include <time.h>
#include <unistd.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/mman.h>

#include "./file.hh"
#include "./kmersConversion.hh"
#include "./analyser.hh"

using namespace std;

template <typename HKMERr>
CLARK<HKMERr>::CLARK(const size_t& 	_kmerLength,
		const char* 		_filesName,
		const char* 		_folderName,
		const size_t&		_weight,
		const std::vector<std::string>& _DSS,
		const ITYPE& 		_minCountT,
		const bool&      	_creatingkmfiles,
		const bool&     	_isLightLoading,
		const bool&		_isSpacedLoading,
		const uint64_t&		_iterKmers,
		const size_t& 		_nbCPU,
		const ITYPE&            _samplingFactor,
		const bool&		_isLongSequence,
		const bool&     	_mmapLoading
		): 
	m_nbCPU(_nbCPU), 
	m_kmerSize(_kmerLength), m_k((uint8_t) _kmerLength), m_weight(_weight),
	m_nbObjects(0),
	m_folder(_folderName), 
	m_minCountTarget(_minCountT), 
	m_powerTable(m_kmerSize),
	m_isFastaFile(true),
	m_useWeight(true),
	m_iterKmers(_iterKmers),
	m_isLightLoading(_isLightLoading),
	m_isSpacedLoading(_isSpacedLoading),
	m_posReads(_nbCPU),
	m_isPaired(false),
	m_isExtended(false),
	m_isLongSequence(_isLongSequence)
{

#ifdef _OPENMP
	omp_set_num_threads(m_nbCPU);
#else
	m_nbCPU = 1;
#endif
	m_posReads.resize(m_nbCPU);
	m_readsLength.resize(m_nbCPU); 
	m_readsEPos.resize(m_nbCPU);
	m_readsSPos.resize(m_nbCPU);
	m_seqENames.resize(m_nbCPU);
	m_seqSNames.resize(m_nbCPU);
	m_targetsBest.resize(m_nbCPU);
	m_resultTargets.resize(m_nbCPU);
	m_ITables.resize(m_nbCPU);
	m_Indexes.resize(m_nbCPU);

	size_t base = 1;
	for(size_t p = 0; p < m_kmerSize ; p++)
	{       
		m_powerTable[p].resize(4);
		m_powerTable[p][0] = 0; m_powerTable[p][1] = base; m_powerTable[p][2] = base*2; m_powerTable[p][3] = base*3;  
		base *= 4;
	}
	m_pTable[0] = m_powerTable[m_kmerSize-1][0]; 
	m_pTable[1] = m_powerTable[m_kmerSize-1][1];
	m_pTable[2] = m_powerTable[m_kmerSize-1][2];
	m_pTable[3] = m_powerTable[m_kmerSize-1][3];

	for(size_t t = 0; t < 256 ; t++)
	{	m_table[t] = -1; m_separators[t] = 0; m_rTable[t] = -1;	m_Letter[t] = -1;}

	// Letter definition (Fastq files)
	for(size_t t = 65; t < 91; t++)
	{	m_Letter[t] = 0;}
	for(size_t t = 97; t < 123; t++)
	{       m_Letter[t] = 0;}
	// Nucleotide definitions (DNA)
	m_table['A'] = 0;	m_table['C'] = 1;       m_table['G'] = 2;	m_table['T'] = 3;
	m_table['a'] = 0;       m_table['c'] = 1;       m_table['g'] = 2;       m_table['t'] = 3;
	// Nucleotide definitions (RNA)
	m_table['U'] = 3; 	m_table['u'] = 3;

	// Other symbol definitions
	m_table['>'] = -2;	m_table['\n']= -10;
	// Reverse - Nucleotide definitions (DNA)
	m_rTable['A'] = 3;      m_rTable['C'] = 2;       m_rTable['G'] = 1;       m_rTable['T'] = 0;
	m_rTable['a'] = 3;      m_rTable['c'] = 2;       m_rTable['g'] = 1;       m_rTable['t'] = 0;
	// Reverse - Nucleotide definitions (RNA)
	m_rTable['U'] = 0;      m_rTable['u'] = 0;

	// Reverse - Other symbol definitions
	m_rTable['>'] = -2;     m_rTable['\n']= -10;

	m_separators[' '] = 1;	m_separators['\t'] = 1; 	m_separators['\n'] = 1;

	for(size_t y = 0; y < _DSS.size(); y++)
	{
		spacedKmer ss(_DSS[y]);	
		m_DSS.push_back(ss);	
	}
	vector<string> filesHT, filesHTC;
	size_t sizeMotherHT = 0;
	m_buildTimeSec = 0.0;
	m_loadTimeSec = 0.0;
	if (!getTargetsData(_filesName, filesHT, filesHTC, _creatingkmfiles, _samplingFactor))
	{
		cerr << "Starting the creation of the database of targets specific " << m_kmerSize << "-mers from input files..." << endl;
		struct timeval buildStart, buildEnd;
		gettimeofday(&buildStart, NULL);
		sizeMotherHT = makeSpecificTargetSets(filesHT, filesHTC);
		gettimeofday(&buildEnd, NULL);
		m_buildTimeSec = (buildEnd.tv_sec - buildStart.tv_sec) + (buildEnd.tv_usec - buildStart.tv_usec) / 1000000.0;
	}
	{
		struct timeval loadStart, loadEnd;
		gettimeofday(&loadStart, NULL);
		loadSpecificTargetSets(filesHT, filesHTC, sizeMotherHT, _samplingFactor, _mmapLoading);
		gettimeofday(&loadEnd, NULL);
		m_loadTimeSec = (loadEnd.tv_sec - loadStart.tv_sec) + (loadEnd.tv_usec - loadStart.tv_usec) / 1000000.0;
	}

	for(size_t t = 0; t < m_nbCPU; t++)
	{
		m_resultTargets[t].resize(m_labels.size() + m_labels_c.size());
		m_ITables[t].resize(m_labels.size() + m_labels_c.size());
		m_Indexes[t].resize(m_labels.size() + m_labels_c.size());
	}

	m_fineProfile = getenv("CLARK_CUD_PROFILE_FINE") != NULL;
	m_tscHz = m_fineProfile ? cudCalibrateTscHz() : 0.0;
	m_matchCycles.assign(m_nbCPU, 0);
	m_hitsCycles.assign(m_nbCPU, 0);
	m_classifyCycles.assign(m_nbCPU, 0);
	m_matchCalls.assign(m_nbCPU, 0);
	m_hitsCalls.assign(m_nbCPU, 0);
	m_classifyCalls.assign(m_nbCPU, 0);
}

	template <typename HKMERr>
CLARK<HKMERr>::~CLARK()
{
	clear();
	delete m_centralHt;
}

	template <typename HKMERr>
void CLARK<HKMERr>::clear()
{
	m_nbObjects = 0;
	m_objectsName.clear();
	m_objectsData.clear();
	m_objectsNorm.clear();
	for(size_t  i = 0; i < m_nbCPU; i++)
	{
		m_targetsBest[i].clear();
		m_seqSNames[i].clear();
		m_seqENames[i].clear();
		m_readsEPos[i].clear();
		m_readsSPos[i].clear();
		m_readsLength[i].clear();

		ITYPE* resultTargets = &m_resultTargets[i].front();
		ITYPE* iTable = &m_ITables[i].front();
		ILBL* idx = &m_Indexes[i].front();
		for(size_t t= 0; t < m_labels.size()+ m_labels_c.size() ;t++)
		{	resultTargets[t] = 0; iTable[t] = 0; idx[t] = 0;	}
	}
	m_ResultsCentral.clear();
	m_resultsFast.clear();
	m_scoresLines.clear();
}

	template <typename HKMERr>
	void CLARK<HKMERr>::createTargetFilesNames(vector<string>& _filesHT, vector<string>& _filesHTC) const
	{
		for(size_t t = 0 ; t < m_labels_c.size(); t++)
		{
			std::string fname;
			if (m_isLightLoading)
			{
				fname = clarkJoinPath(m_folder, m_labels_c[t] + "_k" + clarkSizeToString((size_t) m_kmerSize) + "_light.ht");
			}
			else
			{
				fname = clarkJoinPath(m_folder, m_labels_c[t] + "_k" + clarkSizeToString((size_t) m_kmerSize) + ".ht");
			}
			_filesHTC.push_back(fname);
		}

		for(size_t t = 0 ; t < m_labels.size(); t++)
		{
			std::string fname;
			if (m_isLightLoading)
			{
				fname = clarkJoinPath(m_folder, m_labels[t] + "_k" + clarkSizeToString((size_t) m_kmerSize) + "_light.ht");
			}
			else
			{
				fname = clarkJoinPath(m_folder, m_labels[t] + "_k" + clarkSizeToString((size_t) m_kmerSize) + ".ht");
			}
			_filesHT.push_back(fname);
		}
	}

	template <typename HKMERr>
void CLARK<HKMERr>::run(const char* _filesToObjects, const char* _fileToResults, const size_t& _mode, const ITYPE& _minCountO, const bool& _spectrumAnalysis, const bool& _isExtended, const bool&  _useWeight)
{
	FILE* fd = fopen(_fileToResults, "r");
	m_isPaired = false;
	m_isExtended = _isExtended;
	string mode(_mode == 0 ? "Full": (_mode == 1? "Default" : (_mode == 2? "Express":"Spectrum")));
	if (fd == NULL )
	{
		cerr << "Mode: " << mode<< ",\tProcessing file: " << _filesToObjects << ",\t using "<< m_nbCPU << " CPU." <<  endl;
		CLARK::runSimple(_filesToObjects, _fileToResults, _mode, _minCountO, _spectrumAnalysis, _useWeight);
		return;
	}
	fclose(fd);

	fd = fopen(_filesToObjects, "r");
	string line = "";
	getLineFromFile(fd, line);
	vector<string> ele;
	vector<char> seps;
	seps.push_back(' ');
	seps.push_back('\t');
	seps.push_back(',');
	fclose(fd);

	getElementsFromLine(line, seps, ele);	

	if (line[0] == '>' || line[0] == '@' || ele.size() == 2)
	{
		cerr << "Mode: " << mode<< ",\tProcessing file: " << _filesToObjects << ",\t using "<< m_nbCPU << " CPU." <<  endl;
		CLARK::runSimple(_filesToObjects, _fileToResults, _mode, _minCountO, _spectrumAnalysis, _useWeight);
		return;
	}
	FILE * r_fd = fopen(_fileToResults, "r");
	FILE * o_fd = fopen(_filesToObjects, "r");
	string o_line = "", r_line = "";
	cerr << "Mode: " << mode << " using " << m_nbCPU << " CPU." << endl;
	while (getLineFromFile(o_fd, o_line) && getLineFromFile(r_fd, r_line))
	{
		cerr << "> Processing file: " << o_line.c_str() <<  endl;
		CLARK::runSimple(o_line.c_str(), r_line.c_str(), _mode, _minCountO, _spectrumAnalysis, _useWeight); 
	}
	fclose(r_fd); 
	fclose(o_fd);
	return;
}

        template <typename HKMERr>
void CLARK<HKMERr>::run(const char* _pairedfile1, const char* _pairedfile2, const char* _fileToResults, const size_t& _mode, const ITYPE& _minCountO, const bool& _spectrumAnalysis, const bool& _isExtended, const bool&  _useWeight)
{
        FILE* fd 	= fopen(_fileToResults, "r");
        m_isPaired 	= true;
        m_isExtended 	= _isExtended;
        string mode(_mode == 0 ? "Full": (_mode == 1? "Default" : (_mode == 2? "Express":"Spectrum")));
	string mergedFiles;
        if (fd == NULL )
        {
		// Merge _pairedfile1 + _pairedfile2
		mergedFiles = string(_pairedfile1) + "_ConcatenatedByCLARK.fa";
                mergePairedFiles(_pairedfile1, _pairedfile2, mergedFiles.c_str());
		cerr << "Mode: " << mode<< ",\tProcessing file: " << mergedFiles << ",\t using "<< m_nbCPU << " CPU." <<  endl;
                CLARK::runSimple(mergedFiles.c_str(), _fileToResults, _mode, _minCountO, _spectrumAnalysis, _useWeight);
		// Delete file
		deleteFile(mergedFiles.c_str());
                return;
        }
        fclose(fd);

        fd 		= fopen(_pairedfile1, "r");
        string line 	= "";
        getLineFromFile(fd, line);
        vector<string> ele;
        vector<char> seps;
        seps.push_back(' ');
        seps.push_back('\t');
        seps.push_back(',');
        fclose(fd);

        getElementsFromLine(line, seps, ele);
        if (line[0] == '>' || line[0] == '@' || ele.size() == 2)
        {
		// Merge _pairedfile1 + _pairedfile2 
                mergedFiles = string(_pairedfile1) + "_ConcatenatedByCLARK.fa";
                mergePairedFiles(_pairedfile1, _pairedfile2, mergedFiles.c_str());
		cerr << "Mode: " << mode<< ",\tProcessing file: " << mergedFiles << ",\t using "<< m_nbCPU << " CPU." <<  endl;
                CLARK::runSimple(mergedFiles.c_str(), _fileToResults, _mode, _minCountO, _spectrumAnalysis, _useWeight);
		// Delete file
		deleteFile(mergedFiles.c_str());
                return;
        }
        FILE * r_fd 	= fopen(_fileToResults, "r");
        FILE * o1_fd 	= fopen(_pairedfile1, "r");
	FILE * o2_fd 	= fopen(_pairedfile2, "r");
        string o1_line 	= "", o2_line   = "", r_line = "";
        cerr << "Mode: " << mode << " using " << m_nbCPU << " CPU." << endl;
        while (getLineFromFile(o1_fd, o1_line) && getLineFromFile(o2_fd, o2_line) && getLineFromFile(r_fd, r_line))
        {
		// Merge _pairedfile1 + _pairedfile2 
		mergedFiles = o1_line + "_ConcatenatedByCLARK.fa";
                mergePairedFiles(o1_line.c_str(), o2_line.c_str(), mergedFiles.c_str());

                cerr << "> Processing file: " << mergedFiles <<  endl;
                CLARK::runSimple(mergedFiles.c_str(), r_line.c_str(), _mode, _minCountO, _spectrumAnalysis, _useWeight);
		// Delete file
		deleteFile(mergedFiles.c_str());
        }
        fclose(r_fd);
        fclose(o1_fd);
	fclose(o2_fd);
        return;
}

	template <typename HKMERr>
void CLARK<HKMERr>::runSimple(const char* _fileTofilesname, const char* _fileResult, const size_t& _mode, const ITYPE& _minCountO, const bool& _spectrumAnalysis, const bool& _useWeight)
{
	clear();
	m_isFastaFile		= true;
	m_minCountObject 	= _minCountO;
	m_mode	 		= _mode;
	m_spectrumAnalysis	= _spectrumAnalysis;
	m_useWeight		= _useWeight;

	size_t fileSize = 0;
	int fd;
	uint8_t* map = NULL;
	struct timeval requestStart, requestEnd;

	std::ifstream in(_fileTofilesname, std::ios::binary | std::ios::ate);
	fileSize = in.tellg();
	in.seekg(0,in.beg);
	string fline;
	getline(in,fline);
	in.close();

	if (m_mode == 1 || m_mode == 2)
	{
		fd = open(_fileTofilesname, O_RDONLY);
		if (fd == -1 || fileSize == 0)
		{
			cerr << "Failed to open " << _fileTofilesname << ".\n Please verify you typed the correct filename."<< endl;
			return;
		}
		map = (uint8_t*) mmap(0, fileSize, PROT_READ, MAP_SHARED, fd, 0);
		if ( map == MAP_FAILED )
		{
			close(fd);
			cerr << "Failed to mmapping the file. Please consider using the full mode instead (option: -m 0)." << endl;
			return;
		}
	}
	// Checking file to store result:
	string sfileResult(_fileResult);
	sfileResult += ".csv";
	const char* fileResult = sfileResult.c_str();
	// Try to access and erase content of the file If non-empty
	FILE * _fout = fopen(fileResult,"w");
	if (_fout == NULL)
	{
		cerr << "Failed to create/open file result: "<<fileResult<<".\n Please verify you typed the correct filename." << endl;
		return;
	}

	if (m_mode == 2 && !m_isLightLoading && !m_isSpacedLoading)
	{
		gettimeofday(&requestStart, NULL);
		///////////////////////////////////////////////////////////////////////
		getObjectsDataComputeFast(map, fileSize, _fout);
		///////////////////////////////////////////////////////////////////////
		gettimeofday(&requestEnd, NULL);
		fclose(_fout);
		// Measurement execution time
		printSpeedStats(requestEnd,requestStart,fileResult);

		msync(map, fileSize, MS_SYNC);
		if (munmap(map, fileSize) == -1)
		{	cerr << "Error un-mmapping the file." << endl;}

		close(fd);
		return;
	}
	if (m_mode == 2 && m_isLightLoading)
	{
		gettimeofday(&requestStart, NULL);
		///////////////////////////////////////////////////////////////////////
		getObjectsDataComputeFastLight(map, fileSize, _fout);
		///////////////////////////////////////////////////////////////////////
		gettimeofday(&requestEnd, NULL);
		fclose(_fout);
		// Measurement execution time
		printSpeedStats(requestEnd,requestStart,fileResult);
		printCuDProfile(requestStart,requestEnd);

		msync(map, fileSize, MS_SYNC);
		if (munmap(map, fileSize) == -1)
		{       cerr << "Error un-mmapping the file." << endl;}
		close(fd);

		return;
	}
	if (m_mode == 1 && !m_isSpacedLoading)
	{
		gettimeofday(&requestStart, NULL);
		///////////////////////////////////////////////////////////////////////
		getObjectsDataCompute(map, fileSize, _fout);
		///////////////////////////////////////////////////////////////////////
		gettimeofday(&requestEnd, NULL);
		fclose(_fout);
		// Measurement execution time
		printSpeedStats(requestEnd,requestStart,fileResult);
		printCuDProfile(requestStart,requestEnd);
		printCuDProfileFine();

		msync(map, fileSize, MS_SYNC);
		if (munmap(map, fileSize) == -1)
		{       cerr << "Error un-mmapping the file." << endl;}
		close(fd);
		return;
	}

	fclose(_fout);

	if (m_mode == 2 && m_isSpacedLoading)
	{
		bool isfasta = fline[0] == '>';
		if (!isfasta && fline[0] != '@')
		{
			cerr << "Failed to recognize the format of the file: "<< _fileTofilesname << endl;
			exit(1);
		}
		////////////////////////////////////////////////////////////////////////
		getObjectsDataComputeFastSpaced(isfasta,_fileTofilesname, fileResult);
		////////////////////////////////////////////////////////////////////////
		return;
	}
	if (m_mode == 0 || (m_mode == 1 && m_isSpacedLoading))
	{
		bool isfasta = fline[0] == '>';
		if (!isfasta && fline[0] != '@')
		{
			cerr << "Failed to recognize the format of the file: "<< _fileTofilesname << endl;
			exit(1);
		}
		////////////////////////////////////////////////////////////////////////
		getObjectsDataComputeFull(isfasta,_fileTofilesname, fileResult);
		////////////////////////////////////////////////////////////////////////
		return;
	}
	FILE * sfd = fopen(_fileTofilesname, "r");	
	gettimeofday(&requestStart, NULL);
	///////////////////////////////////////////////////////////////////////
	getObjectsDataSpectrum(sfd);
	loadComputeObjectsSpectrumData();
	printExtendedResultsHeader(fileResult);
	printExtendedResults(fileResult);
	///////////////////////////////////////////////////////////////////////
	gettimeofday(&requestEnd, NULL);
	// Measurement execution time
	printSpeedStats(requestEnd,requestStart,fileResult);

	fclose(sfd);
	return;
}

	template <typename HKMERr>
void CLARK<HKMERr>::loadComputeObjectsSpectrumData()
{
	m_nbObjects = m_objectsName.size();
	m_objectsNorm.resize(m_objectsName.size());
	m_objectsData.resize(m_objectsName.size());
	m_ResultsCentral.resize(5);
	m_scoresLines.resize(m_objectsName.size(),"");

	for(size_t i = 0 ; i < m_ResultsCentral.size(); i++)
	{
		m_ResultsCentral[i].resize(m_objectsName.size(),0); 
	}

	ITYPE size = 0, bacNorm = 0;
	ILBL h;

	HashTop hStore;

	for(size_t t = 0; t < m_objectsName.size(); t++)
	{
		size = 0;
		FILE * sfd = fopen(m_objectsName[t].c_str(), "r");
		if (sfd != NULL)
		{
			if (!m_useWeight)
			{
				string tmp = "";
				ITYPE frq;
				while (getFirstAndSecondElementInLine(sfd, tmp, frq))
				{	if (frq > m_minCountObject)      {size++;}}
				fclose(sfd);
				sfd = fopen(m_objectsName[t].c_str(), "r");
			}
			if (m_spectrumAnalysis)
			{
				fclose(sfd);
				bool _bumpFound = false;
				int minFreq = ((int)m_minCountObject) + 1 , maxFreq = 200000000 ;
				analyser freqExtractor(m_objectsName[t].c_str());
				_bumpFound = freqExtractor.getBumpInterval(minFreq, maxFreq, 1);
				m_objectsData[t].BumpFound = _bumpFound;
				m_objectsData[t].MinCount = minFreq;
				m_objectsData[t].MaxCount = maxFreq;
				sfd = fopen(m_objectsName[t].c_str(), "r");
			}
			string line = "";
			ITYPE val;
			IKMER bKmer;
			bacNorm = 0;
			string lineReduced = "";
			while (getFirstAndSecondElementInLine(sfd, line, val))
			{
				if (val > m_minCountObject)
				{
					bKmer = line;
					if (m_centralHt->queryElement(bKmer, h))
					{	
						//m_ResultsCentral[h][t] += (m_useWeight) ? val: 1;	
						hStore.insert(h,m_useWeight?val:1);
					}
					bacNorm += val;
				}
			}
			fclose(sfd);
			m_objectsNorm[t]  = (m_useWeight) ? bacNorm : size;
			hStore.getBest(m_ResultsCentral[0][t],m_ResultsCentral[1][t]);
			hStore.getSecondBest(m_ResultsCentral[2][t],m_ResultsCentral[3][t]);
			hStore.getTotal(m_ResultsCentral[4][t]);
			hStore.getScoresLine(m_targetsName.size()-1,m_scoresLines[t]);
			hStore.next();
		}
		else
		{
			cerr << "Failed to open " << m_objectsName[t] << endl;
		}
	}
	return;
}

	template <typename HKMERr>
	void CLARK<HKMERr>::getdbName(char *   _dbname, const int& _htID) const
	{
		std::string dbname = getdbNameString(_htID);
		if (dbname.size() >= MXNMLEN)
		{
			std::cerr << "Database path is too long: " << dbname << std::endl;
			exit(1);
		}
		memcpy(_dbname, dbname.c_str(), dbname.size()+1);
	}

	template <typename HKMERr>
	std::string CLARK<HKMERr>::getdbNameString(const int& _htID) const
	{
		size_t sizeHTS =  m_labels.size() + m_labels_c.size();
		std::string basename = "db_central_k" + clarkSizeToString((size_t)m_kmerSize)
			+ "_t" + clarkSizeToString(sizeHTS)
			+ "_s" + clarkSizeToString((size_t)HTSIZE)
			+ "_m" + clarkSizeToString((size_t)m_minCountTarget);

		if (m_isLightLoading)
		{
			return clarkJoinPath(m_folder, basename + "_light_" + clarkSizeToString((size_t)m_iterKmers) + ".tsk");
		}
		if (m_weight == m_kmerSize)
		{
			return clarkJoinPath(m_folder, basename + ".tsk");
		}
		if (_htID <= 0 || (size_t)(_htID - 1) >= m_DSS.size())
		{
			std::cerr << "Invalid spaced-kmer database id: " << _htID << std::endl;
			exit(1);
		}
		return clarkJoinPath(clarkJoinPath(m_folder, m_DSS[_htID-1].getFolder()), basename + "_w" + clarkSizeToString(m_weight) + ".tsk");
	}

	template <typename HKMERr>
void CLARK<HKMERr>::loadSpecificTargetSets(const vector<string>& _filesHT, 
		const vector<string>& 	_filesHTC, 
		const size_t& 		_sizeMotherHT,
		const ITYPE& 		_samplingFactor,	
		const bool&  		_mmapLoading
		)
{
	size_t kmersLoaded = 0;
	ITYPE minCount = m_minCountTarget;
	m_centralHt = new EHashtable<HKMERr, bigElement>(m_weight, m_labels, m_labels_c, m_DSS);
	size_t fileSize = 0;

	if (m_isSpacedLoading /*m_kmerSize != m_weight*/)
	{
		vector<string> names;
		for(size_t u = 0; u < m_DSS.size(); u++)
		{
			names.push_back(getdbNameString(u+1));
		}
		cerr << "Loading databases of spaced k-mers... " << endl;
		m_centralHt->AddDB(names, fileSize, _samplingFactor);
		cerr << "Loading done (" << fileSize / 1000000<<" MB read, with sampling factor " << _samplingFactor << ", " << m_centralHt->Size() << " spaced k-mers stored)" << endl;	
		return;
	}
	std::string cfname = getdbNameString();
	cerr << "Loading database [" << cfname << ".*] ..." << endl;
	if (m_centralHt->Read(cfname.c_str(), fileSize, m_nbCPU, _samplingFactor, _mmapLoading))
	{
		cerr << "Loading done (database size: " << fileSize / 1000000<<" MB read, with sampling factor " << _samplingFactor << ")" << endl;
		return;
	}
	if ( _filesHT.size() + _filesHTC.size() == 0)
	{
		cerr << "Failed to find the database." << endl;
		exit(-1);	
	}
}

template <typename HKMERr>
size_t CLARK<HKMERr>::makeSpecificTargetSets(const vector<string>& _filesHT, const vector<string>& _filesHTC) const
{
	size_t 	nt = 0;
	if (m_isLightLoading)
	{
		EHashtable<HKMERr, lElement> commonKmersHT(m_kmerSize, m_labels, m_labels_c);
		for(size_t t = 0 ; t < m_targetsID.size(); t++)
		{
			FILE* fd = fopen(m_targetsID[t].first.c_str(),"r");
			if (fd == NULL)
			{	cerr << "Failed to open " << m_targetsID[t].first << endl;
				continue;
			}
			ILBL tgt_id;
			commonKmersHT.getTargetID(m_targetsID[t].second, tgt_id);
			char c[MAXRSIZE];
			size_t len = fread(&c, 1, MAXRSIZE, fd), i = 0;
			uint64_t _km_f = 0, _km_r = 0;
			uint8_t cpt = 0;
			uint64_t iter = 0;
			if (len > 0 && c[0] == '>')
			{
				while (len != 0)
				{
					while( i < len)
					{
						if (m_table[c[i]] >= 0)
						{
							nt++;
							_km_r <<= 2;    
							_km_r += 3 - m_table[c[i]];
							if (cpt  == m_kmerSize - 1)
							{
								//[/////////////////////////////////////////////////////////////////////
								if ( iter % m_iterKmers == 0  )
								{       commonKmersHT.addElement(_km_r, m_targetsID[t].second, tgt_id, 1);}
								_km_r = 0; cpt = 0; 
								i++;
								iter++;
								continue;
								//]//////////////////////////////////////////////////////////////////////
							}
							cpt++;  
							i++;
							continue;
						}
						if (m_table[c[i]] == -10)
						{       i++;
							continue;   }
						if (m_table[c[i]] == -1)
						{        nt++; _km_r = 0; cpt = 0; i++;
							continue;  }
						if (m_table[c[i]] == -2)
						{
							_km_r = 0; cpt = 0; 
							while (len != 0)
							{
								while(i < len && c[i] != '\n')
								{       i++;}
								if (i < len)
								{   break;}
								len = fread(&c, 1, MAXRSIZE, fd );
								i = 0;
							}
							i++;
							continue;
						}
						cerr << m_targetsID[t].first << ": " ;
						cerr << "Failed to process sequence -- wrong character found [ASCII]: " << (size_t) c << endl;
					}
					len = fread(&c, 1, MAXRSIZE, fd );
					i = 0;
				}
				fclose(fd);
				cerr << "\r Progress report: (" <<t+1<< "/"<<m_targetsID.size()<<")    ";
				continue;
			}
			if (len > 0 && c[0] == '@')
			{
				while (len != 0)
				{
					while(i < len && c[i] != '\n')
					{       i++;}
					if (i < len)
					{   break;}
					len = fread(&c, 1, MAXRSIZE, fd );
					i = 0;
				}
				i++;
				while (len != 0)
				{
					while( i < len)
					{
						if (m_table[c[i]] >= 0)
						{
							nt++;
							_km_r <<= 2; 
							_km_r += 3 - m_table[c[i]];
							if (cpt  == m_kmerSize - 1)
							{
								//[/////////////////////////////////////////////////////////////////////
								if ( iter % m_iterKmers == 0 )
								{       commonKmersHT.addElement(_km_r, m_targetsID[t].second, tgt_id, 1);}
								_km_r = 0; cpt = 0; 
								i++;
								iter++;
								continue;
								//]//////////////////////////////////////////////////////////////////////
							}
							cpt++;
							i++;
							continue;
						}
						if (m_table[c[i]] == -10)
						{
							_km_r = 0; cpt = 0;
							// Pass third line
							i++;
							while (len != 0)
							{
								while(i < len && c[i] != '\n')
								{       i++;}
								if (i < len)
								{   break;}
								len = fread(&c, 1, MAXRSIZE, fd );
								i = 0;
							}
							i++;
							// Pass fouth line
							while (len != 0)
							{
								while(i < len && c[i] != '\n')
								{       i++;}
								if (i < len)
								{   break;}
								len = fread(&c, 1, MAXRSIZE, fd );
								i = 0;
							}
							i++;
							// Pass first line
							while (len != 0)
							{
								while(i < len && c[i] != '\n')
								{       i++;}
								if (i < len)
								{   break;}
								len = fread(&c, 1, MAXRSIZE, fd );
								i = 0;
							}
							i++;
							continue;
						}
						if (m_table[c[i]] == -1)
						{
							nt++;
							_km_r = 0; cpt = 0; 
							i++;
							continue;
						}
						cerr << m_targetsID[t].first << ": " ;
						cerr << "Failed to process sequence -- wrong character found [ASCII]: " << (size_t) c << endl;
					}
					len = fread(&c, 1, MAXRSIZE, fd );
					i = 0;
				}
				fclose(fd);
				cerr << "\r Progress report: (" <<t+1<< "/"<<m_targetsID.size()<<")              ";
				continue;
			}
			fseek(fd, 0, 0);
			string s_kmer = "";
			ITYPE val = 0;
			uint8_t counter = 0;
			while (getFirstAndSecondElementInLine(fd, s_kmer, val))
			{
				if (s_kmer.size() >= m_kmerSize && counter % m_iterKmers == 0 && val > m_minCountTarget)
				{
					commonKmersHT.addElement(s_kmer, m_targetsID[t].second, (size_t) val);
					counter = 0;
				}
				counter++;
				continue;
			}
			fclose(fd);
			cerr << "\r Progress report: (" <<t+1<< "/"<<m_targetsID.size()<<")              ";
		}
		cerr << nt << " nt read in total." << endl;
		cerr << "Mother Hashtable successfully built. "<<commonKmersHT.Size()<<" " << m_kmerSize << "-mers stored." <<  endl;
		size_t sizeMotherTable = commonKmersHT.Size();

		commonKmersHT.SortAllHashTable(2);
		commonKmersHT.RemoveCommon(m_labels_c, m_minCountTarget);
			std::string cfname = getdbNameString();

			cerr << "Creating light database in disk..." << endl;
			uint64_t nbElement = commonKmersHT.Write(cfname.c_str(),2);
			cerr << nbElement << " " << m_kmerSize << "-mers successfully stored in database." << endl;

		return sizeMotherTable;
	}
	if (_filesHTC.size() + _filesHT.size() == 0)
	{
		EHashtable<HKMERr, lElement> commonKmersHT(m_kmerSize, m_labels, m_labels_c);
		for(size_t t = 0 ; t < m_targetsID.size(); t++)
		{
			FILE* fd = fopen(m_targetsID[t].first.c_str(),"r");
			if (fd == NULL)
			{
				cerr << "Failed to open " << m_targetsID[t].first << endl;
			}
			else
			{
				ILBL tgt_id;
				commonKmersHT.getTargetID(m_targetsID[t].second, tgt_id);
				char c[MAXRSIZE];
				size_t len = fread(&c, 1, MAXRSIZE, fd), i = 0;
				uint64_t _km_f = 0, _km_r = 0;
				bool _isfull = false;
				uint8_t cpt = 0;		
				if (len > 0 && c[0] == '>')
				{
					while (len != 0)
					{
						while( i < len)
						{
							if (m_table[c[i]] >= 0)
							{
								nt++;
								if (_isfull)
								{
									_km_f >>= 2;
									_km_f += m_powerTable[cpt][m_table[c[i]]];

									/////////////////////////////////////////////////////////////////////
									commonKmersHT.addElement(_km_f, m_targetsID[t].second, tgt_id, 1);
									///////////////////////////////////////////////////////////////////////
									i++;
									continue;
								}
								_km_r <<= 2;    
								_km_r += 3 - m_table[c[i]];
								if (cpt  == m_kmerSize -1)
								{       
									_isfull = true;
									//[/////////////////////////////////////////////////////////////////////
									commonKmersHT.addElement(_km_r, m_targetsID[t].second, tgt_id, 1);
									//]//////////////////////////////////////////////////////////////////////
									_km_f = _km_r;
									// The following 6 lines come from Jellyfish source code
									_km_f = ((_km_f >> 2)  & 0x3333333333333333UL) | ((_km_f & 0x3333333333333333UL) << 2);
									_km_f = ((_km_f >> 4)  & 0x0F0F0F0F0F0F0F0FUL) | ((_km_f & 0x0F0F0F0F0F0F0F0FUL) << 4);
									_km_f = ((_km_f >> 8)  & 0x00FF00FF00FF00FFUL) | ((_km_f & 0x00FF00FF00FF00FFUL) << 8);
									_km_f = ((_km_f >> 16) & 0x0000FFFF0000FFFFUL) | ((_km_f & 0x0000FFFF0000FFFFUL) << 16);
									_km_f = ( _km_f >> 32                        ) | ( _km_f                        << 32);
									_km_f = (((uint64_t)-1) - _km_f) >> (64 - (m_k << 1));
									i++;
									continue;
								}
								cpt++;  i++;    
								continue;
							}
							if (m_table[c[i]] == -10)
							{	i++;
								continue;   }
							if (m_table[c[i]] == -1)
							{	nt++; _km_r = 0; cpt = 0; _isfull = false; i++; 
								continue;  }
							if (m_table[c[i]] == -2)
							{
								_km_r = 0; cpt = 0; _isfull = false; 
								while (len != 0)
								{
									while(i < len && c[i] != '\n')  {	i++;}
									if (i < len)		{   break;}
									len = fread(&c, 1, MAXRSIZE, fd );
									i = 0;
								}
								i++;   
								continue;
							}
							cerr << m_targetsID[t].first << ": " ;
							cerr << "Failed to process sequence -- wrong character found [ASCII]: " << (size_t) c << endl;
						}
						len = fread(&c, 1, MAXRSIZE, fd );
						i = 0;
					}
					fclose(fd);
					cerr << "\r Progress report: (" <<t+1<< "/"<<m_targetsID.size()<<")    ";
					continue;	
				}
				if (len > 0 && c[0] == '@')
				{
					i = 0;
					while (len != 0)
					{
						while(i < len && c[i] != '\n')            {	i++;}
						if (i < len)      {   break;}
						len = fread(&c, 1, MAXRSIZE, fd );
						i = 0;
					}
					i++;
					while (len != 0)
					{
						while( i < len)
						{
							if (m_table[c[i]] >= 0)
							{
								nt++;
								if (_isfull)
								{
									_km_f >>= 2; 
									_km_f += m_powerTable[cpt][m_table[c[i]]];
									/////////////////////////////////////////////////////////////////////
									commonKmersHT.addElement(_km_f, m_targetsID[t].second,  tgt_id, 1);
									///////////////////////////////////////////////////////////////////////
									i++;    
									continue;
								}
								_km_r <<= 2;    
								_km_r += 3 - m_table[c[i]];
								if ( cpt  == m_kmerSize - 1)
								{       _isfull = true;
									//[/////////////////////////////////////////////////////////////////////
									commonKmersHT.addElement(_km_r, m_targetsID[t].second, tgt_id, 1);
									//]//////////////////////////////////////////////////////////////////////
									_km_f = _km_r;
									// The following 6 lines come from Jellyfish source code
									_km_f = ((_km_f >> 2)  & 0x3333333333333333UL) | ((_km_f & 0x3333333333333333UL) << 2);
									_km_f = ((_km_f >> 4)  & 0x0F0F0F0F0F0F0F0FUL) | ((_km_f & 0x0F0F0F0F0F0F0F0FUL) << 4);
									_km_f = ((_km_f >> 8)  & 0x00FF00FF00FF00FFUL) | ((_km_f & 0x00FF00FF00FF00FFUL) << 8);
									_km_f = ((_km_f >> 16) & 0x0000FFFF0000FFFFUL) | ((_km_f & 0x0000FFFF0000FFFFUL) << 16);
									_km_f = ( _km_f >> 32                        ) | ( _km_f                        << 32);
									_km_f = (((uint64_t)-1) - _km_f) >> (64 - (m_k << 1));
									i++;
									continue;
								}
								cpt++;  
								i++;    
								continue;
							}
							if (m_table[c[i]] == -10)
							{
								_km_r = 0; cpt = 0; _isfull = false; 
								// Pass third line
								i++;
								while (len != 0)
								{
									while(i < len && c[i] != '\n')		{	i++;}
									if (i < len)   		{   break;}
									len = fread(&c, 1, MAXRSIZE, fd );
									i = 0;
								}
								i++;
								// Pass fouth line
								while (len != 0)
								{
									while(i < len && c[i] != '\n')         {	i++;}
									if (i < len)         {   break;}
									len = fread(&c, 1, MAXRSIZE, fd );
									i = 0;
								}
								i++;
								// Pass first line
								while (len != 0)
								{
									while(i < len && c[i] != '\n')        {		i++;}
									if (i < len){   break;}
									len = fread(&c, 1, MAXRSIZE, fd );
									i = 0;
								}
								i++;    
								continue;
							}
							if (m_table[c[i]] == -1)
							{ 	nt++;   
								_km_r = 0; cpt = 0; _isfull = false;
								i++;    
								continue;
							}
							cerr << m_targetsID[t].first << ": " ;
							cerr << "Failed to process sequence -- wrong character found [ASCII]: " << (size_t) c << endl;
						}
						len = fread(&c, 1, MAXRSIZE, fd );
						i = 0;  
					}
					fclose(fd);
					cerr << "\r Progress report: (" <<t+1<< "/"<<m_targetsID.size()<<")              ";
					continue;
				}
				fseek(fd, 0, 0);
				string s_kmer = "";
				ITYPE val;
				while (getFirstAndSecondElementInLine(fd, s_kmer, val))
				{	if (s_kmer.size() >= m_kmerSize && val > m_minCountTarget)
					{	commonKmersHT.addElement(s_kmer, m_targetsID[t].second, (size_t) val);}
				}
				fclose(fd);
				cerr << "\r Progress report: (" <<t+1<< "/"<<m_targetsID.size()<<")              ";
			}
		}
		cerr << nt << " nt read in total." << endl;
		cerr << "Mother Hashtable successfully built. " << commonKmersHT.Size() << " " << m_kmerSize << "-mers stored." <<  endl;
		size_t sizeMotherTable = commonKmersHT.Size();

		commonKmersHT.SortAllHashTable(2);
		commonKmersHT.RemoveCommon(m_labels_c, m_minCountTarget);
			std::string cfname = getdbNameString();
			cerr << "Creating database in disk..." << endl;
			uint64_t nbElement = commonKmersHT.Write(cfname.c_str(),2);
			cerr << nbElement << " " << m_kmerSize << "-mers successfully stored in database." << endl;

		return sizeMotherTable;
	}
	///////////////////////////////////////////////////////////////////////////////
	EHashtable<HKMERr, Element> commonKmersHT(m_kmerSize, m_labels, m_labels_c);
	for(size_t t = 0 ; t < m_targetsID.size(); t++)
	{
		FILE* fd = fopen(m_targetsID[t].first.c_str(),"r");
		if (fd == NULL)
		{
			cerr << "Failed to open " << m_targetsID[t].first << endl;
		}
		else
		{
			ILBL tgt_id;
			commonKmersHT.getTargetID(m_targetsID[t].second, tgt_id);
			char c[MAXRSIZE];
			size_t len = fread(&c, 1, MAXRSIZE, fd), i = 0;
			uint64_t _km_f = 0, _km_r = 0;
			bool _isfull = false;
			uint8_t cpt = 0;
			if (len > 0 && c[0] == '>')
			{
				while (len != 0)
				{
					while( i < len)
					{
						if (m_table[c[i]] >= 0)
						{
							nt++;
							if (_isfull)
							{
								_km_f >>= 2;
								_km_f += m_powerTable[cpt][m_table[c[i]]];

								/////////////////////////////////////////////////////////////////////
								commonKmersHT.addElement(_km_f, m_targetsID[t].second, tgt_id, 1);
								///////////////////////////////////////////////////////////////////////
								i++;
								continue;
							}
							_km_r <<= 2;    
							_km_r += 3 - m_table[c[i]];
							if (cpt == m_kmerSize - 1)
							{       _isfull = true;
								//[/////////////////////////////////////////////////////////////////////
								commonKmersHT.addElement(_km_r, m_targetsID[t].second, tgt_id, 1);
								//]//////////////////////////////////////////////////////////////////////
								_km_f = _km_r;
								// The following 6 lines come from Jellyfish source code
								_km_f = ((_km_f >> 2)  & 0x3333333333333333UL) | ((_km_f & 0x3333333333333333UL) << 2);
								_km_f = ((_km_f >> 4)  & 0x0F0F0F0F0F0F0F0FUL) | ((_km_f & 0x0F0F0F0F0F0F0F0FUL) << 4);
								_km_f = ((_km_f >> 8)  & 0x00FF00FF00FF00FFUL) | ((_km_f & 0x00FF00FF00FF00FFUL) << 8);
								_km_f = ((_km_f >> 16) & 0x0000FFFF0000FFFFUL) | ((_km_f & 0x0000FFFF0000FFFFUL) << 16);
								_km_f = ( _km_f >> 32                        ) | ( _km_f                        << 32);
								_km_f = (((uint64_t)-1) - _km_f) >> (64 - (m_k << 1));
								i++;
								continue;
							}
							cpt++;  i++;
							continue;
						}
						if (m_table[c[i]] == -10)
						{       i++;
							continue;   }
						if (m_table[c[i]] == -1)
						{        nt++; _km_r = 0; cpt = 0; _isfull = false; i++;
							continue;  }
						if (m_table[c[i]] == -2)
						{
							_km_r = 0; cpt = 0; _isfull = false;
							while (len != 0)
							{
								while(i < len && c[i] != '\n')  {       i++;}
								if (i < len)            {   break;}
								len = fread(&c, 1, MAXRSIZE, fd );
								i = 0;
							}
							i++;
							continue;
						}
						cerr << m_targetsID[t].first << ": " ;
						cerr << "Failed to process sequence -- wrong character found [ASCII]: " << (size_t) c << endl;
					}
					len = fread(&c, 1, MAXRSIZE, fd );
					i = 0;
				}
				fclose(fd);
				cerr << "\r Progress report: (" <<t+1<< "/"<<m_targetsID.size()<<")    ";
				continue;
			}
			if (len > 0 && c[0] == '@')
			{
				while (len != 0)
				{
					while(i < len && c[i] != '\n')            {       i++;}
					if (i < len)      {   break;}
					len = fread(&c, 1, MAXRSIZE, fd );
					i = 0;
				}
				i++;
				while (len != 0)
				{
					while( i < len)
					{
						if (m_table[c[i]] >= 0)
						{
							nt++;
							if (_isfull)
							{
								_km_f >>= 2; 
								_km_f += m_powerTable[cpt][m_table[c[i]]];

								/////////////////////////////////////////////////////////////////////
								commonKmersHT.addElement(_km_f, m_targetsID[t].second, tgt_id, 1);
								///////////////////////////////////////////////////////////////////////
								i++;
								continue;
							}
							_km_r <<= 2;    
							_km_r += 3 - m_table[c[i]];
							if (cpt  == m_kmerSize - 1)
							{       _isfull = true;
								//[/////////////////////////////////////////////////////////////////////
								commonKmersHT.addElement(_km_r, m_targetsID[t].second, tgt_id, 1);
								//]//////////////////////////////////////////////////////////////////////
								_km_f = _km_r;
								// The following 6 lines come from Jellyfish source code
								_km_f = ((_km_f >> 2)  & 0x3333333333333333UL) | ((_km_f & 0x3333333333333333UL) << 2);
								_km_f = ((_km_f >> 4)  & 0x0F0F0F0F0F0F0F0FUL) | ((_km_f & 0x0F0F0F0F0F0F0F0FUL) << 4);
								_km_f = ((_km_f >> 8)  & 0x00FF00FF00FF00FFUL) | ((_km_f & 0x00FF00FF00FF00FFUL) << 8);
								_km_f = ((_km_f >> 16) & 0x0000FFFF0000FFFFUL) | ((_km_f & 0x0000FFFF0000FFFFUL) << 16);
								_km_f = ( _km_f >> 32                        ) | ( _km_f                        << 32);
								_km_f = (((uint64_t)-1) - _km_f) >> (64 - (m_k << 1));

								i++;
								continue;
							}
							cpt++;
							i++;
							continue;
						}
						if (m_table[c[i]] == -10)
						{
							_km_r = 0; cpt = 0; _isfull = false;
							// Pass third line
							i++;
							while (len != 0)
							{
								while(i < len && c[i] != '\n')          {       i++;}
								if (i < len)            {   break;}
								len = fread(&c, 1, MAXRSIZE, fd );
								i = 0;
							}
							i++;
							// Pass fouth line
							while (len != 0)
							{
								while(i < len && c[i] != '\n')         {       i++;}
								if (i < len)         {   break;}
								len = fread(&c, 1, MAXRSIZE, fd );
								i = 0;
							}
							i++;
							// Pass first line
							while (len != 0)
							{
								while(i < len && c[i] != '\n')        {       i++;}
								if (i < len){   break;}
								len = fread(&c, 1, MAXRSIZE, fd );
								i = 0;
							}
							i++;
							continue;
						}
						if (m_table[c[i]] == -1)
						{       nt++;  _km_r = 0; cpt = 0; _isfull = false;
							i++;
							continue;
						}
						cerr << m_targetsID[t].first << ": " ;
						cerr << "Failed to process sequence -- wrong character found [ASCII]: " << (size_t) c << endl;
					}
					len = fread(&c, 1, MAXRSIZE, fd );
					i = 0;
				}
				fclose(fd);
				cerr << "\r Progress report: (" <<t+1<< "/"<<m_targetsID.size()<<")              ";
				continue;
			}
			fseek(fd, 0, 0);
			string s_kmer = "";
			ITYPE val = 0;
			while (getFirstAndSecondElementInLine(fd, s_kmer, val))
			{       if (s_kmer.size() >= m_kmerSize &&  val > m_minCountTarget)
				{       commonKmersHT.addElement(s_kmer, m_targetsID[t].second, (size_t) val);}
			}
			fclose(fd);
			cerr << "\r Progress report: (" <<t+1<< "/"<<m_targetsID.size()<<")              ";
		}
	}
	cerr << nt << " nt read in total." << endl;
	cerr << "Mother Hashtable successfully built. " << commonKmersHT.Size() << " " << m_kmerSize << "-mers stored." <<  endl;
	size_t sizeMotherTable = commonKmersHT.Size();

	commonKmersHT.SaveIntersectionMultiple(_filesHTC, m_labels_c);
	commonKmersHT.SaveMultiple(_filesHT, m_labels);
	commonKmersHT.SortAllHashTable(2);
	commonKmersHT.RemoveCommon(m_labels_c, m_minCountTarget);
	std::string cfname = getdbNameString();
	cerr << "Creating database in disk..." << endl;
	uint64_t nbElement = commonKmersHT.Write(cfname.c_str(),2);
	cerr << nbElement << " " << m_kmerSize << "-mers successfully stored in database." << endl;

	return sizeMotherTable;
	///////////////////////////////////////////////////////////////////////////////
}

	template <typename HKMERr>
void CLARK<HKMERr>::getObjectsDataCompute(const uint8_t * _map, const size_t&  nb, FILE * _fout)
{
	size_t i_r = 0, bigSteps = 1 + nb/m_nbCPU;
	if (_map[0] == '>')
	{
#ifdef _OPENMP
#pragma omp parallel for private(i_r)
#endif
		for (i_r = 0; i_r < m_nbCPU ; i_r++)
		{
			// Variables
			uint64_t _km_f =0, _km_r = 0;
			uint64_t cudT0 = 0, cudT1 = 0, cudT2 = 0, cudT3 = 0;
			bool _isfull = false;
			ILBL h, opt_h = 0, p = 0;
			ITYPE s_best = 0, token = 1;
			ITYPE* resultTargets = &m_resultTargets[i_r].front();
			ITYPE* iTable = &m_ITables[i_r].front();
			size_t iSize = 0, i_c = 0;
			ILBL* idx = &m_Indexes[i_r].front();
			uint16_t capacity;
			size_t readsSPos, readsEPos, readLength, i = bigSteps * i_r, bigMax = bigSteps*(i_r+1);	

			while (_map[i++] != '>')
			{}
			while (true)
			{
				readLength = 1; 
				m_seqSNames[i_r].push_back(i);
				while (i < nb && m_separators[_map[++i]] == 0)
				{}
				m_seqENames[i_r].push_back(i);
				while (i < nb && _map[i++] != '\n')
				{}
				readsSPos = i;
				//readsEPos = i;
				while (i < nb && _map[i]!= '>')
				{
					while (i < nb && _map[i]!= '\n')
					{	i++;	}
					readLength--;
					readsEPos = i++;
				}
				readLength += readsEPos - readsSPos;
				readLength -= m_isPaired ? NBN : 0;
				i_c = readLength < m_kmerSize? readsEPos : readsSPos;
				// Upper-bound for the number of queries to db
				capacity = readLength - m_kmerSize + 1;
				opt_h = 0;
				m_readsLength[i_r].push_back(readLength);
				// Scores the read
				while (i_c < readsEPos)
				{
					if (m_table[_map[i_c]] >= 0)
					{
						if (_isfull)
						{
							_km_f >>= 2;
							_km_f ^= m_pTable[m_table[_map[i_c]]];

							// Query to HashTable (Thread-safe)
							if (m_fineProfile) { cudT0 = cudRdtsc(); }
							if (m_centralHt->queryElement(_km_f, h))
							{
								if (m_fineProfile)
								{
									cudT1 = cudRdtsc();
									m_matchCycles[i_r] += cudT1 - cudT0;
									m_matchCalls[i_r]++;
								}
								updateHits(resultTargets, iTable, idx, iSize, h, token);
								if (m_fineProfile)
								{
									cudT2 = cudRdtsc();
									m_hitsCycles[i_r] += cudT2 - cudT1;
									m_hitsCalls[i_r]++;
								}
								classifyBest(resultTargets[h], h, opt_h, s_best);
								if (m_fineProfile)
								{
									cudT3 = cudRdtsc();
									m_classifyCycles[i_r] += cudT3 - cudT2;
									m_classifyCalls[i_r]++;
								}
								if (resultTargets[h] > capacity)
								{       break;  }
								i_c++;
								continue;
							}
							if (m_fineProfile)
							{
								cudT1 = cudRdtsc();
								m_matchCycles[i_r] += cudT1 - cudT0;
								m_matchCalls[i_r]++;
							}
							capacity--;
							i_c++;
							continue;
						}
						_km_r <<= 2;
						_km_r ^= m_rTable[_map[i_c]];
						if ( p  == m_kmerSize - 1 )
						{
							_isfull = true;
							// Query to HashTable (Thread-safe)
							if (m_fineProfile) { cudT0 = cudRdtsc(); }
							if (m_centralHt->queryElement(_km_r, h))
							{
								if (m_fineProfile)
								{
									cudT1 = cudRdtsc();
									m_matchCycles[i_r] += cudT1 - cudT0;
									m_matchCalls[i_r]++;
								}
								updateHits(resultTargets, iTable, idx, iSize, h, token);
								if (m_fineProfile)
								{
									cudT2 = cudRdtsc();
									m_hitsCycles[i_r] += cudT2 - cudT1;
									m_hitsCalls[i_r]++;
								}
								classifyBest(resultTargets[h], h, opt_h, s_best);
								if (m_fineProfile)
								{
									cudT3 = cudRdtsc();
									m_classifyCycles[i_r] += cudT3 - cudT2;
									m_classifyCalls[i_r]++;
								}
								if (resultTargets[h] > capacity)
								{       break;  }
								i_c++;
								_km_f = _km_r;
								// The following 6 lines come from Jellyfish source code
								_km_f = ((_km_f >> 2)  & 0x3333333333333333UL) | ((_km_f & 0x3333333333333333UL) << 2);
								_km_f = ((_km_f >> 4)  & 0x0F0F0F0F0F0F0F0FUL) | ((_km_f & 0x0F0F0F0F0F0F0F0FUL) << 4);
								_km_f = ((_km_f >> 8)  & 0x00FF00FF00FF00FFUL) | ((_km_f & 0x00FF00FF00FF00FFUL) << 8);
								_km_f = ((_km_f >> 16) & 0x0000FFFF0000FFFFUL) | ((_km_f & 0x0000FFFF0000FFFFUL) << 16);
								_km_f = ( _km_f >> 32                        ) | ( _km_f                        << 32);
								_km_f = (((uint64_t)-1) - _km_f) >> (64 - (m_k << 1));

								continue;
							}
							if (m_fineProfile)
							{
								cudT1 = cudRdtsc();
								m_matchCycles[i_r] += cudT1 - cudT0;
								m_matchCalls[i_r]++;
							}
							capacity--;
							i_c++;
							_km_f = _km_r;
							// The following 6 lines come from Jellyfish source code
							_km_f = ((_km_f >> 2)  & 0x3333333333333333UL) | ((_km_f & 0x3333333333333333UL) << 2);
							_km_f = ((_km_f >> 4)  & 0x0F0F0F0F0F0F0F0FUL) | ((_km_f & 0x0F0F0F0F0F0F0F0FUL) << 4);
							_km_f = ((_km_f >> 8)  & 0x00FF00FF00FF00FFUL) | ((_km_f & 0x00FF00FF00FF00FFUL) << 8);
							_km_f = ((_km_f >> 16) & 0x0000FFFF0000FFFFUL) | ((_km_f & 0x0000FFFF0000FFFFUL) << 16);
							_km_f = ( _km_f >> 32                        ) | ( _km_f                        << 32);
							_km_f = (((uint64_t)-1) - _km_f) >> (64 - (m_k << 1));
							continue;
						}
						p++;
						i_c++;
						continue;
					}
					if (_map[i_c] == '\n')
					{       
						i_c++;  
						continue;
					}
					_km_r = 0; p = 0; _isfull = false; 
					i_c++; 
					capacity--;
				}
				m_targetsBest[i_r].push_back(opt_h);
				iSize = 0;
				token++;	
				s_best = 0;_km_r = 0;  p = 0; _isfull = false;

				if ( (i >= bigMax && _map[i] == '>') || i >= nb)
				{       
					break;
				}
				i++;
			}
		}
	}
	else  if (_map[0] == '@')
	{
		m_posReads[0] = 1;
		size_t pos1, pos2, pos3, pos4, pos5, pos6;
		for(i_r = 1; i_r < m_nbCPU ; i_r++)
		{
			size_t i = bigSteps * i_r;
			while (i < nb && _map[i++] != '\n')
			{	}
			pos1 = i;
			while (i < nb && _map[i++] != '\n')
			{       }
			pos2 = i;
			while (i < nb && _map[i++] != '\n')
			{       }
			pos3 = i;
			while (i < nb && _map[i++] != '\n')
			{       }
			pos4 = i;
			while (i < nb && _map[i++] != '\n')
			{       }
			pos5 = i;
			while (i < nb && _map[i++] != '\n')
			{       }
			pos6 = i;
			if (_map[pos1] == '@')
			{
				i = pos2;
				while (m_Letter[_map[i++]] >= 0)
				{}
				if (i == pos3 && _map[i] == '+')
				{	
					m_posReads[i_r] = pos1+1; 	
					continue;
				}	
			}
			if (_map[pos2] == '@')
			{       
				i = pos3;
				while (m_Letter[_map[i++]] >= 0)
				{}
				if (i == pos4 && _map[i] == '+')
				{       m_posReads[i_r] = pos2+1; 
					continue;
				}
			}
			if (_map[pos3] == '@')
			{
				i = pos4;
				while (m_Letter[_map[i++]] >= 0)
				{}
				if (i == pos5 && _map[i] == '+')
				{       m_posReads[i_r] = pos3+1; 
					continue;
				}
			}
			if (_map[pos4] == '@')
			{
				i = pos5;
				while (m_Letter[_map[i++]] >= 0)
				{}
				if (i == pos6 && _map[i] == '+')
				{       m_posReads[i_r] = pos4+1;
					continue;
				}
			}
		}

#ifdef _OPENMP
#pragma omp parallel for private(i_r)
#endif
		for (i_r = 0; i_r < m_nbCPU ; i_r++)
		{
			// Variables
			uint64_t _km_f =0, _km_r = 0;
			bool _isfull = false;
			ILBL h, opt_h = 0, p = 0;
			ITYPE s_best = 0, token = 1;
			ITYPE* resultTargets = &m_resultTargets[i_r].front();
			ITYPE* iTable = &m_ITables[i_r].front();
			size_t iSize = 0,  i_c = 0;
			ILBL* idx = &m_Indexes[i_r].front();
			uint16_t capacity;
			size_t readsSPos, readsEPos, readLength, iNext = i_r+1 < m_nbCPU ? m_posReads[i_r+1]: nb;
			size_t i = m_posReads[i_r];

			while (true)
			{
				m_seqSNames[i_r].push_back(i);
				while (i < nb && m_separators[_map[++i]] == 0)
				{}
				m_seqENames[i_r].push_back(i);
				while (i < nb && _map[i++] != '\n')
				{}
				readsSPos = i;
				//readsEPos = i;
				while (i < nb && _map[i]!= '\n')
				{       i++;    }
				readsEPos = i++;
				readLength = readsEPos -  readsSPos;
				readLength -= m_isPaired ? NBN : 0;
				// Pass third line
				while (i < nb && _map[i++] != '\n')
				{}
				// Pass fourth line
				while (i < nb && _map[i++] != '\n')
				{}
				//////////////
				i_c = readLength < m_kmerSize? readsEPos : readsSPos;
				// Upper-bound for the number of queries to db
				capacity = readLength - m_kmerSize + 1;
				opt_h = 0;
				m_readsLength[i_r].push_back(readLength);
				// Scores the read
				while (i_c < readsEPos)
				{
					if (m_table[_map[i_c]] >= 0)
					{
						if (_isfull)
						{
							_km_f >>= 2;
							_km_f ^= m_pTable[m_table[_map[i_c]]];
							// Query to HashTable (Thread-safe)
							if (m_centralHt->queryElement(_km_f, h))
							{
								if (iTable[h] != token)
								{
									iTable[h] = token;
									resultTargets[h] = 0;
									idx[iSize++] = h;
								}
								resultTargets[h] +=  2;
								if (resultTargets[h] >= s_best)
								{       opt_h = h+1; s_best = resultTargets[h];  }
								if (resultTargets[h] > capacity)
								{       break;  }
								i_c++;
								continue;
							}
							capacity--;
							i_c++;
							continue;
						}
						_km_r <<= 2;
						_km_r ^= m_rTable[_map[i_c]];
						if ( p  == m_kmerSize - 1 )
						{
							_isfull = true;
							// Query to HashTable (Thread-safe)
							if (m_centralHt->queryElement(_km_r, h))
							{
								if (iTable[h] != token)
								{
									iTable[h] = token;
									resultTargets[h] = 0;
									idx[iSize++] = h;
								}
								resultTargets[h] +=  2;
								if (resultTargets[h] >= s_best)
								{       opt_h = h+1; s_best = resultTargets[h];  }
								if (resultTargets[h] > capacity)
								{       break;  }
								i_c++;
								_km_f = _km_r;// The following 6 lines come from Jellyfish source code
								_km_f = ((_km_f >> 2)  & 0x3333333333333333UL) | ((_km_f & 0x3333333333333333UL) << 2);
								_km_f = ((_km_f >> 4)  & 0x0F0F0F0F0F0F0F0FUL) | ((_km_f & 0x0F0F0F0F0F0F0F0FUL) << 4);
								_km_f = ((_km_f >> 8)  & 0x00FF00FF00FF00FFUL) | ((_km_f & 0x00FF00FF00FF00FFUL) << 8);
								_km_f = ((_km_f >> 16) & 0x0000FFFF0000FFFFUL) | ((_km_f & 0x0000FFFF0000FFFFUL) << 16);
								_km_f = ( _km_f >> 32                        ) | ( _km_f                        << 32);
								_km_f = (((uint64_t)-1) - _km_f) >> (64 - (m_k << 1));
								continue;
							}
							capacity--;
							i_c++;
							_km_f = _km_r;
							// The following 6 lines come from Jellyfish source code
							_km_f = ((_km_f >> 2)  & 0x3333333333333333UL) | ((_km_f & 0x3333333333333333UL) << 2);
							_km_f = ((_km_f >> 4)  & 0x0F0F0F0F0F0F0F0FUL) | ((_km_f & 0x0F0F0F0F0F0F0F0FUL) << 4);
							_km_f = ((_km_f >> 8)  & 0x00FF00FF00FF00FFUL) | ((_km_f & 0x00FF00FF00FF00FFUL) << 8);
							_km_f = ((_km_f >> 16) & 0x0000FFFF0000FFFFUL) | ((_km_f & 0x0000FFFF0000FFFFUL) << 16);
							_km_f = ( _km_f >> 32                        ) | ( _km_f                        << 32);
							_km_f = (((uint64_t)-1) - _km_f) >> (64 - (m_k << 1));
							continue;
						}
						p++;
						i_c++;
						continue;
					}
					_km_r = 0; p = 0; _isfull = false;
					i_c++;
					capacity--;
				}
				m_targetsBest[i_r].push_back(opt_h);
				iSize = 0;
				token++;	
				s_best = 0;_km_r = 0;  p = 0; _isfull = false;

				// Moving on to next read
				if ((++i) >= iNext)
				{       break;}
			}
		}		
	}
	else
	{
		cerr << "Failed to recognize the format of the file." << endl; exit(-1) ;
	}
	gettimeofday(&m_matchEndTime, NULL);
	fprintf(_fout, "Object_ID, Length, Assignment\n");
	size_t i = 0, c = 0;
	for(i_r = 0; i_r < m_nbCPU; i_r++)
	{
		for(i = 0; i < m_seqSNames[i_r].size(); i++)
		{
			for(c = m_seqSNames[i_r][i]; c <  m_seqENames[i_r][i]; c++)
			{	fprintf(_fout, "%c", _map[c]) ;			}
			fprintf(_fout, ",%lu,%s\n", m_readsLength[i_r][i], m_targetsName[m_targetsBest[i_r][i]].c_str() );
		}
		m_nbObjects += m_seqSNames[i_r].size();
	}
	return;

}

	template <typename HKMERr>
void CLARK<HKMERr>::getObjectsDataComputeFastLight(const uint8_t * _map, const size_t&  nb,  FILE * _fout)
{
	size_t i_r = 0, bigSteps = nb/m_nbCPU;
	if (_map[0] == '>')
	{
#ifdef _OPENMP
#pragma omp parallel for private(i_r)
#endif
		for (i_r = 0; i_r < m_nbCPU ; i_r++)
		{
			// Variables
			uint64_t _km_f =0, _km_r = 0;
			bool _isfull = false;
			ILBL h, opt_h = 0, p = 0;
			size_t i_c = 0, readsSPos, readsEPos, readLength;
			size_t i = bigSteps * i_r, bigMax = bigSteps*(i_r+1);

			while (_map[i++] != '>')
			{}
			while (true)
			{
				readLength = 1;
				m_seqSNames[i_r].push_back(i);
				while (i < nb && m_separators[_map[++i]] == 0)
				{}
				m_seqENames[i_r].push_back(i);
				while (i < nb && _map[i++] != '\n')
				{}
				readsSPos = i;
				while (i < nb && _map[i]!= '>')
				{
					while (i < nb && _map[i]!= '\n')
					{       i++;    }
					readLength--;
					readsEPos = i++;
				}
				readLength += readsEPos -  readsSPos;
				readLength -= m_isPaired ? NBN : 0;
				i_c = readLength < m_kmerSize? readsEPos : readsSPos;
				opt_h = 0;
				m_readsLength[i_r].push_back(readLength);
				// Scores the read
				while (i_c < readsEPos)
				{
					if (m_table[_map[i_c]] >= 0)
					{
						if (_isfull)
						{
							_km_f >>= 2;
							_km_f += m_pTable[m_table[_map[i_c]]];
							// Query to HashTable (Thread-safe)
							if (m_centralHt->queryElement(_km_f, h))
							{
								opt_h = h+1; break;
							}
							i_c++;
							continue;
						}
						_km_r <<= 2;
						_km_r ^= m_rTable[_map[i_c]];
						if ( p  == m_kmerSize - 1 )
						{
							_isfull = true;
							// Query to HashTable (Thread-safe)
							if (m_centralHt->queryElement(_km_r, h))
							{
								opt_h = h+1; break;
							}
							i_c++;
							_km_f = _km_r;
							// The following 6 lines come from Jellyfish source code
							_km_f = ((_km_f >> 2)  & 0x3333333333333333UL) | ((_km_f & 0x3333333333333333UL) << 2);
							_km_f = ((_km_f >> 4)  & 0x0F0F0F0F0F0F0F0FUL) | ((_km_f & 0x0F0F0F0F0F0F0F0FUL) << 4);
							_km_f = ((_km_f >> 8)  & 0x00FF00FF00FF00FFUL) | ((_km_f & 0x00FF00FF00FF00FFUL) << 8);
							_km_f = ((_km_f >> 16) & 0x0000FFFF0000FFFFUL) | ((_km_f & 0x0000FFFF0000FFFFUL) << 16);
							_km_f = ( _km_f >> 32                        ) | ( _km_f                        << 32);
							_km_f = (((uint64_t)-1) - _km_f) >> (64 - (m_k << 1));

							continue;
						}
						p++;
						i_c++;
						continue;
					}
					if (_map[i_c] == '\n')
					{
						i_c++;
						continue;
					}
					_km_r = 0; p = 0; _isfull = false;
					i_c++;
				}
				m_targetsBest[i_r].push_back(opt_h);
				_km_r = 0;  p = 0; _isfull = false;

				if ( (i >= bigMax && _map[i] == '>') || i >= nb)
				{       break;}
				i++;
			}
		}
	}
	else  if (_map[0] == '@')
	{
		m_posReads[0] = 1;
		size_t pos1, pos2, pos3, pos4, pos5, pos6;
		for(i_r = 1; i_r < m_nbCPU ; i_r++)
		{
			size_t i = bigSteps * i_r;
			while (i < nb && _map[i++] != '\n')
			{       }
			pos1 = i;
			while (i < nb && _map[i++] != '\n')
			{       }
			pos2 = i;
			while (i < nb && _map[i++] != '\n')
			{       }
			pos3 = i;
			while (i < nb && _map[i++] != '\n')
			{       }
			pos4 = i;
			while (i < nb && _map[i++] != '\n')
			{       }
			pos5 = i;
			while (i < nb && _map[i++] != '\n')
			{       }
			pos6 = i;
			if (_map[pos1] == '@')
			{
				i = pos2;
				while (m_Letter[_map[i++]] >= 0)
				{}
				if (i == pos3 && _map[i] == '+')
				{
					m_posReads[i_r] = pos1+1;
					continue;
				}
			}
			if (_map[pos2] == '@')
			{
				i = pos3;
				while (m_Letter[_map[i++]] >= 0)
				{}
				if (i == pos4 && _map[i] == '+')
				{       m_posReads[i_r] = pos2+1;
					continue;
				}
			}
			if (_map[pos3] == '@')
			{
				i = pos4;
				while (m_Letter[_map[i++]] >= 0)
				{}
				if (i == pos5 && _map[i] == '+')
				{       m_posReads[i_r] = pos3+1;
					continue;
				}
			}
			if (_map[pos4] == '@')
			{
				i = pos5;
				while (m_Letter[_map[i++]] >= 0)
				{}
				if (i == pos6 && _map[i] == '+')
				{       m_posReads[i_r] = pos4+1;
					continue;
				}
			}
		}
#ifdef _OPENMP
#pragma omp parallel for private(i_r)
#endif
		for (i_r = 0; i_r < m_nbCPU ; i_r++)
		{
			// Variables
			uint64_t _km_f =0, _km_r = 0;
			bool _isfull = false;
			ILBL h, opt_h = 0, p = 0;
			size_t i_c = 0, readsSPos, readsEPos, readLength, iNext = i_r+1 < m_nbCPU ? m_posReads[i_r+1]: nb;
			size_t i = m_posReads[i_r];

			while (true)
			{
				m_seqSNames[i_r].push_back(i);
				while (i < nb && m_separators[_map[++i]] == 0)
				{}
				m_seqENames[i_r].push_back(i);
				while (i < nb && _map[i++] != '\n')
				{}
				readsSPos = i;
				while (i < nb && _map[i]!= '\n')
				{       i++;    }
				readsEPos = i++;
				readLength = readsEPos -  readsSPos;
				readLength -= m_isPaired ? NBN : 0;
				// Pass third line
				while (i < nb && _map[i++] != '\n')
				{}
				// Pass fourth line
				while (i < nb && _map[i++] != '\n')
				{}
				//////////////
				i_c = readLength < m_kmerSize? readsEPos : readsSPos;
				opt_h = 0;
				m_readsLength[i_r].push_back(readLength);
				// Scores the read
				while (i_c < readsEPos)
				{
					if (m_table[_map[i_c]] >= 0)
					{
						if (_isfull)
						{
							_km_f >>= 2;
							_km_f += m_pTable[m_table[_map[i_c]]];
							// Query to HashTable (Thread-safe)
							if (m_centralHt->queryElement(_km_f, h))
							{
								opt_h = h+1; break;
							}
							i_c++;
							continue;
						}
						_km_r <<= 2;
						_km_r ^= m_rTable[_map[i_c]];
						if ( p  == m_kmerSize - 1 )
						{
							_isfull = true;
							// Query to HashTable (Thread-safe)
							if (m_centralHt->queryElement(_km_r, h))
							{
								opt_h = h+1; break;
							}
							i_c++;
							_km_f = _km_r;// The following 6 lines come from Jellyfish source code
							_km_f = ((_km_f >> 2)  & 0x3333333333333333UL) | ((_km_f & 0x3333333333333333UL) << 2);
							_km_f = ((_km_f >> 4)  & 0x0F0F0F0F0F0F0F0FUL) | ((_km_f & 0x0F0F0F0F0F0F0F0FUL) << 4);
							_km_f = ((_km_f >> 8)  & 0x00FF00FF00FF00FFUL) | ((_km_f & 0x00FF00FF00FF00FFUL) << 8);
							_km_f = ((_km_f >> 16) & 0x0000FFFF0000FFFFUL) | ((_km_f & 0x0000FFFF0000FFFFUL) << 16);
							_km_f = ( _km_f >> 32                        ) | ( _km_f                        << 32);
							_km_f = (((uint64_t)-1) - _km_f) >> (64 - (m_k << 1));
							continue;
						}
						p++;
						i_c++;
						continue;
					}
					_km_r = 0; p = 0; _isfull = false;
					i_c++;
				}
				m_targetsBest[i_r].push_back(opt_h);
				_km_r = 0;  p = 0; _isfull = false;

				// Moving on to next read
				if (++i >= iNext)
				{       break;}
			}
		}
	}
	else
	{
		cerr << "Failed to recognize the format of the file." << endl; exit(-1) ;
	}
	gettimeofday(&m_matchEndTime, NULL);
	fprintf(_fout, "Object_ID, Length, Assignment\n");
	size_t i = 0, c = 0;
	for(i_r = 0; i_r < m_nbCPU; i_r++)
	{
		for(i = 0; i < m_seqSNames[i_r].size(); i++)
		{
			for(c = m_seqSNames[i_r][i]; c <  m_seqENames[i_r][i]; c++)
			{       fprintf(_fout, "%c", _map[c]) ;                 }
			fprintf(_fout, ",%lu,%s\n", m_readsLength[i_r][i], m_targetsName[m_targetsBest[i_r][i]].c_str() );
		}
		m_nbObjects += m_seqSNames[i_r].size();
	}
	return;
}

	template <typename HKMERr>
void CLARK<HKMERr>::getObjectsDataComputeFast(const uint8_t *      _map,  const size_t&   nb, FILE * _fout)
{
	size_t i_r = 0, bigSteps = 1 + nb/m_nbCPU;
	if (_map[0] == '>')
	{
#ifdef _OPENMP
#pragma omp parallel for private(i_r)
#endif
		for (i_r = 0; i_r < m_nbCPU ; i_r++)
		{
			// Variables
			uint64_t _km_r = 0;
			ILBL h, opt_h = 0, p = 0;
			size_t i_c = 0, readsSPos, readsEPos, readLength;
			size_t i = bigSteps * i_r, bigMax = bigSteps*(i_r+1);

			while (_map[i++] != '>')
			{}
			while (true)
			{
				readLength = 1;
				m_seqSNames[i_r].push_back(i);
				while (i < nb && m_separators[_map[++i]] == 0)
				{}
				m_seqENames[i_r].push_back(i);
				while (i < nb && _map[i++] != '\n')
				{}
				readsSPos = i;
				while (i < nb && _map[i]!= '>')
				{
					while (i < nb && _map[i]!= '\n')
					{       i++;    }
					readLength--;
					readsEPos = i++;
				}
				readLength += readsEPos -  readsSPos;
				readLength -= m_isPaired ? NBN: 0;
				i_c = readLength < m_kmerSize? readsEPos : readsSPos;
				opt_h = 0;
				m_readsLength[i_r].push_back(readLength);
				// Scores the read
				while (i_c < readsEPos)
				{
					if (m_table[_map[i_c]] >= 0)
					{
						_km_r <<= 2;
						_km_r ^= m_rTable[_map[i_c]];
						if ( p  == m_kmerSize - 1 )
						{
							// Query to HashTable (Thread-safe)
							if (m_centralHt->queryElement(_km_r, h))
							{
								opt_h = h+1;
								break;
							}
							i_c++;
							_km_r = 0; p = 0;
							continue;
						}
						p++;
						i_c++;
						continue;
					}
					if (_map[i_c] == '\n')
					{
						i_c++;
						continue;
					}
					_km_r = 0; p = 0;
					i_c++;
				}
				m_targetsBest[i_r].push_back(opt_h);
				_km_r = 0;  p = 0;

				if ((i >= bigMax && _map[i] == '>' )|| i >= nb)
				{       break;}
				i++; 
			}
		}
	}
	else  if (_map[0] == '@')
	{
		m_posReads[0] = 1;
		size_t pos1, pos2, pos3, pos4, pos5, pos6;
		for(i_r = 1; i_r < m_nbCPU ; i_r++)
		{
			size_t i = bigSteps * i_r;
			while (i < nb && _map[i++] != '\n')
			{       }
			pos1 = i;
			while (i < nb && _map[i++] != '\n')
			{       }
			pos2 = i;
			while (i < nb && _map[i++] != '\n')
			{       }
			pos3 = i;
			while (i < nb && _map[i++] != '\n')
			{       }
			pos4 = i;
			while (i < nb && _map[i++] != '\n')
			{       }
			pos5 = i;
			while (i < nb && _map[i++] != '\n')
			{       }
			pos6 = i;
			if (_map[pos1] == '@')
			{
				i = pos2;
				while (m_Letter[_map[i++]] >= 0)
				{}
				if (i == pos3 && _map[i] == '+')
				{
					m_posReads[i_r] = pos1+1;
					continue;
				}
			}
			if (_map[pos2] == '@')
			{
				i = pos3;
				while (m_Letter[_map[i++]] >= 0)
				{}
				if (i == pos4 && _map[i] == '+')
				{       m_posReads[i_r] = pos2+1;
					continue;
				}
			}
			if (_map[pos3] == '@')
			{
				i = pos4;
				while (m_Letter[_map[i++]] >= 0)
				{}
				if (i == pos5 && _map[i] == '+')
				{       m_posReads[i_r] = pos3+1;
					continue;
				}
			}
			if (_map[pos4] == '@')
			{
				i = pos5;
				while (m_Letter[_map[i++]] >= 0)
				{}
				if (i == pos6 && _map[i] == '+')
				{       m_posReads[i_r] = pos4+1;
					continue;
				}
			}
		}

#ifdef _OPENMP
#pragma omp parallel for private(i_r)
#endif
		for (i_r = 0; i_r < m_nbCPU ; i_r++)
		{
			// Variables
			uint64_t _km_r = 0;
			ILBL h, opt_h = 0, p = 0;
			size_t i_c = 0, readsSPos, readsEPos, readLength, iNext = i_r+1 < m_nbCPU ? m_posReads[i_r+1]: nb;
			size_t i = m_posReads[i_r];

			while (true)
			{
				m_seqSNames[i_r].push_back(i);
				while (i < nb && m_separators[_map[++i]] == 0)
				{}
				m_seqENames[i_r].push_back(i);
				while (i < nb && _map[i++] != '\n')
				{}
				readsSPos = i;
				while (i < nb && _map[i]!= '\n')
				{       i++;    }
				readsEPos = i++;
				readLength = readsEPos -  readsSPos;
				readLength -= m_isPaired ? NBN: 0;
				// Pass third line
				while (i < nb && _map[i++] != '\n')
				{}
				// Pass fourth line
				while (i < nb && _map[i++] != '\n')
				{}
				//////////////
				i_c = readLength < m_kmerSize? readsEPos : readsSPos;
				opt_h = 0;
				m_readsLength[i_r].push_back(readLength);
				// Scores the read
				while (i_c < readsEPos)
				{
					if (m_table[_map[i_c]] >= 0)
					{
						_km_r <<= 2;
						_km_r ^= m_rTable[_map[i_c]];
						if ( p  == m_kmerSize - 1 )
						{
							// Query to HashTable (Thread-safe)
							if (m_centralHt->queryElement(_km_r, h))
							{
								opt_h = h+1;   
								break;  
							}
							_km_r = 0; p = 0;
							i_c++;
							continue;
						}
						p++;
						i_c++;
						continue;
					}
					_km_r = 0; p = 0;
					i_c++;
				}
				m_targetsBest[i_r].push_back(opt_h);
				_km_r = 0;  p = 0;

				// Moving on to next read
				if ((++i) >= iNext)
				{       break;}
			}
		}
	}
	else
	{
		cerr << "Failed to recognize the format of the file." << endl; exit(-1) ;
	}
	fprintf(_fout, "Object_ID, Length, Assignment\n");
	size_t i = 0, c = 0;
	for(i_r = 0; i_r < m_nbCPU; i_r++)
	{
		for(i = 0; i < m_seqSNames[i_r].size(); i++)
		{
			for(c = m_seqSNames[i_r][i]; c <  m_seqENames[i_r][i]; c++)
			{       fprintf(_fout, "%c", _map[c]) ;                 }
			fprintf(_fout, ",%lu,%s\n", m_readsLength[i_r][i], m_targetsName[m_targetsBest[i_r][i]].c_str() );
		}
		m_nbObjects += m_seqSNames[i_r].size();
	}
	return;
}


	template <typename HKMERr>
void CLARK<HKMERr>::getObjectsDataComputeFastSpaced(const bool& isfasta, const char* filename, const char* _fileResult)
{
	FileHandler * fdmanager;
	int longFactor = m_isLongSequence ? 10:1;
	if (isfasta)
	{
		fdmanager = new FileHandlerA(filename, m_nbCPU, MAXNBREADS/longFactor);
	}
	else
	{
		fdmanager = new FileHandlerQ(filename, m_nbCPU, MAXNBREADS/longFactor);
	}
	struct timeval requestStart, requestEnd;
	gettimeofday(&requestStart, NULL);

	if (!fdmanager->Open())
	{
		cerr << "(FileHandler) Failed to open "<< filename << endl;
		exit(1);
	}
	size_t eff_nbCPU = fdmanager->GetNbCPU();

#ifdef _OPENMP
	omp_set_num_threads(eff_nbCPU);
#endif

	printExtendedResultsHeader(_fileResult, true);
	std::vector< std::vector<uint8_t> > tReads(eff_nbCPU);

	// Definition of the object's length, based on its type (read or contig/scaffold).      
	size_t length = MAXREADLEN;
	length = m_isLongSequence?MAXLONGLEN:length;

	for(size_t i = 0; i < tReads.size() ; i++)
	{
		tReads[i].resize(length,0);
	}
	std::vector<HashTop> hStore(eff_nbCPU);

	while (fdmanager->Next())
	{
		int i_r = 0;
		m_nbObjects = fdmanager->GetCurrReadsCount();

		if (fdmanager->isStart())
		{
			m_objectsName.resize(m_nbObjects, "");
			m_objectsNorm.resize(m_nbObjects, 0);
			m_ResultsCentral.resize(5);
			if (m_isExtended)
			{       m_scoresLines.resize(m_nbObjects,"");   }
			for(size_t i = 0; i < m_ResultsCentral.size(); i++)
			{
				m_ResultsCentral[i].resize(m_nbObjects,0);
			}
		}
#ifdef _OPENMP
#pragma omp parallel for private(i_r)
#endif
		for(i_r = 0; i_r < eff_nbCPU; i_r++)
		{
			// Variables
			ILBL 		h 	= 0, h1 = 0;
			uint64_t 	i_c 	= 0, rid = 0;
			uint8_t* 	read 	= &tReads[i_r].front();
			uint32_t 	size 	= 0;
			bool 		stat 	= true;

			while (!fdmanager->isOver(i_r))
			{
				rid 	= fdmanager->GetReadID(i_r);
				size 	= 0;
				stat 	= fdmanager->GetRead(i_r, read, size, m_objectsName[rid]);
				m_objectsNorm[rid] = m_isPaired?size-NBN:size;
				size 	= size > 2*m_kmerSize? 2*m_kmerSize: size;
				if ( stat && size >= m_kmerSize)
				{
					i_c = 0;
					// Scores the read
					while (i_c < size - m_kmerSize + 1)
					{
						for(size_t ht = 0; ht < m_DSS.size(); ht++)
						{
							if (m_centralHt->querySpacedElement(read, i_c, size, h, ht))
							{	
								hStore[i_r].insert(h);
							}
						}
						i_c += 2;
					}
				}
				hStore[i_r].getBest(m_ResultsCentral[0][rid], m_ResultsCentral[1][rid]);
				hStore[i_r].getSecondBest(m_ResultsCentral[2][rid], m_ResultsCentral[3][rid]);
				m_ResultsCentral[1][rid] = (m_ResultsCentral[1][rid]>4 && 2*m_ResultsCentral[3][rid]<m_ResultsCentral[1][rid])?m_ResultsCentral[1][rid]:0;
				hStore[i_r].next();
			}
		}
		printExtendedSFResults(_fileResult);
	}
	m_nbObjects = fdmanager->GetReadsCount();
	
	delete fdmanager;

	gettimeofday(&requestEnd, NULL);
	printSpeedStats(requestEnd,requestStart, _fileResult);

	return;
}

	template <typename HKMERr>
void CLARK<HKMERr>::getObjectsDataComputeFull(const bool& isfasta, const char* filename, const char* _fileResult)
{
	FileHandler * fdmanager;
	int longFactor = m_isLongSequence ? 10:1;
	if (isfasta)
	{
		fdmanager = new FileHandlerA(filename, m_nbCPU, MAXNBREADS / longFactor);
	}
	else
	{
		fdmanager = new FileHandlerQ(filename, m_nbCPU, MAXNBREADS / longFactor);
	}

	struct timeval requestStart, requestEnd;
	gettimeofday(&requestStart, NULL);

	if (!fdmanager->Open())
	{
		cerr << "(FileHandler) Failed to open "<< filename << endl;
		exit(1);
	}

	size_t eff_nbCPU = fdmanager->GetNbCPU();

#ifdef _OPENMP
	omp_set_num_threads(eff_nbCPU);
#endif

	printExtendedResultsHeader(_fileResult);
	std::vector< std::vector<uint8_t> > tReads(eff_nbCPU);

	// Definition of the object's length, based on its type (read or contig/scaffold).	
	size_t length = m_isSpacedLoading?MAXREADLEN:MAXSCAFLEN;
	length = m_isLongSequence?MAXLONGLEN:length;

	for(size_t i = 0; i < tReads.size() ; i++)
	{
		tReads[i].resize(length,0);
	}	
	std::vector<HashTop> hStore(eff_nbCPU);

	while (fdmanager->Next())
	{
		int i_r = 0;
		m_nbObjects = fdmanager->GetCurrReadsCount();

		if (fdmanager->isStart())
		{
			m_objectsName.resize(m_nbObjects, "");
			m_objectsNorm.resize(m_nbObjects, 0);
			m_ResultsCentral.resize(5);
			if (m_isExtended)
			{	m_scoresLines.resize(m_nbObjects,"");	}
			for(size_t i = 0; i < m_ResultsCentral.size(); i++)
			{
				m_ResultsCentral[i].resize(m_nbObjects,0);
			}
		}
		if (m_isSpacedLoading)
		{
#ifdef _OPENMP
#pragma omp parallel for private(i_r)
#endif
			for(i_r = 0; i_r < eff_nbCPU; i_r++)
			{
				// Variables
				ILBL 		h 	= 0;
				uint64_t 	i_c 	= 0, rid = 0;
				uint8_t* 	read 	= &tReads[i_r].front();
				uint32_t 	size 	= 0;
				bool 		stat 	= true, qu = false;

				while (!fdmanager->isOver(i_r))
				{
					rid 	= fdmanager->GetReadID(i_r);
					size 	= 0;
					stat 	= fdmanager->GetRead(i_r, read, size, m_objectsName[rid]);
					m_objectsNorm[rid] = m_isPaired?size-NBN:size;
					if ( stat && size >= m_kmerSize)
					{
						i_c = 0;
						// Scores the read
						while (i_c < size - m_kmerSize + 1)
						{
							for(size_t ht = 0; ht < m_DSS.size(); ht++)
							{
								if (m_centralHt->querySpacedElement(read, i_c, size, h, ht))
								{
									hStore[i_r].insert(h);
								}
							}
							i_c++;
						}
					}
					hStore[i_r].getBest(m_ResultsCentral[0][rid], m_ResultsCentral[1][rid]);
					hStore[i_r].getSecondBest(m_ResultsCentral[2][rid], m_ResultsCentral[3][rid]);
					hStore[i_r].getTotal(m_ResultsCentral[4][rid]);
					if (m_isExtended)
					{	hStore[i_r].getScoresLine(m_targetsName.size()-1,m_scoresLines[rid]);	}
					hStore[i_r].next();
				}
			}
			printExtendedSResults(_fileResult);
		}
		else
		{
#ifdef _OPENMP
#pragma omp parallel for private(i_r)
#endif
			for(i_r = 0; i_r < eff_nbCPU; i_r++)
			{
				// Variables
				uint64_t 	_km_f 	= 0, _km_r = 0, i_c = 0;
				ILBL 		h, p 	= 0;
				bool 		isfull 	= false, stat = true;
				uint32_t 	size 	= 0;
				uint64_t 	rid 	= 0;
				uint8_t*	read	= &tReads[i_r].front();

				while (!fdmanager->isOver(i_r))
				{
					rid = fdmanager->GetReadID(i_r);
					size = 0;
					stat = fdmanager->GetRead(i_r, read, size, m_objectsName[rid]);
					m_objectsNorm[rid] = m_isPaired?size-NBN:size;

					if ( stat && size >= m_kmerSize)
					{
						i_c = 0;
						// Scores the read
						while (i_c < size)
						{
							//////////////////////////////////////////////////////////////////////////////
							if (m_table[read[i_c]] >= 0)
							{
								if (isfull)
								{
									_km_f >>= 2;
									_km_f += m_pTable[m_table[read[i_c]]];
									// Query to HashTable (Thread-safe)
									if (m_centralHt->queryElement(_km_f, h))
									{
										hStore[i_r].insert(h);
									}
									i_c++;
									continue;
								}
								_km_r <<= 2;
								_km_r ^= m_rTable[read[i_c]];
								if ( p  == m_kmerSize - 1 )
								{
									isfull = true;
									// Query to HashTable (Thread-safe)
									if (m_centralHt->queryElement(_km_r, h))
									{
										hStore[i_r].insert(h);
									}
									i_c++;
									_km_f = _km_r;
									// The following 6 lines come from Jellyfish source code
									_km_f = ((_km_f >> 2)  & 0x3333333333333333UL) | ((_km_f & 0x3333333333333333UL) << 2);
									_km_f = ((_km_f >> 4)  & 0x0F0F0F0F0F0F0F0FUL) | ((_km_f & 0x0F0F0F0F0F0F0F0FUL) << 4);
									_km_f = ((_km_f >> 8)  & 0x00FF00FF00FF00FFUL) | ((_km_f & 0x00FF00FF00FF00FFUL) << 8);
									_km_f = ((_km_f >> 16) & 0x0000FFFF0000FFFFUL) | ((_km_f & 0x0000FFFF0000FFFFUL) << 16);
									_km_f = ( _km_f >> 32                        ) | ( _km_f                        << 32);
									_km_f = (((uint64_t)-1) - _km_f) >> (64 - (m_k << 1));

									continue;
								}
								p++;
								i_c++;
								continue;
							}
							_km_r = 0; p = 0; isfull = false;
							i_c++;
							///////////////////////////////////////////////////////////////////////////////
						}
						_km_r = 0; p = 0; isfull = false;
					}
					hStore[i_r].getBest(m_ResultsCentral[0][rid], m_ResultsCentral[1][rid]);
					hStore[i_r].getSecondBest(m_ResultsCentral[2][rid], m_ResultsCentral[3][rid]);
					hStore[i_r].getTotal(m_ResultsCentral[4][rid]);
					if (m_isExtended)
					{       hStore[i_r].getScoresLine(m_targetsName.size()-1,m_scoresLines[rid]);  }
					hStore[i_r].next();
				}
			}
			printExtendedResults(_fileResult);
		}
		//printExtendedResults(_fileResult);
	}
	m_nbObjects = fdmanager->GetReadsCount();
	delete fdmanager;

	gettimeofday(&requestEnd, NULL);
	printSpeedStats(requestEnd,requestStart, _fileResult);

	return;
}

	template <typename HKMERr>
bool CLARK<HKMERr>::getObjectsDataSpectrum(FILE * fileToScore)
{
	string subfile = "";
	while(getLineFromFile(fileToScore, subfile))
	{
		m_objectsName.push_back(subfile);
	}
	return true;
}

	template <typename HKMERr>
bool CLARK<HKMERr>::getTargetsData(const char* _filesName, vector<string>& _filesHT, vector<string>& _filesHTC, const bool& _creatingkmfiles, const ITYPE& _samplingfactor)
{
	FILE * meta_f = fopen(_filesName,"r");
	if (meta_f == NULL)
	{
		cerr << "Failed to open targets data in file: " << _filesName << endl;
		exit(-1);
	}

	string subfile = "";
	bool areHTfilespresent = false;
	
	std::map<std::string,size_t> taxidToidx;
	std::map<std::string,size_t>::iterator it;
	size_t cptl = 0;

	std::vector<char> sep;
	sep.push_back(' ');
	sep.push_back('\t');
	sep.push_back('\r');
	std::vector<std::string> ele;

	while(getLineFromFile(meta_f, subfile))
	{
		ele.clear();
		getElementsFromLine(subfile, sep, ele);
		pair< string, string > target;

		target.first = ele[0];

		if (ele.size() > 1)
		{
			target.second = ele[1];
			if (ele[1].size() > 25)
			{
				std::cerr <<  "Please choose a shorter target identifier (i.e., less than 25 characters) for "<< ele[1] << " in line: " << subfile <<std::endl;
				exit(0);
			}
			it = taxidToidx.find(ele[1]); 
			if ( it  == taxidToidx.end())
			{	
				m_labels.push_back(ele[1]); 
				taxidToidx[ele[1]] = cptl++;
			}
		}
		else
		{	cerr << " Missing label for " << ele[0] <<  endl; exit(-1); }

		if (ele.size() > 2)
		{
			std::vector<string>::iterator it = find (m_labels_c.begin(), m_labels_c.end(), ele[2]);
			if ( it  == m_labels_c.end())
			{       m_labels_c.push_back(ele[2]);} 
		}
		m_targetsID.push_back(target);
	}
	fclose(meta_f);
	print(_creatingkmfiles, _samplingfactor);

	if (_creatingkmfiles)
	{
		createTargetFilesNames(_filesHT, _filesHTC);
	}

	size_t sHTS = m_labels_c.size() + m_labels.size();
	if (m_isSpacedLoading)
	{
		if (sHTS > MAXTSSM)
		{
			std::cerr << "The number of targets exceeds the limit (" << MAXTSSM << "). " << std::endl;
			std::cerr << "Please modify the targets definition to satisfy this limit." << std::endl;
			std::cerr << "The program must exit." << std::endl;
			exit(1);
		}
	}
	else
	{
		if (sHTS > MAXTS)
		{
			std::cerr << "The number of targets exceeds the limit (" << MAXTS << "). " << std::endl;
			std::cerr << "Please modify the targets definition to satisfy this limit." << std::endl;
			std::cerr << "The program must exit." << std::endl;
			exit(1);
		}
	}
	m_targetsName.push_back("NA");
	for(size_t t=0 ; t < m_labels.size(); t++)
	{       m_targetsName.push_back(m_labels[t]);       }
	for(size_t t=0 ; t < m_labels_c.size(); t++)
	{       m_targetsName.push_back(m_labels_c[t]);       }

	std::string cfname;
	if (!m_isSpacedLoading)
	{	cfname = getdbNameString();	}
	else
	{	cfname = getdbNameString(1);	}
	
	std::string cfname_ky = cfname + ".ky";
	std::string cfname_lbl = cfname + ".lb";
	std::string cfname_sz = cfname + ".sz";
	FILE * dbfd_ky = fopen(cfname_ky.c_str(), "r");
	FILE * dbfd_lbl = fopen(cfname_lbl.c_str(), "r");
        FILE * dbfd_sz = fopen(cfname_sz.c_str(), "r");

	if (dbfd_sz != NULL && dbfd_ky != NULL && dbfd_lbl != NULL)
	{
		std::ifstream in(cfname_sz.c_str(), std::ios::binary | std::ios::ate);
        	size_t fileSize = in.tellg();

		if (fileSize == HTSIZE || (m_isLightLoading && fileSize == LHTSIZE))
		{
			areHTfilespresent = true;
		}
		else
		{
			cerr << "It seems the database files are incomplete or obsolete. The program will rebuild them." << endl;
			areHTfilespresent = false;
		}
		fclose(dbfd_sz);
		fclose(dbfd_ky);
		fclose(dbfd_lbl);
		in.close();
	}
	else
	{
		cerr << "The program did not find the database files for the provided settings and reference sequences (" << sHTS << " targets). ";
		cerr << "The program will build them." << endl;
	}
	if (!areHTfilespresent && m_isSpacedLoading)
	{
		cerr << "The program cannot find databases of discriminative spaced k-mers." << endl;
		cerr << "Please run the script ``buildSpacedDB.sh'' before running CLARK to compute these databases."<< endl;	
		cerr << "See instructions in http://clark.cs.ucr.edu/Overview/"<< endl;
		exit(1);
	}
	return areHTfilespresent;
}

template <typename HKMERr>
void CLARK<HKMERr>::print(const bool& _creatingkmfiles, const ITYPE& _samplingfactor) const
{
	cerr << "CLARK version " << VERSION << " (Copyright @The Regents of the University of California. All rights reserved.)" << endl;
	if (m_minCountTarget > 0)
	{
		cerr << "Minimum k-mers occurences in Targets is set to " << m_minCountTarget << endl;
	}
	if (_creatingkmfiles)
	{
		cerr << "Creation of targets specific k-mers files requested " << endl;
	}
	if (m_isLightLoading)
	{
		cerr << "Using light database in RAM (" << m_iterKmers << ")" << endl;
	}
}

template <typename HKMERr>
void CLARK<HKMERr>::printSpeedStats(const struct timeval& _requestEnd, const struct timeval& _requestStart, const char* _fileResult) const 
{
	double diffSec = (_requestEnd.tv_sec - _requestStart.tv_sec) + (_requestEnd.tv_usec - _requestStart.tv_usec) / 1000000.0;
	cout <<" - Assignment time: "<<(diffSec * 1e9)<<" ns. Speed: ";
	cout << (size_t) (((double) m_nbObjects)/diffSec)<<" objects/s. ("<< m_nbObjects<<" objects)."<<endl;
	cout <<" - Results stored in " << _fileResult << endl;
}

template <typename HKMERr>
void CLARK<HKMERr>::printCuDProfile(const struct timeval& _requestStart, const struct timeval& _requestEnd) const
{
	if (getenv("CLARK_CUD_PROFILE") == NULL)
	{
		return;
	}
	double matchDiffSec = (m_matchEndTime.tv_sec - _requestStart.tv_sec) + (m_matchEndTime.tv_usec - _requestStart.tv_usec) / 1000000.0;
	double writeDiffSec = (_requestEnd.tv_sec - m_matchEndTime.tv_sec) + (_requestEnd.tv_usec - m_matchEndTime.tv_usec) / 1000000.0;
	// Note: gettimeofday only has microsecond resolution, so these ns
	// values are exact multiples of 1000 -- the unit is nanoseconds, but
	// the underlying resolution is still microseconds.
	cout << "CUD_PROFILE build_ns=" << (m_buildTimeSec * 1e9)
		<< " load_ns=" << (m_loadTimeSec * 1e9)
		<< " match_ns=" << (matchDiffSec * 1e9)
		<< " write_ns=" << (writeDiffSec * 1e9)
		<< " kmer=" << m_kmerSize
		<< " nbObjects=" << m_nbObjects
		<< endl;
}

template <typename HKMERr>
void CLARK<HKMERr>::printCuDProfileFine() const
{
	if (!m_fineProfile)
	{
		return;
	}
	uint64_t matchCycles = 0, hitsCycles = 0, classifyCycles = 0;
	uint64_t matchCalls = 0, hitsCalls = 0, classifyCalls = 0;
	for (size_t t = 0; t < m_nbCPU; t++)
	{
		matchCycles += m_matchCycles[t];
		hitsCycles += m_hitsCycles[t];
		classifyCycles += m_classifyCycles[t];
		matchCalls += m_matchCalls[t];
		hitsCalls += m_hitsCalls[t];
		classifyCalls += m_classifyCalls[t];
	}
	// m_tscHz is cycles-per-second (calibrated at startup on x86; on the
	// non-x86 fallback cudRdtsc() already returns nanoseconds, so m_tscHz
	// is fixed at 1e9 and this division is a no-op).
	double matchNs = (double) matchCycles / m_tscHz * 1e9;
	double hitsNs = (double) hitsCycles / m_tscHz * 1e9;
	double classifyNs = (double) classifyCycles / m_tscHz * 1e9;
	cout << "CUD_PROFILE_FINE unit=ns"
		<< " match_ns=" << matchNs
		<< " hits_ns=" << hitsNs
		<< " classify_ns=" << classifyNs
		<< " match_calls=" << matchCalls
		<< " hits_calls=" << hitsCalls
		<< " classify_calls=" << classifyCalls
		<< " tsc_hz=" << m_tscHz
		<< endl;
}

template <typename HKMERr>
void CLARK<HKMERr>::updateHits(ITYPE* _resultTargets, ITYPE* _iTable, ILBL* _idx, size_t& _iSize, const ILBL& _h, const ITYPE& _token) const
{
	if (_iTable[_h] != _token)
	{
		_iTable[_h] = _token;
		_resultTargets[_h] = 0;
		_idx[_iSize++] = _h;
	}
	_resultTargets[_h] += 2;
}

template <typename HKMERr>
void CLARK<HKMERr>::classifyBest(const ITYPE& _score, const ILBL& _h, ILBL& _opt_h, ITYPE& _s_best) const
{
	if (_score >= _s_best)
	{
		_opt_h = _h + 1;
		_s_best = _score;
	}
}

template <typename HKMERr>
void CLARK<HKMERr>::printExtendedResultsHeader(const char* _fileResult, const bool& _sf) const
{
	ofstream fout(_fileResult);
	if (_sf)
	{
		fout << "Object_ID, Length, Assignment" << endl;
		fout.close();
		return;	
	}
	string header[] = {"Bump_found","Kmers_interval", "Length", "Gamma","1st_assignment", "score1", "2nd_assignment", "score2", "confidence"};
	size_t headerSize = 9;
	if (m_isExtended)
	{
		fout <<"Object_ID";
		for(size_t t = 1 ; t < m_targetsName.size();  t++)
		{
			fout <<"," << m_targetsName[t];
		}
		size_t i_beg = m_spectrumAnalysis ? 0: 2;

		for(size_t t = i_beg ; t < headerSize ;  t++)
		{       fout << "," << header[t]; }
		fout << endl;
		return;
	}
	fout <<"Object_ID";

	for(size_t t = 2 ; t < headerSize ;  t++)
	{      	fout << "," << header[t]; }
	fout << endl;
	fout.close();
	return;
}

template <typename HKMERr>
void CLARK<HKMERr>::printExtendedResults(const char* _fileResult) const
{
	ITYPE best = 0, s_best = 0, indexBest = 0, index_sBest = 0, total = 0;
	double gamma = 0, delta = 0;

	FILE *fout = fopen(_fileResult, "a+");

	if (m_isExtended)
	{
		for(size_t t = 0; t < m_nbObjects; t++)
		{
			indexBest 	= m_ResultsCentral[0][t];
			best 	  	= m_ResultsCentral[1][t];
			index_sBest 	= m_ResultsCentral[2][t];
			s_best 		= m_ResultsCentral[3][t];
			total		= m_ResultsCentral[4][t];

			gamma  = ((double) total)/(((double) m_objectsNorm[t] - m_kmerSize) + 1.0);
			delta = ((double) best + (double) s_best);

			delta = (delta < 0.001) ? 0: ((double) best)/(delta);
			if (m_spectrumAnalysis)
			{
				fprintf(fout,"%s,%s,[%i,%i],",m_objectsName[t].c_str(),(m_objectsData[t].BumpFound ? "Y":"N"),m_objectsData[t].MinCount,m_objectsData[t].MaxCount);
				fprintf(fout,"%s,%u,%g,%s,%u,%s,%u,%g\n",m_scoresLines[t].c_str(),m_objectsNorm[t],gamma,m_targetsName[best==0?0:indexBest+1].c_str(),best,m_targetsName[s_best==0?0:index_sBest+1].c_str(),s_best,delta); 
				continue;
			}
			fprintf(fout,"%s,%s,%u,%g,%s,%u,%s,%u,%g\n",m_objectsName[t].c_str(),m_scoresLines[t].c_str(),m_objectsNorm[t],gamma,m_targetsName[best==0?0:indexBest+1].c_str(),best,m_targetsName[s_best==0?0:index_sBest+1].c_str(),s_best,delta); 
		}
		return;
	}
	for(size_t t = 0; t < m_nbObjects; t++)
	{
		indexBest       = m_ResultsCentral[0][t];
		best            = m_ResultsCentral[1][t];
		index_sBest     = m_ResultsCentral[2][t];
		s_best          = m_ResultsCentral[3][t];
		total           = m_ResultsCentral[4][t];

		gamma  = ((double) total)/(((double) m_objectsNorm[t] - m_kmerSize) + 1.0);
		delta = ((double) best + (double) s_best);
		delta = (delta < 0.001) ? 0.: ((double) best)/(delta);

		fprintf(fout,"%s,%u,%g,%s,%u,%s,%u,%g\n",m_objectsName[t].c_str(),m_objectsNorm[t],gamma,m_targetsName[best==0?0:indexBest+1].c_str(),best,m_targetsName[s_best==0?0:index_sBest+1].c_str(),s_best,delta);
	}
	fclose(fout);
	return;
}

template <typename HKMERr>
void CLARK<HKMERr>::printExtendedSResults(const char* _fileResult) const
{
	ITYPE best = 0, s_best = 0, indexBest = 0, index_sBest = 0, total = 0;
	double gamma = 0, delta = 0;
	string bestAsg = "", sBestAsg = "";
	FILE *fout = fopen(_fileResult, "a+");

	if (m_isExtended)
	{
		for(size_t t = 0; t < m_nbObjects; t++)
		{
			indexBest       = m_ResultsCentral[0][t];
			best            = m_ResultsCentral[1][t];
			index_sBest     = m_ResultsCentral[2][t];
			s_best          = m_ResultsCentral[3][t];
			total           = m_ResultsCentral[4][t];

			gamma  = ((double) total)/(((double) m_objectsNorm[t] - m_kmerSize) + 1.0);
			delta = ((double) best + (double) s_best);
			delta = (delta < 0.001) ? 0: ((double) best)/(delta);

			bestAsg  = m_targetsName[best==0?0:indexBest+1]; 
			sBestAsg = m_targetsName[s_best==0?0:index_sBest+1];
			bestAsg = (delta<MINCFSP || gamma<MINGMSP)?"NA":bestAsg;

			fprintf(fout,"%s,%s,%u,%g,%s,%u,%s,%u,%g\n",m_objectsName[t].c_str(),m_scoresLines[t].c_str(),m_objectsNorm[t],gamma,bestAsg.c_str(),best,sBestAsg.c_str(),s_best,delta);
		}
		return;
	}
	for(size_t t = 0; t < m_nbObjects; t++)
	{
		indexBest       = m_ResultsCentral[0][t];
		best            = m_ResultsCentral[1][t];
		index_sBest     = m_ResultsCentral[2][t];
		s_best          = m_ResultsCentral[3][t];
		total           = m_ResultsCentral[4][t];

		gamma  = ((double) total)/(((double) m_objectsNorm[t] - m_kmerSize) + 1.0);
		delta = ((double) best + (double) s_best);
		delta = (delta < 0.001) ? 0.: ((double) best)/(delta);
		bestAsg  = m_targetsName[best==0?0:indexBest+1]; 
		sBestAsg = m_targetsName[s_best==0?0:index_sBest+1];
		bestAsg = (delta<MINCFSP || gamma<MINGMSP)?"NA":bestAsg;

		fprintf(fout,"%s,%u,%g,%s,%u,%s,%u,%g\n",m_objectsName[t].c_str(),m_objectsNorm[t],gamma,bestAsg.c_str(),best,sBestAsg.c_str(),s_best,delta);
	}
	fclose(fout);
	return;
}

template <typename HKMERr>
void CLARK<HKMERr>::printExtendedSFResults(const char* _fileResult) const
{
	ITYPE best 	= 0, indexBest = 0;
	string bestAsg 	= "";
	FILE *fout = fopen(_fileResult, "a+");

	for(size_t t = 0; t < m_nbObjects; t++)
	{
		indexBest       = m_ResultsCentral[0][t];
		best            = m_ResultsCentral[1][t];
		bestAsg  	= m_targetsName[best==0?0:indexBest+1];
		fprintf(fout,"%s,%u,%s\n",m_objectsName[t].c_str(),m_objectsNorm[t],bestAsg.c_str());
	}
	fclose(fout);
	return;
}
