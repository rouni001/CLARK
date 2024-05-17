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

#ifndef FILEHANDLER_HH
#define FILEHANDLER_HH

#include <vector>
#include <stdint.h>
#include <string>
#include <stdint.h>

#define MAXRDBF 100000
class FileHandler
{
	public:
	FileHandler(const char* filename,const int& _nbCPU, const size_t& _maxNbReads = MAXRDBF);
	~FileHandler();

	bool		 Next();
	bool 		 Close();
	bool 		 Getline(const int& i_cpu, std::string& out, const bool isStoring = true);
	uint64_t	 GetReadsCount() const;	
	size_t		 GetNbCPU() const;
	bool 		 isOver(const int& i_cpu) const;
	uint64_t 	 Size() const;
	uint64_t	 GetCurrReadsCount() const;
	uint64_t         GetReadID(const int& i_cpu) const;
	bool 		 isStart() const;
	bool		 isEnd() const;

	virtual bool     Open();
	virtual bool 	 GetRead(const int& i_cpu, std::string& out, std::string& id);
	virtual bool 	 GetRead(const int& i_cpu, uint8_t* out, uint32_t& size, std::string& id);

	protected:
	virtual bool    SetPositions();
	
	size_t 				nbCPU;
	const char* 			_filename;
	const size_t			maxNbReads;
	uint64_t			_fileSize;
	uint64_t			_nbReads;
	int				_idxFrag;
	FILE *				fd;
	std::vector<uint64_t> 		posReads;
	std::vector<uint64_t>		i_Pos;
	std::vector<uint64_t>		max;
	std::vector<uint64_t>           currNbReads;
	std::vector<uint64_t> 		maxPos;
	uint8_t*			_map;
	std::vector<bool>		i_PosDone;
	std::vector<uint64_t>		_rIndex;
};

#endif
