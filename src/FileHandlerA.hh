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

#ifndef FILEHANDLER_A_HH
#define FILEHANDLER_A_HH

#include <vector>
#include <stdint.h>
#include <string>

#include "./FileHandler.hh"

class FileHandlerA: public FileHandler 
{
	public:
	FileHandlerA(const char* filename,const int& _nbCPU, const size_t& _maxNbReads = MAXRDBF);
	~FileHandlerA();

	bool 	Open();
	bool 	GetRead(const int& i_cpu, std::string& out, std::string& id);
	bool 	GetRead(const int& i_cpu, uint8_t* out, uint32_t& size, std::string& id);
	private:
	bool    SetPositions();
};

#endif
