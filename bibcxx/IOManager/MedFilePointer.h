#ifndef MEDFILEPOINTER_H_
#define MEDFILEPOINTER_H_

/**
 * @file MEdFilePointer.h
 * @brief Fichier entete de la classe MEdFilePointer
 * @author Nicolas Sellenet
 * @section LICENCE
 *   Copyright (C) 1991 - 2026  EDF www.code-aster.org
 *
 *   This file is part of Code_Aster.
 *
 *   Code_Aster is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Code_Aster is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Code_Aster.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "astercxx.h"

#include <filesystem>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#ifdef ASTER_HAVE_MED

enum MedFileAccessType { MedReadOnly, MedReadWrite, MedCreate };

#include "med.h"
/**
 * @class MEdFilePointer
 * @brief Med file pointer interface
 * @author Nicolas Sellenet
 */
class MedFilePointer {
  private:
    /** @brief field name */
    std::string _name;
    /** @brief med field id (after opening) */
    med_idt _fileId = -1;
    /** @brief true if file is open */
    bool _isOpen = false;
    /** @brief true if file is open in parallel */
    bool _parallelOpen = false;
    /** @brief file verion */
    std::array< med_int, 3 > _fileVersion = { 0, 0, 0 };

    /** @brief Check that file version match lib version */
    void checkFileVersion();

  public:
    /** @brief Constructor */
    MedFilePointer() {};

    /** @brief close file */
    int close();

    /** @brief get med file id (-1 if file is not open) */
    med_idt getFileId() const;

    /** @brief true if file is open */
    bool isOpen() const { return _isOpen; };

    /** @brief true if file is open in parallel */
    bool isParallel() const { return _parallelOpen; };

    /** @brief open med file */
    int open( const std::filesystem::path &filename, const MedFileAccessType & );

    /** @brief open med file in parallel */
    int openParallel( const std::filesystem::path &filename, const MedFileAccessType & );
};

#endif
#endif /* MEDFILEPOINTER_H_ */
