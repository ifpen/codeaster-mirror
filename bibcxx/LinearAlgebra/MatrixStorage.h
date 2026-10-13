#ifndef MATRIXSTORAGE_H_
#define MATRIXSTORAGE_H_

/**
 * @file MatrixStorage.h
 * @brief Fichier entete de la classe MatrixStorage
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

#include "DataStructures/DataStructure.h"
#include "MemoryManager/JeveuxCollection.h"
#include "MemoryManager/JeveuxVector.h"

/**
 * @class MatrixStorage
 * @brief Cette classe correspond a un stockage
 * @author Nicolas Sellenet
 */
class MatrixStorage : public DataStructure {
  private:
  public:
    /**
     * @brief Constructeur
     */
    MatrixStorage( const std::string &name ) : DataStructure( name, 19, "STOCKAGE" ) {};
};

/**
 * @class LigneDeCiel
 * @brief Cette classe correspond a un stockage ligne de ciel
 * @author Nicolas Sellenet
 */
class LigneDeCiel : public MatrixStorage {
  private:
    /** @brief Objet Jeveux '.SCBL' */
    JeveuxVectorLong _scbl;
    /** @brief Objet Jeveux '.SCDI' */
    JeveuxVectorLong _scdi;
    /** @brief Objet Jeveux '.SCDE' */
    JeveuxVectorLong _scde;
    /** @brief Objet Jeveux '.SCHC' */
    JeveuxVectorLong _schc;
    /** @brief Objet Jeveux '.SCIB' */
    JeveuxVectorLong _scib;
    /** @brief Objet Jeveux '.M2LC' */
    JeveuxVectorLong _m2lc;
    /** @brief Objet Jeveux '.LC2M' */
    JeveuxVectorLong _lc2m;

  public:
    /**
     * @brief Constructeur
     */
    LigneDeCiel( const std::string &name )
        : MatrixStorage( name ),
          _scbl( JeveuxVectorLong( getName() + ".SCBL" ) ),
          _scdi( JeveuxVectorLong( getName() + ".SCDI" ) ),
          _scde( JeveuxVectorLong( getName() + ".SCDE" ) ),
          _schc( JeveuxVectorLong( getName() + ".SCHC" ) ),
          _scib( JeveuxVectorLong( getName() + ".SCIB" ) ),
          _m2lc( JeveuxVectorLong( getName() + ".M2LC" ) ),
          _lc2m( JeveuxVectorLong( getName() + ".LC2M" ) ) {};

    void updateValuePointers() const {
        _scbl->updateValuePointer();
        _scdi->updateValuePointer();
        _scde->updateValuePointer();
        _schc->updateValuePointer();
        _scib->updateValuePointer();
        _m2lc->updateValuePointer();
        _lc2m->updateValuePointer();
    };
};

/**
 * @typedef LigneDeCielPtr
 * @brief Pointeur intelligent vers un LigneDeCiel
 */
typedef std::shared_ptr< LigneDeCiel > LigneDeCielPtr;

/**
 * @class MorseStorage
 * @brief Cette classe correspond a un stockage ligne de ciel
 * @author Nicolas Sellenet
 */
class MorseStorage : public MatrixStorage {
  private:
    /** @brief Objet Jeveux '.SMDI' */
    JeveuxVectorLong _smdi;
    /** @brief Objet Jeveux '.SMDE' */
    JeveuxVectorLong _smde;
    /** @brief Objet Jeveux '.SMHC' */
    JeveuxVectorShort _smhc;

  public:
    /**
     * @brief Constructeur
     */
    MorseStorage( const std::string &name )
        : MatrixStorage( name ),
          _smdi( JeveuxVectorLong( getName() + ".SMDI" ) ),
          _smde( JeveuxVectorLong( getName() + ".SMDE" ) ),
          _smhc( JeveuxVectorShort( getName() + ".SMHC" ) ) {};

    JeveuxVectorLong getRows() const { return _smdi; };

    JeveuxVectorShort getDiagonalPositions() const { return _smhc; };

    void updateValuePointers() const {
        _smdi->updateValuePointer();
        _smde->updateValuePointer();
        _smhc->updateValuePointer();
    };
};

/**
 * @typedef MorseStoragePtr
 * @brief Pointeur intelligent vers un MorseStorage
 */
typedef std::shared_ptr< MorseStorage > MorseStoragePtr;

/**
 * @class MultFrontStorage
 * @brief Cette classe correspond a un stockage ligne de ciel
 * @author Nicolas Sellenet
 */
class MultFrontStorage : public MatrixStorage {
  private:
    /** @brief Objet Jeveux '.ADNT' */
    JeveuxVectorShort _adnt;
    /** @brief Objet Jeveux '.GLOB' */
    JeveuxVectorShort _glob;
    /** @brief Objet Jeveux '.LOCL' */
    JeveuxVectorShort _locl;
    /** @brief Objet Jeveux '.PNTI' */
    JeveuxVectorShort _pnti;

  public:
    /**
     * @brief Constructeur
     */
    MultFrontStorage( const std::string &name )
        : MatrixStorage( name ),
          _adnt( JeveuxVectorShort( getName() + ".ADNT" ) ),
          _glob( JeveuxVectorShort( getName() + ".GLOB" ) ),
          _locl( JeveuxVectorShort( getName() + ".LOCL" ) ),
          _pnti( JeveuxVectorShort( getName() + ".PNTI" ) ) {};

    void updateValuePointers() const {
        _adnt->updateValuePointer();
        _glob->updateValuePointer();
        _locl->updateValuePointer();
        _pnti->updateValuePointer();
    };
};

/**
 * @typedef MultFrontStoragePtr
 * @brief Pointeur intelligent vers un MultFrontStorage
 */
typedef std::shared_ptr< MultFrontStorage > MultFrontStoragePtr;

#endif /* MATRIXSTORAGE_H_ */
