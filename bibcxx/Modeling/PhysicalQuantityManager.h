#ifndef PHYSICALQUANTITYMANAGER_H_
#define PHYSICALQUANTITYMANAGER_H_

/**
 * @file PhysicalQuantityManager.h
 * @brief Fichier entete de la classe PhysicalQuantityManager
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

#include "MemoryManager/JeveuxCollection.h"

/**
 * @class PhysicalQuantityManager
 * @brief Class to manage "catalogue" interactions
 * @author Nicolas Sellenet
 */
class PhysicalQuantityManager {
  private:
    static NamesMapChar8 _nameOfPhysicalQuantity;
    static JeveuxCollectionChar8 _nameOfCmp;

  public:
    static bool hasQuantityOfName( const std::string );

    static std::string getPhysicalQuantityName( const ASTERINTEGER );

    static std::string getPhysicalQuantityName( const std::string, const std::string );

    static ASTERINTEGER getPhysicalQuantityId( const std::string );

    static ASTERINTEGER getPhysicalQuantityId( const std::string, const std::string );

    static ASTERINTEGER getNumberOfEncodedInteger( const ASTERINTEGER );

    static ASTERINTEGER getNumberOfComponents( const ASTERINTEGER );

    static const VectorString getAllPhysicalQuantityNames();

    static const VectorString getComponentNames( const ASTERINTEGER );
};

#endif /* PHYSICALQUANTITYMANAGER_H_ */
