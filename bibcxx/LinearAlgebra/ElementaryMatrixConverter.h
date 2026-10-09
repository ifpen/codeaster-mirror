#ifndef ELEMENTARYMATRIXCONVERTER_H_
#define ELEMENTARYMATRIXCONVERTER_H_

/**
 * @file ElementaryMatrixConverter.h
 * @brief Elementary matrix transfert
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

#include "aster_fort_calcul.h"

#ifdef ASTER_HAVE_MPI

#include "LinearAlgebra/ElementaryMatrix.h"
#include "LinearAlgebra/ElementaryVector.h"
#include "Modeling/ParallelContactFEDescriptor.h"

ElementaryTermRealPtr transfertToParallelFEDesc( const ElementaryTermRealPtr,
                                                 const ParallelContactFEDescriptorPtr );
ElementaryMatrixDisplacementRealPtr
transfertToParallelFEDesc( const ElementaryMatrixDisplacementRealPtr,
                           const ParallelContactFEDescriptorPtr );

ElementaryVectorDisplacementRealPtr
transfertToParallelFEDesc( const ElementaryVectorDisplacementRealPtr,
                           const ParallelContactFEDescriptorPtr );

#endif /* ASTER_HAVE_MPI */

#endif /* ELEMENTARYMATRIXCONVERTER_H_ */
