/**
 * @file ElementaryMatrixConverter.cxx
 * @brief Implementation de ElementaryMatrixConverter
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

#include "LinearAlgebra/ElementaryMatrixConverter.h"

#ifdef ASTER_HAVE_MPI

ElementaryTermRealPtr transfertToParallelFEDesc( const ElementaryTermRealPtr curElemT,
                                                 const ParallelContactFEDescriptorPtr pFEDesc ) {
    const auto &fEDesc = pFEDesc->getSupportFiniteElementDescriptor();
    const auto &lielC = fEDesc->getListOfGroupsOfElementsExplorer();
    const auto &lielP = pFEDesc->getListOfGroupsOfElementsExplorer();
    const auto &cellMatching = pFEDesc->getCellMatching();
    const auto &lielMatching = pFEDesc->getListOfGroupsOfElementsMatching();
    const auto nbGrelC = fEDesc->getNumberOfGroupsOfElements();
    const auto nbGrelP = pFEDesc->getNumberOfGroupsOfElements();

    // fEDesc->debugPrint();
    // pFEDesc->debugPrint();

    AS_ASSERT( nbGrelP <= nbGrelC );

    ElementaryTermRealPtr termToReturn = std::make_shared< ElementaryTermReal >();
    termToReturn->allocate( pFEDesc, curElemT->getOption(), curElemT->getPhysicalQuantityId(),
                            curElemT->getLocalModeId() );

    if ( nbGrelP <= 0 ) {
        return termToReturn;
    }

    const auto &values = curElemT->getValues();
    values->build();
    values->updateValuePointer();

    auto &newVal = termToReturn->getValues();
    for ( int iC = 0; iC < nbGrelC; ++iC ) {
        const auto iP = lielMatching[iC];
        if ( iP >= 0 ) {
            const auto elemCIdx = iC + 1;
            const auto &curLielP = lielP[iP];
            const auto &cellNbC = lielC[iC].getNumberOfNodes();
            const auto &cellNbP = curLielP.getNumberOfNodes();
            const auto &valueNb = ( *values )[elemCIdx]->size();
            const int matrixTermNbC = valueNb / cellNbC;
            AS_ASSERT( matrixTermNbC * cellNbC == valueNb );
            const int valueNbP = matrixTermNbC * cellNbP;
            // std::cout << "ET: " << iC << ", " << iP << ", " << valueNb << ", " << matrixTermNbC
            //           << std::endl;
            // std::cout << cellNbC << ", " << cellNbP << std::endl;

            AS_ASSERT( valueNbP > 0 );
            const auto elemPIdx = iP + 1;
            newVal->allocateObject( elemPIdx, valueNbP );
            ( *values )[elemCIdx]->updateValuePointer();
            ( *newVal )[elemPIdx]->updateValuePointer();
            const auto &cellM = cellMatching[iC];
            int curPos = 0;
            for ( const auto &cellId : curLielP ) {
                const auto &cellIdC = cellM[curPos];
                for ( int iMatTerm = 0; iMatTerm < matrixTermNbC; ++iMatTerm ) {
                    const auto pos1 = curPos * matrixTermNbC + iMatTerm;
                    const auto pos2 = cellIdC * matrixTermNbC + iMatTerm;
                    ( *( *newVal )[elemPIdx] )[pos1] = ( *( *values )[elemCIdx] )[pos2];
                }
                ++curPos;
            }
        }
    }

    return termToReturn;
}

ElementaryMatrixDisplacementRealPtr
transfertToParallelFEDesc( const ElementaryMatrixDisplacementRealPtr eMIn,
                           const ParallelContactFEDescriptorPtr pFEDesc ) {
    auto eM = std::make_shared< ElementaryMatrixDisplacementReal >();
    auto elemTerms = eMIn->getElementaryTerms();
    const auto &fEDesc = pFEDesc->getSupportFiniteElementDescriptor();
    const auto &lielC = fEDesc->getListOfGroupsOfElementsExplorer();

#ifdef ASTER_DEBUG_CXX
    for ( const auto &colObj : lielC ) {
        for ( const auto &val : colObj ) {
            if ( val >= 0 )
                throw std::runtime_error( "Inconsistent situation" );
        }
    }
#endif

    for ( const auto &curElemT : elemTerms ) {
        const auto curFEDesc = curElemT->getFiniteElementDescriptor();
        if ( curFEDesc != fEDesc ) {
            throw std::runtime_error(
                "SupportFiniteElementDescriptor of "
                "ParallelContactFEDescriptor must be the same than elementary matrix "
                "FiniteElementDescriptor" );
        }

        eM->addElementaryTerm( transfertToParallelFEDesc( curElemT, pFEDesc ) );
    }
    return eM;
};

ElementaryVectorDisplacementRealPtr
transfertToParallelFEDesc( const ElementaryVectorDisplacementRealPtr eVIn,
                           const ParallelContactFEDescriptorPtr pFEDesc ) {
    auto eV = std::make_shared< ElementaryVectorDisplacementReal >( eVIn->getModel() );
    auto elemTerms = eVIn->getElementaryTerms();
    const auto &fEDesc = pFEDesc->getSupportFiniteElementDescriptor();
    const auto &lielC = fEDesc->getListOfGroupsOfElementsExplorer();

#ifdef ASTER_DEBUG_CXX
    for ( const auto &colObj : lielC ) {
        for ( const auto &val : colObj ) {
            if ( val >= 0 )
                throw std::runtime_error( "Inconsistent situation" );
        }
    }
#endif

    for ( const auto &curElemT : elemTerms ) {
        const auto curFEDesc = curElemT->getFiniteElementDescriptor();
        if ( curFEDesc != fEDesc ) {
            throw std::runtime_error(
                "SupportFiniteElementDescriptor of "
                "ParallelContactFEDescriptor must be the same than elementary matrix "
                "FiniteElementDescriptor" );
        }

        eV->addElementaryTerm( transfertToParallelFEDesc( curElemT, pFEDesc ) );
    }
    return eV;
};

#endif /* ASTER_HAVE_MPI */
