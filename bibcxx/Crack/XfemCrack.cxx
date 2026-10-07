/**
 * @file XfemCrack.cxx
 * @brief Implementation de Material
 * @author Nicolas Tardieu
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

#include "Crack/XfemCrack.h"

#include "Crack/CrackShape.h"
#include "Supervis/CommandSyntax.h"
#include "Supervis/ResultNaming.h"

XfemCrack::XfemCrack( const std::string name, MeshPtr mesh )
    : DataStructure( name, 8, "FISS_XFEM" ),
      ListOfTables( name ),
      _mesh( mesh ),
      _auxiliaryGrid( MeshPtr() ),
      _existingCrackWithGrid( XfemCrackPtr() ),
      _discontinuityType( "Crack" ),
      _crackLipsEntity(),
      _crackTipEntity(),
      _normalLevelSetFunction( FunctionPtr() ),
      _tangentialLevelSetFunction( FunctionPtr() ),
      _crackShape( CrackShapePtr() ),
      _enrichedCells(),
      _discontinuousField( "DEPL" ),
      _enrichmentType( std::string() ),
      _enrichmentRadiusZone( 0 ),
      _enrichedLayersNumber( 0 ),

      _tangentialLevelSetField( std::make_shared< FieldOnNodesReal >( getName() + ".LTNO      " ) ),
      _normalLevelSetField( std::make_shared< FieldOnNodesReal >( getName() + ".LNNO      " ) ),
      _tangentialLevelSetGradient(
          std::make_shared< FieldOnNodesReal >( getName() + ".GRLTNO    " ) ),
      _normalLevelSetGradient( std::make_shared< FieldOnNodesReal >( getName() + ".GRLNNO    " ) ),
      _localBasis( std::make_shared< FieldOnNodesReal >( getName() + ".BASLOC    " ) ),
      _crackTipCoords( JeveuxVectorReal( getName() + ".FONDFISS" ) ),
      _crackTipBasis( JeveuxVectorReal( getName() + ".BASEFOND" ) ),
      _crackTipMultiplicity( JeveuxVectorLong( getName() + ".FONDMULT" ) ),
      _crackTipCharacteristics( JeveuxVectorReal( getName() + ".CARAFOND" ) ),
      _elementSize( JeveuxVectorReal( getName() + ".FOND.TAILLE_R" ) ),
      _enrichedNodes( JeveuxVectorLong( getName() + ".GROUP_NO_ENRI" ) ),
      _crackTipCells( JeveuxVectorLong( getName() + ".MAILFISS.CTIP" ) ),
      _heavisideCells( JeveuxVectorLong( getName() + ".MAILFISS.HEAV" ) ),
      _crackTipAndHeavisideCells( JeveuxVectorLong( getName() + ".MAILFISS.HECT" ) ),
      _nodeStatusField( std::make_shared< FieldOnNodesLong >( getName() + ".STNO" ) ),
      _crackTipCellGroup( JeveuxVectorLong( getName() + ".MAILFISS.MAFOND" ) ),
      _enrichedCellGroup( JeveuxVectorLong( getName() + ".GROUP_MA_ENRI" ) ),
      _meshName( JeveuxVectorChar8( getName() + ".MAILLAGE" ) ),
      _info( JeveuxVectorChar16( getName() + ".INFO" ) ),
      _crackTipNodeFacesField( JeveuxVectorLong( getName() + ".NOFACPTFON" ) ) {
    _tangentialLevelSetField->setDescription(
        std::make_shared< EquationNumbering >( getName() + ".LTNO.NUMEQ" ) );
    _normalLevelSetField->setDescription( _tangentialLevelSetField->getDescription() );
    _nodeStatusField->setDescription(
        std::make_shared< EquationNumbering >( getName() + ".STNO.NUMEQ" ) );
    _tangentialLevelSetGradient->setDescription(
        std::make_shared< EquationNumbering >( getName() + ".GRLT.NUMEQ" ) );
    _normalLevelSetGradient->setDescription( _tangentialLevelSetGradient->getDescription() );
    _localBasis->setDescription(
        std::make_shared< EquationNumbering >( getName() + ".BASL.NUMEQ" ) );
};

XfemCrack::XfemCrack( MeshPtr mesh ) : XfemCrack( ResultNaming::getNewResultName(), mesh ) {};

const JeveuxVectorReal XfemCrack::getCrackTipCoords() {
    _crackTipCoords->updateValuePointer();
    return _crackTipCoords;
};

const JeveuxVectorReal XfemCrack::getCrackTipBasis() {
    _crackTipBasis->updateValuePointer();
    return _crackTipBasis;
};

const JeveuxVectorLong XfemCrack::getCrackTipMultiplicity() {
    _crackTipMultiplicity->updateValuePointer();
    return _crackTipMultiplicity;
};

const std::string XfemCrack::getTipType() {
    _info->updateValuePointer();
    return ( *_info )[2].toString();
};

const JeveuxVectorLong XfemCrack::getCrackTipNodeFacesField() {
    _crackTipNodeFacesField->updateValuePointer();
    return _crackTipNodeFacesField;
};

const JeveuxVectorReal XfemCrack::getCrackFrontRadius() { return _elementSize; };

bool XfemCrack::build() {
    CommandSyntax cmdSt( "DEFI_FISS_XFEM" );
    cmdSt.setResult( ResultNaming::getCurrentName(), "FISS_XFEM" );

    SyntaxMapContainer dict;

    dict.container["MAILLAGE"] = _mesh->getName();

    if ( _auxiliaryGrid ) {
        dict.container["MAILLAGE_GRILLE"] = _auxiliaryGrid->getName();
    }
    if ( _existingCrackWithGrid ) {
        dict.container["FISS_GRILLE"] = _existingCrackWithGrid->getName();
    }

    if ( _discontinuityType == "Crack" )
        dict.container["TYPE_DISCONTINUITE"] = "FISSURE";
    if ( _discontinuityType == "Interface" )
        dict.container["TYPE_DISCONTINUITE"] = "INTERFACE";
    if ( _discontinuityType == "Cohesive" )
        dict.container["TYPE_DISCONTINUITE"] = "COHESIVE";

    dict.container["CHAM_DISCONTINUITE"] = _discontinuousField;

    SyntaxMapContainer dict2;
    if ( _crackLipsEntity.size() != 0 ) {
        dict2.container["GROUP_MA_FISS"] = _crackLipsEntity;
        if ( _crackTipEntity.size() != 0 ) {
            dict2.container["GROUP_MA_FOND"] = _crackTipEntity;
        }
    }

    if ( _normalLevelSetFunction ) {
        dict2.container["FONC_LN"] = _normalLevelSetFunction->getName();
        if ( _tangentialLevelSetFunction ) {
            dict2.container["FONC_LT"] = _tangentialLevelSetFunction->getName();
        }
    }

    if ( _crackShape->getShape() != Shape::NoShape ) {
        if ( _crackShape->getShape() == Shape::Ellipse ) {
            dict2.container["DEMI_GRAND_AXE"] = _crackShape->getSemiMajorAxis();
            dict2.container["DEMI_PETIT_AXE"] = _crackShape->getSemiMinorAxis();
            dict2.container["CENTRE"] = _crackShape->getCenter();
            dict2.container["VECT_X"] = _crackShape->getVectX();
            dict2.container["VECT_Y"] = _crackShape->getVectY();
            dict2.container["COTE_FISS"] = _crackShape->getCrackSide();
        } else if ( _crackShape->getShape() == Shape::Square ) {
            dict2.container["DEMI_GRAND_AXE"] = _crackShape->getSemiMajorAxis();
            dict2.container["DEMI_PETIT_AXE"] = _crackShape->getSemiMinorAxis();
            dict2.container["RAYON_CONGE"] = _crackShape->getFilletRadius();
            dict2.container["CENTRE"] = _crackShape->getCenter();
            dict2.container["VECT_X"] = _crackShape->getVectX();
            dict2.container["VECT_Y"] = _crackShape->getVectY();
            dict2.container["COTE_FISS"] = _crackShape->getCrackSide();
        } else if ( _crackShape->getShape() == Shape::Cylinder ) {
            dict2.container["DEMI_GRAND_AXE"] = _crackShape->getSemiMajorAxis();
            dict2.container["DEMI_PETIT_AXE"] = _crackShape->getSemiMinorAxis();
            dict2.container["CENTRE"] = _crackShape->getCenter();
            dict2.container["VECT_X"] = _crackShape->getVectX();
            dict2.container["VECT_Y"] = _crackShape->getVectY();
        } else if ( _crackShape->getShape() == Shape::Notch ) {
            dict2.container["DEMI_LONGUEUR"] = _crackShape->getHalfLength();
            dict2.container["RAYON_CONGE"] = _crackShape->getFilletRadius();
            dict2.container["CENTRE"] = _crackShape->getCenter();
            dict2.container["VECT_X"] = _crackShape->getVectX();
            dict2.container["VECT_Y"] = _crackShape->getVectY();
        } else if ( _crackShape->getShape() == Shape::HalfPlane ) {
            dict2.container["PFON"] = _crackShape->getEndPoint();
            dict2.container["DTAN"] = _crackShape->getTangent();
            dict2.container["NORMALE"] = _crackShape->getNormal();
        } else if ( _crackShape->getShape() == Shape::Segment ) {
            dict2.container["PFON_ORIG"] = _crackShape->getStartingPoint();
            dict2.container["PFON_EXTR"] = _crackShape->getEndPoint();
        } else if ( _crackShape->getShape() == Shape::HalfLine ) {
            dict2.container["PFON"] = _crackShape->getStartingPoint();
            dict2.container["DTAN"] = _crackShape->getTangent();
        } else if ( _crackShape->getShape() == Shape::Line ) {
            dict2.container["POINT"] = _crackShape->getStartingPoint();
            dict2.container["DTAN"] = _crackShape->getTangent();
        }
        dict2.container["FORM_FISS"] = _crackShape->getShapeName();
    }

    if ( _normalLevelSetField ) {
        dict2.container["CHAM_NO_LSN"] = _normalLevelSetField->getName();
        if ( _tangentialLevelSetField ) {
            dict2.container["CHAM_NO_LST"] = _tangentialLevelSetField->getName();
        }
    }

    ListSyntaxMapContainer crackDefinition;
    crackDefinition.push_back( dict2 );
    dict.container["DEFI_FISS"] = crackDefinition;

    if ( _enrichedCells.size() != 0 )
        dict.container["GROUP_MA_ENRI"] = _enrichedCells;

    if ( _enrichmentType == "GEOMETRIQUE" ) {
        dict.container["TYPE_FOND_ENRI"] = "GEOMETRIQUE";
        if ( _enrichmentRadiusZone ) {
            dict.container["RAYON_ENRI"] = _enrichmentRadiusZone;
        } else {
            dict.container["NB_COUCHES"] = _enrichedLayersNumber;
        }
    }

    if ( _junctingCracks.size() != 0 ) {
        SyntaxMapContainer dict3;
        VectorString junctingCracksNames;
        for ( std::vector< XfemCrackPtr >::iterator crack = _junctingCracks.begin();
              crack != _junctingCracks.end(); ++crack ) {
            junctingCracksNames.push_back( ( *( *crack ) ).getName() );
        }
        dict3.container["FISSURE"] = junctingCracksNames;
        dict3.container["POINT"] = _pointForJunctingCracks;
        ListSyntaxMapContainer junctionDefinition;
        junctionDefinition.push_back( dict3 );
    }

    cmdSt.define( dict );

    ASTERINTEGER op = 41;
    CALL_EXECOP( &op );

    return update_tables();
};

ModelPtr XfemCrack::enrichModelWithXfem( ModelPtr &baseModel ) {
    CommandSyntax cmdSt( "MODI_MODELE_XFEM" );

    SyntaxMapContainer dict;

    dict.container["MODELE_IN"] = baseModel->getName();
    if ( baseModel->isEmpty() )
        throw std::runtime_error( "The Model must be built first" );
    dict.container["FISSURE"] = this->getName();

    // Set some default kwd values (TODO : support other values for the kwd)
    dict.container["DECOUPE_FACETTE"] = "DEFAUT";
    dict.container["CONTACT"] = "SANS";
    dict.container["PRETRAITEMENTS"] = "AUTO";

    cmdSt.define( dict );

    // Create model and get its name
    ModelPtr newModelPtr = std::make_shared< Model >( baseModel->getMesh() );
    cmdSt.setResult( newModelPtr->getName(), "MODELE" );

    // Call  OP00113
    ASTERINTEGER op = 113;
    CALL_EXECOP( &op );

    return newModelPtr;
};
