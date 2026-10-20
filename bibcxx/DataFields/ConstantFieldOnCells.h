#ifndef CONSTANTFIELDONCELLS_H_
#define CONSTANTFIELDONCELLS_H_

/**
 * @file ConstantFieldOnCells.h
 * @brief Fichier entete de la classe ConstantFieldOnCells
 * @author Natacha Bereux
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
#include "aster_fort_ds.h"
#include "aster_fort_utils.h"

#include "DataFields/DataField.h"
#include "MemoryManager/JeveuxCollection.h"
#include "MemoryManager/JeveuxVector.h"
#include "Meshes/BaseMesh.h"
#include "Meshes/MeshEntities.h"
#include "Modeling/FiniteElementDescriptor.h"
#include "Modeling/PhysicalQuantityManager.h"
#include "Supervis/ResultNaming.h"

#include <stdexcept>
#include <string>

#include <assert.h>

/**
 * @class ConstantFieldOnZone Constant Field Zone
 * @author Natacha Bereux
 */
class ConstantFieldOnZone {
  public:
    enum LocalizationType {
        AllMesh,
        AllDelayedCells,
        OnGroupOfCells,
        ListOfCells,
        ListOfDelayedCells
    };

  private:
    BaseMeshPtr _mesh;
    FiniteElementDescriptorPtr _ligrel;
    const LocalizationType _localisation;
    GroupOfCellsPtr _grp;
    VectorLong _indexes;

  public:
    ConstantFieldOnZone( BaseMeshPtr mesh )
        : _mesh( mesh ), _localisation( AllMesh ), _grp( new GroupOfCells( "" ) ) {};

    ConstantFieldOnZone( FiniteElementDescriptorPtr ligrel )
        : _ligrel( ligrel ), _localisation( AllDelayedCells ), _grp( new GroupOfCells( "" ) ) {};

    ConstantFieldOnZone( BaseMeshPtr mesh, GroupOfCellsPtr grp )
        : _mesh( mesh ), _localisation( OnGroupOfCells ), _grp( grp ) {};

    ConstantFieldOnZone( BaseMeshPtr mesh, const VectorLong &indexes )
        : _mesh( mesh ),
          _localisation( ListOfCells ),
          _grp( new GroupOfCells( "" ) ),
          _indexes( indexes ) {};

    ConstantFieldOnZone( FiniteElementDescriptorPtr ligrel, const VectorLong &indexes )
        : _ligrel( ligrel ), _localisation( ListOfDelayedCells ), _indexes( indexes ) {};

    BaseMeshPtr getMesh() const {
        if ( _localisation != AllMesh and _localisation != OnGroupOfCells and
             _localisation != ListOfCells )
            throw std::runtime_error( "Zone not on a mesh" );
        return _mesh;
    };

    const FiniteElementDescriptorPtr &getFiniteElementDescriptor() const {
        if ( _localisation != AllDelayedCells and _localisation != ListOfDelayedCells )
            throw std::runtime_error( "Zone not on a FiniteElementDescriptor" );
        return _ligrel;
    };

    LocalizationType getLocalizationType() const { return _localisation; };

    GroupOfCellsPtr getGroup() const { return _grp; };

    const VectorLong &getListOfCells() const { return _indexes; };
};

/**
 * @class ConstantFieldValues Constant Field values
 * @author Natacha Bereux
 */
template < class ValueType >
class ConstantFieldValues {
  private:
    VectorString _components;
    std::vector< ValueType > _values;

  public:
    ConstantFieldValues( const VectorString &comp, const std::vector< ValueType > &val )
        : _components( comp ), _values( val ) {};

    const VectorString &getComponents() const { return _components; };

    const std::vector< ValueType > &getValues() const { return _values; };
};

/**
 * @class ConstantFieldOnCells Constant Field on Mesh template
 * @brief Cette classe permet de definir une carte (champ défini sur les mailles)
 * @author Natacha Bereux
 */
template < class ValueType >
class ConstantFieldOnCells : public DataField {
  private:
    /** @brief Vecteur Jeveux '.NOMA' */
    JeveuxVectorChar8 _meshName;
    /** @brief Vecteur Jeveux '.DESC' */
    JeveuxVectorLong _descriptor;
    /** @brief Vecteur Jeveux '.NOLI' */
    JeveuxVectorChar24 _nameOfLigrels;
    /** @brief Collection  '.LIMA' */
    JeveuxCollectionLong _listOfMeshCells;
    /** @brief Vecteur Jeveux '.VALE' */
    JeveuxVector< ValueType > _values;
    /** @brief Maillage sous-jacent */
    const BaseMeshPtr _mesh;
    /** @brief Ligrel */
    FiniteElementDescriptorPtr _FEDesc;
    /** @brief Objet temporaire '.NCMP' */
    JeveuxVectorChar8 _componentNames;
    /** @brief Objet temporaire '.VALV' */
    JeveuxVector< ValueType > _valuesTmp;

  private:
    void fortranAddValues( const ASTERINTEGER &code, const std::string &grp,
                           const std::string &mode, const ASTERINTEGER &nma,
                           const JeveuxVectorLong &limanu, const JeveuxVectorChar8 &component,
                           JeveuxVector< ValueType > &values ) {
        if ( ( code == -1 || code == -3 ) && !_FEDesc )
            throw std::runtime_error(
                "Build of ConstantFieldOnCells impossible, FiniteElementDescriptor is missing" );
        _componentNames->updateValuePointer();
        _valuesTmp->updateValuePointer();
        const ASTERINTEGER taille = _componentNames->size();

        const ASTERINTEGER tVerif1 = component->size();
        const ASTERINTEGER tVerif2 = values->size();
        if ( tVerif1 > taille || tVerif2 > taille || tVerif1 != tVerif2 )
            throw std::runtime_error( "Unconsistent size" );

        for ( int position = 0; position < tVerif1; ++position ) {
            ( *_componentNames )[position] = ( *component )[position];
            ( *_valuesTmp )[position] = ( *values )[position];
        }

        const std::string limano( " " );

        CALLO_NOCARTC( getName(), &code, &tVerif1, grp, mode, &nma, limano, &( *limanu )[0],
                       _FEDesc->getName() );
    };

    void fortranAddValues( const ASTERINTEGER &code, const std::string &grp,
                           const std::string &mode, const ASTERINTEGER &nma,
                           const JeveuxVectorLong &limanu, const VectorString &component,
                           const std::vector< ValueType > &values ) {
        if ( ( code == -1 || code == -3 ) && !_FEDesc )
            throw std::runtime_error(
                "Build of ConstantFieldOnCells impossible, FiniteElementDescriptor is missing" );
        if ( !_componentNames.exists() )
            _componentNames->allocate( 30 );
        _componentNames->updateValuePointer();
        if ( !_valuesTmp.exists() )
            _valuesTmp->allocate( 30 );
        _valuesTmp->updateValuePointer();
        const ASTERINTEGER taille = _componentNames->size();

        const ASTERINTEGER tVerif1 = component.size();
        const ASTERINTEGER tVerif2 = values.size();
        if ( tVerif1 > taille || tVerif2 > taille || tVerif1 != tVerif2 )
            throw std::runtime_error( "Unconsistent size" );

        for ( int position = 0; position < tVerif1; ++position ) {
            ( *_componentNames )[position] = component[position];
            ( *_valuesTmp )[position] = values[position];
        }
        std::string feDescName( " " );
        if ( _FEDesc != nullptr )
            feDescName = _FEDesc->getName();

        const std::string limano( " " );

        CALLO_NOCARTC( getName(), &code, &tVerif1, grp, mode, &nma, limano, &( *limanu )[0],
                       feDescName );
    };

    void fortranAllocate( const std::string base, const std::string quantity ) {

        CALLO_ALCART( base, getName(), _mesh->getName(), quantity );
    };

  public:
    /**
     * @typedef ConstantFieldOnBaseMeshPtr
     * @brief Pointeur intelligent vers un ConstantFieldOnCells
     */
    typedef std::shared_ptr< ConstantFieldOnCells > ConstantFieldOnBaseMeshPtr;

    /**
     * @brief Constructeur
     * @param name Nom Jeveux de la carte
     * @param mesh Maillage
     */
    ConstantFieldOnCells( const std::string &name, const BaseMeshPtr &mesh )
        : DataField( name, "CARTE" ),
          _meshName( JeveuxVectorChar8( getName() + ".NOMA" ) ),
          _descriptor( JeveuxVectorLong( getName() + ".DESC" ) ),
          _nameOfLigrels( JeveuxVectorChar24( getName() + ".NOLI" ) ),
          _listOfMeshCells( JeveuxCollectionLong( getName() + ".LIMA" ) ),
          _values( JeveuxVector< ValueType >( getName() + ".VALE" ) ),
          _mesh( mesh ),
          _FEDesc( FiniteElementDescriptorPtr() ),
          _componentNames( getName() + ".NCMP" ),
          _valuesTmp( getName() + ".VALV" ) {};

    /**
     * @brief Constructeur
     * @param name Nom Jeveux de la carte
     * @param ligrel Ligrel support
     */
    ConstantFieldOnCells( std::string name, const FiniteElementDescriptorPtr &ligrel )
        : ConstantFieldOnCells( name, ligrel->getMesh() ) {
        _FEDesc = ligrel;
    };

    /**
     * @brief Constructeur
     * @param mesh Maillage
     * @param name Nom Jeveux de la carte
     */
    ConstantFieldOnCells( const BaseMeshPtr &mesh )
        : ConstantFieldOnCells( ResultNaming::getNewResultName(), mesh ) {};

    /**
     * @brief Constructeur
     * @param ligrel Ligrel support
     * @param name Nom Jeveux de la carte
     */
    ConstantFieldOnCells( const FiniteElementDescriptorPtr &ligrel )
        : ConstantFieldOnCells( ResultNaming::getNewResultName(), ligrel ) {};

    /**
     * @brief Constructeur
     * @param ligrel Ligrel support
     * @param name Nom Jeveux de la carte
     */
    ConstantFieldOnCells( const std::string &name, const ConstantFieldOnCells &toCopy )
        : ConstantFieldOnCells( name, toCopy.getMesh() ) {
        *( _meshName ) = *( toCopy._meshName );
        *( _descriptor ) = *( toCopy._descriptor );
        *( _values ) = *( toCopy._values );
        *( _nameOfLigrels ) = *( toCopy._nameOfLigrels );
        *( _listOfMeshCells ) = *( toCopy._listOfMeshCells );
        *( _componentNames ) = *( toCopy._componentNames );
        *( _valuesTmp ) = *( toCopy._valuesTmp );
        _FEDesc = toCopy._FEDesc;

        updateValuePointers();
    };

    typedef std::shared_ptr< ConstantFieldOnCells< ValueType > > ConstantFieldOnCellsValueTypePtr;

    /**
     * @brief Destructeur
     */
    ~ConstantFieldOnCells() {};

    bool exists() const { return _meshName.exists() && _descriptor.exists() && _values.exists(); };

    /**
     * @brief Allocation de la carte
     * @return true si l'allocation s'est bien deroulee, false sinon
     */
    void allocate( const std::string componant ) {
        if ( _mesh.use_count() == 0 || _mesh->isEmpty() ) {
            AS_ABORT( "Mesh is empty" );
        }

        std::string strJeveuxBase( "G" );
        fortranAllocate( strJeveuxBase, componant );
    };

    /**
     * @brief Allocation de la carte
     * @return true si l'allocation s'est bien deroulee, false sinon
     */
    void allocate( const ConstantFieldOnCellsValueTypePtr &model ) {
        auto componant = model->getPhysicalQuantityName();
        std::string strJeveuxBase( "G" );
        fortranAllocate( strJeveuxBase, componant );
    };

    /**
     * @brief Deallocate the ConstantField
     */
    void deallocate() {
        _meshName->deallocate();
        _descriptor->deallocate();
        _nameOfLigrels->deallocate();
        _listOfMeshCells->deallocate();
        _values->deallocate();
        _componentNames->deallocate();
        _valuesTmp->deallocate();
    };

    /**
     * @brief Get support physical quantity
     */
    std::string getPhysicalQuantityName() const {
        _descriptor->updateValuePointer();
        ASTERINTEGER gdeur = ( *_descriptor )[0];
        return PhysicalQuantityManager::getPhysicalQuantityName( gdeur );
    };

    /**
     * @brief Get mesh
     */
    BaseMeshPtr getMesh() const { return _mesh; };

    /**
     * @brief Get values of a zone
     */
    ConstantFieldValues< ValueType > getValues( const int &position ) const {
        if ( position >= size() )
            throw std::runtime_error( "Out of ConstantFieldOnCells bound" );
        _values->updateValuePointer();

        ASTERINTEGER nbZoneMax = ( *_descriptor )[1];
        ASTERINTEGER gdeur = ( *_descriptor )[0];
        const auto name1 = PhysicalQuantityManager::getPhysicalQuantityName( gdeur );
        ASTERINTEGER nec = PhysicalQuantityManager::getNumberOfEncodedInteger( gdeur );
        ASTERINTEGER ndim = PhysicalQuantityManager::getNumberOfComponents( gdeur );
        const auto &compNames = PhysicalQuantityManager::getComponentNames( gdeur );
        const ASTERINTEGER nbCmpMax = compNames.size();
        VectorString cmpToReturn;
        cmpToReturn.reserve( 30 * nec );
        std::vector< ValueType > valToReturn;
        valToReturn.reserve( 30 * nec );
        for ( int i = 0; i < nec; ++i ) {
            ASTERINTEGER encodedInt = ( *_descriptor )[3 + 2 * nbZoneMax + position * nec + i];
            VectorLong vecOfComp( ndim, -1 );
            CALL_ISDECO( &encodedInt, vecOfComp.data(), &ndim );
            ASTERINTEGER pos = 0;
            for ( const auto &val : vecOfComp ) {
                if ( val == 1 ) {
                    cmpToReturn.push_back( compNames[pos + i * 30] );
                    const ASTERINTEGER posInVale = pos + i * 30 + nbCmpMax * position;
                    valToReturn.push_back( ( *_values )[posInVale] );
                }
                ++pos;
            }
        }
        return ConstantFieldValues< ValueType >( cmpToReturn, valToReturn );
    };

    /**
     * @brief Get values of all zones
     */
    std::vector< ConstantFieldValues< ValueType > > getValues() const {
        _values->updateValuePointer();
        _descriptor->updateValuePointer();
        ASTERINTEGER nbZoneMax = ( *_descriptor )[1];
        ASTERINTEGER gdeur = ( *_descriptor )[0];
        auto size = ( *_descriptor )[2];
        ASTERINTEGER nec = PhysicalQuantityManager::getNumberOfEncodedInteger( gdeur );
        ASTERINTEGER ndim = PhysicalQuantityManager::getNumberOfComponents( gdeur );
        const auto &compNames = PhysicalQuantityManager::getComponentNames( gdeur );
        const ASTERINTEGER nbCmpMax = compNames.size();

        std::vector< ConstantFieldValues< ValueType > > vectorOfConstantFieldValues;
        vectorOfConstantFieldValues.reserve( size );

        for ( int position = 0; position < size; ++position ) {
            VectorString cmpToReturn;
            cmpToReturn.reserve( 30 * nec );
            std::vector< ValueType > valToReturn;
            valToReturn.reserve( 30 * nec );
            for ( int i = 0; i < nec; ++i ) {
                ASTERINTEGER encodedInt = ( *_descriptor )[3 + 2 * nbZoneMax + position * nec + i];
                VectorLong vecOfComp( ndim, -1 );
                CALL_ISDECO( &encodedInt, vecOfComp.data(), &ndim );
                ASTERINTEGER pos = 0;
                for ( const auto &val : vecOfComp ) {
                    if ( val == 1 ) {
                        cmpToReturn.push_back( compNames[pos + i * 30] );
                        const ASTERINTEGER posInVale = pos + i * 30 + nbCmpMax * position;
                        valToReturn.push_back( ( *_values )[posInVale] );
                    }
                    ++pos;
                }
            }
            vectorOfConstantFieldValues.push_back(
                ConstantFieldValues< ValueType >( cmpToReturn, valToReturn ) );
        }
        return vectorOfConstantFieldValues;
    };

    /**
     * @brief Get zone description
     */
    ConstantFieldOnZone getZoneDescription( const int &position ) const {
        if ( position >= size() )
            throw std::runtime_error( "Out of ConstantFieldOnCells bound" );

        ASTERINTEGER code = ( *_descriptor )[3 + 2 * position];
        if ( code == 1 )
            return ConstantFieldOnZone( _mesh );
        else if ( code == -1 )
            return ConstantFieldOnZone( _FEDesc );
        else if ( code == 2 ) {
            const auto numGrp = ( *_descriptor )[4 + 2 * position];
            const auto &map = _mesh->getGroupsOfCellsMap();
            const auto name = map->getStringFromIndex( numGrp );
            return ConstantFieldOnZone( _mesh, GroupOfCellsPtr( new GroupOfCells( name ) ) );
        } else if ( code == 3 ) {
            const auto numGrp = ( *_descriptor )[4 + 2 * position];
            const auto &object = ( *_listOfMeshCells )[numGrp];
            return ConstantFieldOnZone( _mesh, object->toVector() );
        } else if ( code == -3 ) {
            const auto numGrp = ( *_descriptor )[4 + 2 * position];
            const auto &object = ( *_listOfMeshCells )[numGrp];
            return ConstantFieldOnZone( _FEDesc, object->toVector() );
        } else
            throw std::runtime_error( "Error in ConstantFieldOnCells" );
    };

    bool setValueOnCells( const VectorLong cells, const VectorString cmp,
                          const std::vector< ValueType > values ) {
        return setValueOnZone( ConstantFieldOnZone( _mesh, cells ),
                               ConstantFieldValues< ValueType >( cmp, values ) );
    };

    /**
     * @brief Fixer une valeur sur tout le maillage
     * @param component JeveuxVectorChar8 contenant le nom des composantes à fixer
     * @param values JeveuxVector< ValueType > contenant les valeurs
     * @return renvoit true si l'ajout s'est bien deroulee, false sinon
     */
    bool setValueOnMesh( const JeveuxVectorChar8 &component,
                         const JeveuxVector< ValueType > &values ) {
        if ( !_mesh || _mesh->isEmpty() ) {
            AS_ABORT( "Mesh is empty" );
        }

        const ASTERINTEGER code = 1;
        const std::string grp( " " );
        const std::string mode( " " );
        const ASTERINTEGER nbMa = 0;
        JeveuxVectorLong limanu( "empty" );
        limanu->allocate( 1 );
        fortranAddValues( code, grp, mode, nbMa, limanu, component, values );
        return true;
    };

    /**
     * @brief Fixer une valeur sur un groupe de mailles
     * @param component JeveuxVectorChar8 contenant le nom des composantes à fixer
     * @param values JeveuxVector< ValueType > contenant les valeurs
     * @param grp Groupe de mailles
     * @return renvoit true si l'ajout s'est bien deroulee, false sinon
     */
    bool setValueOnListOfDelayedCells( const JeveuxVectorChar8 &component,
                                       const JeveuxVector< ValueType > &values,
                                       const VectorLong &grp ) {
        if ( !_mesh || _mesh->isEmpty() ) {
            AS_ABORT( "Mesh is empty" );
        }

        const ASTERINTEGER code = -3;
        const std::string grp2( " " );
        const std::string mode( "NUM" );
        const ASTERINTEGER nbMa = 0;
        JeveuxVectorLong limanu( "&&TEMPORARY" );
        limanu->allocate( grp.size() );
        for ( ASTERINTEGER pos = 0; pos < grp.size(); ++pos )
            ( *limanu )[pos] = grp[pos];
        fortranAddValues( code, grp2, mode, nbMa, limanu, component, values );
        return true;
    };

    /**
     * @brief Fixer une valeur sur un groupe de mailles
     * @param component JeveuxVectorChar8 contenant le nom des composantes à fixer
     * @param values JeveuxVector< ValueType > contenant les valeurs
     * @param grp Groupe de mailles
     * @return renvoit true si l'ajout s'est bien deroulee, false sinon
     */
    bool setValueOnGroupOfCells( const JeveuxVectorChar8 &component,
                                 const JeveuxVector< ValueType > &values,
                                 const GroupOfCells &grp ) {
        if ( _mesh.use_count() == 0 || _mesh->isEmpty() ) {
            AS_ABORT( "Mesh is empty" );
        }
        if ( !_mesh->hasGroupOfCells( grp.getName() ) )
            throw std::runtime_error( "Group " + grp.getName() + " not in mesh" );

        const ASTERINTEGER code = 2;
        const std::string mode( " " );
        const ASTERINTEGER nbMa = 0;
        JeveuxVectorLong limanu( "empty" );
        limanu->allocate( 1 );
        fortranAddValues( code, grp.getName(), mode, nbMa, limanu, component, values );
        return true;
    };

    /**
     * @brief Fixer une valeur sur un groupe de mailles
     * @param zone Zone sur laquelle on alloue la carte
     * @param values Valeur a allouer
     * @return renvoit true si l'ajout s'est bien deroulee, false sinon
     */
    bool setValueOnZone( const ConstantFieldOnZone &zone,
                         const ConstantFieldValues< ValueType > &values ) {
        if ( _mesh.use_count() == 0 || _mesh->isEmpty() ) {
            AS_ABORT( "Mesh is empty" );
        }

        ASTERINTEGER code = 0;
        std::string grp( " " );
        std::string mode( " " );
        ASTERINTEGER nbMa = 0;
        JeveuxVectorLong limanu( "&&TEMPORARY" );
        if ( zone.getLocalizationType() == ConstantFieldOnZone::AllMesh ) {
            code = 1;
            limanu->allocate( 1 );
        } else if ( zone.getLocalizationType() == ConstantFieldOnZone::AllDelayedCells ) {
            code = -1;
            limanu->allocate( 1 );
        } else if ( zone.getLocalizationType() == ConstantFieldOnZone::OnGroupOfCells ) {
            code = 2;
            grp = zone.getGroup()->getName();
            limanu->allocate( 1 );
        } else if ( zone.getLocalizationType() == ConstantFieldOnZone::ListOfCells ) {
            code = 3;
            mode = "NUM";
            const auto &vecTmp = zone.getListOfCells();
            nbMa = vecTmp.size();
            // si taille 1 on doit pouvoir préallouer
            limanu->allocate( nbMa );
            for ( ASTERINTEGER pos = 0; pos < nbMa; ++pos )
                ( *limanu )[pos] = vecTmp[pos];
        } else if ( zone.getLocalizationType() == ConstantFieldOnZone::ListOfDelayedCells ) {
            code = -3;
            mode = "NUM";
            const auto &vecTmp = zone.getListOfCells();
            nbMa = vecTmp.size();
            // si taille 1 on doit pouvoir préallouer
            limanu->allocate( nbMa );
            for ( ASTERINTEGER pos = 0; pos < nbMa; ++pos )
                ( *limanu )[pos] = vecTmp[pos];
        }
        fortranAddValues( code, grp, mode, nbMa, limanu, values.getComponents(),
                          values.getValues() );
        return true;
    };

    /**
     * @brief Get number of zone in ConstantFieldOnCells
     */
    ASTERINTEGER size() const {
        _descriptor->updateValuePointer();
        return ( *_descriptor )[2];
    };

    /**
     * @brief Mise a jour des pointeurs Jeveux
     * @return true si la mise a jour s'est bien deroulee, false sinon
     */
    void updateValuePointers() {
        _meshName->updateValuePointer();
        _descriptor->updateValuePointer();
        _values->updateValuePointer();
        // Les deux elements suivants sont facultatifs
        _listOfMeshCells->updateValuePointer();
        if ( _nameOfLigrels.exists() )
            _nameOfLigrels->updateValuePointer();
    };

    bool build( bool force = false ) {
        if ( !_listOfMeshCells->isBuilt() || force ) {
            _listOfMeshCells->build();
            return true;
        }

        return false;
    };
};

/** @typedef ConstantFieldOnCellsReal Class d'une carte de double */
typedef ConstantFieldOnCells< ASTERDOUBLE > ConstantFieldOnCellsReal;
/** @typedef ConstantFieldOnCellsLong Class d'une carte de long */
typedef ConstantFieldOnCells< ASTERINTEGER > ConstantFieldOnCellsLong;
/** @typedef ConstantFieldOnCellsComplex Class d'une carte de complexe */
typedef ConstantFieldOnCells< ASTERCOMPLEX > ConstantFieldOnCellsComplex;
/** @typedef ConstantFieldOnCellsChar8 Class d'une carte de char*8 */
typedef ConstantFieldOnCells< JeveuxChar8 > ConstantFieldOnCellsChar8;
/** @typedef ConstantFieldOnCellsChar16 Class d'une carte de char*16 */
typedef ConstantFieldOnCells< JeveuxChar16 > ConstantFieldOnCellsChar16;
/** @typedef ConstantFieldOnCellsChar24 Class d'une carte de char*16 */
typedef ConstantFieldOnCells< JeveuxChar24 > ConstantFieldOnCellsChar24;

/**
 * @typedef ConstantFieldOnBaseMeshPtrReal
 * @brief   Definition d'une carte de double
 */
typedef std::shared_ptr< ConstantFieldOnCellsReal > ConstantFieldOnCellsRealPtr;

/**
 * @typedef ConstantFieldOnCellsLongPtr
 * @brief   Definition d'une carte de double
 */
typedef std::shared_ptr< ConstantFieldOnCellsLong > ConstantFieldOnCellsLongPtr;

/**
 * @typedef ConstantFieldOnBaseMeshPtrComplex
 * @brief   Definition d'une carte de complexe
 */
typedef std::shared_ptr< ConstantFieldOnCellsComplex > ConstantFieldOnCellsComplexPtr;

/**
 * @typedef ConstantFieldOnBaseMeshPtrChar8 Definition d'une carte de char[8]
 * @brief Pointeur intelligent vers un ConstantFieldOnCells
 */
typedef std::shared_ptr< ConstantFieldOnCellsChar8 > ConstantFieldOnCellsChar8Ptr;

/**
 * @typedef ConstantFieldOnBaseMeshPtrChar16 Definition d'une carte de char[16]
 * @brief Pointeur intelligent vers un ConstantFieldOnCells
 */
typedef std::shared_ptr< ConstantFieldOnCellsChar16 > ConstantFieldOnCellsChar16Ptr;

/**
 * @typedef ConstantFieldOnBaseMeshPtrChar16 Definition d'une carte de char[24]
 * @brief Pointeur intelligent vers un ConstantFieldOnCells
 */
typedef std::shared_ptr< ConstantFieldOnCellsChar24 > ConstantFieldOnCellsChar24Ptr;

#endif /* CONSTANTFIELDONCELLS_H_ */
