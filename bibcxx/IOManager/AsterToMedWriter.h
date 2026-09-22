#ifndef ASTERTOMEDWRITER_H_
#define ASTERTOMEDWRITER_H_

/**
 * @file AsterToMedWriter.h
 * @brief Fichier entete de la classe AsterToMedWriter
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

#include "IOManager/MedFileReader.h"
#include "Meshes/ConnectionMesh.h"
#include "Meshes/Mesh.h"
#include "Meshes/ParallelMesh.h"
#include "Results/Result.h"

#include <filesystem>

/**
 * @class AsterToMedWriter
 * @brief writer Med file interface
 * @author Nicolas Sellenet
 */
class AsterToMedWriter {
  private:
#ifdef ASTER_HAVE_MED
    enum entityType { Cells, Nodes };
    void _createGroups( const BaseMesh &, MedMeshPtr, std::vector< med_int > &, const VectorLong &,
                        entityType, bool );
#ifdef ASTER_HAVE_MPI
    void _buildFilterInformations(
        const VectorLong &nodeVector, const VectorLong &cellVector, JeveuxVectorLong cellTypeVec,
        std::array< med_int, 3 > &cumNodeNb,
        std::map< med_geometry_type, std::array< med_int, 3 > > &mapMedTypeCellsByProc,
        VectorOfVectorsLong & );
#endif
    void _createMedGlobalNumbering( VectorLong &globNum, const VectorLong &innerNodes, int startNum,
                                    int localNodeNumber );

    bool _printMeshFromList( const BaseMesh &, const std::filesystem::path &filename,
                             const VectorLong &nodeList, const VectorLong &cellList, bool local,
                             const std::string &meshName, const std::string &mode = "a" );

    void _sortCellsByType( const VectorLong &cellVector, JeveuxVectorLong cellTypeVec,
                           VectorOfVectorsLong &cellIdByType ) const;
#endif

  public:
    /**
     * @typedef AsterToMedWriterPtr
     * @brief Pointeur intelligent vers un AsterToMedWriter
     */
    typedef std::shared_ptr< AsterToMedWriter > AsterToMedWriterPtr;

    /** @brief Constructor */
    AsterToMedWriter() {};

    ~AsterToMedWriter() {};

#ifdef ASTER_HAVE_MED
    /** @brief print Mesh object */
    bool printMesh( const Mesh &, const std::filesystem::path &filename, bool local = true,
                    const std::string &meshName = "", const std::string &mode = "a" );
#ifdef ASTER_HAVE_MPI
    /** @brief print ParallelMesh object */
    bool printMesh( const ParallelMesh &, const std::filesystem::path &filename, bool local = true,
                    const std::string &meshName = "", const std::string &mode = "a" );

    /** @brief print ConnectionMesh object */
    bool printMesh( const ConnectionMesh &, const std::filesystem::path &filename,
                    bool local = true, const std::string &meshName = "",
                    const std::string &mode = "a" );
#endif
    /** @brief print MeshPtr object */
    bool printMesh( const MeshPtr &mesh, const std::filesystem::path &filename, bool local = true,
                    const std::string &meshName = "", const std::string &mode = "a" );
#ifdef ASTER_HAVE_MPI
    /** @brief print ParallelMeshPtr object */
    bool printMesh( const ParallelMeshPtr &mesh, const std::filesystem::path &filename,
                    bool local = true, const std::string &meshName = "",
                    const std::string &mode = "a" );

    /** @brief print ConnectionMeshPtr object */
    bool printMesh( const ConnectionMeshPtr &mesh, const std::filesystem::path &filename,
                    bool local = true, const std::string &meshName = "",
                    const std::string &mode = "a" );
#endif
    /** @brief print ResultPtr object */
    bool printResult( const ResultPtr &, const std::filesystem::path &filename, bool local = true,
                      const std::string &mode = "a" );
#endif
};

/**
 * @typedef AsterToMedWriterPtr
 * @brief Pointeur intelligent vers un AsterToMedWriter
 */
typedef std::shared_ptr< AsterToMedWriter > AsterToMedWriterPtr;

#endif /* ASTERTOMEDWRITER_H_ */
