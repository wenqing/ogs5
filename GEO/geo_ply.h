/**
 * \copyright
 * Copyright (c) 2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/**************************************************************************
   GeoLib - Object: Polyline
   Task:
   Programing:
   07/2003 OK/CC/TK/WW GEOLib1
   07/2005 CC/OK  GEOLib2 Design
**************************************************************************/
#ifndef gl_ply_INC
#define gl_ply_INC
#include "geo_lin.h"
#include "geo_pnt.h"
#include <fstream>
#include <string>
#define PLY_FILE_EXTENSION ".ply"
// enum PLY_TYPES {GEO,MSH,IC,BC,ST};

class CGLPolyline
{
private:
    friend class Surface;
    std::vector<double> sbuffer;
    std::vector<int> ibuffer;

    // TF 10/2010 - made attributes private
    // properties
    size_t id;  // CC
    std::string name;
    int type;
    int data_type;
    int mat_group;
    //	std::string ply_data;//CC9999
    //	double minDis;
    //	double mesh_density;//CC9999 ply density
    //	bool closed;
    //	double min_plg_Dis;
    // components
    bool computeline;
    //	std::string ply_file_name;
    std::vector<int> OrderedPoint;

    /// is the epsilon value set external or is it on default value
    bool _set_eps;

protected:
    std::vector<CGLLine*> line_vector;

public:
    // constructor
    CGLPolyline(void);
    CGLPolyline(const std::string&);  // OK
    // destructor
    ~CGLPolyline(void);

    /// is the epsilon value set external or is it on default value
    /// used in mesh node search algorithms
    bool isSetEps() const { return _set_eps; }
    const std::string& getName() const;
    void setName(const std::string& nname);
    size_t getID() const;
    void setID(size_t nid);
    // properties
    size_t getType() const { return type; }
    /**
     * if the data is read from the gli-file, data_type is 0
     * if the data is read from an external file or is generated by OpenGeoSys,
     * data_type is 1
     * @return the data type
     */
    int getDataType() const { return data_type; }
    void setDataType(int ndata_type) { data_type = ndata_type; }
    int getMatGroup() const { return mat_group; }
    double epsilon;

    // components
    const std::vector<CGLLine*>& getLineVector() const { return line_vector; }
    std::vector<CGLLine*>& getLineVector() { return line_vector; }
    const std::vector<int>& getOrderedPoints() const { return OrderedPoint; }
    std::vector<CGLPoint*> point_vector;

    std::vector<double>& getSBuffer() { return sbuffer; }
    std::vector<int>& getIBuffer() { return ibuffer; }
    // I/O
    std::ios::pos_type Read(std::ifstream&);  // TF , const std::string &);//CC
    void Write(char* file_name);
    // display
    //    void AssignColor();//CC
    // method
    void ComputeLines();
    bool PointExists(CGLPoint* point, CGLPoint* point1);
    void AddPoint(CGLPoint* m_point);
    CGLPoint* CenterPoint(void);
    // point vector
    void WritePointVector(const std::string&);  // CC
    void ReadPointVector(const std::string&);   // CC
    void SortPointVectorByDistance();
    // write tecplot file
    void WriteTecplot(const std::string&);  // CC
    // Meshing
    std::vector<long> msh_nodes_vector;
    std::vector<double*> msh_coor_vector;

    void GetPointOrderByDistance();
    void SetPointOrderByDistance(CGLPoint*);  // OK
//	void CalcMinimumPointDistance(); //OK
#ifdef RFW_FRACTURE
    double CalcPolylineLength();  // RFW
#endif
};

extern std::vector<CGLPolyline*> polyline_vector;          // CC
extern std::vector<CGLPolyline*> GetPolylineVector(void);  // CC
// Access
extern CGLPolyline* GEOGetPLYByName(const std::string&);
extern CGLPolyline* GEOGetPLYById(long);  // CC
// methods
extern void GEOPolylineGLI2GEO(FILE* geo_file);
// extern void GEOUnselectPLY(); //OK
// Remove
extern void GEORemoveAllPolylines();     // CC
extern void GEORemovePolyline(long);     // CC 03/06
extern void GEORemovePLY(CGLPolyline*);  // OK
// I/O
extern void GEOReadPolylines(const std::string& file_name_path_base);
extern void GEOWritePolylines(char* file_name);  // CC
// RF
extern void InterpolationAlongPolyline(CGLPolyline* plyL,
                                       std::vector<double>& bcNodalValue);

class CColumn : public CGLPolyline  // OK
{
private:
    CGLLine* m_lin;

public:
    ~CColumn();
    void deleteLines();
    double geo_area;
    double center_point[3];
};

class CSoilProfile : public CGLPolyline  // YD
{
public:
    CSoilProfile();
    ~CSoilProfile();
    long profile_type;
    std::string soil_name;
    std::vector<long> soil_type;
    std::vector<double> soil_layer_thickness;
};

extern void COLDeleteLines();
extern void COLDelete();
extern CColumn* COLGet(int);
extern CColumn* COLGet(const std::string&);

extern std::vector<CColumn*> column_vector;
extern std::vector<CSoilProfile*> profile_vector;  // YD
#endif
