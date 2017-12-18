#ifndef gl_lin_INC
#define gl_lin_INC

#include "geo_pnt.h" 

/*---------------------------------------------------------------*/
class CGLLine {
private:
	bool marked; // For mesh
	friend class Surface; //WW

public:
	CGLLine(void);
	~CGLLine(void);
	CGLPoint *m_point1;
	CGLPoint *m_point2;
	std::string name;
	double value;
	long mesh_index, gli_line_id; // mesh_index: unique index for mesh
	long point1, point2;
	int orientation;
	double epsilon;
	long *msh_nodes;
	long no_msh_nodes;
	int mat_group;
	int display_mode;
	//MSH
	std::vector<double*> nodes_coor_vector;
	//Method
	CGLLine* GEOGetLine(long);
	CGLLine *CheckLineOutPut();
	CGLLine *Exists();
	// void SetRFIPointsClose();//MSH + geomathlib
	//void CreateMSHLines(void);//MSH
};

extern std::vector<CGLLine*> GEOLIB_GetGLILines_Vector(void);
extern std::vector<CGLLine*> gli_lines_vector;
extern std::vector<CGLLine*> gli_file_lines_vector;
extern void Clear_LineVector();

#endif
