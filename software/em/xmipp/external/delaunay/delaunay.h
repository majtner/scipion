#ifndef DELAUNAY_H
#define DELAUNAY_H

#include "defines.h"
#include "dcel.h"
#include "graph.h"
#include "point.h"



/**************************************************************************
* Data types definition
**************************************************************************/
struct Delaunay_T
{
	struct Graph_T 	graph;
	struct DCEL_T	*dcel;
};


/**************************************************************************
* Public functions declaration
**************************************************************************/
int		init_Delaunay( struct Delaunay_T *delaunay, int nPoints);
void	delete_Delaunay( struct Delaunay_T *delaunay);
void	insert_Point( struct Delaunay_T *delaunay, double x, double y);
void	create_Delaunay_Triangulation( struct Delaunay_T *delaunay);
int 	get_Face_Points( struct Delaunay_T *delaunay, int face_ID,
												struct Point_T *p,
												struct Point_T *q,
												struct Point_T *r);

int 	initialize_Delaunay(struct Delaunay_T *delaunay, struct DCEL_T *dcel);
void 	finalize_Delaunay(struct Delaunay_T *delaunay);

void	incremental_Delaunay(struct Delaunay_T *delaunay);
void 	build_Delaunay_From_Triangulation(struct DCEL_T *dcel);
void    purge_Delaunay( struct Delaunay_T *delaunay);

int		is_Interior_To_Node(struct DCEL_T *dcel, struct Graph_T *graph,
														struct Point_T *p,
														int node_Index);
int  	select_Closest( struct Delaunay_T *delaunay, int index);
void    select_Two_Closest( struct Delaunay_T *delaunay, int *first, int *second);

enum Turn_T return_Turn( struct DCEL_T *dcel, struct Point_T *p,
												int source_ID,
												int dest_ID);

#endif

