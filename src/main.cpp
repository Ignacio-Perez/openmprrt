/*  OPENMP-RRT
    Copyright (C) 2019  Ignacio Perez-Hurtado
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.*/
    

#include <iostream>
#include <ctime>
#include <string>
#include <pgm.h>
#include <rrt.hpp>
#include <cstring>
#include <omp.h>
#include <cmath>



#define MAP1 "../maps/map.pgm"
#define MAP2 "../maps/office.pgm"
#define MAP3 "../maps/labyrinth.pgm"
#define MAP4 "../maps/ccia_h.pgm"

#define RESOLUTION 0.05
#define RADIUS 0.2
#define NODES 4096

#define INF 99999
#define DELTA 0.15

struct Point
{
	double x;
	double y;
};

struct XYD
{
	double x;
	double y;
	double d;
};

struct Map
{
	PGM* pgm;
	double resolution;
	Point init; //initial robot position
};


struct Tree
{
	Point *nodes;
	Point *parents;
	unsigned size;
};


XYD xyd_min2(XYD a, XYD b)
{
	return a.d < b.d ? a : b;
}

#pragma omp declare reduction(xyd_min : XYD : omp_out=xyd_min2(omp_out,omp_in))\
		initializer(omp_priv={0,0,INF})



float rnd()
{
	return (float)rand()/(float)(RAND_MAX);
}



void rrt_star_algorithm(const Map& map, unsigned size, Tree& t)
{
	t.nodes = new Point[size];
	t.parents = new Point[size];
	t.size = size;
	double* distances = new double[size];
	double* cost = new double[size];
	
	t.nodes[0].x = map.init.x;
	t.nodes[0].y = map.init.y;
	cost[0] = 0;
	
	double p = map.pgm->width * map.resolution;
	double q = map.pgm->height * map.resolution;
	
	unsigned index=1;
	
	
	while(index<size) {
		Point random = {p*rnd(), q*rnd()};
		
		#pragma omp parallel for
		for (unsigned i=0;i<index;i++) {
			distances[i] = (t.nodes[i].x - random.x) * (t.nodes[i].x - random.x)  + (t.nodes[i].y - random.y) * (t.nodes[i].y - random.y);
		}
		
		// compute minimun distance and nearest point
		XYD value = {0,0,INF};
		#pragma omp parallel for reduction(xyd_min:value)
		for (unsigned i=0;i<index;i++) {
			XYD new_value = {t.nodes[i].x,t.nodes[i].y,distances[i]};
			value = xyd_min2(value,new_value);
		}
		double d = sqrt(value.d);
		Point nearest = {value.x,value.y};
		
		t.nodes[index] = {nearest.x + DELTA * (random.x - nearest.x)/d, nearest.y + DELTA * (random.y - nearest.y)/d};
			
		int x0 = t.nodes[index].x / map.resolution;
		int y0 = t.nodes[index].y / map.resolution;
		int x1 = nearest.x / map.resolution;
		int y1 = nearest.y / map.resolution;
		
		if (detect_obstacle(map.pgm, x0,  y0, x1,y1, 250)) {
			continue;
		}
		
		value = {0,0,INF};
		#pragma omp parallel for reduction(xyd_min:value)
		for (unsigned i =0; i<index;i++) {
			int x1 = t.nodes[i].x / map.resolution;
			int y1 = t.nodes[i].y / map.resolution;
			XYD new_value = {t.nodes[i].x,t.nodes[i].y,cost[i] + (t.nodes[index].x - t.nodes[i].x) * (t.nodes[index].x - t.nodes[i].x) + 
				(t.nodes[index].y - t.nodes[i].y) * (t.nodes[index].y - t.nodes[i].y) + detect_obstacle(map.pgm, x0,  y0, x1,y1, 250)*INF };
			value = xyd_min2(value,new_value);
		}
		
		t.parents[index] = {value.x,value.y};
		cost[index] = value.d;
		
		#pragma omp parallel for 
		for (unsigned i=0;i<index;i++) {
			int x1 = t.nodes[i].x / map.resolution;
			int y1 = t.nodes[i].y / map.resolution;
			double c = cost[index] + (t.nodes[index].x - t.nodes[i].x )*(t.nodes[index].x - t.nodes[i].x) + (t.nodes[index].y - t.nodes[i].y )*(t.nodes[index].y - t.nodes[i].y);
			if (c < cost[i] && !detect_obstacle(map.pgm,x0,y0,x1,y1,250)) {
				cost[i] = c;
				t.parents[i] = t.nodes[index];
			}
		}
		
		index++;
	}
	
	
}


void rrt_algorithm(const Map& map, unsigned size, Tree& t)
{
	t.nodes = new Point[size];
	t.parents = new Point[size];
	t.size = size;
	
	double* distances = new double[size];
	
	t.nodes[0].x = map.init.x;
	t.nodes[0].y = map.init.y;
	
	double p = map.pgm->width * map.resolution;
	double q = map.pgm->height * map.resolution;
	
	unsigned index=1;
	
	
	while(index<size) {
		Point random = {p*rnd(), q*rnd()};
		
		#pragma omp parallel for
		for (unsigned i=0;i<index;i++) {
			distances[i] = (t.nodes[i].x - random.x) * (t.nodes[i].x - random.x)  + (t.nodes[i].y - random.y) * (t.nodes[i].y - random.y);
		}
		
		// compute minimun distance and nearest point
		XYD value = {0,0,INF};
		#pragma omp parallel for reduction(xyd_min:value)
		for (unsigned i=0;i<index;i++) {
			XYD new_value = {t.nodes[i].x,t.nodes[i].y,distances[i]};
			value = xyd_min2(value,new_value);
		}
		double d = sqrt(value.d);
		Point nearest = {value.x,value.y};
		
		t.nodes[index] = {nearest.x + DELTA * (random.x - nearest.x)/d, nearest.y + DELTA * (random.y - nearest.y)/d};
		t.parents[index] = nearest;
	
		int x0 = t.nodes[index].x / map.resolution;
		int y0 = t.nodes[index].y / map.resolution;
		int x1 = t.parents[index].x / map.resolution;
		int y1 = t.parents[index].y / map.resolution;
		
		if (!detect_obstacle(map.pgm, x0,  y0, x1,y1, 250)) {
			index++;
		}
		
	}
	
}



int main(int argc, char* argv[])
{
		
	int selection=1;
	int algorithm = RRT_ALGORITHM;
	int seed = time(NULL);
	
	if (argc==1) {
		printf("\nFormat: %s [m] [a] [s]\n",argv[0]);
		printf("Where: \n");
		printf("\n[m] is the map index\n");
		printf("\t1 = %s\n",MAP1);
		printf("\t2 = %s\n",MAP2);
		printf("\t3 = %s\n",MAP3);
		printf("\t4 = %s\n",MAP4);
		printf("\n[a] is the algorithm index\n");
		printf("\t1 = RRT ALGORITHM\n");
		printf("\t2 = RRT-STAR ALGORITHM\n");
		printf("\n[s] is the random seed\n");
		return 0;
	}
		
	if (argc>1) {
		selection = atoi(argv[1]);	
	}
	
	if (argc>2) {
		algorithm = atoi(argv[2])-1;
	}

	if (argc>3) {
		seed = atoi(argv[3]);
	}

	if (algorithm!= RRT_ALGORITHM && algorithm != RRT_STAR_ALGORITHM) {
		printf("Invalid algorithm\n");
		return 0;
	}
	
	srand(seed);
	
	Map map;
	PGM* pgm;
	map.resolution = RESOLUTION;
	
	switch(selection) {
		case 1:
			map.pgm =load_pgm(MAP1);
			pgm = load_pgm(MAP1);
			map.init.x = 8;
			map.init.y = 10;
		break;
		case 2:
			map.pgm = load_pgm(MAP2);
			pgm = load_pgm(MAP2);
			map.init.x = 32;
			map.init.y = 9.3;
		break;
		case 3:
			map.pgm = load_pgm(MAP3);
			pgm = load_pgm(MAP3);
			map.init.x = 21.5;
			map.init.y = 21.5;
		break;
		case 4:
			map.pgm = load_pgm(MAP4);
			pgm = load_pgm(MAP4);
			map.init.x = 5.25;
			map.init.y = 30.54;
		break;
		default:
			printf("Invalid map\n");
			return 0;
	}
	
	inflate_obstacles(map.pgm, RADIUS/RESOLUTION);
	//remove_inner_obstacles(map.pgm);
	
	Tree tree;
	
	if (algorithm == RRT_ALGORITHM) {
		rrt_algorithm(map,NODES,tree);
	} else {
		rrt_star_algorithm(map,NODES,tree);
	}
	
	
	for (int i=1;i<NODES;i++) {
		int x0 = (int)round(tree.nodes[i].x / RESOLUTION);
		int y0 = (int)round(tree.nodes[i].y / RESOLUTION);
		int x1 = (int)round(tree.parents[i].x / RESOLUTION);
		int y1 = (int)round(tree.parents[i].y  / RESOLUTION);
		draw_line(pgm,x0,y0,x1,y1,0);
	}
	strcpy(pgm->file,"output.pgm");
	save_pgm(pgm);
	
	
	return 0;
	
}    
