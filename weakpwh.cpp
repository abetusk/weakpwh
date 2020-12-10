/*
    Copyright (C) 2013 Abram Connelly

    This file is part of weakpwh v0.1.

    weakpwh is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    weakpwh is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with weakpwh.  If not, see <http://www.gnu.org/licenses/>.

*/


// Bare bones command line utility to convert a polygon with holes
// into a single path (a weak polygon with holes).
//
// Assumes integer input.
//
// to compile:
//   g++ -g sweep/advancing_front.cc sweep/cdt.cc sweep/sweep.cc sweep/sweep_context.cc common/shapes.cc weakpwh.cpp -o weakpwh
//
// sample usage:
//   cat testdata/test4.gp | ./weakpwh
//

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <unistd.h>
#include <string.h>
#include <getopt.h>

#include <errno.h>

#include <algorithm>
#include <map>

#include "poly2tri.h"
#include "clipper.hpp"

#define WEAKPWH_VER "0.2.0"

bool gDebugFlag = false;

long long gMulFactor = 1000;

double gPremul = 1000000000.0;

typedef signed long long cInt;
typedef unsigned long long cUInt;


int gFloatRead = 0;
int gIgnoreColinear = 1;


cInt dtocint( double d ) {
  if (d < 0.0) return (unsigned long long)(d-0.5);
  return (signed long long)(d+0.5);
}

//---

class iPnt {
  public:
  long long int X, Y;
  iPnt() { X=0; Y=0; }
  iPnt(long long int a, long long int b) { X = a; Y = b; }
} ;

struct iPnt_cmp {
  bool operator() (const iPnt &lhs, const iPnt &rhs) const {
    if ( lhs.X == rhs.X ) return lhs.Y < rhs.Y ;
    return lhs.X < rhs.X;
  }
};

class PointInfo {
  public:
  int pointInd, polygonInd;
  int dir;
  bool visited, is_outer_boundary, is_clockwise;
} ;

typedef std::map< iPnt , PointInfo, iPnt_cmp > PointInfoMap;
typedef std::map< iPnt, std::vector< iPnt >, iPnt_cmp > PointAdjacency;

//---

class iEdge {
  public:
    iPnt v[2];
};

struct iEdge_cmp {
  bool operator() (const iEdge &lhs, const iEdge &rhs) const {
    if ( ( lhs.v[0].X == rhs.v[0].X ) && 
         ( lhs.v[0].Y == rhs.v[0].Y ) ) {
      if ( lhs.v[1].X == rhs.v[1].X ) return lhs.v[1].Y < rhs.v[1].Y;
      return lhs.v[1].X < rhs.v[1].X;
    }
    
    if ( lhs.v[0].X == rhs.v[0].X ) return lhs.v[0].Y < rhs.v[0].Y ;
    return lhs.v[0].X < rhs.v[0].X;
  }
};

typedef std::map< iEdge, int, iEdge_cmp > EdgeMap;

//---

class SimplePolyTreeNode {
  public:
    int boundaryIndex;
    std::vector<int> holeIndex;
} ;

class SimplePolyTree {
  public:
    std::vector<SimplePolyTreeNode> node;
    std::vector<SimplePolyTree *>children;
} ;

//---

bool consistency_check( std::vector< std::vector<iPnt> > &pwh, PointAdjacency &edge, PointInfoMap &pointInfo) {
  int i, j, k;
  iPnt xy;
  PointInfo pi;

  for (i=0; i<pwh.size(); i++) {
    for (j=0; j<pwh[i].size(); j++) {
      xy.X = pwh[i][j].X;
      xy.Y = pwh[i][j].Y;

      if ( pointInfo.find(xy) == pointInfo.end() ) {
        fprintf(stderr, "ERROR: point (%lli,%lli) in pwh[%i][%i] but not found in pointInfo\n",
            xy.X, xy.Y, i, j);
        exit(2);
      }

      pi = pointInfo[xy];

      if (pi.polygonInd != i ) {
        fprintf(stderr, "ERROR: point (%lli,%lli) in pwh[%i][%i] does not match polygonInd (%i)\n",
            xy.X, xy.Y, i, j, pi.polygonInd);
        exit(2);
      }

      if (pi.pointInd != j ) {
        fprintf(stderr, "ERROR: point (%lli,%lli) in pwh[%i][%i] does not match pointInd (%i)\n",
            xy.X, xy.Y, i, j, pi.pointInd);
        exit(2);
      }

      if ((i==0) && (!pi.is_outer_boundary)) {
        fprintf(stderr, "ERROR: point (%lli,%lli) in pwh[%i][%i] should be boundary but isn't\n",
            xy.X, xy.Y, i, j );
        exit(2);
      }

      if ((i!=0) && (pi.is_outer_boundary)) {
        fprintf(stderr, "ERROR: point (%lli,%lli) in pwh[%i][%i] should be not boundary but is\n",
            xy.X, xy.Y, i, j );
        exit(2);
      }

    }
  }

  for (i=0; i<pwh.size(); i++) {
    for (j=0; j<pwh[i].size(); j++) {
      xy.X = pwh[i][j].X;
      xy.Y = pwh[i][j].Y;
      pi = pointInfo[xy];
    }
  }

  PointAdjacency::iterator e_it;
  for (e_it = edge.begin(); e_it != edge.end(); ++e_it ) {
    xy = e_it->first;
    std::vector< iPnt > v = e_it->second;

    if ( pointInfo.find(xy) == pointInfo.end() ) {
      fprintf(stderr, "ERROR: adjacency lookup fail: (%lli,%lli) not found in PointInfoMap\n",
          xy.X, xy.Y);
      exit(2);
    }

    for (i=0; i<v.size(); i++) {

      if ( pointInfo.find(v[i]) == pointInfo.end() ) {
        fprintf(stderr, "ERROR: adjacency lookup fail: edge (%lli,%lli)->(%lli,%lli) [%i] not found in PointInfoMap\n",
            xy.X, xy.Y, v[i].X, v[i].Y, i);
        exit(2);
      }

    }
  }

  return true;
}

//---

void addEdgesFromTri( PointAdjacency &edge, p2t::Triangle *tri, EdgeMap &edgeMap ){
  iEdge uv;
  int a,b,i, perm[][2] = { {0,1}, {1,0}, {0,2}, {2,0}, {1,2}, {2,1} };

  for (i=0; i<6; i++) {
    a = perm[i][0];
    b = perm[i][1];
    uv.v[0].X = dtocint( tri->GetPoint(a)->x );
    uv.v[0].Y = dtocint( tri->GetPoint(a)->y );
    uv.v[1].X = dtocint( tri->GetPoint(b)->x );
    uv.v[1].Y = dtocint( tri->GetPoint(b)->y );

    if ( edgeMap.find( uv ) != edgeMap.end() ) { continue; }

    edge[ uv.v[0] ].push_back( uv.v[1] );
    edgeMap[ uv ] = 1;
  }

}

void debugPrintEdge( PointAdjacency &edge, PointInfoMap &bmap ) {
  PointAdjacency::iterator it;
  iPnt u, p;
  std::vector< iPnt > *v;

  printf("#debug edge print\n");

  for (it = edge.begin(); it != edge.end(); ++it ) {
    printf("#\n");

    u = it->first;
    v = &(it->second);

    for (int i=0; i<v->size(); i++) {
      p = v->at(i);
      printf("#[%i] u(%lli, %lli)(%i) -> v(%lli %lli)(%i)\n", i, 
          u.X, u.Y, (int)(bmap[ u ].is_outer_boundary ),
          p.X, p.Y, (int)(bmap[ p ].is_outer_boundary ) );
    }

  }
}

typedef std::vector< iPnt > Path;
typedef std::vector< std::vector<iPnt> > Paths;

void ConvertPWHToWeakPath_r( Path &outPath,
                             Paths &paths,
                             PointAdjacency &edge,
                             PointInfoMap &pointInfoMap,
                             std::map< int, bool> &visitedPoly,
                             int polygonInd,
                             int pointInd ) {

  int i, k, n_nei;
  int end_ind, m;
  PointInfo *pi_u, *pi_v;


  m = paths[polygonInd].size();
  end_ind = (pointInd + m - 1) % m;

  visitedPoly[ polygonInd ] = true;

  for ( i=pointInd; i != end_ind; i = ((i+1)%m) ) {

    outPath.push_back( iPnt( paths[polygonInd][i].X, paths[polygonInd][i].Y ) );

    pi_u = &(pointInfoMap[ paths[polygonInd][i] ]);
    pi_u->visited = true;

    n_nei = edge[ paths[polygonInd][i] ].size();

    for (k=0; k<n_nei; k++) {

      pi_v = &(pointInfoMap[ edge[ paths[polygonInd][i] ][k] ]);
      if ( visitedPoly[pi_v->polygonInd] || (pi_v->polygonInd == pi_u->polygonInd) ) { 
        continue; 
      }

      ConvertPWHToWeakPath_r( outPath, paths, edge, pointInfoMap, visitedPoly, pi_v->polygonInd, pi_v->pointInd );

      outPath.push_back( iPnt( paths[polygonInd][i].X, paths[polygonInd][i].Y ) );
    }

  }

  outPath.push_back( iPnt( paths[polygonInd][end_ind].X, paths[polygonInd][end_ind].Y ) );
  outPath.push_back( iPnt( paths[polygonInd][pointInd].X, paths[polygonInd][pointInd].Y ) );

}


// Go through each vertex of each polygon.  If the polygon has already
// been visited, skip over it.  If the 
//
int ConvertPolygonWithHolesToWeakPath( Path &outPath, 
                                       Paths &paths, 
                                       PointAdjacency &edge, 
                                       PointInfoMap &pointInfoMap ) {

  int i, j, k, n, m, n_nei;
  iPnt u, v, w;
  PointInfo *pi_u, *pi_v;

  std::map< int, bool > visitedPoly;

  outPath.clear();

  n = paths.size();
  for (i=0; i<n; i++)  {
    m = paths[i].size();
    if (m < 3) { return -1; }
    visitedPoly[ i ] = false;
  }


  for (i=0; i<n; i++) {

    if (visitedPoly[i]) { continue; }

    visitedPoly[ i ] = true;

    m = paths[i].size();
    for (j=0; j<m; j++) {

      outPath.push_back( iPnt( paths[i][j].X, paths[i][j].Y ) );
      
      pi_u = &(pointInfoMap[ paths[i][j] ]);
      pi_u->visited = true;

      n_nei = edge[ paths[i][j] ].size();

      for (k=0; k<n_nei; k++) {
        pi_v = &(pointInfoMap[ edge[ paths[i][j] ][k] ]);

        // If we've already visited the polygon contour (or it's this polygon contour), skip
        //
        if ( visitedPoly[pi_v->polygonInd] || (pi_v->polygonInd == pi_u->polygonInd) ) { 
          continue; 
        }

        // Otherwise, start the recursion
        //
        ConvertPWHToWeakPath_r( outPath, paths, edge, pointInfoMap, visitedPoly, pi_v->polygonInd, pi_v->pointInd );

        // Push our return point back
        //
        outPath.push_back( iPnt( paths[i][j].X, paths[i][j].Y ) );

      }

    }

    // Make sure to tie off the point from which we started
    //
    outPath.push_back( iPnt( paths[i][0].X, paths[i][0].Y ) );

  }

  return 0;

}

struct option gLongOption[] =
{
  {"input", required_argument,0,'i'},
  {"output", required_argument,0,'o'},
  {"premul", required_argument,0,'F'},
  {"version", no_argument, 0, 'V'},
  {"verbose", no_argument, 0, 'v'},
  {"help", no_argument, 0, 'h'},
  {0,0,0,0}
};

char gOptionDescription[][1024] =
{
  "input",
  "output",
  "premul",
  "version",
  "verbose",
  "help",
  "n/a"
};

void show_version(FILE *fp) {
  fprintf(fp, "version %s\n", WEAKPWH_VER);
}

void show_help(FILE *fp) {
  int i, j, k;
  int len;


  fprintf(fp, "\n");
  fprintf(fp, "Convert a polygon with holes to a weakly simple polygon\n");
  show_version(fp);
  fprintf(fp, "\n");

  fprintf(fp, "usage:\n");
  fprintf(fp, "\n");
  fprintf(fp, "  weakpwh [-i inputfile] [-h] [-v] [-V] [-o outfile] [inputfile]\n");
  fprintf(fp, "\n");

  for (i=0; gLongOption[i].name; i++) {
    len = strlen(gLongOption[i].name);

    fprintf(fp, "  -%c, --%s", gLongOption[i].val, gLongOption[i].name);
    if (gLongOption[i].has_arg) {
      fprintf(fp, " %s", gLongOption[i].name);
      len = 2*len + 3;
    }
    else {
      len = len + 2;
    }
    for (j=0; j<(32-len); j++) {
      fprintf(fp, " ");
    }

    fprintf(fp, "%s\n", gOptionDescription[i]);
  }
  fprintf(fp, "\n");

}

bool gVerboseFlag = false;
char *gInputFilename = NULL;
char *gOutputFilename = NULL;

FILE *gInputFp = NULL;
FILE *gOutputFp = NULL;

void process_command_line_options( int argc, char **argv ) {
  extern char *optarg;
  extern int optind;
  int option_index;

  char ch;

  while ((ch = getopt_long(argc, argv, "hi:o:vVF:", gLongOption, &option_index)) > 0) switch(ch) {
    case 'h':
      show_help(stdout);
      exit(0);
      break;
    case 'F':
      gPremul = strtod(optarg, NULL);
      gFloatRead = 1;
      break;
    case 'V':
      show_version(stdout);
      exit(0);
      break;
    case 'v':
      gVerboseFlag = true;
      break;
    case 'i':
      gInputFilename = strdup(optarg);
      break;
    case 'o':
      gOutputFilename = strdup(optarg);
      break;
    case 0:
    default:
      fprintf(stderr, "bad option\n");
      exit(2);
      break;
  }

  if (!gInputFilename) {
    if (optind < argc) {
      gInputFilename = strdup(argv[optind]);
    } else {
      fprintf(stderr, "provide input file\n");
      show_help(stderr);
      exit(2);
    }
  }

  if (gInputFilename[0] == '-') {
    gInputFp = stdin;
  } else {
    gInputFp = fopen(gInputFilename, "r");
    if (!gInputFp) {
      perror(gInputFilename);
      exit(errno);
    }
  }
  free(gInputFilename);


  if ((!gOutputFilename) || (gOutputFilename[0] == '-')) {
    gOutputFp = stdout;
  } else {
    gOutputFp = fopen(gOutputFilename, "w");
    if (!gOutputFp) {
      perror(gOutputFilename);
      exit(errno);
    }
  }

}

void print_poly_nodes( ClipperLib::PolyNode *node, int depth ) {
  int i, j, k;

  printf("# [%i] isHole %i\n", depth, node->IsHole() ? 1 : 0 );
  for (i=0; i<node->Contour.size(); i++) {
    printf("%lli %lli\n", node->Contour[i].X, node->Contour[i].Y );
  }
  printf("\n");

  for (i=0; i<node->Childs.size(); i++) {
    print_poly_nodes( node->Childs[i], depth+1 );
  }

}

int constructPolygonsWithHoles( ClipperLib::PolyNode *node,
                                std::vector< std::vector< std::vector< iPnt > > > &pwhs,
                                int cur_pwh_ind,
                                int next_ind ) {
  int i, j, k, n, ind;
  std::vector< iPnt > path;
  std::vector< std::vector< iPnt > > pwh;

  for (i=0; i<node->Contour.size(); i++) {
    path.push_back( iPnt( node->Contour[i].X, node->Contour[i].Y ) );
  }

  ind = cur_pwh_ind;

  if (!node->IsHole()) {
    ind = next_ind;
    next_ind++;
    pwhs.push_back( pwh );
  }

  pwhs[ ind ].push_back( path );

  for (i=0; i<node->Childs.size(); i++) {
    next_ind = constructPolygonsWithHoles( node->Childs[i], pwhs, ind, next_ind );
  }

  return next_ind;

}


int clip_polygon_is_A_simply_inside_B( ClipperLib::Path &A, ClipperLib::Path &B ) {
  ClipperLib::Clipper clipAmB, clipBmA;
  ClipperLib::Paths solnAmB, solnBmA;

  ClipperLib::Path x, y;

  x = A;
  y = B;

  if ( ClipperLib::Area( x ) < 0 ) { ClipperLib::ReversePath( x ); }
  if ( ClipperLib::Area( y ) < 0 ) { ClipperLib::ReversePath( y ); }

  clipAmB.AddPath( x, ClipperLib::ptSubject, true );
  clipAmB.AddPath( y, ClipperLib::ptClip, true );

  clipAmB.Execute( ClipperLib::ctDifference,
                   solnAmB,
                   ClipperLib::pftNonZero,
                   ClipperLib::pftNonZero  );

  clipBmA.AddPath( y, ClipperLib::ptSubject, true );
  clipBmA.AddPath( x, ClipperLib::ptClip, true );

  clipBmA.Execute( ClipperLib::ctDifference,
                   solnBmA,
                   ClipperLib::pftNonZero,
                   ClipperLib::pftNonZero  );

  if ( (solnAmB.size() == 0) && (solnBmA.size() > 0) ) {
    return 1;
  }
  return 0;

}


void emitWeakPWH( std::vector< std::vector< iPnt > > pwh ) {
  int i, j, k;
  int polyname = 0;
  int line_no=0;
  iPnt xy;
  PointInfoMap pointInfoMap;
  EdgeMap edgeMap;
  std::map< int, bool > path_visited_map;

  PointAdjacency edge;
  PointInfo pi;
  iEdge uv;
  std::vector< iPnt > ipoly;

  std::vector<p2t::Point *> outer_boundary;
  std::vector< std::vector<p2t::Point *> > holes;
  std::vector<p2t::Point *> poly;

  std::vector<p2t::Point *> garbage_poly;

  double x_d, y_d;


  for (i=0; i<pwh.size(); i++) {
    for (j=0; j<pwh[i].size(); j++) {

      xy.X = pwh[i][j].X;
      xy.Y = pwh[i][j].Y;

      pi.pointInd = j;
      pi.polygonInd = i;
      pi.is_outer_boundary = ((i==0) ? true : false );
      pi.is_clockwise = ((i==0) ? true : false );
      pi.dir = ((i==0) ? 1 : 1 );

      if ( pointInfoMap.find( xy ) != pointInfoMap.end() ) {
        fprintf(stderr, "DUPLICATE POINT FOUND! line %i, exiting\n", line_no); exit(2);
      }

      pointInfoMap[xy] = pi;

    }
  }

  for (i=0; i<pwh[0].size(); i++) {
    outer_boundary.push_back( new p2t::Point( pwh[0][i].X, pwh[0][i].Y ) );
    garbage_poly.push_back( outer_boundary[ outer_boundary.size()-1 ] );
  }

  for (i=1; i<pwh.size(); i++) {
    poly.clear();
    for (j=0; j<pwh[i].size(); j++) {
      poly.push_back( new p2t::Point( pwh[i][j].X, pwh[i][j].Y ) );
      garbage_poly.push_back( poly[ poly.size()-1 ] );
    }
    holes.push_back(poly);
  }

  if (outer_boundary.size() < 3) {
    fprintf(stderr, "invalid number of points for outer boundary (%lu), exiting\n", outer_boundary.size()); exit(2); 
  }

  for (i=0; i<holes.size(); i++) {
    if (holes[i].size() < 3) {
      fprintf(stderr ,"hole %i has less than 3 points (%lu), exiting\n", i, holes[i].size() );
      exit(2);
    }
  }


  p2t::CDT *cdt = new p2t::CDT(outer_boundary);
  for (i=0; i<holes.size(); i++ ) {
    cdt->AddHole(holes[i]);
  }

  cdt->Triangulate();

  std::vector< p2t::Triangle * > tris = cdt->GetTriangles();

  for (i=0; i<tris.size(); i++) {
    p2t::Triangle * tri = tris[i];
    addEdgesFromTri( edge, tri, edgeMap );
  }

  if (!consistency_check( pwh, edge, pointInfoMap)) {
  } else { }


  Path opath;

  ConvertPolygonWithHolesToWeakPath( opath, pwh, edge, pointInfoMap );

  int n = opath.size();

  if (gFloatRead) {
    for (i=0; i<n; i++) {
      x_d = (double) (opath[i].X/gMulFactor);
      y_d = (double) (opath[i].Y/gMulFactor);
      x_d /= gPremul;
      y_d /= gPremul;
      fprintf(gOutputFp, "%f %f\n", x_d, y_d);
    }
  }
  else {
    for (i=0; i<n; i++) {
      fprintf(gOutputFp, "%lli %lli\n", opath[i].X/gMulFactor, opath[i].Y/gMulFactor );
    }
  }

  delete cdt;
  for (i=0; i<garbage_poly.size(); i++) {
    delete garbage_poly[i];
  }

}

int main(int argc, char **argv) {
  int i, j, k;
  long long int X, Y;
  int line_no=-1;
  bool first = true;
  bool is_outer_boundary = true;


  std::vector<p2t::Point *> outer_boundary;
  std::vector< std::vector<p2t::Point *> > boundaries;
  std::vector< std::vector<p2t::Point *> > holes;
  std::vector<p2t::Point *> poly;

  int polyname = 0;
  double X_d, Y_d;

  PointInfoMap pointInfoMap;
  EdgeMap edgeMap;
  std::map< int, bool > path_visited_map;

  PointAdjacency edge;

  std::vector< std::vector<iPnt> > pwh;

  iPnt xy;
  PointInfo pi;

  iEdge uv;

  std::vector< iPnt > ipoly;

  int debug_var = 0;

  //char buf[1024];
  std::string buf;
  int ch;
  int process_input = 1;

  ClipperLib::Path cur_clip_path;
  ClipperLib::Paths clip_paths;
  ClipperLib::Paths clip_paths_res;

  std::vector<ClipperLib::Path> clip_boundaries;
  std::vector<ClipperLib::Path> clip_holes;

  ClipperLib::Path tpath;

  process_command_line_options( argc, argv );

  //while (fgets(buf, 1024, gInputFp)) {
  while (process_input) {

    ch=fgetc(gInputFp);
    if (ch == EOF) {
      buf += '\n';
      process_input = 0;
    }
    else {
      buf += ch;
    }

    if ((ch != '\n') && (ch != EOF)) { continue; }

    // we've reached a line for processing
    //

    line_no++;

    // if the buffer is empty or it's a comment,
    // skip
    //
    if ( (buf.size()==0) ||
         (buf[0] == '#') ) {

      buf.clear();

      // hacky way to handle special case if file
      // ends with a comment and no new line
      //
      if (ch == EOF) {
        buf += '\n';
      }
      else {
        continue;
      }
    }

    // an empty line (or the end-of-file)
    // means we've finished a polygon and we should
    // add it to the current polygon
    //
    if (buf[0] != '\n') {

      if (gFloatRead) {
        k = sscanf(buf.c_str(), "%lf %lf", &X_d, &Y_d);
      }
      else {
        k = sscanf(buf.c_str(), "%lli %lli", &X, &Y);
      }

      if (k!=2) {
        fprintf(stderr, "invalid read on line %i, exiting\n", line_no);
        exit(2);
      }

      if (gFloatRead) {
        X = (long long int)(gPremul * X_d);
        Y = (long long int)(gPremul * Y_d);
      }

      cur_clip_path.push_back( ClipperLib::IntPoint(X*gMulFactor, Y*gMulFactor) );
    }


    // if we've tied of the polygon with a newline (or eof),
    // save the polygon
    //
    if ( (buf[0] == '\n') ||
         (ch == EOF) ) {

      if (cur_clip_path.size() == 0) {
        buf.clear();
        continue;
      }

      poly.clear();
      ipoly.clear();

      clip_paths.clear();
      clip_paths_res.clear();

      printf("### %i\n", (int)cur_clip_path.size());

      clip_paths.push_back( cur_clip_path );

      ClipperLib::SimplifyPolygons( clip_paths, clip_paths_res, ClipperLib::pftNonZero );

      if (clip_paths_res.size() == 0) {

        if (gIgnoreColinear) {
          cur_clip_path.clear();

          buf.clear();
          continue;
        }

        fprintf(stderr, "Resulting polygon invalid (colinear?).  line_no %i\n", line_no);
        exit(2);
      }

      // copy over to p2t library
      //
      for (i=0; i<clip_paths_res[0].size(); i++) {
        xy.X = clip_paths_res[0][i].X;
        xy.Y = clip_paths_res[0][i].Y;
        poly.push_back( new p2t::Point( xy.X, xy.Y ) );

        ipoly.push_back( xy );
      }

      if (ClipperLib::Area( cur_clip_path ) < 0.0) {
        boundaries.push_back( poly );
        clip_boundaries.push_back( cur_clip_path );
      } else {
        holes.push_back( poly );
        clip_holes.push_back( cur_clip_path );
      }
      pwh.push_back(ipoly);

      cur_clip_path.clear();

      path_visited_map[polyname] = false;

      polyname ++;
      is_outer_boundary = false;

      buf.clear();
      continue;
    }

    buf.clear();

  }

  // Construct a PolyTree that has the heirarchy of polygons with holes
  //
  ClipperLib::PolyTree pt;
  ClipperLib::Clipper clip;

  if (clip_boundaries.size()==0) {
    fprintf(stderr, "no polygons to process?\n");
    exit(-1);
  }

  for (i=0; i<clip_boundaries.size(); i++) {

    // If we have non-trivial boundaries/holes, we need
    // to separate them out.  For example
    //
    //  ---
    //  | |
    //  ---
    //    |
    //  ---
    //  | |
    //  ---
    //
    ClipperLib::Paths   dirty_polys, clean_polys;
    ClipperLib::Path    clean_poly;
    ClipperLib::Clipper clean_clip;

    dirty_polys.push_back( clip_boundaries[i] );
    clean_clip.AddPaths( dirty_polys, ClipperLib::ptSubject, true );
    clean_clip.Execute( ClipperLib::ctUnion, clean_polys, ClipperLib::pftNegative, ClipperLib::pftNegative );

    for (j=0; j<clean_polys.size(); j++) {
      clip.AddPath( clean_polys[j], ClipperLib::ptSubject, true );
    }

  }

  for (i=0; i<clip_holes.size(); i++) {

    // Do the same as above with the boundaries as with holes.
    //
    ClipperLib::Paths   dirty_polys, clean_polys;
    ClipperLib::Path    clean_poly;
    ClipperLib::Clipper clean_clip;

    ClipperLib::ReversePath( clip_holes[i] );
    dirty_polys.push_back( clip_holes[i] );
    clean_clip.AddPaths( dirty_polys, ClipperLib::ptSubject, true );
    clean_clip.Execute( ClipperLib::ctUnion, clean_polys, ClipperLib::pftNegative, ClipperLib::pftNegative );

    for (j=0; j<clean_polys.size(); j++) {
      ClipperLib::ReversePath( clean_polys[j] );
      clip.AddPath( clean_polys[j], ClipperLib::ptSubject, true );
    }
  }


  // Now do the big operation..
  //
  clip.Execute( ClipperLib::ctUnion,
                pt,
                ClipperLib::pftNonZero,
                ClipperLib::pftNonZero );


  std::vector< std::vector< std::vector< iPnt > > > pwhs;
  std::vector< std::vector< iPnt > > dummy_pwh;

  constructPolygonsWithHoles( pt.AllNodes[0], pwhs, -1, 0 );

  for (i=0; i<pwhs.size(); i++) {
    fprintf( gOutputFp, "# CPWH %i\n", i);
    emitWeakPWH( pwhs[i] );
    fprintf( gOutputFp, "\n" );
  }

  if (gOutputFp != stdout) {
    fclose(gOutputFp);
  }

}
