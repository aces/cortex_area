
/*
   CORTEX AREA

   Compute the area of each patch of triangular faces defining zones
   on the cortex.

   cortex_area -surface cortex_file_name -zone zone_file_name -output out_file

   Values: cortex_file_name = file name for a cortical surface, where the 
                              input_signal is diffused. (D. MacDonalds' .obj 
                              file format)
           zone_file_name = list of zone id's, by node
           out_file = file name for output (if none, use stdout)

   HISTORY:   VERSION 1.0  May, 2005 Initial implementation (Claude Lepage)

   COPYRIGHT: 
              Permission to use, copy, modify, and distribute this
              software and its documentation for any purpose and without
              fee is hereby granted, provided that the above copyright
              notice appear in all copies.  The author and McGill University
              make no representations about the suitability of this
              software for any purpose.  It is provided "as is" without
              express or implied warranty.
*/

#include <stdio.h>

#include "volume_io/internal_volume_io.h"
#include "bicpl.h"

// Prototypes of functions in this file.

private void usage( char * );
private Status read_surface_obj( STRING, int *, Point *[],
                                 Vector *[], int *[], int **[] );
private Status read_zone_ids( STRING, int, int [] );
private Status write_zones( STRING, int, Real areas[] );
private Status compute_zones( int, Point *, int *, int **, int *,
                              int *, Real *[] );
private Status get_surface_neighbours( polygons_struct *, int *[],
                                       int ** [] );

// Main program.

int main( int argc, char * argv[] ) {

  char   * cortex_file_name;   // file name for cortex file
  char   * zone_file_name;     // file name for cortex zones
  char   * output_file_name;   // file name for the output
  int      n_points;           // number of grid points
  Point  * coords;             // coordinates
  Vector * normals;            // normal vectors
  int      n_zones;            // number of zones
  int    * zone_ids;           // zone ids of nodes
  Real   * areas;              // area of each zone
  int    * n_ngh;              // node neighbours (inverse connectivity)
  int   ** ngh;
  STRING   arg;

  // Defaults
  cortex_file_name = NULL;
  zone_file_name = NULL;
  output_file_name = NULL;

  // Parse the command line arguments for the options.

  initialize_argument_processing( argc, argv );

  while( get_string_argument( NULL, &arg ) ) {

    if( equal_strings( arg, "-surface" ) ) {
      if( !get_string_argument( NULL, &cortex_file_name ) ) {
        print_error( "Error in -surface arguments.\n" );
        usage( argv[0] );
        return( 1 );
      }
    } else if( equal_strings( arg, "-zone" ) ) {
      if( !get_string_argument( NULL, &zone_file_name ) ) {
        print_error( "Error in -zone arguments.\n" );
        usage( argv[0] );
        return( 1 );
      }
    } else if( equal_strings( arg, "-output" ) ) {
      if( !get_string_argument( NULL, &output_file_name ) ) {
        print_error( "Error in -output arguments.\n" );
        usage( argv[0] );
        return( 1 );
      }
    } else {
      usage( argv[0] );
      return( 1 );
    }
  }

  // Check validity of the arguments.

  if( cortex_file_name == NULL ) {
    print_error( "Must supply a name for the surface file.\n" );
    usage( argv[0] );
    return( 1 );
  }
  if( zone_file_name == NULL ) {
    print_error( "Must supply a name for the cortex zone file.\n" );
    usage( argv[0] );
    return( 1 );
  }

  // Read the surface file.

  if( read_surface_obj( cortex_file_name, &n_points,
                        &coords, &normals, &n_ngh, &ngh ) != OK ) {
    return 1;
  }
  FREE( normals );   // not needed

  // Read the signal file.

  ALLOC( zone_ids, n_points );
  if( read_zone_ids( zone_file_name, n_points, zone_ids ) != OK ) {
    return 1;
  }

  if( compute_zones( n_points, coords, n_ngh, ngh, zone_ids, 
                     &n_zones, &areas ) != OK ) {
    return 1;
  }

  if( write_zones( output_file_name, n_zones, areas ) != OK ) {
    return 1;
  }

  FREE( areas );
  FREE( zone_ids );
  FREE( coords );
  FREE( n_ngh );
  FREE( ngh[0] );   // this is ngh_array
  FREE( ngh );
}


// -------------------------------------------------------------------
// Help message on how to use this module.
//
private void usage( char * executable_name ) {

  STRING  usage_format = "\
Usage: %s -surface file.obj -zone file.txt -output file.out\n\n";

  print_error( usage_format, executable_name );
}


// -------------------------------------------------------------------
// Load the cortical surface.
//
// filename: name of the .obj file
// n_points: the number of the vertices
// points: (x,y,z) coordinates
// normals: normal vectors
// n_neighbours: number of vertices around each node
// neighbours: the set of ordered triangle consisting of the vertices
//
private Status read_surface_obj( STRING filename,
                                 int * n_points,
                                 Point * points[],
                                 Vector * normals[],
                                 int * n_neighbours[],
                                 int ** neighbours[] ) {

  int               i, n_objects;
  object_struct  ** object_list;
  polygons_struct * surface;
  File_formats      format;
  STRING            expanded;

  expanded = expand_filename( filename );   // why?????

  int err = input_graphics_file( expanded, &format, &n_objects,
                                 &object_list );

  if( err != OK ) {
    print_error( "Error reading file %s\n", expanded );
    return( ERROR );
  }

  if( n_objects != 1 || 
      ( n_objects == 1 && get_object_type(object_list[0]) != POLYGONS ) ) {
    print_error( "Error in contents of file %s\n", expanded );
    return( ERROR );
  }

  delete_string( expanded );

  surface = get_polygons_ptr( object_list[0] );

  // Make a copy of the coordinates and the normals, since
  // delete_object_list will destroy them.

  *n_points = surface->n_points;
  ALLOC( *points, surface->n_points );
  ALLOC( *normals, surface->n_points );
  for( i = 0; i < *n_points; i++ ) {
    (*points)[i].coords[0] = surface->points[i].coords[0];
    (*points)[i].coords[1] = surface->points[i].coords[1];
    (*points)[i].coords[2] = surface->points[i].coords[2];
    (*normals)[i].coords[0] = surface->normals[i].coords[0];
    (*normals)[i].coords[1] = surface->normals[i].coords[1];
    (*normals)[i].coords[2] = surface->normals[i].coords[2];
  }

  get_surface_neighbours( surface, n_neighbours, neighbours );

  delete_object_list( n_objects, object_list );

  return( OK );
}

// -------------------------------------------------------------------
// Read the zone ids. These should be integers, but somehow, they
// are floats.
//
private Status read_zone_ids( STRING filename,
                              int n_points,
                              int zone_ids[] ) {

  int    i;
  FILE * fp;
  double val;

  fp = fopen( filename, "r" );
  if( fp == NULL ) {
    print_error( "Error opening file %s\n", filename );
    return( ERROR );
  }

  for( i = 0; i < n_points; i++ ) {
    if( fscanf( fp, "%lf", &val ) != 1 ) break;
    zone_ids[i] = (int)val;
  }

  fclose( fp );

  if( i != n_points ) return( ERROR );

  return( OK );
}


// -------------------------------------------------------------------
// Write the area of each cortex zone.
//
private Status write_zones( STRING filename, 
                            int n_zones,
                            Real areas[] ) {

  int    i;
  Real   total_area = 0.0;

  if( filename ) {

    // Save to a file
    FILE * fp;
    fp = fopen( filename, "w" );
    if( fp == NULL ) {
      print_error( "Error opening file %s\n", filename );
      return( ERROR );
    }

    for( i = 0; i < n_zones; i++ ) {
      fprintf( fp, "%i  %f\n", i, areas[i] );
      total_area += areas[i];
    }
    fprintf( fp, "Total %f\n", total_area );

    fclose( fp );

  } else {

    // Write to the screen.

    for( i = 0; i < n_zones; i++ ) {
      printf( "%i  %f\n", i, areas[i] );
      total_area += areas[i];
    }
    printf( "Total %f\n", total_area );

  }

  return( OK );
}


// -------------------------------------------------------------------
// Calculation of the area of each cortex zone.
// 
// n_points: the number of the vertices
// xyz: (x,y,z) coordinates
// n_ngh: number of vertices around each node
// ngh: the set of ordered triangle consisting of the vertices
// zone_ids: id of cortex zones, defined at the nodes
// z_zones: number of zones
// areas: area of each zone
//
private Status compute_zones( int n_points,
                              Point * xyz,
                              int * n_ngh,
                              int ** ngh,
                              int * zone_ids,
                              int * n_zones,
                              Real * areas[] ) {

  int     i, j, nn, zid;

  // Find the maximum zone id for this cortex surface.

  *n_zones = 0;
  for( i = 0; i < n_points; i++ ) {
    if( zone_ids[i] > *n_zones ) *n_zones = zone_ids[i];
  }
  (*n_zones)++;  // count 0

  // Allocation of the memory to store the area for the zones.

  Real * z_area;
  ALLOC( z_area, *n_zones );
  *areas = z_area;
  for( zid = 0; zid < *n_zones; zid++ ) {
    z_area[zid] = 0;
  }

  // Loop over each node to compute an area associated to it.

  for( i = 0; i < n_points; i++ ) {
    zid = zone_ids[i];

    for( nn = 0; nn < n_ngh[i]; nn++ ) {

      int n1 = ngh[i][nn];
      int n2 = ngh[i][(nn+1)%n_ngh[i]];

      Real v1[3], v2[3], vcent[3], vcross[3];

      for( j = 0; j < 3; j++ ) {
        v1[j] = 0.5 * ( xyz[n1].coords[j] - xyz[i].coords[j] );
        v2[j] = 0.5 * ( xyz[n2].coords[j] - xyz[i].coords[j] );
        vcent[j] = ( xyz[i].coords[j] + xyz[n1].coords[j] + 
                     xyz[n2].coords[j] ) / 3.0 - xyz[i].coords[j];
      }

      vcross[0] = v1[1]*vcent[2] - v1[2]*vcent[1];
      vcross[1] = v1[2]*vcent[0] - v1[0]*vcent[2];
      vcross[2] = v1[0]*vcent[1] - v1[1]*vcent[0];
      z_area[zid] += sqrt( vcross[0]*vcross[0] + vcross[1]*vcross[1] +
                           vcross[2]*vcross[2] );

      vcross[0] = v2[1]*vcent[2] - v2[2]*vcent[1];
      vcross[1] = v2[2]*vcent[0] - v2[0]*vcent[2];
      vcross[2] = v2[0]*vcent[1] - v2[1]*vcent[0];
      z_area[zid] += sqrt( vcross[0]*vcross[0] + vcross[1]*vcross[1] +
                           vcross[2]*vcross[2] );
    }
  }

  // Must scale by 2.0 (since vcross = 2 * area).
  for( zid = 0; zid < *n_zones; zid++ ) {
    z_area[zid] *= 0.5;
  }

  return( OK );

}

// -------------------------------------------------------------------
// Construct the edges around each node. The edges are sorted to
// make an ordered closed loop to reconstruct the triangles.
// 
private Status get_surface_neighbours( polygons_struct * surface,
                                       int * n_neighbours_return[],
                                       int ** neighbours_return[] ) {

  int    i, j, k, jj;
  int  * tri;
  int  * n_ngh;
  int ** ngh;
  int  * ngh_array;

  // Check if all polygons are triangles.

  if( 3 * surface->n_items != surface->end_indices[surface->n_items-1] ) {
    printf( "Surface must contain only triangular polygons.\n" );
    return ERROR;
  }

  // Check if the node numbering starts at 0 or 1.

  int min_idx, max_idx;

  min_idx = 100*surface->n_points;  // anything big
  max_idx = 0;                      // anything small

  for( i = 0; i < 3*surface->n_items; i++ ) {
    if( surface->indices[i] < min_idx ) min_idx = surface->indices[i];
    if( surface->indices[i] > max_idx ) max_idx = surface->indices[i];
  }

  // Shift numbering to start at zero, for array indexing. Note 
  // that we don't care if surface->indices array is modified.

  if( min_idx != 0 ) {
    for( i = 0; i < 3*surface->n_items; i++ ) {
      surface->indices[i] -= min_idx;
    }
  }

  // Count number of triangles attached to each node.

  ALLOC( n_ngh, surface->n_points );
  ALLOC( ngh, surface->n_points );
  ALLOC( ngh_array, 3*surface->n_items );

  for( i = 0; i < surface->n_points; i++ ) {
    n_ngh[i] = 0;
  }

  for( i = 0; i < 3*surface->n_items; i++ ) {
    n_ngh[surface->indices[i]]++;
    ngh_array[i] = -1;
  }

  int max_ngh = 0;
  int sum_ngh = 0;
  for( i = 0; i < surface->n_points; i++ ) {
    ngh[i] = &(ngh_array[sum_ngh]);
    sum_ngh += n_ngh[i];
    max_ngh = MAX( max_ngh, n_ngh[i] );
  }

  // At first, store the indices of the triangles in the neighbours.
  for( i = 0; i < surface->n_items; i++ ) {
    for( j = 0; j < 3; j++ ) {
      jj = surface->indices[3*i+j];
      for( k = 0; k < n_ngh[jj]; k++ ) {
        if( ngh[jj][k] == -1 ) {
          ngh[jj][k] = i;
          break;
        }
      }
    }
  }

  // Now create a sort closed loop of the node neighbours.
  // This is needed by the parametric=0 FEM algorithm.
  //
  //         1 ----- 2
  //          /\   /\
  //         /  \ /  \
  //       0 ----P---- 3
  //         \  / \  /
  //          \/   \/
  //         5 ----- 4
  //

  int * tmp;
  ALLOC( tmp, 2*max_ngh );

  for( i = 0; i < surface->n_points; i++ ) {
    for( k = 0; k < n_ngh[i]; k++ ) {
      tri = &(surface->indices[3*ngh[i][k]]);
      for( j = 0; j < 3; j++ ) {
        if( tri[j] == i ) break;
      }
      tmp[2*k+0] = tri[(j+1)%3];
      tmp[2*k+1] = tri[(j+2)%3];
    }
  
    ngh[i][0] = tmp[0];
    ngh[i][1] = tmp[1];
    for( k = 2; k < n_ngh[i]; k++ ) {
      for( j = 1; j < n_ngh[i]; j++ ) {
        if( tmp[2*j] == ngh[i][k-1] || tmp[2*j+1] == ngh[i][k-1] ) {
          if( tmp[2*j] == ngh[i][k-1] ) {
            ngh[i][k] = tmp[2*j+1];
          } else {
            ngh[i][k] = tmp[2*j];
          }
          tmp[2*j] = -1;
          tmp[2*j+1] = -1;
          break;
        }
      }
    }
  }

  *n_neighbours_return = n_ngh;
  *neighbours_return = ngh;

  FREE( tmp );

}

