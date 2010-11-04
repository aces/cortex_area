
/*
   CORTEX VOLUME

   Compute the volume of each prism between the white and gray surfaces of
   the cortex.

   cortex_volume [-zero] -white white_surface.obj -gray gray_surface.obj [-output out_file]

   Options: -zero = set negative volumes to zero (crossing/inverted surfaces)
   Values: white_surface.obj = file name for the white cortical surface
           gray_surface.obj = file name for the gray cortical surface
           out_file = file name for output (if none, use stdout)

   HISTORY:   VERSION 1.3  Nov, 2010 Initial implementation (Claude Lepage)

   COPYRIGHT: Copyright Alan C. Evans
              Professor of Neurology
              McGill University
*/

#include <stdio.h>

#include "volume_io/internal_volume_io.h"
#include "bicpl.h"

// Prototypes of functions in this file.

private void usage( char * );
private Status read_surface_obj( STRING, int *, Point *[],
                                 Vector *[], int *[], int **[] );
private Status write_vertex( STRING, int, Real areas[] );
private Status compute_vertex_volume( int, Point *, Point *,
                                      int *, int **, Real [], int );
private Real prism_volume( Point *, Point *, Point *, Point *, Point *, Point * );
private Status get_surface_neighbours( polygons_struct *, int *[],
                                       int ** [] );

// Main program.

int main( int argc, char * argv[] ) {

  char   * white_file_name;    // file name for cortex file
  char   * gray_file_name;     // file name for cortex file
  char   * output_file_name;   // file name for the output
  int      wn_points;          // number of grid points (white)
  Point  * wcoords;            // coordinates (white)
  int    * wn_ngh;             // node neighbours (inverse connectivity) (white)
  int   ** wngh;
  int      gn_points;          // number of grid points (gray)
  Point  * gcoords;            // coordinates (gray)
  int    * gn_ngh;             // node neighbours (inverse connectivity) (gray)
  int   ** gngh;
  Vector * normals;            // normal vectors

  Real   * vert_volume;        // volume at each vertex
  STRING   arg;
  int      status;
  int      zero;               // set -ve volume to zero

  // Defaults
  zero = 0;
  white_file_name = NULL;
  gray_file_name = NULL;
  output_file_name = NULL;

  // Parse the command line arguments for the options.

  initialize_argument_processing( argc, argv );

  while( get_string_argument( NULL, &arg ) ) {

    if( equal_strings( arg, "-white" ) ) {
      if( !get_string_argument( NULL, &white_file_name ) ) {
        print_error( "Error in -white argument.\n" );
        usage( argv[0] );
        return( 1 );
      }
    } else if( equal_strings( arg, "-gray" ) ) {
      if( !get_string_argument( NULL, &gray_file_name ) ) {
        print_error( "Error in -gray argument.\n" );
        usage( argv[0] );
        return( 1 );
      }
    } else if( equal_strings( arg, "-output" ) ) {
      if( !get_string_argument( NULL, &output_file_name ) ) {
        print_error( "Error in -output argument.\n" );
        usage( argv[0] );
        return( 1 );
      }
    } else if( equal_strings( arg, "-zero" ) ) {
      zero = 1;
    } else {
      usage( argv[0] );
      return( 1 );
    }
  }

  // Check validity of the arguments.

  if( white_file_name == NULL ) {
    print_error( "Must supply a name for the white surface file.\n" );
    usage( argv[0] );
    return( 1 );
  }

  if( gray_file_name == NULL ) {
    print_error( "Must supply a name for the gray surface file.\n" );
    usage( argv[0] );
    return( 1 );
  }

  // Read the surface file.

  if( read_surface_obj( white_file_name, &wn_points,
                        &wcoords, &normals, &wn_ngh, &wngh ) != OK ) {
    return 1;
  }
  FREE( normals );   // not needed

  if( read_surface_obj( gray_file_name, &gn_points,
                        &gcoords, &normals, &gn_ngh, &gngh ) != OK ) {
    FREE( wcoords );
    FREE( wn_ngh );
    FREE( wngh[0] );   // this is ngh_array
    FREE( wngh );
    return 1;
  }
  FREE( normals );   // not needed
  FREE( gn_ngh );
  FREE( gngh[0] );   // this is ngh_array
  FREE( gngh );

  // The two surfaces must have the same number of points.
  if( wn_points != gn_points ) {
    FREE( wcoords );
    FREE( wn_ngh );
    FREE( wngh[0] );   // this is ngh_array
    FREE( wngh );
    FREE( gcoords );
    return 1;
  }

  ALLOC( vert_volume, wn_points );
  status = compute_vertex_volume( wn_points, wcoords, gcoords, wn_ngh, 
                                  wngh, vert_volume, zero );

  FREE( wcoords );
  FREE( wn_ngh );
  FREE( wngh[0] );   // this is ngh_array
  FREE( wngh );
  FREE( gcoords );

  if( status != OK ) return 1;

  // Simply write out the volume at each vertex.
  status =  write_vertex( output_file_name, wn_points, vert_volume );

  FREE( vert_volume );

  return( status );
}


// -------------------------------------------------------------------
// Help message on how to use this module.
//
private void usage( char * executable_name ) {

  STRING  usage_format = "\
Usage: %s [-zero] -white white_surface.obj -gray gray_surface.obj [-output file.out]\n\n\
Options: -zero = set negative volumes to zero (crossing/inverted surfaces)\n\n\
Copyright Alan C. Evans\nProfessor of Neurology\nMcGill University\n";

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
// Write the area of each vertex.
//
private Status write_vertex( STRING filename, int n_points,
                             Real volumes[] ) {

  int    i;
  Real   total_volume = 0.0;

  if( filename ) {

    // Save to a file
    FILE * fp;
    fp = fopen( filename, "w" );
    if( fp == NULL ) {
      print_error( "Error opening file %s\n", filename );
      return( ERROR );
    }

    for( i = 0; i < n_points; i++ ) {
      fprintf( fp, "%f\n", volumes[i] );
      total_volume += volumes[i];
    }

    fclose( fp );

  } else {

    // Write to the screen.

    for( i = 0; i < n_points; i++ ) {
      printf( "%f\n", volumes[i] );
      total_volume += volumes[i];
    }

  }
  printf( "Total %f\n", total_volume );

  return( OK );
}


// -------------------------------------------------------------------
// Calculation of the area of each vertex.
// 
// n_points: the number of the vertices
// xyz: (x,y,z) coordinates
// n_ngh: number of vertices around each node
// ngh: the set of ordered triangle consisting of the vertices
// vert_area: area of each vertex
//
private Status compute_vertex_volume( int n_points,
                                      Point * wxyz, Point * gxyz,
                                      int * n_ngh,
                                      int ** ngh,
                                      Real vert_volume[],
                                      int zero ) {

  int     i, j, nn;
  Real    vol;

  // Loop over each node to compute an area associated to it.

  for( i = 0; i < n_points; i++ ) {
    vert_volume[i] = 0.0;
  }

  for( i = 0; i < n_points; i++ ) {

    for( nn = 0; nn < n_ngh[i]; nn++ ) {

      int n1 = ngh[i][nn];
      int n2 = ngh[i][(nn+1)%n_ngh[i]];

      if( i < n1 && i < n2 ) {

        vol = prism_volume( &(wxyz[i]), &(wxyz[n1]), &(wxyz[n2]),
                            &(gxyz[i]), &(gxyz[n1]), &(gxyz[n2]) );
        if( zero && ( vol < 0.0 ) ) vol = 0.0;

        vert_volume[i] += vol;
        vert_volume[n1] += vol;
        vert_volume[n2] += vol;
      }

    }
  }

  for( i = 0; i < n_points; i++ ) {
    vert_volume[i] /= 3.0;
  }

  return( OK );

}

// -------------------------------------------------------------------
// Compute the volume of a prisms with vertices (v0,v1,v2,v3,v4,v5).
//
//              5-----------4        N0 = ( 1 - xi - eta ) * ( 1 - zeta ) / 2
//              |\         /|        N1 = xi * ( 1 - zeta ) / 2
//              | \       / |        N2 = eta * ( 1 - zeta ) / 2
//              |  \     /  |        N3 = ( 1 - xi - eta ) * ( 1 + zeta ) / 2
//              |   \   /   |        N4 = xi * ( 1 + zeta ) / 2
//              |    \ /    |        N5 = eta * ( 1 + zeta ) / 2
//              |     3     |
//              |     |     |
//              |     |     |
//              2-----|-----1
//               \    |    /
//                \   |   /
//                 \  |  /
//                  \ | /
//                   \|/
//                    0
// 
private Real prism_volume( Point * v0, Point * v1, Point * v2,
                           Point * v3, Point * v4, Point * v5 ) {

  int  i, j, nz, ntri;
  Real N0_xi, N1_xi, N3_xi, N4_xi;
  Real N0_eta, N2_eta, N3_eta, N5_eta;
  Real N0_zeta, N1_zeta, N2_zeta, N3_zeta, N4_zeta, N5_zeta;
  Real dx_dxi, dx_deta, dx_dzeta, dy_dxi, dy_deta, dy_dzeta, 
       dz_dxi, dz_deta, dz_dzeta;
  Real xi[10], eta[10], zeta[10], wgt_tri[10], wgt_lin[10];
  Real det;

  // Reference: http://www.cs.rpi.edu/~flaherje/pdf/fea6.pdf
  // Integration in zeta:
  //                    zeta        wgt
  // n=1, degree=1:       1          2
  // n=2, degree=3:    sqrt(1/3)     1
  //                  -sqrt(1/3)     1
  // n=3, degree=5:   0.0           8/9
  //                  sqrt(15)/5    5/9
  //                 -sqrt(15)/5    5/9
  // n=4, degree=7:   sqrt(525-70*sqrt(30))/35    (18+sqrt(30))/36
  //                 -sqrt(525-70*sqrt(30))/35    (18+sqrt(30))/36
  //                  sqrt(525+70*sqrt(30))/35    (18-sqrt(30))/36
  //                 -sqrt(525+70*sqrt(30))/35    (18-sqrt(30))/36
  // Integration over triangle (symmetric):
  //                     1-xi-eta         xi          eta        wgt
  // n=1, degree=1     0.333333333   0.33333333  0.33333333  1.0000000
  // n=3, degree=2     0.666666666   0.16666666  0.16666666  0.3333333
  //                   0.166666666   0.66666666  0.16666666  0.3333333
  //                   0.166666666   0.16666666  0.66666666  0.3333333
  // n=4, degree=4     0.333333333   0.333333333  0.333333333 -0.5625000
  //                   0.600000000   0.200000000  0.200000000  0.5208333
  //                   0.200000000   0.600000000  0.200000000  0.5208333
  //                   0.200000000   0.200000000  0.600000000  0.5208333

  // ntri = 1 is enough to integrate "det" exactly (a polynomial
  // of max degree 1 in xi, eta).
  ntri = 1;
  if( ntri == 1 ) {
    xi[0] = 0.333333333333;    eta[0] = 0.333333333333;
    wgt_tri[0] = 1.000000000000;
  } else if( ntri == 3 ) {
    xi[0] = 0.166666666667;    eta[0] = 0.166666666667;
    xi[1] = 0.666666666667;    eta[1] = 0.166666666667;
    xi[2] = 0.166666666667;    eta[2] = 0.666666666667;
    wgt_tri[0] = 0.333333333333;
    wgt_tri[1] = 0.333333333333;
    wgt_tri[2] = 0.333333333333;
  } else if( ntri == 4 ) {
    xi[0] = 0.333333333333;    eta[0] = 0.333333333333;
    xi[1] = 0.200000000000;    eta[1] = 0.200000000000;
    xi[2] = 0.600000000000;    eta[2] = 0.200000000000;
    xi[3] = 0.200000000000;    eta[3] = 0.600000000000;
    wgt_tri[0] = -0.562500000000;
    wgt_tri[1] =  0.520833333333;
    wgt_tri[2] =  0.520833333333;
    wgt_tri[3] =  0.520833333333;
  }

  // nz = 2 is enough to integrate "det" exactly (a polynomial
  // of max degree 3 in zeta).
  nz = 2;
  if( nz == 1 ) {
    zeta[0] =  0.000000000000;   wgt_lin[0] = 1.000000000000;
  } else if( nz == 2 ) {
    zeta[0] =  0.577350269190;   wgt_lin[0] = 0.500000000000;
    zeta[1] = -0.577350269190;   wgt_lin[1] = 0.500000000000;
  } else if( nz == 3 ) {
    zeta[0] =  0.000000000000;   wgt_lin[0] = 0.444444444444;
    zeta[1] =  0.774596669240;   wgt_lin[1] = 0.277777777778;
    zeta[2] = -0.774596669240;   wgt_lin[2] = 0.277777777778;
  } else if( nz == 4 ) {
    zeta[0] =  0.339981043570;   wgt_lin[0] = 0.326072577430;
    zeta[1] = -0.339981043570;   wgt_lin[1] = 0.326072577430;
    zeta[2] =  0.861136311600;   wgt_lin[2] = 0.173927422570;
    zeta[3] = -0.861136311600;   wgt_lin[3] = 0.173927422570;
  }

  det = 0.0;
  for( i = 0; i < ntri; i++ ) {
    for( j = 0; j < nz; j++ ) {

      // xi-derivatives (factor of 2):
      N0_xi = -( 1.0 - zeta[j] );
      N1_xi =  ( 1.0 - zeta[j] );
      N3_xi = -( 1.0 + zeta[j] );
      N4_xi =  ( 1.0 + zeta[j] );

      // eta-derivatives (factor of 2):
      N0_eta = -( 1.0 - zeta[j] );
      N2_eta = ( 1.0 - zeta[j] );
      N3_eta = -( 1.0 + zeta[j] );
      N5_eta = ( 1.0 + zeta[j] );

      // zeta-derivatives (factor of 2):
      N0_zeta = -( 1.0 - xi[i] - eta[i] );
      N1_zeta = -xi[i];
      N2_zeta = -eta[i];
      N3_zeta = ( 1.0 - xi[i] - eta[i] );
      N4_zeta = xi[i];
      N5_zeta = eta[i];
 
      dx_dxi = N0_xi * v0->coords[0] + N1_xi * v1->coords[0] + 
               N3_xi * v3->coords[0] + N4_xi * v4->coords[0];
      dy_dxi = N0_xi * v0->coords[1] + N1_xi * v1->coords[1] + 
               N3_xi * v3->coords[1] + N4_xi * v4->coords[1];
      dz_dxi = N0_xi * v0->coords[2] + N1_xi * v1->coords[2] + 
               N3_xi * v3->coords[2] + N4_xi * v4->coords[2];
  
      dx_deta = N0_eta * v0->coords[0] + N2_eta * v2->coords[0] + 
                N3_eta * v3->coords[0] + N5_eta * v5->coords[0];
      dy_deta = N0_eta * v0->coords[1] + N2_eta * v2->coords[1] +
                N3_eta * v3->coords[1] + N5_eta * v5->coords[1];
      dz_deta = N0_eta * v0->coords[2] + N2_eta * v2->coords[2] + 
                N3_eta * v3->coords[2] + N5_eta * v5->coords[2];
  
      dx_dzeta = N0_zeta * v0->coords[0] + N1_zeta * v1->coords[0] + 
                 N2_zeta * v2->coords[0] + N3_zeta * v3->coords[0] + 
                 N4_zeta * v4->coords[0] + N5_zeta * v5->coords[0];
      dy_dzeta = N0_zeta * v0->coords[1] + N1_zeta * v1->coords[1] + 
	         N2_zeta * v2->coords[1] + N3_zeta * v3->coords[1] + 
                 N4_zeta * v4->coords[1] + N5_zeta * v5->coords[1];
      dz_dzeta = N0_zeta * v0->coords[2] + N1_zeta * v1->coords[2] + 
	         N2_zeta * v2->coords[2] + N3_zeta * v3->coords[2] + 
	         N4_zeta * v4->coords[2] + N5_zeta * v5->coords[2];

      det += ( dx_dxi * ( dy_deta * dz_dzeta - dy_dzeta * dz_deta ) +
               dx_deta * ( dy_dzeta * dz_dxi - dy_dxi * dz_dzeta ) +
               dx_dzeta * ( dy_dxi * dz_deta - dy_deta * dz_dxi ) ) * 
              wgt_tri[i] * wgt_lin[j];
    }
  }

  det *= 0.125;  // factor of 1/8 in 3D

  return( det );

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

