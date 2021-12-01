/*
 * fem element type.h
 *
 *  Created on: March 28, 2013
 *      Author: Christoph Paulus
 */

#include "utils.h"

// defines
#ifndef __FEM_ELEMENTTYPE_H
#define __FEM_ELEMENTTYPE_H

#define DIM 3

/*
class fem element type:
- basis for the classes like the class linear tetrahedra

could be the basis for other classes that inherit from this class, for example quadratic tetra hedra
*/

class FEM_Elementtype
{
public:
	// functions
	virtual ~FEM_Elementtype() {};
};

#endif
