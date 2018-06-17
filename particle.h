#ifndef	_PARTICLE_H
#define	_PARTICLE_H
#include "common.h"
// all of the particle attributes are normalized.
class Particle
{
  public:
    int         ID;          // species number
    Vector3     grid;        // 3d-vector partice position
    Vector3     uu;	     // ux=vx*gamma
    Scalar      gamma;       // Lorenta factor
    Scalar      q;           // charge (e)
    Scalar      m;           // mass   (me)
    Vector3	e_part;	     // electric field experienced by this particle.
    Vector3	b_part;
    Vector3	a_part;
#ifdef QED_BLOCK
    Scalar      eta;
    Scalar	optical_depth;
    int 	radn;
    Scalar	radt;
    Scalar	rad_px;
    Scalar	rad_py;
    Scalar	rad_pz;
#endif
  
    Particle(Particle *_part) {
	ID        = _part->ID;
	grid = _part->grid;
	uu   = _part->uu;
	gamma     =   _part->gamma;
	q         =   _part->q;
	m         =   _part->m;
	e_part    =   _part->e_part;
	b_part    =   _part->b_part;
	a_part    =   _part->a_part;
#ifdef QED_BLOCK
    eta       = _part->eta;
	optical_depth = _part->optical_depth;
  	radn      = _part->radn;
	radt      = _part->radt;  	
	rad_px    = _part->rad_px;  	
	rad_py    = _part->rad_py;  	
	rad_pz    = _part->rad_pz;  	
#endif
    };

    Particle& operator = (Particle _part) {
        ID        = _part.ID;
        grid = _part.grid;
        uu   = _part.uu;
        gamma     =   _part.gamma;
        q         =   _part.q;
        m         =   _part.m;
        e_part    =   _part.e_part;
        b_part    =   _part.b_part;
        a_part    =   _part.a_part;
#ifdef QED_BLOCK
        eta = _part.eta;
	    optical_depth = _part.optical_depth;
  	    radn = _part.radn;
	    radt = _part.radt;  	
	    rad_px = _part.rad_px;  	
	    rad_py = _part.rad_py;  	
	    rad_pz = _part.rad_pz;  	
#endif
        return *this;
    };

    Particle() {
    };
};
#endif //_PARTICLE_H
