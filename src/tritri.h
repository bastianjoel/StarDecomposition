/* Triangle/triangle intersection test routine,
 * by Tomas Moller, 1997.
 * See article "A Fast Triangle-Triangle Intersection Test",
 * Journal of Graphics Tools, 2(2), 1997
 *
 * int tri_tri_intersect(float V0[3],float V1[3],float V2[3],
 *                         float U0[3],float U1[3],float U2[3])
 *
 * parameters: vertices of triangle 1: V0,V1,V2
 *             vertices of triangle 2: U0,U1,U2
 * result    : returns 1 if the triangles intersect, otherwise 0
 *
 */

#pragma once
#include "vectorq.h"
#include <math.h>

/* sort so that a<=b */
#define SORT(a,b)       \
             if(a>b)    \
             {          \
               mpq_class c; \
               c=a;     \
               a=b;     \
               b=c;     \
             }

#define ISECT(VV0,VV1,VV2,D0,D1,D2,isect0,isect1) \
              isect0=VV0+(VV1-VV0)*D0/(D0-D1);    \
              isect1=VV0+(VV2-VV0)*D0/(D0-D2);


#define COMPUTE_INTERVALS(VV0,VV1,VV2,D0,D1,D2,D0D1,D0D2,isect0,isect1) \
  if(D0D1>0.0f)                                         \
  {                                                     \
    /* here we know that D0D2<=0.0 */                   \
    /* that is D0, D1 are on the same side, D2 on the other or on the plane */ \
    ISECT(VV2,VV0,VV1,D2,D0,D1,isect0,isect1);          \
  }                                                     \
  else if(D0D2>0.0f)                                    \
  {                                                     \
    /* here we know that d0d1<=0.0 */                   \
    ISECT(VV1,VV0,VV2,D1,D0,D2,isect0,isect1);          \
  }                                                     \
  else if(D1*D2>0.0f || D0!=0.0f)                       \
  {                                                     \
    /* here we know that d0d1<=0.0 or that D0!=0.0 */   \
    ISECT(VV0,VV1,VV2,D0,D1,D2,isect0,isect1);          \
  }                                                     \
  else if(D1!=0.0f)                                     \
  {                                                     \
    ISECT(VV1,VV0,VV2,D1,D0,D2,isect0,isect1);          \
  }                                                     \
  else if(D2!=0.0f)                                     \
  {                                                     \
    ISECT(VV2,VV0,VV1,D2,D0,D1,isect0,isect1);          \
  }                                                     \
  else                                                  \
  {                                                     \
    /* triangles are coplanar */                        \
    return coplanar_tri_tri(N1,V0,V1,V2,U0,U1,U2);      \
  }



/* this edge to edge test is based on Franlin Antonio's gem:
   "Faster Line Segment Intersection", in Graphics Gems III,
   pp. 199-202 */ 
#define EDGE_EDGE_TEST(V0,U0,U1)                      \
  Bx=U0[i0]-U1[i0];                                   \
  By=U0[i1]-U1[i1];                                   \
  Cx=V0[i0]-U0[i0];                                   \
  Cy=V0[i1]-U0[i1];                                   \
  f=Ay*Bx-Ax*By;                                      \
  d=By*Cx-Bx*Cy;                                      \
  if((f>0 && d>=0 && d<=f) || (f<0 && d<=0 && d>=f))  \
  {                                                   \
    e=Ax*Cy-Ay*Cx;                                    \
    if(f>0)                                           \
    {                                                 \
      if(e>=0 && e<=f) return 1;                      \
    }                                                 \
    else                                              \
    {                                                 \
      if(e<=0 && e>=f) return 1;                      \
    }                                                 \
  }                                

#define EDGE_AGAINST_TRI_EDGES(V0,V1,U0,U1,U2) \
{                                              \
  mpq_class Ax,Ay,Bx,By,Cx,Cy,e,d,f;           \
  Ax=V1[i0]-V0[i0];                            \
  Ay=V1[i1]-V0[i1];                            \
  /* test edge U0,U1 against V0,V1 */          \
  EDGE_EDGE_TEST(V0,U0,U1);                    \
  /* test edge U1,U2 against V0,V1 */          \
  EDGE_EDGE_TEST(V0,U1,U2);                    \
  /* test edge U2,U1 against V0,V1 */          \
  EDGE_EDGE_TEST(V0,U2,U0);                    \
}

#define POINT_IN_TRI(V0,U0,U1,U2)           \
{                                           \
  mpq_class a,b,c,d0,d1,d2;                 \
  /* is T1 completly inside T2? */          \
  /* check if V0 is inside tri(U0,U1,U2) */ \
  a=U1[i1]-U0[i1];                          \
  b=-(U1[i0]-U0[i0]);                       \
  c=-a*U0[i0]-b*U0[i1];                     \
  d0=a*V0[i0]+b*V0[i1]+c;                   \
                                            \
  a=U2[i1]-U1[i1];                          \
  b=-(U2[i0]-U1[i0]);                       \
  c=-a*U1[i0]-b*U1[i1];                     \
  d1=a*V0[i0]+b*V0[i1]+c;                   \
                                            \
  a=U0[i1]-U2[i1];                          \
  b=-(U0[i0]-U2[i0]);                       \
  c=-a*U2[i0]-b*U2[i1];                     \
  d2=a*V0[i0]+b*V0[i1]+c;                   \
  if(d0*d1>0.0)                             \
  {                                         \
    if(d0*d2>0.0) return 1;                 \
  }                                         \
}

inline bool coplanar_tri_tri(const Vector3q& N, const Vector3q& V0, const Vector3q& V1, const Vector3q& V2,
                     const Vector3q& U0, const Vector3q& U1, const Vector3q& U2)
{
  Vector3q A;
  short i0,i1;
  /* first project onto an axis-aligned plane, that maximizes the area */
  /* of the triangles, compute indices: i0,i1. */
  A[0] = abs(N[0]);
  A[1] = abs(N[1]);
  A[2] = abs(N[2]);
  if(A[0] > A[1]) {
    if(A[0] > A[2]) {
      i0 = 1;      /* A[0] is greatest */
      i1 = 2;
    } else {
      i0 = 0;      /* A[2] is greatest */
      i1 = 1;
    }
  } else { /* A[0]<=A[1] */
    if(A[2] > A[1]) {
      i0 = 0;      /* A[2] is greatest */
      i1 = 1;                                           
    } else {
      i0 = 0;      /* A[1] is greatest */
      i1 = 2;
    }
  }               
                
  /* test all edges of triangle 1 against the edges of triangle 2 */
  EDGE_AGAINST_TRI_EDGES(V0, V1, U0, U1, U2);
  EDGE_AGAINST_TRI_EDGES(V1, V2, U0, U1, U2);
  EDGE_AGAINST_TRI_EDGES(V2, V0, U0, U1, U2);
                
  /* finally, test if tri1 is totally contained in tri2 or vice versa */
  POINT_IN_TRI(V0,U0,U1,U2);
  POINT_IN_TRI(U0,V0,V1,V2);

  return 0;
}


inline bool tri_tri_intersect(const Vector3q& V0, const Vector3q& V1, const Vector3q& V2, const Vector3q& N1, const Vector3q& U0, const Vector3q& U1, const Vector3q& U2, const Vector3q& N2)
{
  mpq_class d1, d2;
  mpq_class du0, du1, du2, dv0, dv1, dv2;
  Vector3q D;
  mpq_class isect1[2], isect2[2];
  mpq_class du0du1, du0du2, dv0dv1, dv0dv2;
  short index;
  mpq_class vp0,vp1,vp2;
  mpq_class up0,up1,up2;
  mpq_class b,c,max;

  /* compute plane equation of triangle(V0,V1,V2) */
  d1 = -N1.dot(V0);
  /* plane equation 1: N1.X+d1=0 */

  /* put U0,U1,U2 into plane equation 1 to compute signed distances to the plane*/
  du0 = N1.dot(U0) + d1;
  du1 = N1.dot(U1) + d1;
  du2 = N1.dot(U2) + d1;

  du0du1 = du0 * du1;
  du0du2 = du0 * du2;

  if(du0du1 > 0.0f && du0du2 > 0.0f) /* same sign on all of them + not equal 0 ? */
    return 0;                    /* no intersection occurs */

  /* compute plane of triangle (U0,U1,U2) */
  d2 = -N2.dot(U0);
  /* plane equation 2: N2.X+d2=0 */

  /* put V0,V1,V2 into plane equation 2 */
  dv0 = N2.dot(V0) + d2;
  dv1 = N2.dot(V1) + d2;
  dv2 = N2.dot(V2) + d2;

  dv0dv1 = dv0 * dv1;
  dv0dv2 = dv0 * dv2;
        
  if(dv0dv1 > 0.0f && dv0dv2 > 0.0f) /* same sign on all of them + not equal 0 ? */
    return 0;                    /* no intersection occurs */

  /* compute direction of intersection line */
  D = N1.cross(N2);

  /* compute and index to the largest component of D */
  index = 0;
  max = abs(D[0]);
  b = abs(D[1]);
  c = abs(D[2]);
  if(b > max) max=b, index=1;
  if(c > max) max=c, index=2;

  /* this is the simplified projection onto L*/
  vp0 = V0[index];
  vp1 = V1[index];
  vp2 = V2[index];

  up0 = U0[index];
  up1 = U1[index];
  up2 = U2[index];

  /* compute interval for triangle 1 */
  COMPUTE_INTERVALS(vp0,vp1,vp2,dv0,dv1,dv2,dv0dv1,dv0dv2,isect1[0],isect1[1]);

  /* compute interval for triangle 2 */
  COMPUTE_INTERVALS(up0,up1,up2,du0,du1,du2,du0du1,du0du2,isect2[0],isect2[1]);

  SORT(isect1[0],isect1[1]);
  SORT(isect2[0],isect2[1]);

  if(isect1[1]<isect2[0] || isect2[1]<isect1[0]) return 0;
  return 1;
}

inline bool tri_tri_intersect(const Vector3q& V0, const Vector3q& V1, const Vector3q& V2, const Vector3q& U0, const Vector3q& U1, const Vector3q& U2)
{
  Vector3q E1, E2;
  Vector3q N1, N2;
  E1 = V1 - V0;
  E2 = V2 - V0;
  N1 = E1.cross(E2);

  E1 = U1 - U0;
  E2 = U2 - U0;
  N2 = E1.cross(E2);
  return tri_tri_intersect(V0, V1, V2, N1, U0, U1, U2, N2);
}

