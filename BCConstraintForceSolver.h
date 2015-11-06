/* 
 This file is a part of BetterContactPlugin.
 
 Author: Shin'ichiro Nakaoka
 Author: Ryo Kikuuwe
 
 Copyright (c) 2007-2015 Shin'ichiro Nakaoka
 Copyright (c) 2014-2015 Ryo Kikuuwe
 Copyright (c) 2007-2015 National Institute of Advanced Industrial
                         Science and Technology (AIST)
 Copyright (c) 2014-2015 Kyushu University

 BetterContactPlugin is a plugin for better simulation of frictional contacts.
 
 BetterContactPlugin is a free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.
 BetterContactPlugin is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with BetterContactPlugin; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
 
 Contact: Ryo Kikuuwe, kikuuwe@ieee.org
*/
/*************************************************

ORIGINAL FILE: src/Body/ConstraintForceSolver.h

**************************************************/


#ifndef CNOID_BCCONSTRAINT_FORCE_SOLVER_H
#define CNOID_BCCONSTRAINT_FORCE_SOLVER_H

#include <cnoid/CollisionSeq>
//#include "exportdecl.h"

namespace cnoid
{
class Link;
class BCCFSImpl;
class WorldBase;
class CollisionDetector;
typedef boost::shared_ptr<CollisionDetector> CollisionDetectorPtr;

class CNOID_EXPORT BCConstraintForceSolver
{
    BCCFSImpl* impl;
		
public:
    BCConstraintForceSolver(WorldBase& world);
    ~BCConstraintForceSolver();
		
    void setCollisionDetector(CollisionDetectorPtr detector);
    CollisionDetectorPtr collisionDetector();

    void setFriction(double staticFriction, double slipFliction);
    double staticFriction() const;
    double slipFriction() const;

    void setContactCullingDistance(double thresh);
    double contactCullingDistance() const;
    
    void setContactCullingDepth(double depth);
    double contactCullingDepth();
    
    void setCoefficientOfRestitution(double epsilon);
    double coefficientOfRestitution() const;

    void setGaussSeidelErrorCriterion(double e);
    double gaussSeidelErrorCriterion();

    void setGaussSeidelMaxNumIterations(int n);
    int gaussSeidelMaxNumIterations();

    void setContactDepthCorrection(double depth, double velocityRatio);
    double contactCorrectionDepth();
    double contactCorrectionVelocityRatio();

    void set2Dmode(bool on);
    void enableConstraintForceOutput(bool on);


    void initialize(void);
    void solve();
    void clearExternalForces();

    CollisionLinkPairListPtr getCollisions();

#ifdef ENABLE_SIMULATION_PROFILING
    double getCollisionTime();
#endif

    double penaltyKpCoef();
    double penaltyKvCoef();
    double penaltySizeRatio();
    int    solverID();
    void setPenaltyKpCoef(double aKpCoef);
    void setPenaltyKvCoef(double aKvCoef);
    void setPenaltySizeRatio(double aSizeRatio);
    void setSolverID     (int arg);
};

};

#endif
