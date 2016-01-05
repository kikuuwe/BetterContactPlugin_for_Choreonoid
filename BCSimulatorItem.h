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

ORIGINAL FILE: src/BodyPlugin/AISTSimulatorItem.h

**************************************************/


#ifndef CNOID_BC_SIMULATOR_ITEM_H
#define CNOID_BC_SIMULATOR_ITEM_H

#include <cnoid/SimulatorItem>
#include <cnoid/EigenTypes>
//#include "exportdecl.h"

namespace cnoid {

class BCSimulatorItemImpl;
        
class CNOID_EXPORT BCSimulatorItem : public SimulatorItem
{
public:
    static void initializeClass(ExtensionManager* ext);
        
    BCSimulatorItem();
    BCSimulatorItem(const BCSimulatorItem& org);
    virtual ~BCSimulatorItem();

    enum DynamicsMode    { FORWARD_DYNAMICS = 0, HG_DYNAMICS, KINEMATICS, N_DYNAMICS_MODES };
    enum IntegrationMode { EULER_INTEGRATION = 0, RUNGE_KUTTA_INTEGRATION, N_INTEGRATION_MODES };
/*BC*/ enum SolverMode      { SLV_GAUSS_SEIDEL = 0, SLV_SICONOS, SLV_QMR, N_SOLVER_MODES };

    void setDynamicsMode(int mode);
    void setIntegrationMode(int mode);
/*BC*/ void setSolverMode(int mode);
    void setGravity(const Vector3& gravity);
    void setStaticFriction(double value);
    void setSlipFriction(double value);
    void setContactCullingDistance(double value);        
    void setContactCullingDepth(double value);        
    void setErrorCriterion(double value);        
    void setMaxNumIterations(int value);
    void setContactCorrectionDepth(double value);
    void setContactCorrectionVelocityRatio(double value);
    void setEpsilon(double epsilon);
    void set2Dmode(bool on);
    void setKinematicWalkingEnabled(bool on); 

    virtual void setForcedBodyPosition(BodyItem* bodyItem, const Position& T);
    virtual void clearForcedBodyPositions();
    
protected:
    virtual SimulationBody* createSimulationBody(Body* orgBody);
    virtual ControllerItem* createBodyMotionController(BodyItem* bodyItem, BodyMotionItem* bodyMotionItem);
    virtual bool initializeSimulation(const std::vector<SimulationBody*>& simBodies);
    virtual bool stepSimulation(const std::vector<SimulationBody*>& activeSimBodies);
    virtual void finalizeSimulation();
    virtual CollisionLinkPairListPtr getCollisions();
        
    virtual Item* doDuplicate() const;
    virtual void doPutProperties(PutPropertyFunction& putProperty);
    virtual bool store(Archive& archive);
    virtual bool restore(const Archive& archive);
#ifdef ENABLE_SIMULATION_PROFILING
    virtual void getProfilingNames(std::vector<std::string>& profilingNames);
    virtual void getProfilingTimes(std::vector<double>& profilingTimes);
#endif

private:
    BCSimulatorItemImpl* impl;
    friend class BCSimulatorItemImpl;
};

typedef ref_ptr<BCSimulatorItem> BCSimulatorItemPtr;
}

#endif
